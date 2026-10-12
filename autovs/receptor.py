"""受体准备与质量门：结构获取 → 受体 PDBQT → 对接盒 → 共晶配体重对接 RMSD。

质量门（Gate）是全库筛选的放行条件：
  共晶配体重对接 frame RMSD < gate.rmsd_threshold 才允许进入批量对接。
  RMSD 计算沿用 GeneDrug-AutoMD 的 receptor_frame_rmsd 方法学
  （rdkit CalcRMS，对称性感知、不做叠合）。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import Config
from .utils import ensure_dir, get_logger, run_cmd, write_json

log = get_logger("receptor")


class ReceptorError(RuntimeError):
    pass


# ---------------- 结构获取 ----------------
def fetch_pdb(pdb_id: str, out_dir: Path, timeout: int = 120) -> Path:
    """下载 PDB 结构到本地（仅此步骤联网；之后全部离线）。"""
    import requests

    out = ensure_dir(out_dir) / f"{pdb_id.lower()}.pdb"
    if out.exists():
        return out
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    resp = requests.get(url, timeout=timeout)
    if resp.status_code != 200:
        raise ReceptorError(f"PDB 下载失败（HTTP {resp.status_code}）: {url}")
    out.write_bytes(resp.content)
    log.info("PDB 结构已本地化: %s", out)
    return out


# ---------------- 共晶配体提取 ----------------
def extract_reference_ligand(pdb_path: Path, resname: str, chain: str, out_dir: Path) -> Path:
    """从 PDB 提取共晶配体结合构象，用 RCSB CCD 模板补键级，输出 bound SDF。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    serials: set[str] = set()
    atom_lines: list[str] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == resname and line[21:22].strip() == chain \
                and (len(line) < 17 or line[16:17] in {" ", "A"}):
            atom_lines.append(line)
            serials.add(line[6:11].strip())
    if not atom_lines:
        raise ReceptorError(f"PDB 中未找到共晶配体 {resname}（chain {chain}）")
    conect = [l for l in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines()
              if l.startswith("CONECT") and (set(l[6:].split()) & serials)]
    bound = Chem.MolFromPDBBlock("\n".join(atom_lines + conect + ["END", ""]),
                                 removeHs=False, sanitize=False)
    if bound is None:
        raise ReceptorError("RDKit 无法读取共晶配体坐标")

    # CCD 模板补键级（与 automd.prepare_reference_ligand 同法）
    # 硬约束：模板不可用必须报错——静默跳过会让距离猜键的垃圾键级进入对接（环变可旋转、
    # 芳香性丢失、加氢坐标出原点），实测导致质量门 RMSD 全废
    import requests
    ccd = ensure_dir(out_dir) / f"{resname}_ccd.sdf"
    if not ccd.exists():
        try:
            resp = requests.get(f"https://files.rcsb.org/ligands/download/{resname}.sdf", timeout=60)
            if resp.status_code == 200:
                ccd.write_bytes(resp.content)
        except requests.RequestException:
            pass
    if not ccd.exists():
        raise ReceptorError(
            f"CCD 模板 {resname} 不可用（下载失败且无本地缓存）。"
            f"请手动放置 {resname}.sdf 至 {ccd} 后重试")
    template = next((m for m in Chem.SDMolSupplier(str(ccd), removeHs=False) if m), None)
    if template is None:
        raise ReceptorError(f"CCD 模板文件无效: {ccd}")
    try:
        bound = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(template), Chem.RemoveHs(bound))
        Chem.SanitizeMol(bound)
    except Exception as exc:
        raise ReceptorError(f"CCD 模板键级分配失败: {exc}") from exc

    out = out_dir / "reference_ligand_bound.sdf"
    bound_h = Chem.AddHs(bound, addCoords=True)  # Meeko 要求显式氢
    # RDKit 隐性坑：AddHs(addCoords=True) 偶发个别 H 落在坐标原点，
    # 幽灵氢进入 pdbqt 会破坏 Vina 内坐标（实测致质量门 RMSD 全废），此处强制重排
    conf_h = bound_h.GetConformer()
    for idx in range(bound_h.GetNumAtoms()):
        p = conf_h.GetAtomPosition(idx)
        if abs(p.x) + abs(p.y) + abs(p.z) < 1e-6:
            nbr = bound_h.GetAtomWithIdx(idx).GetNeighbors()[0].GetIdx()
            Chem.rdmolops.SetTerminalAtomCoords(bound_h, idx, nbr)
            log.warning("共晶配体加氢坐标异常(原点)，已重排原子 %d", idx)
    bound_h.SetProp("_Name", f"{resname}_bound")
    writer = Chem.SDWriter(str(out))
    writer.write(bound_h)
    writer.close()
    log.info("共晶配体已提取: %s（%d 重原子）", out, bound.GetNumHeavyAtoms())
    return out


# ---------------- 对接盒 ----------------
def compute_box(ref_sdf: Path, box_cfg) -> tuple[list[float], list[float], str]:
    """由共晶配体坐标计算对接盒（span + padding，限幅 [min_box, max_box]）。"""
    from rdkit import Chem

    if box_cfg.mode == "manual":
        return list(map(float, box_cfg.center)), list(map(float, box_cfg.size)), "manual"

    mol = next((m for m in Chem.SDMolSupplier(str(ref_sdf), removeHs=True) if m), None)
    if mol is None:
        raise ReceptorError(f"无法读取配体坐标: {ref_sdf}")
    conf = mol.GetConformer()
    xs, ys, zs = [], [], []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(p.x); ys.append(p.y); zs.append(p.z)
    center = [(min(c) + max(c)) / 2.0 for c in (xs, ys, zs)]
    span = [max(c) - min(c) for c in (xs, ys, zs)]
    size = [max(float(box_cfg.min_box), min(float(box_cfg.max_box), s + 2.0 * float(box_cfg.padding)))
            for s in span]
    return center, size, f"co_crystal:{ref_sdf.stem}"


# ---------------- 受体准备 ----------------
def clean_receptor_pdb(source: Path, out: Path) -> Path:
    """受体预处理：仅保留蛋白 ATOM 记录（去配体/辅因子/离子/结晶添加剂/水），altloc 归一。

    共晶配体由 extract_reference_ligand 单独从原始 PDB 提取，故此处全部丢弃 HETATM。
    """
    kept: list[str] = []
    for line in source.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ATOM  "):
            if line[16:17] not in {" ", "A"}:  # 丢弃非首选 altloc
                continue
            line = line[:16] + " " + line[17:]  # altloc 位归一
            kept.append(line)
    out.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return out


def prepare_receptor(cfg: Config) -> Path:
    """受体准备：结构本地化 → 蛋白化清洗 → receptor.pdbqt（产物缓存）。"""
    assert cfg.target is not None
    work = ensure_dir(cfg.paths.temp_dir / "receptor" / cfg.target.name)
    pdb = cfg.target.structure_path or fetch_pdb(cfg.target.pdb_id, work,
                                                 timeout=cfg.runtime.request_timeout)
    clean_pdb = work / "receptor_clean.pdb"
    if not clean_pdb.exists():
        clean_receptor_pdb(pdb, clean_pdb)

    rec_pdbqt = work / "receptor.pdbqt"
    if not rec_pdbqt.exists():
        base_cmd = [str(cfg.tools.mk_prepare_receptor), "-i", str(clean_pdb),
                    "-o", str(work / "receptor"), "-p", "--default_altloc", "A"]
        try:
            run_cmd(base_cmd, cwd=work, log=work / "meeko_receptor.log")
        except Exception:
            # 兜底：自动剔除无法匹配模板的残基（风险已记录在日志中）
            run_cmd(base_cmd + ["--allow_bad_res"], cwd=work, log=work / "meeko_receptor.log")
        if not rec_pdbqt.exists():
            alternatives = list(work.glob("receptor*.pdbqt"))
            if alternatives:
                alternatives[0].rename(rec_pdbqt)
        if not rec_pdbqt.exists():
            raise ReceptorError("mk_prepare_receptor 未生成 receptor.pdbqt")
    log.info("受体已准备: %s", rec_pdbqt)
    return rec_pdbqt


def prepare_reference(cfg: Config) -> Path:
    """提取共晶配体结合构象（用于盒子定义与质量门）。"""
    assert cfg.target is not None and cfg.target.reference_ligand
    work = ensure_dir(cfg.paths.temp_dir / "receptor" / cfg.target.name)
    pdb = cfg.target.structure_path or fetch_pdb(cfg.target.pdb_id, work,
                                                 timeout=cfg.runtime.request_timeout)
    ref = cfg.target.reference_ligand
    out = work / "reference_ligand_bound.sdf"
    if out.exists():
        return out
    return extract_reference_ligand(pdb, str(ref.get("resname")), str(ref.get("chain", "A")), work)


# ---------------- 质量门：共晶配体重对接 ----------------
def frame_rmsd(reference_sdf: Path, pose_sdf: Path) -> float:
    """对称性感知的受体坐标系 RMSD（不叠合），与 automd.receptor_frame_rmsd 同法。"""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    reference = next((m for m in Chem.SDMolSupplier(str(reference_sdf), removeHs=True) if m), None)
    pose = next((m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=True) if m), None)
    if reference is None or pose is None:
        raise ReceptorError("RMSD 计算失败：结构不可读")
    if reference.GetNumAtoms() != pose.GetNumAtoms():
        raise ReceptorError(
            f"重原子数不一致（晶体 {reference.GetNumAtoms()} vs 重对接 {pose.GetNumAtoms()}），"
            f"请核对参考配体与配体状态设置"
        )
    return float(rdMolAlign.CalcRMS(pose, reference, maxMatches=10000))


def evaluate_redocking(
    ref_sdf: Path, poses_sdf: Path, scores: list[float],
    rmsd_threshold: float, score_window: float,
) -> dict[str, Any]:
    """逐 pose 计算受体坐标系 RMSD，并按打分窗口判据给出质量门结论。

    判据（行业通用 redocking 标准）：打分 ≤ 最优 + score_window 的 pose 中，
    存在 RMSD < rmsd_threshold 即认为对接设置能复现晶体学结合模式。
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    ref = next((m for m in Chem.SDMolSupplier(str(ref_sdf), removeHs=True) if m), None)
    poses = [m for m in Chem.SDMolSupplier(str(poses_sdf), removeHs=True) if m]
    if ref is None or not poses:
        raise ReceptorError("RMSD 评估失败：结构不可读")
    rows = []
    for i, (pose, score) in enumerate(zip(poses, scores), 1):
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            raise ReceptorError("重原子数不一致，请核对参考配体与配体状态设置")
        rows.append({"mode": i, "score_kcal_mol": round(float(score), 3),
                     "rmsd": round(float(rdMolAlign.CalcRMS(pose, ref, maxMatches=10000)), 3)})
    best_score = min(r["score_kcal_mol"] for r in rows)
    in_window = [r for r in rows if r["score_kcal_mol"] <= best_score + float(score_window)]
    min_rmsd_in_window = min(r["rmsd"] for r in in_window)
    return {
        "per_mode": rows,
        "top_pose_rmsd": rows[0]["rmsd"],
        "min_rmsd_in_window": min_rmsd_in_window,
        "score_window": score_window,
        "passed": min_rmsd_in_window < float(rmsd_threshold),
    }


def run_gate(cfg: Config) -> dict[str, Any]:
    """执行质量门：参考配体重对接并评估晶体姿势复现，产出 gate_result.json。"""
    assert cfg.target is not None and cfg.gate is not None and cfg.docking is not None
    from .docking import dock_single

    work = ensure_dir(cfg.paths.temp_dir / "receptor" / cfg.target.name)
    rec_pdbqt = prepare_receptor(cfg)
    ref_sdf = prepare_reference(cfg)
    center, size, source = compute_box(ref_sdf, cfg.target.box)

    # 参考配体转 PDBQT（与库内配体同一条准备链）
    ref_pdbqt = work / "reference_ligand.pdbqt"
    if not ref_pdbqt.exists():
        run_cmd([str(cfg.tools.mk_prepare_ligand), "-i", str(ref_sdf), "-o", str(ref_pdbqt)],
                cwd=work, log=work / "meeko_ref_ligand.log")

    gate_dir = ensure_dir(work / "gate_redocking")
    result = dock_single(
        receptor_pdbqt=rec_pdbqt, ligand_pdbqt=ref_pdbqt,
        center=center, size=size, out_dir=gate_dir,
        exhaustiveness=int(cfg.gate.exhaustiveness), num_modes=int(cfg.gate.num_modes),
        seed=int(cfg.docking.redock_seed), tools=cfg.tools, vina_threads=cfg.runtime.vina_threads,
    )
    ev = evaluate_redocking(
        ref_sdf, Path(result["pose_sdf"]), result["scores"],
        float(cfg.gate.rmsd_threshold), float(cfg.gate.score_window),
    )
    payload = {
        "target": cfg.target.name,
        "box": {"center": center, "size": size, "source": source},
        "reference_ligand": str(ref_sdf),
        "top_score_kcal_mol": result["top_score"],
        "redocking_rmsd": ev["min_rmsd_in_window"],
        "top_pose_rmsd": ev["top_pose_rmsd"],
        "rmsd_threshold": cfg.gate.rmsd_threshold,
        "score_window": cfg.gate.score_window,
        "per_mode": ev["per_mode"],
        "passed": ev["passed"],
    }
    write_json(work / "gate_result.json", payload)
    log.info("质量门结果: 打分窗口内最小 RMSD=%.3f Å（阈值 %.1f，窗口 %.1f）-> %s",
             ev["min_rmsd_in_window"], cfg.gate.rmsd_threshold, cfg.gate.score_window,
             "PASS，允许全库筛选" if ev["passed"] else "FAIL，禁止进入批量对接")
    return payload
