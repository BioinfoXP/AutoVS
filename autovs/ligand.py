"""配体准备：标准化（去盐）→ 质子化 → 3D 构象 → PDBQT。

方法学与 GeneDrug-AutoMD 的 normalize_ligand / enumerate_ligand_states 一致：
  - 去盐：rdkit MolStandardize.FragmentParent
  - 质子化：dimorphite-dl（pH 7.4 ± 容差），取前 max_variants 个状态
  - 3D：ETKDGv3 生成 + MMFF94 优化（失败回退 UFF）
  - PDBQT：Meeko mk_prepare_ligand.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .utils import ensure_dir, get_logger, run_cmd

log = get_logger("ligand")


class LigandError(RuntimeError):
    pass


# ---------------- 3D 构象 ----------------
def prepare_3d_molecule(mol, seed: int = 20261010):
    """生成单个低能 3D 构象（ETKDG + MMFF94，回退 UFF）。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise LigandError("3D 构象嵌入失败")
    try:
        result = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        if result != 0:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    return mol


# ---------------- 状态枚举（质子化） ----------------
def enumerate_states(smiles: str, prep_cfg, seed: int = 20261010) -> list[dict[str, Any]]:
    """按 pH 枚举质子化状态，返回 [{state_id, smiles, sdf_mol}]。"""
    from dimorphite_dl import protonate_smiles
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise LigandError(f"SMILES 不可读: {smiles}")
    if prep_cfg.fragment_parent:
        mol = rdMolStandardize.FragmentParent(mol)
    base = Chem.MolToSmiles(mol, isomericSmiles=True)

    variants = list(protonate_smiles(
        base,
        ph_min=float(prep_cfg.ph) - float(prep_cfg.ph_tolerance),
        ph_max=float(prep_cfg.ph) + float(prep_cfg.ph_tolerance),
        precision=float(prep_cfg.protonation_precision),
        max_variants=max(int(prep_cfg.max_variants) * 4, 4),
    ) or [base])

    states: list[dict[str, Any]] = []
    seen: set[str] = set()
    for smi in variants:
        pmol = Chem.MolFromSmiles(smi)
        if pmol is None:
            continue
        if int(prep_cfg.max_tautomers_per_state) <= 1:
            candidates = [pmol]
        else:
            candidates = list(rdMolStandardize.TautomerEnumerator().Enumerate(pmol))[
                : int(prep_cfg.max_tautomers_per_state)
            ]
        for cand in candidates:
            canonical = Chem.MolToSmiles(cand, isomericSmiles=True)
            if canonical in seen:
                continue
            seen.add(canonical)
            states.append({
                "state_id": f"s{len(states) + 1}",
                "smiles": canonical,
                "formal_charge": int(Chem.GetFormalCharge(cand)),
                "mol": cand,
            })
            if len(states) >= int(prep_cfg.max_variants):
                return states
    if not states:
        raise LigandError(f"未能生成可读的质子化状态: {smiles}")
    return states


# ---------------- 单条目产物（SDF + PDBQT） ----------------
def build_entry(
    entry_id: str,
    smiles: str,
    out_dir: Path,
    prep_cfg,
    tools_cfg,
    seed: int = 20261010,
) -> dict[str, Any]:
    """为单个化合物生成 3D SDF 与 PDBQT；返回状态字典（供索引写入）。

    out_dir 结构：out_dir/<shard>/<entry_id>.sdf|.pdbqt
    失败时返回 {entry_id, status: 'failed', error: ...}，不抛异常（失败隔离）。
    """
    from rdkit import Chem

    sdf_path = out_dir / f"{entry_id}.sdf"
    pdbqt_path = out_dir / f"{entry_id}.pdbqt"
    if prep_cfg.skip_existing and sdf_path.exists() and pdbqt_path.exists():
        return {"entry_id": entry_id, "status": "cached", "sdf_path": str(sdf_path), "pdbqt_path": str(pdbqt_path)}

    try:
        states = enumerate_states(smiles, prep_cfg, seed=seed)
        state = states[0]
        mol3d = prepare_3d_molecule(state["mol"], seed=seed)
        mol3d.SetProp("_Name", entry_id)
        mol3d.SetProp("entry_id", entry_id)
        mol3d.SetProp("smiles", state["smiles"])
        ensure_dir(out_dir)
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol3d)
        writer.close()

        run_cmd([str(tools_cfg.mk_prepare_ligand), "-i", str(sdf_path), "-o", str(pdbqt_path)],
                cwd=out_dir, log=out_dir / "prepare_ligand.log")
        if not pdbqt_path.exists():
            raise LigandError("mk_prepare_ligand 未生成 PDBQT")

        return {
            "entry_id": entry_id,
            "status": "ok",
            "state_id": state["state_id"],
            "smiles": state["smiles"],
            "formal_charge": state["formal_charge"],
            "sdf_path": str(sdf_path.resolve()),
            "pdbqt_path": str(pdbqt_path.resolve()),
        }
    except Exception as exc:  # 失败隔离：单分子失败不中断批次
        return {"entry_id": entry_id, "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
