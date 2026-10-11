"""Vina 批量对接引擎：单条对接 + 并行批量 + 断点续跑 + 失败隔离。

约定：
  - 每个条目独立目录 results/<entry_id>/result.json，产物存在且配置指纹一致即跳过；
  - 打分解析自 Vina 输出 PDBQT 的 REMARK VINA RESULT 行（与 automd.vina_scores 同法）；
  - 单条失败不中断批次，错误写入 result.json 的 error 字段。
"""

from __future__ import annotations

import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from .config import Config
from .utils import CommandError, config_hash, ensure_dir, get_logger, read_json, run_cmd, write_json

log = get_logger("docking")


class DockingError(RuntimeError):
    pass


# ---------------- 打分解析 ----------------
def vina_scores(pose_pdbqt: Path) -> list[float]:
    """解析 Vina 输出 PDBQT 中各 mode 的结合能。"""
    import re

    scores: list[float] = []
    for line in Path(pose_pdbqt).read_text(encoding="utf-8", errors="ignore").splitlines():
        m = re.match(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", line.strip())
        if m:
            scores.append(float(m.group(1)))
    return scores


# ---------------- 单条对接 ----------------
def dock_single(
    *,
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str | Path,
    center: list[float],
    size: list[float],
    out_dir: str | Path,
    exhaustiveness: int,
    num_modes: int,
    seed: int,
    tools,
    vina_threads: int = 1,
    log_path: Path | None = None,
    timeout_s: int | None = None,
) -> dict[str, Any]:
    """执行一次 Vina 对接，返回 {top_score, scores, pose_pdbqt, pose_sdf, elapsed_s}。"""
    out_dir = ensure_dir(Path(out_dir))
    poses_pdbqt = out_dir / "poses.pdbqt"
    pose_sdf = out_dir / "pose.sdf"

    cmd = [
        str(tools.vina),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cpu", str(max(1, int(vina_threads))),
        "--seed", str(seed),
        "--out", str(poses_pdbqt),
    ]
    t0 = time.time()
    try:
        run_cmd(cmd, cwd=out_dir, log=log_path or (out_dir / "vina.log"),
                timeout=int(timeout_s) if timeout_s and int(timeout_s) > 0 else None)
    except subprocess.TimeoutExpired as exc:
        raise DockingError(f"Vina 超时（>{timeout_s}s），疑似超大/高柔性配体，已终止: {ligand_pdbqt}") from exc
    scores = vina_scores(poses_pdbqt)
    if not scores:
        raise DockingError(f"Vina 未返回打分: {out_dir}")

    # 最优 pose 转 SDF（供 RMSD / 可视化）
    run_cmd([str(tools.mk_export), str(poses_pdbqt), "-s", str(pose_sdf)],
            cwd=out_dir, log=log_path or (out_dir / "export.log"))
    return {
        "top_score": scores[0],
        "scores": scores,
        "pose_pdbqt": str(poses_pdbqt.resolve()),
        "pose_sdf": str(pose_sdf.resolve()),
        "elapsed_s": round(time.time() - t0, 1),
    }


# ---------------- 批量对接 ----------------
def _dock_entry_job(args: tuple) -> dict[str, Any]:
    """线程池 worker：单条目对接（子进程调用 Vina，GIL 不受限）。"""
    (entry_id, ligand_pdbqt, receptor_pdbqt, center, size, out_dir,
     exh, num_modes, seed, tools_dict, threads, cfg_hash, timeout_s) = args
    from .config import ToolsConfig
    tools = ToolsConfig(**tools_dict)
    result_json = Path(out_dir) / "result.json"

    # 断点续跑：产物存在且配置指纹一致 -> 跳过
    if result_json.exists():
        try:
            old = read_json(result_json)
            if old.get("config_hash") == cfg_hash and old.get("status") in {"ok", "cached"}:
                return {**old, "status": "cached"}
        except Exception:
            pass

    try:
        res = dock_single(
            receptor_pdbqt=receptor_pdbqt, ligand_pdbqt=ligand_pdbqt,
            center=center, size=size, out_dir=out_dir,
            exhaustiveness=exh, num_modes=num_modes, seed=seed,
            tools=tools, vina_threads=threads, timeout_s=timeout_s,
        )
        payload = {"entry_id": entry_id, "status": "ok", "config_hash": cfg_hash, **res}
    except (CommandError, DockingError, subprocess.TimeoutExpired) as exc:
        payload = {"entry_id": entry_id, "status": "failed", "config_hash": cfg_hash,
                   "error": f"{type(exc).__name__}: {exc}"}
    write_json(result_json, payload)
    return payload


def run_batch(
    cfg: Config,
    entries: list[Any],
    receptor_pdbqt: Path,
    center: list[float],
    size: list[float],
    out_root: Path,
) -> list[dict[str, Any]]:
    """批量对接：并行 + 断点续跑 + 失败隔离，返回全部结果字典。"""
    assert cfg.docking is not None
    out_root = ensure_dir(out_root)
    cfg_hash = config_hash({
        "exhaustiveness": cfg.docking.exhaustiveness,
        "num_modes": cfg.docking.num_modes,
        "seed": cfg.docking.redock_seed,
        "center": center, "size": size,
    })

    jobs = []
    results: list[dict[str, Any]] = []
    n_filtered = 0
    assert cfg.docking is not None
    max_rot = int(getattr(cfg.docking, "max_rotatable_bonds", 0) or 0)
    max_mw = float(getattr(cfg.docking, "max_mol_weight", 0) or 0)
    timeout_s = int(getattr(cfg.docking, "ligand_timeout_s", 0) or 0)
    for entry in entries:
        if not entry.pdbqt_path or not Path(entry.pdbqt_path).exists():
            jobs.append(None)
            continue
        # 离群大分子预过滤：MW/可旋转键超限的直接跳过（提交前拦截，
        # 避免 Vina 按天卡死与接触面积虚高假阳性；阈值可在 target.yaml 调整）
        if max_rot or max_mw:
            try:
                rot = int(float(entry.row.get("rotatable_bonds", 0) or 0))
                mw = float(entry.row.get("mol_weight", 0) or 0)
            except ValueError:
                rot, mw = 0, 0.0
            reasons = []
            if max_rot and rot > max_rot:
                reasons.append(f"rotatable_bonds={rot}>{max_rot}")
            if max_mw and mw > max_mw:
                reasons.append(f"MW={mw:.0f}>{max_mw:.0f}")
            if reasons:
                n_filtered += 1
                payload = {"entry_id": entry.entry_id, "status": "skipped", "config_hash": cfg_hash,
                           "reason": "; ".join(reasons), "name": entry.name}
                write_json(Path(out_root) / entry.entry_id / "result.json", payload)
                results.append(payload)
                continue
        jobs.append((
            entry.entry_id, entry.pdbqt_path, str(receptor_pdbqt),
            center, size, str(out_root / entry.entry_id),
            int(cfg.docking.exhaustiveness), int(cfg.docking.num_modes),
            int(cfg.docking.redock_seed) + hash(entry.entry_id) % 10000,
            cfg.tools.as_dict(), cfg.runtime.vina_threads, cfg_hash, timeout_s,
        ))

    pending = [j for j in jobs if j is not None]
    skipped = len(jobs) - len(pending)
    if skipped:
        log.warning("跳过 %d 个缺少 PDBQT 的条目", skipped)
    if n_filtered:
        log.warning("预过滤 %d 个离群大分子（超 MW≤%s Da / 可旋转键≤%s 阈值，见各 result.json）",
                    n_filtered, max_mw or "∞", max_rot or "∞")

    workers = max(1, int(cfg.runtime.workers))
    log.info("批量对接启动: %d 个条目（workers=%d）", len(pending), workers)
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_dock_entry_job, job) for job in pending]
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 20 == 0 or done == len(pending):
                log.info("对接进度: %d / %d", done, len(pending))

    n_ok = sum(1 for r in results if r["status"] in {"ok", "cached"})
    n_fail = sum(1 for r in results if r["status"] == "failed")
    n_skip = sum(1 for r in results if r["status"] == "skipped")
    log.info("批量对接完成: 成功 %d（其中缓存复用 %d），失败 %d，预过滤跳过 %d",
             n_ok, sum(1 for r in results if r["status"] == "cached"), n_fail, n_skip)
    return results
