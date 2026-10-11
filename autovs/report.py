"""任务报告生成：排行榜 + 质量门 + 可视化产物 -> 自包含 Markdown 报告。

输出 output_figure/<靶点>/<批次>/report.md（与图谱同目录，相对链接直接可读，
整包拷走后报告仍可离线浏览）。dock 汇总-可视化后自动生成，也可独立执行
`autovs report --target ... --tag ...`。
"""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

from .config import Config
from .utils import get_logger

log = get_logger("report")


def generate_report(cfg: Config, tag: str, top_show: int = 20) -> Path:
    """生成 Markdown 任务报告，返回报告路径。"""
    assert cfg.target is not None and cfg.docking is not None
    tgt = cfg.target.name
    data_dir = cfg.paths.output_data_dir / tgt
    fig_dir = cfg.paths.output_figure_dir / tgt
    task_dir = fig_dir / tag
    ranking_csv = data_dir / f"{tag}_ranking.csv"
    if not ranking_csv.exists():
        raise FileNotFoundError(f"排行榜不存在: {ranking_csv}（先执行 autovs dock）")
    with ranking_csv.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    threshold = float(cfg.docking.score_threshold)
    n_hits = sum(1 for r in rows if r["top_score_kcal_mol"] and float(r["top_score_kcal_mol"]) <= threshold)
    best = rows[0] if rows else None

    # ---- 全量排行 CSV 副本：报告目录自包含（整包拷走即含 rank 排行）----
    rank_copy = task_dir / ranking_csv.name
    shutil.copyfile(ranking_csv, rank_copy)

    # ---- 汇总统计（summary.json 由 aggregate 写出，含失败/预过滤计数）----
    summary_path = data_dir / f"{tag}_summary.json"
    summary: dict[str, Any] = {}
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))

    # ---- 质量门（gate_result.json 由 autovs gate 写出）----
    gate_path = cfg.paths.temp_dir / "receptor" / tgt / "gate_result.json"
    gate: dict[str, Any] = {}
    if gate_path.exists():
        gate = json.loads(gate_path.read_text(encoding="utf-8"))

    # ---- 可视化对照图（MANIFEST 由 autovs visualize 写出；缺失则该节自动跳过）----
    manifest_path = task_dir / "MANIFEST.csv"
    manifest: list[dict[str, str]] = []
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as fh:
            manifest = list(csv.DictReader(fh))

    L: list[str] = []
    ap = L.append
    ap(f"# AutoVS 虚拟筛选任务报告：{tgt}（批次 {tag}）")
    ap("")
    ap(f"> 生成时间：{datetime.now().strftime('%Y-%m-%d %H:%M')} ｜ 框架：AutoVS（AutoDock Vina 批量虚拟筛选）")
    ap("")

    # ---- 1. 任务概览 ----
    ap("## 1. 任务概览")
    ap("")
    ap("| 项目 | 值 |")
    ap("|---|---|")
    ap(f"| 靶点 | {tgt} |")
    ap(f"| 结构 | PDB {cfg.target.pdb_id or '—'}（{cfg.target.structure_path or '—'}）|")
    box = (gate.get("box") if gate else None) or {"center": cfg.target.box.center, "size": cfg.target.box.size}
    box_ctr = ", ".join(f"{v:.1f}" for v in box["center"])
    box_sz = ", ".join(f"{v:.0f}" for v in box["size"])
    ap(f"| 对接口袋 | 中心 ({box_ctr}) Å，盒子 {box_sz} Å{('，' + box['source']) if box.get('source') else ''} |")
    ap(f"| 对接参数 | exhaustiveness={cfg.docking.exhaustiveness}, num_modes={cfg.docking.num_modes}, seed={cfg.docking.redock_seed} |")
    ap(f"| 筛选库 | FDA 已上市药物（ChEMBL max_phase=4 本地化索引）|")
    if gate:
        ap(f"| 质量门（共晶重对接）| **{'PASS' if gate.get('passed') else 'FAIL'}**：打分窗口 {gate.get('score_window')} kcal/mol 内最小 RMSD = "
           f"**{gate.get('redocking_rmsd')} Å**（阈值 < {gate.get('rmsd_threshold')} Å；top-1 pose RMSD {gate.get('top_pose_rmsd')} Å）|")
    ap(f"| 命中阈值 | {threshold} kcal/mol |")
    ap("")

    # ---- 2. 结果统计 ----
    ap("## 2. 结果统计")
    ap("")
    ap(f"- 实际对接：**{len(rows)}** 个化合物")
    ap(f"- 阈值命中（≤ {threshold} kcal/mol）：**{n_hits}** 个（{len(rows) and n_hits / len(rows) * 100:.1f}%）")
    if best:
        ap(f"- 最佳打分：**{best['top_score_kcal_mol']} kcal/mol**（{best['name'] or best['entry_id']}，rank {best['rank']}）")
    if summary.get("n_skipped_outlier"):
        ap(f"- 离群预过滤：{summary['n_skipped_outlier']} 个（MW ≤ {cfg.docking.max_mol_weight} Da 且可旋转键 ≤ {cfg.docking.max_rotatable_bonds} 之外的大分子，防止 Vina 卡死与虚高假阳性）")
    if summary.get("n_failed"):
        ap(f"- 对接失败：{summary['n_failed']} 个（已失败隔离，详见各 result.json 的 error 字段）")
    ap(f"- 全量排行榜：[`{ranking_csv.name}`]({ranking_csv.name})（{len(rows)} 条按 rank 排序，与图谱同包）")
    ap("")
    dist_rel = f"../{tag}_score_distribution.pdf"
    if (fig_dir / f"{tag}_score_distribution.pdf").exists():
        ap(f"![打分分布]({dist_rel})")
        ap("")
        ap("*全库打分分布（PDF 矢量图），红色虚线为命中阈值。*")
        ap("")

    # ---- 3. Top 排行 ----
    ap(f"## 3. Top {min(top_show, len(rows))} 排行")
    ap("")
    ap("| rank | 化合物 | entry_id | 打分 (kcal/mol) | 中位打分 | MW | logP | 首次批准 |")
    ap("|---:|---|---|---:|---:|---:|---:|---|")
    for r in rows[:top_show]:
        ap(f"| {r['rank']} | {r['name'] or '—'} | {r['entry_id']} | {r['top_score_kcal_mol']} | "
           f"{r['median_score_kcal_mol']} | {r['mol_weight']} | {r['logp']} | {r['first_approval'] or '—'} |")
    ap("")

    # ---- 4. Top 命中结合模式 ----
    if manifest:
        ap("## 4. Top 命中结合模式（3D 口袋 | 2D 相互作用 对照图）")
        ap("")
        for m in manifest:
            pair = m.get("fig_pair", "")
            pair_name = Path(pair).name if pair else ""
            ap(f"### rank {m['rank']} · {m['name'] or m['entry_id']}（{m.get('score_kcal_mol', '')} kcal/mol）")
            ap("")
            if pair_name and (task_dir / pair_name).exists():
                ap(f"![rank {m['rank']} 3D|2D 对照]({pair_name})")
                ap("")
            inter = m.get("interactions", "")
            if inter:
                ap(f"**相互作用**：{inter}")
                ap("")

    # ---- 5. 方法与解读注意 ----
    ap("## 5. 方法与解读注意事项")
    ap("")
    ap("- 打分来自 AutoDock Vina 半经验力场，仅用于富集排序；**绝对值不等于结合自由能**，命中需 MD/MMPBSA 或实验复核。")
    ap("- Vina 对大环、多芳环、高柔性分子存在系统性打分偏好；Top 命中中该类骨架占比偏高时，建议按分子量/柔性分箱后再比较。")
    ap("- 2D 相互作用为几何判据推断（氢键 H···A ≤2.6 Å 且 D-H···A >110°；疏水 C···C ≤4.0 Å；盐桥 ≤4.0 Å；π 堆积环心 ≤5.5 Å；阳离子-π ≤5.0 Å），供快速检视，不构成定量结论。")
    ap("- 3D 图为正交投影口袋特写：天蓝半透明 cartoon、配体与口袋残基 stick、绿色虚线为氢键；标签经像素级标定，与残基一一对应。")
    ap(f"- 完整排行榜见本包 `{ranking_csv.name}`（{len(rows)} 条）；图谱索引见 `MANIFEST.csv`。")
    ap("")

    out = task_dir / "report.md"
    out.write_text("\n".join(L), encoding="utf-8")
    log.info("任务报告已生成: %s", out)
    return out
