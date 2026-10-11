"""结果汇总：对接结果 × 库索引 -> 排行榜 CSV + 打分分布图 + Top N 摘要。"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from .config import Config
from .library import LocalLibrary
from .utils import ensure_dir, get_logger, write_json

log = get_logger("aggregate")

RANK_FIELDS = [
    "rank", "entry_id", "lib_id", "name", "smiles", "inchi_key",
    "top_score_kcal_mol", "median_score_kcal_mol", "n_modes",
    "mol_weight", "logp", "hbd", "hba", "rotatable_bonds",
    "first_approval", "atc", "pose_sdf",
]


def aggregate_results(cfg: Config, results: list[dict[str, Any]], tag: str = "screen") -> dict[str, Any]:
    """汇总对接结果：过滤 -> 排序 -> 排行榜 CSV + 打分分布 PDF。"""
    assert cfg.docking is not None
    out_data = ensure_dir(cfg.paths.output_data_dir / cfg.target.name if cfg.target else cfg.paths.output_data_dir)
    out_fig = ensure_dir(cfg.paths.output_figure_dir / cfg.target.name if cfg.target else cfg.paths.output_figure_dir)

    lib = LocalLibrary(cfg.library.local.index_csv) if cfg.library else None
    rows = []
    for res in results:
        if res.get("status") not in {"ok", "cached"}:
            continue
        row: dict[str, Any] = {
            "entry_id": res["entry_id"],
            "top_score_kcal_mol": round(float(res["top_score"]), 2),
            "median_score_kcal_mol": round(sorted(res["scores"])[len(res["scores"]) // 2], 2),
            "n_modes": len(res["scores"]),
            "pose_sdf": res.get("pose_sdf", ""),
        }
        if lib is not None:
            entry = lib.get(res["entry_id"])
            row.update({
                "lib_id": entry.lib_id, "name": entry.name, "smiles": entry.smiles,
                "inchi_key": entry.inchi_key,
                "mol_weight": entry.row.get("mol_weight", ""), "logp": entry.row.get("logp", ""),
                "hbd": entry.row.get("hbd", ""), "hba": entry.row.get("hba", ""),
                "rotatable_bonds": entry.row.get("rotatable_bonds", ""),
                "first_approval": entry.row.get("first_approval", ""), "atc": entry.row.get("atc", ""),
            })
        rows.append(row)

    rows.sort(key=lambda r: r["top_score_kcal_mol"])
    threshold = float(cfg.docking.score_threshold)
    hits = [r for r in rows if r["top_score_kcal_mol"] <= threshold]

    rank_csv = out_data / f"{tag}_ranking.csv"
    with rank_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=RANK_FIELDS)
        writer.writeheader()
        for i, row in enumerate(rows, 1):
            writer.writerow({k: row.get(k, "") for k in RANK_FIELDS} | {"rank": i})

    fig_path = out_fig / f"{tag}_score_distribution.pdf"
    _plot_scores(rows, threshold, fig_path)

    summary = {
        "screen": tag,
        "n_docked": len(rows),
        "n_hits_below_threshold": len(hits),
        "n_failed": sum(1 for r in results if r.get("status") == "failed"),
        "n_skipped_outlier": sum(1 for r in results if r.get("status") == "skipped"),
        "score_threshold_kcal_mol": threshold,
        "best_score_kcal_mol": rows[0]["top_score_kcal_mol"] if rows else None,
        "top_n": [dict(rank=i, entry_id=r["entry_id"], name=r["name"],
                       top_score_kcal_mol=r["top_score_kcal_mol"])
                  for i, r in enumerate(rows[: int(cfg.docking.top_n)], 1)],
        "ranking_csv": str(rank_csv),
        "score_figure": str(fig_path),
    }
    write_json(out_data / f"{tag}_summary.json", summary)
    log.info("汇总完成: 对接 %d，阈值命中 %d -> %s", len(rows), len(hits), rank_csv)
    return summary


def _plot_scores(rows: list[dict[str, Any]], threshold: float, out_pdf: Path) -> None:
    """打分分布直方图（PDF 矢量输出）。"""
    if not rows:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scores = [r["top_score_kcal_mol"] for r in rows]
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.hist(scores, bins=40, color="#4c78a8", edgecolor="white", linewidth=0.4)
    ax.axvline(threshold, color="#e45756", linestyle="--", linewidth=1.2,
               label=f"threshold = {threshold} kcal/mol")
    ax.set_xlabel("Vina best score (kcal/mol)")
    ax.set_ylabel("Count")
    ax.set_title("FDA library docking score distribution", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
