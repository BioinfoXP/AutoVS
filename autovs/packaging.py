"""任务打包：把 output_figure/<靶点>/<批次>/ 任务自包含目录压缩为 zip 交付包。

zip 位于任务目录同级（output_figure/<靶点>/ 下），命名 <靶点>_<批次>_results_<日期>.zip；
包内保留 <批次>/ 一级目录，解压即还原任务目录（报告/排行副本/MANIFEST/图谱）。
dock 全链路末尾自动调用，也可 `autovs package` 手动执行。
"""

from __future__ import annotations

import shutil
import time
import zipfile
from pathlib import Path
from typing import Any

from .config import Config
from .utils import ensure_dir, get_logger

log = get_logger("package")


def package_task(cfg: Config, tag: str, out: Path | None = None) -> dict[str, Any]:
    """打包任务目录为 zip，返回 {zip, n_files, size_mb, task_dir}。

    out 为 None 时输出到 output_figure/<靶点>/<靶点>_<批次>_results_<YYYYMMDD>.zip。
    重复打包覆盖旧包（幂等，重跑可视化/报告后再打包即得到最新交付包）。
    """
    tag = tag or "screen"
    task_dir = cfg.paths.output_figure_dir / cfg.target.name / tag
    if not task_dir.is_dir():
        raise FileNotFoundError(f"任务目录不存在: {task_dir}（先执行 dock/visualize/report）")
    date = time.strftime("%Y%m%d")
    zip_path = Path(out) if out else \
        ensure_dir(cfg.paths.output_figure_dir / cfg.target.name) / f"{cfg.target.name}_{tag}_results_{date}.zip"
    zip_path = Path(str(zip_path))
    if zip_path.exists():
        zip_path.unlink()
    base = shutil.make_archive(str(zip_path.with_suffix("")), "zip",
                               root_dir=str(task_dir.parent), base_dir=task_dir.name)
    with zipfile.ZipFile(base) as zf:
        n_files = sum(1 for info in zf.infolist() if not info.is_dir())
    size_mb = Path(base).stat().st_size / 1e6
    log.info("任务打包完成: %s（%d 文件, %.1f MB）", base, n_files, size_mb)
    return {"zip": base, "n_files": n_files, "size_mb": round(size_mb, 1), "task_dir": str(task_dir)}
