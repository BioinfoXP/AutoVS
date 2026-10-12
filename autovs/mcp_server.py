"""AutoVS MCP Server：把虚拟筛选全流程暴露为 MCP 工具，供 AI 客户端（Trae/Claude/Cursor）驱动。

启动（stdio 传输，写进客户端 mcp 配置）：
    python -m autovs.mcp_server           # 或 autovs mcp

工具面（全流程）：
    doctor            环境/工具/配置体检
    library_stats     FDA 库规模统计
    library_search    按名称/适应症搜库
    prepare_receptor  受体清洗与加氢准备
    run_gate          共晶重对接质量门（PASS 才允许全库筛选）
    submit_dock       提交批量对接（后台进程，立即返回 job_id）
    dock_status       查询对接进度（支持 server 重启后按 target+tag 追查）
    query_ranking     查询排行榜（top_n / 阈值过滤）
    visualize         3D 口袋图 + 2D 相互作用图 + 对照图（Top N 或指定化合物）
    generate_report   Markdown 任务报告
    package           任务目录 zip 打包（自包含交付包，dock 链路末尾自动执行）

长任务设计：对接为小时级任务，submit_dock 以子进程启动 `autovs dock`，
job 映射持久化到 temp/mcp_jobs.json（server 重启不丢）；dock_status 直接数
result.json 给出真实进度，与提交方式无关。
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from .config import load_config
from .utils import get_logger

log = get_logger("mcp")
mcp = FastMCP("autovs")

BASE = Path(__file__).resolve().parent.parent
DEFAULT_FRAMEWORK = BASE / "config" / "framework.yaml"
DEFAULT_LIBRARY = BASE / "config" / "library.yaml"
PY = sys.executable                     # 当前解释器（autovs conda 环境）
_JOBS: dict[str, dict] = {}             # job_id -> {target_yaml, tag, proc, log}
JOB_FILE = BASE / "temp" / "mcp_jobs.json"


def _cfg(target_yaml: str):
    return load_config(DEFAULT_FRAMEWORK, DEFAULT_LIBRARY, target_yaml)


def _save_jobs() -> None:
    JOB_FILE.parent.mkdir(parents=True, exist_ok=True)
    JOB_FILE.write_text(json.dumps(
        {k: {"target_yaml": v["target_yaml"], "tag": v["tag"],
             "pid": v["proc"].pid if v.get("proc") else v.get("pid"),
             "log": str(v["log"])} for k, v in _JOBS.items()}, indent=2), encoding="utf-8")


def _load_jobs() -> None:
    if JOB_FILE.exists() and not _JOBS:
        for k, v in json.loads(JOB_FILE.read_text(encoding="utf-8")).items():
            _JOBS[k] = {**v, "proc": None}   # 重启后 proc 丢失，进度靠数文件


def _count_done(cfg, tag: str) -> tuple[int, int]:
    """（已完成 result.json 数, 库总数）——进度与提交方式解耦。"""
    from .library import LocalLibrary
    done = len(list((cfg.paths.temp_dir / "docking" / cfg.target.name / tag).glob("*/result.json"))) \
        if (cfg.paths.temp_dir / "docking" / cfg.target.name / tag).exists() else 0
    lib = LocalLibrary(cfg.library.local.index_csv)
    return done, len(list(lib.iter_entries()))


# ---------------------------------------------------------------- doctor / library
@mcp.tool()
def doctor() -> str:
    """环境与工具体检：Python、AutoDock Vina、PyMOL、Meeko、FDA 库索引是否可用。"""
    import subprocess as sp
    cfg = load_config(DEFAULT_FRAMEWORK, DEFAULT_LIBRARY, None)
    lines = [f"python: {sys.executable}"]
    for name, path in [("vina", cfg.tools.vina), ("pymol", cfg.tools.pymol),
                       ("mk_prepare_receptor", cfg.tools.mk_prepare_receptor)]:
        ok = path and Path(str(path)).exists()
        lines.append(f"{name}: {'OK ' + str(path) if ok else 'MISSING'}")
    try:
        from .library import LocalLibrary
        lib = LocalLibrary(cfg.library.local.index_csv)
        lines.append(f"fda_library: {len(list(lib.iter_entries()))} entries")
    except Exception as exc:                       # noqa: BLE001
        lines.append(f"fda_library: ERROR {exc}")
    r = sp.run([str(Path(cfg.tools.vina)), "--version"], capture_output=True, text=True) if cfg.tools.vina and Path(str(cfg.tools.vina)).exists() else None
    if r:
        lines.append(f"vina_version: {r.stdout.strip().splitlines()[0] if r.stdout else '?'}")
    return "\n".join(lines)


@mcp.tool()
def library_stats() -> str:
    """FDA 本地库统计：化合物总数、属性分布概览。"""
    from .library import LocalLibrary
    cfg = load_config(DEFAULT_FRAMEWORK, DEFAULT_LIBRARY, None)
    lib = LocalLibrary(cfg.library.local.index_csv)
    entries = list(lib.iter_entries())
    mws = [float(e.row["mol_weight"]) for e in entries if e.row.get("mol_weight")]
    return json.dumps({"total": len(entries),
                       "mol_weight_mean": round(sum(mws) / len(mws), 1) if mws else None}, ensure_ascii=False)


@mcp.tool()
def library_search(query: str, limit: int = 10) -> str:
    """按药名/entry_id 子串搜索 FDA 本地库，返回化合物列表（名称/ID/分子量/logP）。"""
    from .library import LocalLibrary
    cfg = load_config(DEFAULT_FRAMEWORK, DEFAULT_LIBRARY, None)
    lib = LocalLibrary(cfg.library.local.index_csv)
    q = query.strip().lower()
    hits = [e for e in lib.iter_entries()
            if q in (e.name or "").lower() or q in e.entry_id.lower()][: max(1, min(int(limit), 50))]
    return json.dumps([{"entry_id": e.entry_id, "name": e.name,
                        "mol_weight": e.row.get("mol_weight"), "logp": e.row.get("logp")} for e in hits],
                      ensure_ascii=False)


# ---------------------------------------------------------------- receptor / gate
@mcp.tool()
def prepare_receptor(target_yaml: str) -> str:
    """受体准备：清洗（去配体/水/altloc 归一）+ 加氢，输出供对接与相互作用检测的受体文件。"""
    from .receptor import prepare_receptor as prep
    cfg = _cfg(target_yaml)
    pdb = prep(cfg)
    return f"受体已就绪: {pdb}"


@mcp.tool()
def run_gate(target_yaml: str) -> str:
    """共晶配体重对接质量门：打分窗口内最小 RMSD < 2 Å 才 PASS；FAIL 时禁止全库筛选。"""
    from .receptor import run_gate as gate
    cfg = _cfg(target_yaml)
    result = gate(cfg)
    return json.dumps(result, ensure_ascii=False, default=str)


# ---------------------------------------------------------------- dock（长任务）
@mcp.tool()
def submit_dock(target_yaml: str, tag: str, limit: int | None = None) -> str:
    """提交批量对接为后台任务（断点续跑，自动跳过已完成化合物），立即返回 job_id。
    MW/可旋转键超限的离群大分子提交前自动预过滤（status=skipped），单配体超时/失败自动隔离不阻塞批次。
    limit 用于小规模冒烟。完成后自动汇总排行 + 可视化 + 报告（含排行 CSV 副本）+ zip 打包。"""
    _load_jobs()
    job_id = time.strftime("%Y%m%d_%H%M%S")
    log_path = BASE / "temp" / f"mcp_dock_{job_id}.log"
    cmd = [PY, "-m", "autovs.cli", "dock", "--target", target_yaml, "--tag", tag]
    if limit:
        cmd += ["--limit", str(int(limit))]
    proc = subprocess.Popen(cmd, cwd=str(BASE), stdout=open(log_path, "w"), stderr=subprocess.STDOUT)
    _JOBS[job_id] = {"target_yaml": target_yaml, "tag": tag, "proc": proc, "log": log_path}
    _save_jobs()
    return json.dumps({"job_id": job_id, "pid": proc.pid, "log": str(log_path),
                       "hint": "用 dock_status(job_id) 查进度"}, ensure_ascii=False)


@mcp.tool()
def dock_status(job_id: str = "", target_yaml: str = "", tag: str = "") -> str:
    """查询对接进度：完成数/总数、进程是否存活、最新日志行。job_id 或 target_yaml+tag 二选一。"""
    _load_jobs()
    info = _JOBS.get(job_id) if job_id else None
    if not info and job_id:
        return f"未知 job_id: {job_id}（历史任务可用 target_yaml+tag 查）"
    ty = info["target_yaml"] if info else target_yaml
    tg = info["tag"] if info else tag
    cfg = _cfg(ty)
    done, total = _count_done(cfg, tg)
    alive = None
    if info:
        if info.get("proc"):
            alive = info["proc"].poll() is None
        elif info.get("pid"):
            try:
                Path(f"/proc/{info['pid']}")     # linux: 重启后用 pid 探活
                alive = True
            except OSError:
                alive = False
    last_log = ""
    lp = Path(info["log"]) if info else None
    if lp and lp.exists():
        last_log = lp.read_text(encoding="utf-8").strip().splitlines()[-1][:200] if lp.read_text(encoding="utf-8").strip() else ""
    return json.dumps({"job_id": job_id or f"{cfg.target.name}/{tg}", "done": done, "total": total,
                       "running": alive, "last_log": last_log,
                       "ranking_ready": (cfg.paths.output_data_dir / cfg.target.name / f"{tg}_ranking.csv").exists()},
                      ensure_ascii=False)


# ---------------------------------------------------------------- ranking / visualize / report
@mcp.tool()
def query_ranking(target_yaml: str, tag: str, top_n: int = 20, below_threshold_only: bool = False) -> str:
    """查询排行榜：默认 Top N（rank/化合物/打分/分子量/logP/首次批准），可只看阈值命中。"""
    cfg = _cfg(target_yaml)
    path = cfg.paths.output_data_dir / cfg.target.name / f"{tag}_ranking.csv"
    if not path.exists():
        return f"排行榜不存在: {path}（先 submit_dock）"
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if below_threshold_only:
        th = float(cfg.docking.score_threshold)
        rows = [r for r in rows if r["top_score_kcal_mol"] and float(r["top_score_kcal_mol"]) <= th]
    keys = ["rank", "entry_id", "name", "top_score_kcal_mol", "mol_weight", "logp", "first_approval"]
    return json.dumps([{k: r.get(k) for k in keys} for r in rows[: max(1, min(int(top_n), 100))]],
                      ensure_ascii=False)


@mcp.tool()
def visualize(target_yaml: str, tag: str, top_n: int = 10, entries: str = "") -> str:
    """生成 Top N（或 entries 指定的化合物，逗号分隔药名/ID）的 3D 口袋图、2D 相互作用图与 3D|2D 对照图。"""
    from .visualize import run_visualize
    cfg = _cfg(target_yaml)
    ent = [e for e in entries.split(",") if e.strip()] or None
    summary = run_visualize(cfg, tag, top=int(top_n), entries=ent)
    out = cfg.paths.output_figure_dir / cfg.target.name / tag
    return json.dumps({"n_poses": summary["n_poses"], "output_dir": str(out),
                       "manifest": str(out / "MANIFEST.csv")}, ensure_ascii=False)


@mcp.tool()
def generate_report(target_yaml: str, tag: str, top_show: int = 20) -> str:
    """生成 Markdown 任务报告（概览/质量门/统计/Top 排行/图谱对照/解读注意），
    并把全量 rank 排行 CSV 副本输出到报告同目录（任务包自包含）。返回报告路径。"""
    from .report import generate_report as gen
    cfg = _cfg(target_yaml)
    path = gen(cfg, tag, top_show=int(top_show))
    csv_copy = path.parent / f"{tag}_ranking.csv"
    return f"任务报告: {path}（排行副本: {csv_copy}）"


@mcp.tool()
def package(target_yaml: str, tag: str, out: str = "") -> str:
    """把任务自包含目录（报告/排行副本/图谱/MANIFEST）压缩为 zip 交付包。
    默认输出 output_figure/<靶点>/<靶点>_<批次>_results_<日期>.zip；dock 完成链路末尾自动打包，本工具用于补打或自定义路径。"""
    from .packaging import package_task
    cfg = _cfg(target_yaml)
    r = package_task(cfg, tag, out=Path(out) if out else None)
    return json.dumps(r, ensure_ascii=False)


def main() -> None:
    """stdio 传输启动（AI 客户端以子进程方式拉起本 server）。"""
    mcp.run()


if __name__ == "__main__":
    main()
