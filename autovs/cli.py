"""AutoVS 命令行入口。

常用命令：
  autovs doctor                     环境 / 工具 / 配置体检
  autovs library fetch              FDA 库联网拉取（唯一联网步骤）
  autovs library build              FDA 库本地化构建（纯离线，可断点续跑）
  autovs library stats              本地库统计
  autovs library verify-offline     本地库离线可用性验证（禁网执行）
  autovs receptor prepare           受体准备（PDB 本地化 + PDBQT）
  autovs gate                       质量门：共晶配体重对接 RMSD
  autovs dock --limit 50            批量对接（可限量冒烟）+ 自动汇总
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import Config, ConfigError, load_config
from .utils import get_logger, setup_logging

log = get_logger("cli")

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FRAMEWORK = ROOT / "config" / "framework.yaml"
DEFAULT_LIBRARY = ROOT / "config" / "library.yaml"


def _load(args) -> Config:
    return load_config(args.framework, args.library, getattr(args, "target", None))


# ---------------- doctor ----------------
def cmd_doctor(args) -> int:
    cfg = _load(args)
    print(f"AutoVS {cfg.framework.get('version', '?')} 环境体检")
    print(f"  root        : {cfg.paths.root}")
    ok = True
    for name, path in cfg.tools.as_dict().items():
        exists = path and Path(path).exists()
        print(f"  tool {name:<22}: {'OK ' if exists else 'MISSING'} {path}")
        ok = ok and bool(exists)
    for d in (cfg.paths.library_dir, cfg.paths.temp_dir,
              cfg.paths.output_data_dir, cfg.paths.output_figure_dir):
        exists = d.exists()
        print(f"  dir  {str(d.name):<22}: {'OK ' if exists else 'NEW'} {d}")
    print(f"  workers     : {cfg.runtime.workers}")
    print("结论:", "全部就绪" if ok else "存在缺失工具，请检查 tools 配置")
    return 0 if ok else 1


# ---------------- library ----------------
def cmd_library(args) -> int:
    from . import library as libmod

    cfg = _load(args)
    if args.action == "fetch":
        libmod.library_fetch(cfg)
    elif args.action == "build":
        libmod.library_build(cfg, limit=args.limit)
    elif args.action == "stats":
        stats = libmod.LocalLibrary(cfg.library.local.index_csv).stats()
        for k, v in stats.items():
            print(f"  {k}: {v}")
    elif args.action == "verify-offline":
        return _verify_offline(cfg)
    return 0


def _verify_offline(cfg: Config) -> int:
    """禁网验证本地库：索引可读、产物文件齐全、抽样条目可直接用于对接。"""
    import socket

    from .library import LocalLibrary

    def _blocked(*a, **k):  # 禁网哨兵：任何网络调用立即失败
        raise RuntimeError("离线验证期间禁止网络访问")

    socket.socket = _blocked  # type: ignore[assignment]
    socket.create_connection = _blocked  # type: ignore[assignment]

    lib = LocalLibrary(cfg.library.local.index_csv)
    stats = lib.stats()
    print(f"[offline] 索引读取 OK: {stats['total']} 条，可对接 {stats['dockable']} 条")
    missing = 0
    checked = 0
    for entry in lib.sample(min(20, stats["dockable"]), seed=cfg.seed):
        checked += 1
        if not (entry.pdbqt_path and Path(entry.pdbqt_path).exists()
                and entry.sdf_path and Path(entry.sdf_path).exists()):
            missing += 1
    print(f"[offline] 抽样 {checked} 个条目产物完整性: 缺失 {missing}")
    print("[offline] 结论:", "本地库功能完整，可离线调用" if missing == 0 else "存在缺失产物")
    return 0 if missing == 0 else 1


# ---------------- receptor / gate ----------------
def cmd_receptor(args) -> int:
    from .receptor import prepare_reference, prepare_receptor

    cfg = _load(args)
    rec = prepare_receptor(cfg)
    print(f"受体 PDBQT: {rec}")
    if cfg.target and cfg.target.reference_ligand:
        ref = prepare_reference(cfg)
        print(f"参考配体  : {ref}")
    return 0


def cmd_gate(args) -> int:
    from .receptor import run_gate

    cfg = _load(args)
    result = run_gate(cfg)
    print(f"redocking RMSD : {result['redocking_rmsd']} Å（阈值 {result['rmsd_threshold']}）")
    print(f"对接盒         : center={result['box']['center']} size={result['box']['size']}")
    print("质量门         :", "PASS，可进入全库筛选" if result["passed"] else "FAIL，禁止进入批量对接")
    return 0 if result["passed"] else 2


# ---------------- dock ----------------
def cmd_dock(args) -> int:
    from .aggregate import aggregate_results
    from .docking import run_batch
    from .library import LocalLibrary
    from .receptor import compute_box, prepare_reference, prepare_receptor
    from .utils import read_json

    cfg = _load(args)
    work = cfg.paths.temp_dir / "receptor" / cfg.target.name
    gate_file = work / "gate_result.json"
    if not args.skip_gate_check:
        if not gate_file.exists():
            log.error("未找到质量门结果，请先运行: autovs gate --target %s", args.target)
            return 2
        gate = read_json(gate_file)
        if not gate.get("passed") and not gate.get("waived"):
            log.error("质量门未通过（RMSD=%s Å），禁止批量对接", gate.get("redocking_rmsd"))
            return 2
        if not gate.get("passed") and gate.get("waived"):
            log.warning("质量门 FAIL 已豁免放行：%s", gate.get("waive_reason", "（未注明原因）"))
        center, size = gate["box"]["center"], gate["box"]["size"]
    else:
        ref = prepare_reference(cfg)
        center, size, _ = compute_box(ref, cfg.target.box)

    rec = prepare_receptor(cfg)
    lib = LocalLibrary(cfg.library.local.index_csv)
    entries = lib.sample(args.limit, seed=cfg.seed) if args.limit else list(lib.iter_entries())
    log.info("本次对接条目数: %d", len(entries))

    results = run_batch(cfg, entries, rec, center, size,
                        cfg.paths.temp_dir / "docking" / cfg.target.name / (args.tag or "screen"))
    summary = aggregate_results(cfg, results, tag=args.tag or "screen")
    print(f"对接完成: {summary['n_docked']} 条，命中（≤{summary['score_threshold_kcal_mol']}）"
          f"{summary['n_hits_below_threshold']} 条"
          f"，失败 {summary.get('n_failed', 0)}，离群预过滤 {summary.get('n_skipped_outlier', 0)}")
    print(f"排行榜  : {summary['ranking_csv']}")

    # 自动可视化：Top 命中 3D 渲染 + 2D 相互作用图（失败不影响主流程）
    try:
        from .visualize import run_visualize
        viz = run_visualize(cfg, args.tag or "screen", top=min(int(cfg.docking.top_n), 10))
        print(f"可视化  : {viz['n_poses']} 个姿势 -> {cfg.paths.output_figure_dir / cfg.target.name / (args.tag or 'screen')}")
    except Exception as exc:
        log.warning("自动可视化跳过（%s: %s）；可稍后手动执行 autovs visualize", type(exc).__name__, exc)

    # 自动任务报告：Markdown（失败不影响主流程）
    try:
        from .report import generate_report
        rpt = generate_report(cfg, args.tag or "screen")
        print(f"任务报告: {rpt}")
    except Exception as exc:
        log.warning("自动报告跳过（%s: %s）；可稍后手动执行 autovs report", type(exc).__name__, exc)

    # 自动打包：任务自包含目录 -> zip 交付包（失败不影响主流程）
    try:
        from .packaging import package_task
        pkg = package_task(cfg, args.tag or "screen")
        print(f"任务打包: {pkg['zip']}（{pkg['n_files']} 文件, {pkg['size_mb']} MB）")
    except Exception as exc:
        log.warning("自动打包跳过（%s: %s）；可稍后手动执行 autovs package", type(exc).__name__, exc)
    return 0


# ---------------- visualize ----------------
def cmd_visualize(args) -> int:
    from .visualize import run_visualize

    cfg = _load(args)
    entries = [e for e in (args.entry or "").split(",") if e.strip()]
    viz = run_visualize(cfg, args.tag, top=int(args.top), entries=entries or None)
    for p in viz["3d_png"]:
        print(f"  3D: {p}")
    for item in viz["2d_items"]:
        print(f"  2D: rank {item['rank']} {item['name']} -> {item['pdf']}")
    return 0


# ---------------- package ----------------
def cmd_package(args) -> int:
    from .packaging import package_task

    cfg = _load(args)
    r = package_task(cfg, args.tag, out=Path(args.out) if args.out else None)
    print(f"任务打包: {r['zip']}（{r['n_files']} 文件, {r['size_mb']} MB）")
    return 0


# ---------------- report ----------------
def cmd_report(args) -> int:
    from .report import generate_report

    cfg = _load(args)
    rpt = generate_report(cfg, args.tag, top_show=int(args.top))
    print(f"任务报告: {rpt}")
    return 0


# ---------------- mcp ----------------
def cmd_mcp(args) -> int:
    from .mcp_server import main as mcp_main

    mcp_main()          # stdio 传输，客户端断开即退出
    return 0


# ---------------- main ----------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="autovs", description="AutoVS：FDA 药物库本地化批量虚拟筛选框架")
    parser.add_argument("--framework", default=str(DEFAULT_FRAMEWORK), help="framework.yaml 路径")
    parser.add_argument("--library", default=str(DEFAULT_LIBRARY), help="library.yaml 路径")
    parser.add_argument("--log-level", default=None, help="日志级别（默认取配置）")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("doctor", help="环境 / 工具 / 配置体检")

    p_lib = sub.add_parser("library", help="FDA 库本地化")
    p_lib.add_argument("action", choices=["fetch", "build", "stats", "verify-offline"])
    p_lib.add_argument("--limit", type=int, default=None, help="build 限量（调试用）")

    p_rec = sub.add_parser("receptor", help="受体准备")
    p_rec.add_argument("action", choices=["prepare"])
    p_rec.add_argument("--target", required=True, help="targets/*.yaml 路径")

    p_gate = sub.add_parser("gate", help="质量门（共晶配体重对接 RMSD）")
    p_gate.add_argument("--target", required=True, help="targets/*.yaml 路径")

    p_dock = sub.add_parser("dock", help="批量对接 + 汇总")
    p_dock.add_argument("--target", required=True, help="targets/*.yaml 路径")
    p_dock.add_argument("--limit", type=int, default=None, help="限量对接（冒烟测试）")
    p_dock.add_argument("--tag", default=None, help="结果批次标签")
    p_dock.add_argument("--skip-gate-check", action="store_true", help="跳过质量门检查（不推荐）")

    p_viz = sub.add_parser("visualize", help="Top 命中姿势 3D 渲染 + 2D 相互作用图")
    p_viz.add_argument("--target", required=True, help="targets/*.yaml 路径")
    p_viz.add_argument("--tag", default="screen", help="结果批次标签（与 dock 一致）")
    p_viz.add_argument("--top", type=int, default=10, help="取 Top N 姿势")
    p_viz.add_argument("--entry", default=None,
                       help="指定化合物可视化：entry_id 或药名子串，逗号分隔（给出时忽略 --top）")

    p_rpt = sub.add_parser("report", help="生成 Markdown 任务报告（排行/质量门/图谱索引）")
    p_rpt.add_argument("--target", required=True, help="targets/*.yaml 路径")
    p_rpt.add_argument("--tag", default="screen", help="结果批次标签（与 dock 一致）")
    p_rpt.add_argument("--top", type=int, default=20, help="报告中展示的 Top N 行数")

    p_pkg = sub.add_parser("package", help="任务自包含目录打包为 zip 交付包")
    p_pkg.add_argument("--target", required=True, help="targets/*.yaml 路径")
    p_pkg.add_argument("--tag", default="screen", help="结果批次标签（与 dock 一致）")
    p_pkg.add_argument("--out", default=None, help="zip 输出路径（默认 output_figure/<靶点>/<靶点>_<批次>_results_<日期>.zip）")

    sub.add_parser("mcp", help="启动 MCP server（stdio，供 AI 客户端接入虚拟筛选全流程）")

    args = parser.parse_args(argv)
    cfg_preview = None
    try:
        cfg_preview = load_config(args.framework, args.library, None)
        setup_logging(args.log_level or cfg_preview.runtime.log_level)
    except ConfigError as exc:
        print(f"配置错误: {exc}", file=sys.stderr)
        return 1

    handlers = {
        "doctor": cmd_doctor,
        "library": cmd_library,
        "receptor": cmd_receptor,
        "gate": cmd_gate,
        "dock": cmd_dock,
        "visualize": cmd_visualize,
        "report": cmd_report,
        "package": cmd_package,
        "mcp": cmd_mcp,
    }
    try:
        return handlers[args.command](args)
    except ConfigError as exc:
        print(f"配置错误: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        log.error("%s: %s", type(exc).__name__, exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
