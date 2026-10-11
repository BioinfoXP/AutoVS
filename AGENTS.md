# AutoVS 项目上下文备份（AI 交接文档）

> 本文档供 AI 助手/协作者快速恢复项目上下文。覆盖：环境、路径、MCP、GitHub、运行命令、当前状态、领域硬约束。更新日期：2026-10-11。

## 1. 项目定位

AutoVS：靶点结构 → FDA 全库虚拟筛选 → 出版级 3D/2D 可视化 → Markdown 任务报告的**全自动化流水线**（断点续跑、失败隔离、离群大分子自动防护）。宿主目录 `/home/xiap/bio_project/for_others/AutoMD/`，仓库根 `AutoVS/`。

姊妹目录：`GeneDrug-AutoMD-Skill/`（GROMACS MD 技能包，备将来命中分子动力学验证用，与 AutoVS 无代码耦合）。

## 2. 环境配置（关键！）

| 项 | 值 |
|---|---|
| Python 解释器 | `/home/xiap/opt/gene-drug-automd/env/bin/python`（conda 环境 `gene-drug-automd`） |
| CLI 入口 | `/home/xiap/opt/gene-drug-automd/env/bin/autovs` |
| 系统默认环境 | **没有** gene-drug-automd——直接 `python`/`pytest` 会失败，必须用上表绝对路径或先 `conda activate gene-drug-automd` |
| 新机复现 | `conda env create -f environment.yml`（conda-forge + bioconda，含 vina/meeko/openbabel/fpocket/rdkit）；或 pip：`python3.11 -m venv .venv && pip install -e ".[dev,mcp]"`（外部工具自备并登记 config） |
| 本体安装 | editable 模式（`pip install -e .`），改代码即生效，无需重装 |
| 外部工具 | vina/meeko(mk_prepare_receptor·mk_prepare_ligand)/obabel/fpocket/pymol，路径登记在 `config/framework.yaml` 的 `tools:` 段；体检命令 `autovs doctor` |

## 3. 路径速查（相对 `AutoVS/`）

```
autovs/                  # 包本体：config.py(docker三层配置) library.py fetching.py
                         #   receptor.py gating.py docking.py aggregate.py
                         #   visualize.py(3D+2D渲染) report.py mcp_server.py cli.py utils.py
config/framework.yaml    # 全局框架配置（工具路径/并行度/目录）
config/libraries/*.yaml  # 数据库配置（fda 默认：ChEMBL max_phase=4）
config/targets/*.yaml    # 靶点配置（PDB/口袋/对接参数），如 hsp90b1_7ull.yaml
data/fda_library/index/fda_library_index.csv   # FDA 库本地索引（3,311 条）
temp/                    # 缓存：receptor/<靶点>/gate_result.json、docking/<靶点>/<批次>/、mcp_jobs.json
output_data/<靶点>/      # full_ranking.csv / full_summary.json（对接产物）
output_figure/<靶点>/<批次>/   # 任务自包含目录，整包可拷走：
                         #   report.md + full_ranking.csv(全量排行副本) + MANIFEST.csv
                         #   + pose_R{rank}_{药名}_{3d|2d|pair}.png + top_overview.png
tests/                   # pytest，15 项；运行：绝对路径 python -m pytest tests/ -q
docs/                    # framework.md（架构）mcp_guide.md（MCP 使用）
```

## 4. MCP

- 启动：`autovs mcp` 或 `python -m autovs.mcp_server`（stdio 传输）
- 依赖 pin：`mcp>=1.30,<2`（2.x 改名 MCPServer 不兼容）
- 10 工具：doctor / library_stats / library_search / prepare_receptor / run_gate / submit_dock / dock_status / query_ranking / visualize / generate_report
- 长任务模式：submit_dock 子进程提交 + job 映射持久化 `temp/mcp_jobs.json`（server 重启可追查）+ dock_status 读 result.json 真实进度
- 客户端接入 JSON 见 README「MCP 配置」节；依赖与筛选运行环境相互独立

## 5. GitHub 仓库

| 项 | 值 |
|---|---|
| 仓库 | `BioinfoXP/AutoVS`（branch main） |
| remote | `ssh://git@ssh.github.com:443/BioinfoXP/AutoVS.git` |
| 推送通道 | **必须走 443**——本机 22 端口连 github.com 超时（`ssh -p 443 -T git@ssh.github.com` 认证 OK）；勿改回 22 端口 |
| 提交规范 | 中文标题 + 正文要点（feat/fix 前缀），如 `23472ba` |

## 6. 常用命令

```bash
PY=/home/xiap/opt/gene-drug-automd/env/bin
$PY/autovs doctor                                              # 环境体检
$PY/autovs prepare-receptor --target config/targets/hsp90b1_7ull.yaml
$PY/autovs gate --target config/targets/hsp90b1_7ull.yaml      # 共晶配体 redocking 质量门
$PY/autovs dock --target config/targets/hsp90b1_7ull.yaml --tag full    # 断点续跑
$PY/autovs visualize --target ... --tag full --top 10          # 可独立执行；--entry 指定分子
$PY/autovs report --target ... --tag full --top 20
$PY/python -m pytest tests/ -q                                 # 15 项回归
```

## 7. 当前状态快照（2026-10-11）

- **HSP90B1（GRP94，PDB 7ULL 链A）全库筛选已完成**：3,196 ok + 6 失败隔离 + 80 离群预过滤 = 3,282；命中（≤-7.0）1,686；最佳 MIDOSTAURIN -10.8 kcal/mol；质量门 PASS（窗口内 RMSD 1.172 Å）
- 全部产物在 `output_data/HSP90B1/` 与 `output_figure/HSP90B1/full/`
- 框架健康：doctor 全绿，测试 15/15，GitHub main 已同步（`23472ba`）

## 8. 领域硬约束与关键经验（违反即返工）

1. **负控先行**：蛋白-蛋白对接前必须跑 lysozyme × myoglobin 负控；ipTM < 0.16 为本项目确认的假阳性基线
2. **质量门**：全库筛选前必须过共晶配体 redocking 质量门；判据是"打分窗口(1.0 kcal/mol)内最小 RMSD < 2 Å"，**不能用 top-1 pose RMSD**（Vina 排序近简并）
3. **离群防护三参数**（target.yaml docking 节，0=关闭）：`max_rotatable_bonds=30`、`max_mol_weight=1000`（提交前预过滤，糖肽类大分子会让 Vina 卡死数天且靠接触面积刷虚高分）、`ligand_timeout_s=1800`（单配体超时兜底）
4. **受体清洗**：仅保留 ATOM 记录（去全部 HETATM），共晶配体单独提取；altloc 归一为 A 且传 `--default_altloc A`
5. **pdbfixer 坑**：`PDBFile.writeFile` 必须传 `keepIds=True`，否则残基重编号导致 2D/3D 编号体系分裂
6. **PyMOL 坑**：脚本模式 `cmd.orient/zoom` 需 `animate=0`；颜色名只有 `forest` 无 `forestgreen`；measurement 对象名必须带下划线（hb_0 非 hb0）；get_view 前 9 数需转置才是真实旋转矩阵；`clip('near',d)` 是相对移动不可用于标签深度
7. **ColabFold**：1.6.3 起 `--rank-by` 改名 `--rank`
8. **ZINC 不可脚本化下载**（反爬验证页），FDA 库用 EBI ChEMBL REST（max_phase=4）
9. **MANIFEST.csv 为合并式索引**：`--entry` 增量可视化按 entry_id 更新行，禁止整体覆盖
10. **2D 图不显示距离标注**、同残基同类型多线合并；3D 径向邻近标注、虚线引导线全画面预算 3 条（详见 visualize.py 注释与 project_memory）
