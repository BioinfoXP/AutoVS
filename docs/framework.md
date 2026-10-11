# AutoVS 框架搭建文档

> 版本 v0.1.0 ｜ 2026-10-10 ｜ 位置：`/home/xiap/bio_project/for_others/AutoMD/AutoVS`

AutoVS（Automated Virtual Screening）是 FDA 已上市药物库本地化批量虚拟筛选框架：从库资源本地化、受体准备、质量门验证，到 Vina 批量对接与结果汇总的全链路标准化流水线。命中分子可继续走 GeneDrug-AutoMD 的 MD / MM-PBSA 验证闭环。

---

## 1. 技术栈选型

| 层 | 选型 | 理由 |
|---|---|---|
| 语言 | Python 3.11 | 与现有 `gene-drug-automd` 环境一致，化学信息学生态最全 |
| 配置 | YAML + 标准库 dataclass | 无额外依赖、类型化、可校验；配置即单一事实源 |
| 分子处理 | RDKit 2025.09 | 3D 构象（ETKDGv3 + MMFF94）、标准化（MolStandardize）、RMSD |
| 质子化 | dimorphite-dl 2.1 | pH 7.4 生理状态枚举，与 automd 方法学一致 |
| 对接 | AutoDock Vina | CPU 并行成熟、可复现（固定 seed）、批量吞吐高 |
| 配体/受体准备 | Meeko 0.8 | mk_prepare_ligand / mk_prepare_receptor / mk_export |
| 并行 | ThreadPoolExecutor（Vina 子进程）+ ProcessPoolExecutor（RDKit 构建） | CPU 密集型两层并行，36 核机器 workers=30 |
| 依赖管理 | pip（editable 安装） + environment.yml | 开发迭代快；环境可整体复现 |
| 测试 | pytest | 13 项单元/集成测试，全程离线可跑 |

**不引入**：Web 框架、数据库服务、消息队列——批量 VS 是本地批处理场景，YAML + CSV + 文件系统即最简可靠方案。

## 2. 目录结构规范

```
AutoVS/
├── pyproject.toml          # 依赖声明 + CLI 入口（autovs）
├── environment.yml         # conda 环境复现（与 gene-drug-automd 对齐）
├── config/                 # ★ 配置体系（单一配置源）
│   ├── framework.yaml      #   框架级：路径 / 工具 / 运行时
│   ├── library.yaml        #   库级：数据源 / 本地化参数
│   └── targets/            #   项目级：靶点 / 盒子 / 质量门 / 对接参数
│       ├── template.yaml   #   新靶点模板（复制即用）
│       └── demo_egfr_1m17.yaml
├── autovs/                 # ★ 核心包（模块划分见 §3）
├── tests/                  # pytest 测试套件
├── docs/                   # 框架文档 / 本地化说明 / 验证报告
├── data/fda_library/       # FDA 库本地化资源（raw → prepared → index）
├── output_data/            # 结果表格（排行榜 CSV、summary JSON）
├── output_figure/          # 结果图（PDF 矢量）
└── temp/                   # 缓存与中间产物（受体、对接 pose，可删除重建）
```

布局纪律（与个人分析项目约定一致）：
- **代码进 `autovs/` 与 `tests/`**，表格进 `output_data/`，图进 `output_figure/`，可重建的中间产物进 `temp/`；
- **配置全部绝对路径**（不拼接 `file.path` 式相对路径），改机器只改 `config/framework.yaml`；
- 所有产物路径镜像项目层级（如 `output_data/<靶点名>/`）。

## 3. 核心模块划分

| 模块 | 职责 | 关键 API |
|---|---|---|
| `autovs/config.py` | 三层 YAML 加载 + dataclass 类型化 + 严格校验 | `load_config()` |
| `autovs/library.py` | FDA 库本地化：可插拔数据源 → 本地索引 → 离线访问 | `library_fetch()` / `library_build()` / `LocalLibrary` |
| `autovs/ligand.py` | 单分子准备链：去盐 → 质子化 → 3D → PDBQT | `build_entry()` / `enumerate_states()` |
| `autovs/receptor.py` | 受体清洗/准备、对接盒、质量门（redocking RMSD） | `prepare_receptor()` / `run_gate()` / `evaluate_redocking()` |
| `autovs/docking.py` | Vina 批量引擎：并行 + 断点续跑 + 失败隔离 + 离群预过滤/单配体超时 | `dock_single()` / `run_batch()` |
| `autovs/aggregate.py` | 汇总排序、阈值过滤、打分分布图（含失败/预过滤计数） | `aggregate_results()` |
| `autovs/visualize.py` | Top 命中自动可视化：PyMOL 3D 渲染 + DS 风格 2D 相互作用图 | `run_visualize()` / `detect_interactions()` |
| `autovs/report.py` | Markdown 任务报告（概览/质量门/统计/排行/图谱对照）+ 全量排行 CSV 副本 | `generate_report()` |
| `autovs/mcp_server.py` | MCP server：全流程 10 工具暴露给 AI 客户端（stdio） | `autovs mcp` / `python -m autovs.mcp_server` |
| `autovs/cli.py` | 统一命令行入口 | `autovs <command>` |
| `autovs/utils.py` | 日志 / 子进程 / JSON / 配置指纹 | `run_cmd()` / `config_hash()` |

**MCP 接入示例**（客户端 mcp 配置 JSON，如 Trae/Claude Desktop/Cursor）：

```json
{
  "mcpServers": {
    "autovs": {
      "command": "/home/xiap/opt/gene-drug-automd/env/bin/python",
      "args": ["-m", "autovs.mcp_server"],
      "cwd": "/home/xiap/bio_project/for_others/AutoMD/AutoVS"
    }
  }
}
```

MCP 工具面：`doctor` / `library_stats` / `library_search` / `prepare_receptor` / `run_gate` / `submit_dock`（后台提交，返回 job_id）/ `dock_status`（进度与进程探活，server 重启后按 target+tag 追查）/ `query_ranking` / `visualize`（支持 entries 指定化合物）/ `generate_report`。依赖 `mcp<2`（1.x SDK）。

### 3.1 数据流

```
library fetch (联网一次) ──► raw/source_records.json
                                │
library build (纯离线)  ──► prepared/<shard>/*.{sdf,pdbqt} + index/*.csv
                                │
receptor prepare        ──► temp/receptor/<靶点>/receptor.pdbqt
                                │
gate (质量门)           ──► gate_result.json  ──不过──► 硬停止
                                │ PASS
dock --limit N / 全库   ──► temp/docking/<靶点>/<批次>/result.json × N（断点续跑）
                                │
aggregate               ──► output_data/*_ranking.csv + output_figure/*.pdf
                                │
visualize（dock 自动触发）──► output_figure/<靶点>/<批次>/   # 任务自包含目录，整包可拷走
                                ├─  pose_R{rank}_{药名}_3d.png       # 3D 口袋特写
                                ├─  pose_R{rank}_{药名}_2d.{pdf,png} # 2D 相互作用图
                                ├─  pose_R{rank}_{药名}_pair.png     # 3D|2D 横拼对照图
                                ├─  top_overview.png / MANIFEST.csv / <批次>_visualize_summary.json
report（dock 自动触发）  ──► 同目录 report.md + <批次>_ranking.csv 副本（全量排行，报告自包含）
```

可视化模块要点：

- **3D（PyMOL 批处理）**：`pymol -cq` 独立进程渲染 Top N 叠加总览 + 逐姿势特写（天蓝透明 cartoon、紫碳配体、非疏水相互作用残基黄碳侧链 + 黑体残基标签、绿色氢键虚线）+ 可交互 `.pse` 会话。
- **2D（DS 风格）**：自研几何判据（氢键 H···A ≤2.6 Å 且 D-H···A >110°；疏水 C···C ≤4.0 Å；盐桥 ≤4.0 Å；π 堆积环心 ≤5.5 Å；阳离子-π ≤5.0 Å）→ 配体键线式居中、残基椭圆环绕节点（氢键深绿/疏水浅绿/π 粉紫/盐桥红）、彩色短虚线连到配体具体原子，输出 PDF + PNG。
- **3D↔2D 对应关系**：同一姿势的 3D/2D/对照图共用文件前缀 `pose_R{rank}_{药名}`；`_pair.png` 将两图等高横拼，残基一一对应可直接互读；MANIFEST.csv 汇总 rank/化合物/打分/各图路径/相互作用清单。
- **指定化合物可视化**：`autovs visualize --target ... --tag ... --entry "tadalafil,FDA00466"`（entry_id 精确或药名子串，大小写不敏感，逗号分隔；给出时忽略 --top），产物落入同一任务目录、同一命名规则。
- **受体加氢**：pdbfixer + OpenMM（pH 7.4，`keepIds=True` 保留原始残基编号，保证 2D/3D 编号一致）。
- **自动化**：`autovs dock` 汇总后自动调用 `run_visualize()` 与 `generate_report()`（失败仅告警不影响主流程）；也可独立执行 `autovs visualize --target ... --tag ... --top N`、`autovs report --target ... --tag ...`。报告内附全量 rank 排行 CSV 副本与相对链接，任务目录整包拷走即含排行。

### 3.2 关键设计决策

1. **质量门（Gate）是硬性放行条件**：共晶配体重对接打分窗口（默认 1.0 kcal/mol）内最小 RMSD < 2.0 Å 才允许全库筛选。判据采用行业通用 redocking 标准——晶体学姿势须在 top-ranked poses 中被复现，而非苛求 top-1 pose（Vina 打分排序存在近简并歧义）。
2. **断点续跑以配置指纹为准**：每个条目 `result.json` 记录 `config_hash`（exhaustiveness/num_modes/seed/盒子），指纹不一致自动重算，一致则跳过。
3. **失败隔离**：单分子准备或对接失败只写 `error` 字段，批次不中断，事后统一统计重试。
4. **离群大分子自动防护**（docking 三参数，可在 target.yaml 调整，0=关闭）：`max_rotatable_bonds=30` / `max_mol_weight=1000` 提交前预过滤（MW/柔性超限的糖肽类大分子 Vina 会卡死数天且打分靠接触面积虚高，属假阳性来源），写 `status=skipped + reason`；`ligand_timeout_s=1800` 对每个 vina 子进程限时，超时转 DockingError 走失败隔离。summary.json 记录 `n_failed` / `n_skipped_outlier`。
5. **可插拔数据源**：`chembl`（默认）/ `zinc` / `local` 三源适配，网络层仅存在于 fetch 阶段。

## 4. 依赖管理方案

环境配置与 MCP 配置相互独立，均可独立复现：

- **conda 方式（推荐，外部工具一并装好）**：`conda env create -f environment.yml && conda activate autovs`（conda-forge + bioconda，含 vina/meeko/openbabel/fpocket/rdkit）；或复用已就绪环境 `conda activate "$HOME/opt/gene-drug-automd/env"`；
- **pip 方式（裸环境，外部工具自备）**：`python3.11 -m venv .venv && source .venv/bin/activate && pip install -e ".[dev,mcp]"`；vina/meeko/obabel/fpocket/pymol 需自行安装并登记到 `config/framework.yaml` 的 `tools:`；
- **框架本体**：`pip install -e .`（editable，改代码即生效），CLI 注册为 `autovs`；
- **外部二进制**（vina/meeko/fpocket）路径统一登记在 `config/framework.yaml` 的 `tools:`，`autovs doctor` 校验；
- **MCP 依赖单独剥离**：extras `mcp`（`pip install "mcp>=1.30,<2"`），仅 AI 客户端驱动需要，纯命令行使用可不装；
- 版本约束只锁下限（关键库），避免与 conda-forge 滚动更新冲突。

## 5. 基础配置体系

三层配置，职责分离，均可独立变更：

| 文件 | 内容 | 变更频率 |
|---|---|---|
| `framework.yaml` | 路径、工具绝对路径、workers、seed、日志 | 换机器时 |
| `library.yaml` | 数据源（chembl/zinc/local）、pH、质子化状态数、分片大小 | 换库/换标准时 |
| `targets/<名>.yaml` | PDB、链、共晶配体、对接盒、质量门阈值、对接参数 | 每个新靶点 |

约定：
- 任何命令的参数差异都体现在配置文件里，不写死在代码；
- 配置校验失败直接报 `ConfigError` 并指出缺失字段（fail fast）；
- 靶点配置从 `config/targets/template.yaml` 复制，改完即用。

## 6. 开发约束（标准化）

1. 新功能进对应模块，禁止在 `cli.py` 堆业务逻辑；
2. 外部命令统一走 `utils.run_cmd()`（带日志、可追溯）；
3. 产物路径由配置推导，代码中不出现裸相对路径；
4. 缓存产物必须带有效性判据（文件存在 + `config_hash`），跳过逻辑可 `force` 覆盖；
5. 每个公开函数写中文 docstring，说明方法学来源；
6. 提交前 `pytest tests/ -q` 全绿 + `autovs doctor` 全绿。

## 7. 快速上手（开发人员视角）

```bash
conda activate "$HOME/opt/gene-drug-automd/env"
cd /home/xiap/bio_project/for_others/AutoMD/AutoVS

autovs doctor                                        # 环境体检
autovs library fetch && autovs library build         # FDA 库本地化（一次性）
autovs receptor prepare --target config/targets/demo_egfr_1m17.yaml
autovs gate       --target config/targets/demo_egfr_1m17.yaml   # 不过则禁止筛选
autovs dock       --target config/targets/demo_egfr_1m17.yaml --limit 50   # 冒烟
autovs dock       --target config/targets/demo_egfr_1m17.yaml              # 全库
autovs library stats && autovs library verify-offline
```

新靶点开发 = 复制 `config/targets/template.yaml` → 填 PDB/共晶配体 → `receptor prepare` → `gate` → `dock`。无需写任何代码。
