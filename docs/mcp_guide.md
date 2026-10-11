# AutoVS MCP 指南：配置与使用

MCP（Model Context Protocol）是 AI 客户端（Trae / Claude Desktop / Cursor 等）与外部工具之间的标准接口。把 AutoVS 接入 MCP 后，你可以直接用自然语言让 AI 帮你跑虚拟筛选——AI 会自己调用体检、质量门、对接、可视化、报告等工具，你只看结果。

本指南分三部分：**配置接入 → 工具详解 → 实战话术**。

---

## 1. 配置接入

### 1.1 前置条件（一次性）

> MCP 依赖与筛选运行环境（conda/pip，见 README「环境配置」）**相互独立**：先把 AutoVS 本体装好并跑通 `autovs doctor`，再叠加这一步。

```bash
# 在 AutoVS 的 Python 环境里装 MCP SDK（注意必须 <2，2.x 改名不兼容）
pip install "mcp>=1.30,<2"        # 安装本体时也可带上：pip install -e ".[mcp]"

# 验证：能列出 10 个工具即 OK
python -c "import autovs.mcp_server; print('OK')"
```

### 1.2 通用配置模板

MCP server 以 **stdio 子进程**方式运行：AI 客户端负责拉起 `python -m autovs.mcp_server` 并通过标准输入输出通信。你只需要在客户端的 MCP 配置文件里加一段 JSON：

```json
{
  "mcpServers": {
    "autovs": {
      "command": "<python 绝对路径>",
      "args": ["-m", "autovs.mcp_server"],
      "cwd": "<AutoVS 仓库绝对路径>"
    }
  }
}
```

三个值按你的机器替换：

| 字段 | 填什么 | 怎么查 |
|---|---|---|
| `command` | 装 AutoVS 依赖的那个 Python 解释器**绝对路径** | `conda activate <env>` 后 `which python`；或 `autovs doctor` 输出的第一行 |
| `args` | 固定 `["-m", "autovs.mcp_server"]` | 不用改 |
| `cwd` | AutoVS 仓库根目录（`config/` 所在处） | `cd` 到仓库后 `pwd` |

> **为什么 `cwd` 重要**：server 默认从 `config/framework.yaml` 等相对位置加载三层配置，`cwd` 必须指向仓库根。命令行也可以带 `--framework/--library` 参数，但 MCP 模式下用 `cwd` 最省事。

### 1.3 各客户端配置位置

**Trae**：设置 → MCP → 添加 MCP Server，或直接编辑 `~/.trae-cn/mcp.json`（内容即上面模板）。

**Claude Desktop**：

- macOS：`~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows：`%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "autovs": {
      "command": "/home/you/miniconda3/envs/autovs/bin/python",
      "args": ["-m", "autovs.mcp_server"],
      "cwd": "/home/you/AutoVS"
    }
  }
}
```

保存后**完全退出并重启** Claude Desktop（菜单栏图标 → Quit，不是关窗口）。

**Cursor**：`~/.cursor/mcp.json`，格式同模板；或在设置 → MCP → Add new global MCP server 里粘贴。

### 1.4 验证接入成功

重启客户端后，新建对话，看工具/锤子图标是否出现 **10 个 autovs 工具**（doctor、library_stats、library_search、prepare_receptor、run_gate、submit_dock、dock_status、query_ranking、visualize、generate_report）。也可以直接问 AI：

> "调用 doctor 工具体检一下环境"

AI 返回 vina/pymol/库条数即接入成功。

### 1.5 排错：server 没起来

| 现象 | 排查 |
|---|---|
| 工具列表为空 / 无 autovs | 手动跑 `cd <cwd> && <command> -m autovs.mcp_server`——应挂起等待输入（stdio 正常现象），报错则看堆栈：多为 `ModuleNotFoundError: mcp`（装错环境/没装） |
| `Exception: FastMCP not found` 或 API 报错 | 装到了 mcp 2.x：`pip install "mcp>=1.30,<2"` |
| 工具在但调用报 FileNotFoundError | `cwd` 不对，或 `config/framework.yaml` 里工具绝对路径失效——先 `autovs doctor` |
| 对话中出现多个重复工具 | 同一配置粘贴了两次，删掉多余的 |

---

## 2. 十个工具详解

按虚拟筛选的执行顺序介绍。**你通常不需要记参数**——AI 会从对话上下文里填，但知道每个工具做什么有助于你下指令。

| # | 工具 | 干什么 | 关键参数 | 典型时机 |
|---|---|---|---|---|
| 1 | `doctor` | 环境体检：vina/pymol/mk_prepare 是否就绪、库有多少条 | 无 | 一开始；换机器后 |
| 2 | `library_stats` | 库规模统计（总数、平均分子量） | 无 | 想确认库 |
| 3 | `library_search` | 按药名/ID 子串搜库 | `query`, `limit` | "库里有没有伊马替尼？" |
| 4 | `prepare_receptor` | 受体清洗（去配体/水/altloc 归一）+ 加氢 | `target_yaml` | 质量门前（也可跳过，gate 会自动准备） |
| 5 | `run_gate` | 质量门：共晶配体重对接，窗口内最小 RMSD < 2 Å 才 PASS | `target_yaml` | **全库筛选前必须 PASS** |
| 6 | `submit_dock` | **后台提交**批量对接（断点续跑；MW/可旋转键超限大分子自动预过滤，超时/失败自动隔离），立即返回 job_id | `target_yaml`, `tag`, `limit`(可选冒烟) | 正式筛选 |
| 7 | `dock_status` | 查真实进度：已完成数/总数/进程存活/最后日志 | `job_id` 或 `target_yaml+tag` | 对接中轮询；**server 重启后也可按靶点+批次追查** |
| 8 | `query_ranking` | 查排行榜 Top N（可只看阈值命中） | `target_yaml`, `tag`, `top_n`, `below_threshold_only` | 对接完/中途看结果 |
| 9 | `visualize` | 生成 3D 口袋图 + 2D 相互作用图 + 3D\|2D 对照图 | `target_yaml`, `tag`, `top_n` 或 `entries` | 深入看具体分子 |
| 10 | `generate_report` | 生成 Markdown 任务报告，并在报告同目录附全量 rank 排行 CSV 副本 | `target_yaml`, `tag`, `top_show` | 收尾交付 |

细节说明：

- **长任务模式**：对接 3,282 个药要几小时，所以 `submit_dock` 只负责"提交"（返回 job_id 和日志路径），真正的进度用 `dock_status` 查——它直接数磁盘上的 `result.json`，和 server 进程生死无关，重启 AI 客户端也不丢任务。
- **`entries` 参数**（visualize）：逗号分隔药名子串或 entry_id，如 `"tadalafil,gefitinib"`；给了就忽略 `top_n`。
- **产物位置**：所有工具产出都落在 `output_data/<靶点>/` 与 `output_figure/<靶点>/<批次>/`（详见 [README 输出布局](../README.md#输出布局)），MCP 只是"遥控器"，产物仍是普通文件。

---

## 3. 实战话术

配置好后，你这样说，AI 这样做：

### 场景 A：从零跑一个新靶点

> **你**："用 AutoVS 检查一下环境，然后给 PDB 1M17（EGFR，链 A，共晶配体 AQ4）建个靶点配置，跑质量门，PASS 了就全库筛选，批次叫 screen。"

AI 的动作链：`doctor` → 创建 `config/targets/egfr_1m17.yaml`（按仓库模板）→ `run_gate` → PASS 后 `submit_dock(tag="screen")` → 用 `dock_status` 定期向你汇报进度。

### 场景 B：查进度（对接中另开对话）

> **你**："EGFR 那个 screen 批次跑到哪了？"

AI 调 `dock_status(target_yaml=..., tag="screen")`，回你"已完成 1520/3282，进程存活，最近日志 …"。

### 场景 C：看结果

> **你**："现在 Top 20 是哪些？只看低于 -7 的。"
> **你**："把第一名和 erlotinib 的结合模式图画出来。"

AI 调 `query_ranking(top_n=20, below_threshold_only=true)` 和 `visualize(entries="RISDIPLAM,erlotinib")`，并告诉你图在 `output_figure/EGFR_demo/screen/pose_R01_RISDIPLAM_pair.png`。

### 场景 D：出交付物

> **你**："生成完整任务报告，Top 20 展示。"

AI 调 `generate_report`，返回 `report.md` 路径——整包 `output_figure/<靶点>/<批次>/` 拷走即可分发。

---

## 4. 与命令行的关系

MCP 工具与 CLI 是**同一套核心函数**的两个入口，产物完全一致、混用无冲突：

| 你想做的 | CLI | MCP |
|---|---|---|
| 后台跑全库 | `autovs dock ...`（自己 nohup） | `submit_dock` + `dock_status` |
| 快速单次操作 | `autovs visualize/report/...` | 对话式调用 |

习惯命令行的可以在终端干活，AI 对话里查进度、改需求更方便——两者可同时进行。
