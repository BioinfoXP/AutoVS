# AutoVS 从零上手教程

本教程带你从裸机开始，一步步搭好环境、跑通第一个虚拟筛选、看懂每一个输出。预计 30 分钟内完成第一次全流程体验（含 20 个化合物的冒烟对接）。

> 已有环境？直接跳到[第 3 节](#3-五分钟跑通第一个筛选egfr-demo)。

---

## 1. AutoVS 是什么、什么时候用它

一句话：**把 3,282 个 FDA 已上市药物挨个对接进你的靶蛋白口袋，按打分排序，自动出版级的结合模式图和任务报告。**

典型场景：

- 老药新用（drug repurposing）：给目标靶点找已有的临床药物
- 快速评估某靶点口袋的"可药性"：看全库打分分布
- 为后续 MD / FREE 能计算 / 实验验证提供候选清单

**它不是**：高精度结合自由能预测工具。Vina 打分只用于富集排序，Top 命中必须经 MD/实验复核（详见 [方法论要点](../README.md#方法论要点解读结果前必读)）。

---

## 2. 环境搭建（一次性）

### 2.1 创建 conda 环境（只装 pip 装不了的 CLI 二进制）

```bash
conda env create -f environment.yml   # 环境名 autovs：python + vina/openbabel/fpocket/pymol
conda activate autovs
```

> 注意：对接引擎用 conda-forge 的 **`vina`** 包（environment.yml 已声明）。bioconda 的
> `autodock-vina` 官方已只剩 1.1.2（2011 年古董版），不要使用。

### 2.2 安装 Python 依赖

上一步 environment.yml 的 pip 段已自动执行 `pip install -e ".[dev,mcp]"`；如单独执行：

```bash
cd AutoVS        # 本仓库根目录
pip install -e ".[dev,mcp]"
```

这一步会装上：RDKit（化学信息学）、Meeko（PDBQT 准备）、pdbfixer/OpenMM（受体加氢）、scipy/gemmi（meeko 0.8 隐性依赖）、matplotlib（分布图）、mcp（AI 客户端接入，可选）。

### 2.3 安装外部工具

AutoVS 通过**绝对路径**调用外部工具，全部路径集中在 `config/framework.yaml` 的 `tools:` 段。

| 工具 | 用途 | 安装方式 |
|---|---|---|
| AutoDock Vina | 对接引擎 | `conda install -c conda-forge vina` 或[官方发布页](https://github.com/ccsb-scripps/AutoDock-Vina/releases) |
| Meeko 脚本 | PDB/SDF → PDBQT | 随 `pip install meeko` 生成 `mk_prepare_receptor.py` 等可执行脚本 |
| Open Babel | 格式转换（可选） | `conda install -c conda-forge openbabel` |
| PyMOL | 3D 渲染（可选） | `conda install -c schrodinger pymol`（开源版 `conda install -c conda-forge pymol-open-source`） |
| fpocket | 口袋检测（可选） | `conda install -c conda-forge fpocket` |

> **PyMOL 缺失不会阻塞流程**——自动跳过 3D 渲染，仍出 2D 相互作用图与报告。

### 2.4 把工具路径写进配置

编辑 `config/framework.yaml`，将 `tools:` 段改成你的实际路径（`which vina` 等命令可查）：

```yaml
tools:
  python: /home/you/miniconda3/envs/autovs/bin/python
  vina: /home/you/miniconda3/envs/autovs/bin/vina
  mk_prepare_receptor: /home/you/miniconda3/envs/autovs/bin/mk_prepare_receptor.py
  mk_prepare_ligand: /home/you/miniconda3/envs/autovs/bin/mk_prepare_ligand.py
  mk_export: /home/you/miniconda3/envs/autovs/bin/mk_export.py
  obabel: /home/you/miniconda3/envs/autovs/bin/obabel
  pymol: /home/you/miniconda3/envs/autovs/bin/pymol
```

同时把 `paths:` 段的 `root` 改为你的 AutoVS 克隆路径。

### 2.5 验证：环境体检

```bash
autovs doctor
```

预期输出（全部 OK 即就绪）：

```
python: /home/you/miniconda3/envs/autovs/bin/python
vina: OK /home/you/miniconda3/envs/autovs/bin/vina
pymol: OK /home/you/miniconda3/envs/autovs/bin/pymol
mk_prepare_receptor: OK /home/you/miniconda3/envs/autovs/bin/mk_prepare_receptor.py
fda_library: 3282 entries
vina_version: AutoDock Vina v1.2.5
```

### 2.6 验证：离线可用性

```bash
autovs library verify-offline     # 库索引与 prepared 结构可离线读取
autovs library stats              # {"total": 3282, ...}
```

---

## 3. 五分钟跑通第一个筛选（EGFR demo）

仓库自带一个配好的示范靶点 `config/targets/demo_egfr_1m17.yaml`（EGFR 激酶域 × ERLOCOXIB 共晶，1M17）。直接冒烟：

```bash
# 只对接前 20 个化合物，验证整条链路（受体→对接→汇总→可视化→报告）
autovs dock --target config/targets/demo_egfr_1m17.yaml --tag smoke --limit 20
```

你会看到四段进度：

```
[1/4] 受体准备 ...        temp/receptor/EGFR_1M17/receptor_h.pdb
[2/4] 质量门  ... PASS    （窗口内最小 RMSD x.xx Å < 2.0）——不 PASS 会直接终止！
[3/4] 批量对接 ...        temp/docking/EGFR_1M17/smoke/FDA00001/...
[4/4] 汇总 + 可视化 + 报告 ...
```

完成后打开（直接双击或 `open`）：

```
output_figure/EGFR_1M17/smoke/report.md        ← 任务报告（从这里开始看）
output_figure/EGFR_1M17/smoke/pose_R01_*_pair.png  ← 最佳命中的结合模式对照图
output_data/EGFR_1M17/smoke_ranking.csv        ← 排行榜
```

**看到 report.md 就代表全流程通了。** 接下来换成你自己的靶点（第 4 节），或直接全量筛选（第 5 节）。

> 仓库自带的 demo 报告在 `output_figure/EGFR_demo/smoke/report.md`，可直接当"参考答案"对照。

---

## 4. 换成你自己的靶点

### 4.1 找结构

去 [RCSB PDB](https://www.rcsb.org/) 搜你的靶蛋白，选一个**共晶结构**（有结合配体的），记下三样东西：

1. **PDB ID**（如 `1M17`）
2. **链**（如 `A`）
3. **共晶配体三字母码**（PDB 页面右侧 "Ligands" 处，如 `AQ4`）

### 4.2 写靶点配置

```bash
cp config/targets/template.yaml config/targets/my_target.yaml
```

最少只需改 4 行：

```yaml
target:
  name: MYTARGET            # ① 靶点名（决定输出目录名）
  pdb_id: 1M17              # ② PDB 编号（自动下载并缓存）
  chains: [A]               # ③ 对接哪条链
  reference_ligand:
    resname: AQ4            # ④ 共晶配体三字母码
```

盒子不用手算：`box.mode: reference_ligand` 会让程序以共晶配体质心为中心、外扩 `padding: 5.0` Å 自动生成对接盒（尺寸钳制在 `min_box`~`max_box`）。如果你的配体不是口袋探针（例如是别构脂质），才需要手动填 `center`/`size`。

### 4.3 质量门：先证明口袋"可复现"

```bash
autovs receptor --target config/targets/my_target.yaml   # 清洗 + 加氢
autovs gate     --target config/targets/my_target.yaml   # 共晶配体重对接
```

质量门逻辑：把**晶体里已有的配体拿掉再重新对接回去**，如果 Vina 能复现晶体姿势，说明这套盒子+参数对你的靶点有效，之后的全库筛选结果才可信。

- **PASS**（窗口内最小 RMSD < 2 Å）→ 放行
- **FAIL** → 程序拒绝全库筛选。排查顺序：共晶配体名对不对 → 链对不对 → 口袋是不是太刚/太柔 → 换 `gate.exhaustiveness: 64` 再试

> 为什么判据是"打分窗口内最小 RMSD"而不是 top-1 pose？Vina 打分存在近简并歧义，晶体姿势可能排在 mode 4/5 但打分只差 0.3 kcal/mol。按 top-1 判会错杀好口袋（这是本项目实测踩过的坑）。

---

## 5. 全库筛选

```bash
# 建议先 50 个冒烟（约 1-2 分钟），确认无误再全量
autovs dock --target config/targets/my_target.yaml --tag screen --limit 50
autovs dock --target config/targets/my_target.yaml --tag screen          # 全量 3,282 个
```

要点：

- **断点续跑**：每个化合物独立落盘 `result.json`。中途 Ctrl+C / 断电 / 崩溃，**重跑同一条命令自动跳过已完成的**。
- **并行度**：`framework.yaml` 的 `runtime.workers`（默认 30）。改大前确认 CPU 核数。
- **进度查看**：对接中另开终端

  ```bash
  ls temp/docking/MYTARGET/screen/*/result.json | wc -l    # 已完成数
  ```
- **耗时参考**：36 核机器全库约 2-4 小时（exhaustiveness=8）。

跑完自动生成排行榜、Top 10 可视化与报告，无需额外操作。

---

## 6. 看懂每一个输出

### 6.1 排行榜 `output_data/<靶点>/<批次>_ranking.csv`

| 列 | 含义 |
|---|---|
| `rank` | 按 `top_score` 升序的名次 |
| `entry_id / name` | 化合物 ChEMBL ID / 药名 |
| `top_score_kcal_mol` | 9 个 pose 中最好的 Vina 打分（越负越好） |
| `median_score_kcal_mol` | 9 pose 中位数——中位数也好才算稳，最好值好但中位数差可能是"侥幸姿势" |
| `mol_weight / logp` | 分子量 / 脂水分配系数（LogP） |
| `first_approval` | 首次批准年份（老药安全性数据更充分） |
| `pose_sdf` | 最佳姿势 3D 坐标文件路径 |

### 6.2 打分分布图 `<批次>_score_distribution.pdf`

全库打分直方图 + 红色阈值虚线（默认 -7.0）。**看整体形状比看单个值重要**：正态偏右说明口袋区分度好；一堆 -11 的异常尖峰提示大环/多芳环偏好。

### 6.3 结合模式对照图（核心产物）

每张 `_pair.png` 左边 3D、右边 2D，**残基一一对应**：

- **3D 半边**：天蓝半透明 cartoon = 蛋白骨架；紫/黄 stick = 配体与口袋残基；绿色虚线 = 氢键；黑体文字 = 残基标签（下划线连到残基位置）
- **2D 半边**：中心键线式 = 配体；环绕椭圆 = 相互作用残基——**深绿实心 = 氢键、浅绿 = 疏水、粉紫 = π 堆积/阳离子-π、红 = 盐桥**；虚线连到配体具体参与原子

真实例子（demo 靶点 EGFR 1M17 冒烟筛选 Top1，RISDIPLAM）：

![RISDIPLAM 3D|2D 对照图](../output_figure/EGFR_demo/smoke/pose_R01_RISDIPLAM_pair.png)

同一张图里可以同时读：3D 空间摆放（LYS114/ASP149 在口袋哪个方向）+ 2D 作用类型（哪些是氢键哪些是疏水）。

### 6.4 任务报告 `<批次>/report.md`

五章结构：任务概览（含质量门 RMSD）→ 结果统计（命中率）→ Top 20 排行表 → Top 命中对照图 → 方法论注意事项。**报告是与图谱同目录的自包含文件，整个 `<批次>/` 文件夹拷给同事就能看。**

---

## 7. 对指定化合物深入分析

排行榜里看到感兴趣的药，或想验证某个已知结合药：

```bash
# 按药名子串或 entry_id，逗号分隔任意多个（大小写不敏感）
autovs visualize --target config/targets/my_target.yaml --tag screen \
                 --entry "imatinib,FDA09265,gefitinib"
```

产物落在同一 `<批次>/` 目录、同一命名规则（`pose_R{它的rank}_{药名}_*`），并自动合并进 `MANIFEST.csv` 与重跑 `autovs report` 即可纳入报告。

---

## 8. 让 AI 帮你跑（MCP 接入）

把下面配置加进 AI 客户端（Trae / Claude Desktop / Cursor）的 MCP 设置：

```json
{
  "mcpServers": {
    "autovs": {
      "command": "<你的 python 绝对路径>",
      "args": ["-m", "autovs.mcp_server"],
      "cwd": "<AutoVS 绝对路径>"
    }
  }
}
```

> 各客户端的配置文件位置、11 个工具逐个详解、更多实战话术与排错表，见 **[MCP 指南](mcp_guide.md)**。

重启客户端后即可自然语言驱动，例如：

> "用 AutoVS 给 1M17 这个靶点（EGFR，链 A，共晶配体 AQ4）跑一遍虚拟筛选，先跑质量门，PASS 后全量筛，完成后给我 Top 10 的对照图和报告。"

AI 会依次调用 `doctor → run_gate → submit_dock → dock_status(轮询) → visualize → generate_report`。对接这种小时级任务走后台提交模式，AI 查进度用的是"数结果文件"的真实进度，server 重启也不丢任务。

---

## 9. FAQ 与排错

| 现象 | 原因与解决 |
|---|---|
| `gate` FAIL | 共晶配体名/链写错（最常见）；口袋定义不适合；`gate.exhaustiveness` 提到 64 重试 |
| `mk_prepare_receptor` 报 altloc 错误 | 结构有替代位标残基——框架已自动归一为 A 链；若仍报错，检查自定义结构的 ALTLOC |
| 某化合物对接失败 | 该化合物 3D 构象生成失败，排行榜里会缺行；`temp/docking/<靶点>/<批次>/<ID>/` 里有错误日志 |
| 对接比上次快很多且秒完成 | 正常，断点续跑在跳过已完成化合物 |
| 换了配置后结果"没变化" | 配置指纹校验生效，配置变更后中间产物自动失效重算；确认你改的配置在 `--framework/--library` 指向的文件里 |
| 3D 图没生成只有 2D | PyMOL 未安装或 `tools.pymol` 路径不对；不影响其余产物 |
| `--entry` 找不到化合物 | 匹配的是排行榜里的 `entry_id` 精确匹配或药名**子串**；用 `autovs library stats` 确认库里存在 |
| 打分普遍很平（±0.5 内） | 口袋太小或盒子太大稀释了区分度；检查盒子是否钳到了 `max_box` |

---

## 10. 下一步

- 读 [框架架构](framework.md) 了解模块划分与数据流
- 读 [FDA 库本地化](fda_library_localization.md) 了解数据从哪来、怎么重建
- 读 [验证测试报告](validation_report.md) 了解质量门与断点续跑的实测数据
- 读 [README 方法论要点](../README.md#方法论要点解读结果前必读)，理解打分边界，别把 Vina 分当自由能
