# FDA 库本地化改造说明

> 版本 v1 ｜ 2026-10-10 ｜ 库：`fda_approved_drugs`（data/fda_library/）

## 1. 本地化目标与范围

| 要求 | 改造落点 | 状态 |
|---|---|---|
| 资源文件本地化 | 全量化合物记录 + 3D SDF + PDBQT + 索引 CSV 落盘本地 | ✅ |
| 接口请求本地化适配 | 网络层收敛到唯一 `fetch` 步骤（可插拔数据源适配器），运行期零远程调用 | ✅ |
| 本地环境兼容性处理 | 路径绝对化、工具链显式登记、失败隔离、断点续跑 | ✅ |
| 本地化配置接入 | `config/library.yaml` 单一配置源，数据源/参数全部可配 | ✅ |
| 功能完整、无异常报错 | `autovs library stats` / `verify-offline` 验证 + pytest 13 项 | ✅（见验证报告） |

## 2. 数据源决策（本地化适配过程）

### 2.1 首选源与实际采用源

| 候选源 | 探测结果 | 结论 |
|---|---|---|
| ZINC15 FDA 子集（`zinc15.docking.org/substances/subsets/fda/`，约 1,379 化合物） | 站点启用反爬验证（返回 "Verification Required" HTML），curl/UA 均无法脚本化下载 | **保留为可选源**（`type: zinc`），浏览器下载后可用 `type: local` 导入 |
| EBI ChEMBL REST（`max_phase=4` 已上市小分子） | 官方开放 API，可脚本化、可分页、字段丰富 | **采用为默认源** |
| DrugBank | 需机构账号 | 不采用 |

**采用 ChEMBL max_phase=4 的说明**：该集合为全球主要监管机构（含 FDA、EMA 等）批准的已上市小分子药物，是 FDA 药物库的**超集**（3,475 条记录，含药物名、ATC 分类、首次批准年份），用于药物重定位虚拟筛选覆盖度更高；严格 FDA-only 规模约 2,000+，属其子集。数据源切换只需改 `library.yaml` 的 `source.type`，无需改代码。

### 2.2 数据本地化流水线

```
ChEMBL REST（4 页 × 1000）
   │  library fetch  ←—— 唯一联网步骤
   ▼
raw/source_records.json（3,311 条含 SMILES）+ raw/manifest.json（来源指纹/时间戳/条数）
   │  library build  ←—— 纯离线
   ▼
去重（InChIKey）→ 去盐（FragmentParent）→ 质子化（dimorphite-dl pH 7.4）
   → 3D 构象（ETKDGv3 + MMFF94）→ PDBQT（Meeko）→ 分片存储
   ▼
prepared/shard_XXX/FDAxxxxx_s1.{sdf,pdbqt}   （3,282 条可对接）
index/fda_library_index.csv                  （本地索引，21 字段）
```

## 3. 本地索引设计（离线访问接口）

`data/fda_library/index/fda_library_index.csv`，一行一个可对接条目：

| 字段组 | 字段 | 用途 |
|---|---|---|
| 标识 | entry_id, lib_id, state_id, source, source_id, name | 唯一键 / 溯源 / 药名 |
| 化学 | smiles, inchi_key, formal_charge, heavy_atoms | 结构标识 |
| 类药性 | mol_weight, logp, hbd, hba, rotatable_bonds | Lipinski 过滤（标记不删除） |
| 药物元数据 | first_approval, atc | 首次批准年份 / ATC 适应症分类 |
| 产物 | status, sdf_path, pdbqt_path, error | 本地文件定位与质量审计 |

访问接口（`autovs/library.py::LocalLibrary`）：`iter_entries()` / `get(entry_id)` / `sample(n, seed)` / `stats()`，**全部纯本地文件读写**。

## 4. 本地环境兼容性处理

1. **绝对路径**：库资源、工具、输出全部绝对路径登记于配置，换机器只改 `config/framework.yaml`；
2. **工具链显式化**：vina / mk_prepare_ligand / mk_export 等二进制路径进配置，`autovs doctor` 逐项校验；
3. **失败隔离**：29 条边缘分子（无机盐 NaF、含硒/金属氨基酸、复杂嵌合体构象嵌入失败）准备失败，仅标记 `status=failed` 并记录 `error`，不影响其余 3,282 条；
4. **断点续跑**：`prepare.skip_existing=true`，重建时已有 SDF+PDBQT 直接跳过，中断可续；
5. **确定性**：全局 seed 写入配置，构象生成与抽样可复现；
6. **离线验证**：`autovs library verify-offline` 在禁网（socket 哨兵）条件下验证索引读取与产物完整性。

## 5. 性能与规模实测

| 指标 | 数值 |
|---|---|
| 原始记录（ChEMBL max_phase=4 含 SMILES） | 3,311 |
| 去重后条目 | 3,311（InChIKey 无重复） |
| 可对接条目（SDF+PDBQT 齐全） | 3,282（99.1%） |
| 失败隔离条目 | 29（0.9%：无机盐/含金属/构象嵌入失败） |
| 本地化构建耗时（30 workers） | ~20 min |
| 本地资源占用（SDF+PDBQT） | ~40 MB |
| 离线查询延迟（索引加载 3,311 行） | <0.1 s |

## 6. 本地化后调用方式（零网络）

```bash
autovs library stats              # 统计（纯本地）
autovs library verify-offline     # 禁网自检
autovs dock --target ...          # 批量对接直接消费本地 PDBQT
```

代码内消费：

```python
from autovs.library import LocalLibrary
lib = LocalLibrary("/home/xiap/bio_project/for_others/AutoMD/AutoVS/data/fda_library/index/fda_library_index.csv")
entry = lib.get("FDA00001_s1")     # 药名、SMILES、类药性、PDBQT 路径
```

## 7. 维护与升级

- **更新库版本**：`autovs library fetch && autovs library build`（增量跳过已有产物）；
- **切换 ZINC FDA 子集**：浏览器下载 `.smi` → `library.yaml` 设 `type: local` + `local_file` → 重新 build；
- **放宽/收紧质子化**：改 `library.yaml` 的 `prepare.ph / max_variants` 后重建（注意断点缓存需清理 `prepared/`）。
