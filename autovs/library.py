"""FDA 已上市药物库——本地化核心模块。

设计原则（本地化适配）：
  1. 网络访问只存在于 fetch 阶段（可插拔数据源适配器）；
     fetch 产物落盘 raw/source_records.json + manifest.json。
  2. build / query / dock 全部基于本地文件，运行期零网络依赖。
  3. 本地索引 = 单一 CSV（fda_library_index.csv），一行一个可对接条目。
"""

from __future__ import annotations

import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

from .config import Config
from .ligand import build_entry
from .utils import config_hash, ensure_dir, get_logger, read_json, write_json

log = get_logger("library")


class LibraryError(RuntimeError):
    pass


# =====================================================================
# 1. 数据源适配器（唯一联网层）
# =====================================================================
def _http_get(url: str, *, timeout: int, retries: int, params: dict | None = None) -> bytes:
    """带重试的 HTTP GET（requests），供数据源适配器使用。"""
    import requests

    last_exc: Exception | None = None
    for attempt in range(1, int(retries) + 1):
        try:
            resp = requests.get(url, params=params, timeout=timeout,
                                headers={"User-Agent": "AutoVS/0.1 (localization pipeline)"})
            ctype = resp.headers.get("Content-Type", "")
            if resp.status_code == 200 and "html" not in ctype.lower():
                return resp.content
            if resp.status_code == 200 and "html" in ctype.lower():
                # 反爬验证页（如 ZINC）：返回 HTML 而非数据
                raise LibraryError(
                    f"数据源返回 HTML 而非数据（疑似反爬/验证页）: {url}\n"
                    f"请改用 browser 下载后以 source.type=local 导入，或切换 chembl 源。"
                )
            raise LibraryError(f"HTTP {resp.status_code}: {url}")
        except Exception as exc:  # 网络抖动重试
            last_exc = exc
            log.warning("下载失败（第 %d/%d 次）: %s", attempt, retries, exc)
            time.sleep(2 * attempt)
    raise LibraryError(f"下载失败（已重试 {retries} 次）: {url}\n{last_exc}")


def fetch_chembl(source_cfg, timeout: int, retries: int) -> list[dict[str, Any]]:
    """ChEMBL REST 分页拉取 max_phase=4 已上市小分子（FDA/EMA 等主要监管批准）。"""
    records: list[dict[str, Any]] = []
    offset = 0
    page_size = int(source_cfg.page_size)
    while True:
        url = f"{source_cfg.base_url}/molecule.json"
        params = {
            "max_phase": source_cfg.max_phase,
            "molecule_type": source_cfg.molecule_type,
            "limit": page_size,
            "offset": offset,
        }
        blob = _http_get(url, timeout=timeout, retries=retries, params=params)
        page = json.loads(blob)
        mols = page.get("molecules", [])
        for m in mols:
            structures = m.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            if not smiles:
                continue
            records.append({
                "source": "chembl",
                "source_id": m.get("molecule_chembl_id"),
                "name": m.get("pref_name") or "",
                "smiles": smiles,
                "inchi_key": structures.get("standard_inchi_key") or "",
                "first_approval": m.get("first_approval"),
                "atc": ";".join(m.get("atc_classifications") or []),
            })
        total = int(page.get("page_meta", {}).get("total_count", 0))
        offset += len(mols)
        log.info("ChEMBL 拉取进度: %d / %d", offset, total)
        if offset >= total or not mols:
            break
    return records


def fetch_zinc(source_cfg, timeout: int, retries: int) -> list[dict[str, Any]]:
    """ZINC FDA 子集 SMILES 下载（站点若启用反爬验证会给出明确指引）。"""
    blob = _http_get(source_cfg.zinc_url, timeout=timeout, retries=retries)
    records = []
    for line in blob.decode("utf-8", errors="ignore").splitlines():
        parts = line.split()
        if not parts:
            continue
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else ""
        records.append({"source": "zinc", "source_id": name or "", "name": name,
                        "smiles": smiles, "inchi_key": "", "first_approval": None, "atc": ""})
    return records


def fetch_local(source_cfg) -> list[dict[str, Any]]:
    """导入本地 SDF / SMI 文件（浏览器下载等手动渠道的统一入口）。"""
    from rdkit import Chem

    path = Path(str(source_cfg.local_file)).expanduser().resolve()
    if not path.exists():
        raise LibraryError(f"本地库文件不存在: {path}")
    records: list[dict[str, Any]] = []
    if path.suffix.lower() in {".sdf", ".sd", ".mol"}:
        for i, mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=False), 1):
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"local_{i:05d}"
            records.append({"source": "local", "source_id": name, "name": name, "smiles": smiles,
                            "inchi_key": "", "first_approval": None, "atc": ""})
    else:  # .smi / .smiles / .txt / .csv 首两列 = SMILES 名称
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            parts = line.replace(",", " ").split()
            if not parts:
                continue
            records.append({"source": "local", "source_id": parts[1] if len(parts) > 1 else "",
                            "name": parts[1] if len(parts) > 1 else "", "smiles": parts[0],
                            "inchi_key": "", "first_approval": None, "atc": ""})
    return records


def library_fetch(cfg: Config) -> Path:
    """执行数据源拉取（唯一联网步骤），落盘 raw/source_records.json。"""
    assert cfg.library is not None
    src = cfg.library.source
    ensure_dir(cfg.library.local.raw_json.parent)
    log.info("数据源: %s（联网拉取，仅此步骤依赖远程接口）", src.type)
    if src.type == "chembl":
        records = fetch_chembl(src, cfg.runtime.request_timeout, cfg.runtime.request_retries)
    elif src.type == "zinc":
        records = fetch_zinc(src, cfg.runtime.request_timeout, cfg.runtime.request_retries)
    elif src.type == "local":
        records = fetch_local(src)
    else:
        raise LibraryError(f"未知数据源: {src.type}")

    write_json(cfg.library.local.raw_json, records)
    write_json(cfg.library.local.manifest, {
        "library": cfg.library.name,
        "version": cfg.library.version,
        "source_type": src.type,
        "source_url": src.base_url or src.zinc_url or src.local_file,
        "record_count": len(records),
        "raw_json": str(cfg.library.local.raw_json),
        "fetched_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "source_config_hash": config_hash(src.__dict__),
    })
    log.info("本地化原始库完成: %d 条记录 -> %s", len(records), cfg.library.local.raw_json)
    return cfg.library.local.raw_json


# =====================================================================
# 2. 本地化构建（纯离线）：去重 → 3D → 质子化 → PDBQT → 索引
# =====================================================================
INDEX_FIELDS = [
    "entry_id", "lib_id", "state_id", "source", "source_id", "name",
    "smiles", "inchi_key", "formal_charge", "mol_weight", "logp",
    "hbd", "hba", "rotatable_bonds", "heavy_atoms",
    "first_approval", "atc", "status", "sdf_path", "pdbqt_path", "error",
]


def _props_from_smiles(smiles: str) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mol_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
    }


def _build_worker(args: tuple) -> dict[str, Any]:
    """进程池 worker：单条目构建（3D + 质子化 + PDBQT）。"""
    entry_id, smiles, out_dir, prep_dict, tools_dict, seed = args
    from .config import LibraryPrepConfig, ToolsConfig
    return build_entry(
        entry_id, smiles, Path(out_dir),
        LibraryPrepConfig(**prep_dict), ToolsConfig(**tools_dict), seed=seed,
    )


def library_build(cfg: Config, limit: int | None = None) -> Path:
    """构建本地化 FDA 库（纯离线）：生成 prepared/<shard>/ 与本地索引 CSV。"""
    assert cfg.library is not None
    if not cfg.library.local.raw_json.exists():
        raise LibraryError("缺少原始库文件，请先执行: autovs library fetch")
    records = read_json(cfg.library.local.raw_json)
    if limit:
        records = records[: int(limit)]

    # 1) 去重（InChIKey 优先，退化为规范 SMILES）
    seen: set[str] = set()
    unique: list[dict[str, Any]] = []
    for rec in records:
        key = (rec.get("inchi_key") or "").strip() or rec["smiles"]
        if key in seen:
            continue
        seen.add(key)
        unique.append(rec)
    log.info("原始 %d 条 -> 去重后 %d 条", len(records), len(unique))

    # 2) 任务分片
    shard_size = int(cfg.library.prepare.shard_size)
    prepared_dir = ensure_dir(cfg.library.local.prepared_dir)
    jobs = []
    for i, rec in enumerate(unique):
        lib_id = f"FDA{i + 1:05d}"
        entry_id = f"{lib_id}_s1"
        shard = f"shard_{i // shard_size:03d}"
        jobs.append((entry_id, rec["smiles"], str(prepared_dir / shard),
                     cfg.library.prepare.__dict__, cfg.tools.as_dict(), cfg.seed))

    # 3) 并行构建（失败隔离：单分子失败不中断）
    results: dict[str, dict[str, Any]] = {}
    workers = max(1, int(cfg.runtime.workers))
    log.info("开始并行准备 %d 个化合物（workers=%d）", len(jobs), workers)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_build_worker, job): job[0] for job in jobs}
        done = 0
        for fut in as_completed(futures):
            res = fut.result()
            results[res["entry_id"]] = res
            done += 1
            if done % 200 == 0 or done == len(jobs):
                log.info("配体准备进度: %d / %d", done, len(jobs))

    # 4) 写本地索引 CSV
    ensure_dir(cfg.library.local.index_csv.parent)
    with cfg.library.local.index_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=INDEX_FIELDS)
        writer.writeheader()
        for i, rec in enumerate(unique):
            lib_id = f"FDA{i + 1:05d}"
            entry_id = f"{lib_id}_s1"
            res = results.get(entry_id, {"status": "missing"})
            row = {
                "entry_id": entry_id, "lib_id": lib_id, "state_id": "s1",
                "source": rec.get("source", ""), "source_id": rec.get("source_id", ""),
                "name": rec.get("name", ""), "smiles": res.get("smiles") or rec["smiles"],
                "inchi_key": rec.get("inchi_key", ""),
                "formal_charge": res.get("formal_charge", ""),
                "first_approval": rec.get("first_approval") or "",
                "atc": rec.get("atc", ""),
                "status": res.get("status", "missing"),
                "sdf_path": res.get("sdf_path", ""), "pdbqt_path": res.get("pdbqt_path", ""),
                "error": res.get("error", ""),
            }
            row.update(_props_from_smiles(rec["smiles"]))
            writer.writerow({k: row.get(k, "") for k in INDEX_FIELDS})

    n_ok = sum(1 for r in results.values() if r.get("status") in {"ok", "cached"})
    log.info("本地化构建完成: 成功 %d / %d，索引 -> %s", n_ok, len(unique), cfg.library.local.index_csv)
    return cfg.library.local.index_csv


# =====================================================================
# 3. 本地库访问（100% 离线）
# =====================================================================
@dataclass
class LibraryEntry:
    entry_id: str
    lib_id: str
    name: str
    smiles: str
    inchi_key: str
    status: str
    pdbqt_path: str
    sdf_path: str
    row: dict[str, str]


class LocalLibrary:
    """本地 FDA 库离线访问接口（仅读本地 CSV 与产物文件）。"""

    def __init__(self, index_csv: Path):
        self.index_csv = Path(index_csv)
        if not self.index_csv.exists():
            raise LibraryError(f"本地索引不存在: {self.index_csv}（先执行 autovs library build）")
        with self.index_csv.open(encoding="utf-8") as fh:
            self._rows = list(csv.DictReader(fh))

    def __len__(self) -> int:
        return len(self._rows)

    def iter_entries(self, status: str = "ok") -> Iterator[LibraryEntry]:
        """遍历条目；status='ok' 时同时接受 'cached'，'all' 时不过滤。"""
        for row in self._rows:
            row_status = row.get("status", "")
            if status == "all" or row_status == status or (status == "ok" and row_status == "cached"):
                yield self._to_entry(row)

    def _to_entry(self, row: dict[str, str]) -> LibraryEntry:
        return LibraryEntry(
            entry_id=row["entry_id"], lib_id=row["lib_id"], name=row.get("name", ""),
            smiles=row.get("smiles", ""), inchi_key=row.get("inchi_key", ""),
            status=row.get("status", ""), pdbqt_path=row.get("pdbqt_path", ""),
            sdf_path=row.get("sdf_path", ""), row=row,
        )

    def get(self, entry_id: str) -> LibraryEntry:
        for row in self._rows:
            if row["entry_id"] == entry_id:
                return self._to_entry(row)
        raise LibraryError(f"条目不存在: {entry_id}")

    def sample(self, n: int, seed: int = 20261010) -> list[LibraryEntry]:
        """确定性抽样 n 个可用于对接的条目。"""
        import random

        usable = [r for r in self._rows if r.get("status") in {"ok", "cached"} and r.get("pdbqt_path")]
        rng = random.Random(int(seed))
        rng.shuffle(usable)
        return [self._to_entry(r) for r in usable[: int(n)]]

    def stats(self) -> dict[str, Any]:
        total = len(self._rows)
        ok = sum(1 for r in self._rows if r.get("status") in {"ok", "cached"})
        failed = sum(1 for r in self._rows if r.get("status") == "failed")
        sources: dict[str, int] = {}
        for r in self._rows:
            sources[r.get("source", "?")] = sources.get(r.get("source", "?"), 0) + 1
        return {"total": total, "dockable": ok, "failed": failed, "by_source": sources,
                "index_csv": str(self.index_csv)}
