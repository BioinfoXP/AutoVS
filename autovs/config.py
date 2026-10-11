"""AutoVS 配置体系：YAML 加载 + 类型化 dataclass + 严格校验。

配置分三层：
  1. framework.yaml —— 路径 / 外部工具 / 运行时（框架级，一次配置）
  2. library.yaml   —— FDA 库数据源与本地化参数（库级）
  3. targets/*.yaml —— 靶点 / 对接盒 / 质量门 / 对接参数（项目级）
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .utils import expand


class ConfigError(ValueError):
    """配置错误：字段缺失 / 类型错误 / 路径无效。"""


def _load_yaml(path: Path) -> dict[str, Any]:
    path = expand(path)
    if not path.exists():
        raise ConfigError(f"配置文件不存在: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ConfigError(f"配置文件格式错误（应为 YAML 映射）: {path}")
    return data


def _require(section: dict[str, Any], key: str, where: str) -> Any:
    if key not in section or section[key] is None:
        raise ConfigError(f"配置缺少必填字段 {where}.{key}")
    return section[key]


# ---------------- 类型化配置 ----------------
@dataclass
class PathsConfig:
    root: Path
    library_dir: Path
    temp_dir: Path
    output_data_dir: Path
    output_figure_dir: Path


@dataclass
class ToolsConfig:
    python: Path
    vina: Path
    mk_prepare_receptor: Path
    mk_prepare_ligand: Path
    mk_export: Path
    obabel: Path
    fpocket: Path | None = None
    pymol: Path | None = None   # 可选：3D 姿势渲染（缺失时跳过 3D 环节）

    def as_dict(self) -> dict[str, str | None]:
        return {k: (str(v) if v else None) for k, v in self.__dict__.items()}

    def missing(self) -> list[str]:
        return [k for k, v in self.__dict__.items() if v and not Path(v).exists()]


@dataclass
class RuntimeConfig:
    workers: int = 30
    log_level: str = "INFO"
    vina_threads: int = 1
    request_timeout: int = 120
    request_retries: int = 3


@dataclass
class LibrarySourceConfig:
    type: str
    base_url: str | None = None
    max_phase: int = 4
    molecule_type: str = "Small molecule"
    page_size: int = 1000
    zinc_url: str | None = None
    local_file: str | None = None


@dataclass
class LibraryPrepConfig:
    ph: float = 7.4
    ph_tolerance: float = 1.0
    protonation_precision: float = 1.0
    max_variants: int = 1
    max_tautomers_per_state: int = 1
    fragment_parent: bool = True
    embed: str = "etkdg"
    optimize: str = "mmff"
    shard_size: int = 200
    skip_existing: bool = True


@dataclass
class LibraryLocalConfig:
    raw_json: Path
    manifest: Path
    prepared_dir: Path
    index_csv: Path


@dataclass
class LibraryConfig:
    name: str
    version: str
    source: LibrarySourceConfig
    prepare: LibraryPrepConfig
    local: LibraryLocalConfig


@dataclass
class BoxConfig:
    mode: str = "reference_ligand"
    padding: float = 5.0
    min_box: float = 22.0
    max_box: float = 30.0
    center: list[float] | None = None
    size: list[float] | None = None


@dataclass
class TargetConfig:
    name: str
    pdb_id: str | None
    structure_path: Path | None
    chains: list[str]
    reference_ligand: dict[str, Any] | None
    box: BoxConfig


@dataclass
class GateConfig:
    rmsd_threshold: float = 2.0
    exhaustiveness: int = 32
    num_modes: int = 9
    score_window: float = 1.0   # 打分窗口（kcal/mol）：窗口内任一 pose 复现晶体姿势即通过


@dataclass
class DockingParams:
    exhaustiveness: int = 8
    num_modes: int = 9
    score_threshold: float = -7.0
    top_n: int = 50
    redock_seed: int = 20261010


@dataclass
class Config:
    framework: dict[str, Any]
    paths: PathsConfig
    tools: ToolsConfig
    runtime: RuntimeConfig
    library: LibraryConfig
    target: TargetConfig | None = None
    gate: GateConfig | None = None
    docking: DockingParams | None = None

    @property
    def seed(self) -> int:
        return int(self.framework.get("seed", 20261010))


# ---------------- 加载入口 ----------------
def load_config(
    framework_yaml: str | Path,
    library_yaml: str | Path | None = None,
    target_yaml: str | Path | None = None,
) -> Config:
    """加载并校验三层配置，返回类型化 Config。"""
    fw = _load_yaml(framework_yaml)

    paths_raw = _require(fw, "paths", "framework.yaml")
    paths = PathsConfig(**{k: expand(_require(paths_raw, k, "paths")) for k in
                           ("root", "library_dir", "temp_dir", "output_data_dir", "output_figure_dir")})

    tools_raw = _require(fw, "tools", "framework.yaml")
    tools = ToolsConfig(**{k: (expand(v) if v else None) for k, v in tools_raw.items()})

    runtime_raw = fw.get("runtime", {}) or {}
    runtime = RuntimeConfig(**{k: v for k, v in runtime_raw.items() if k in RuntimeConfig.__dataclass_fields__})

    # 库配置（可选：部分命令不需要）
    library: LibraryConfig | None = None
    if library_yaml:
        lb = _load_yaml(library_yaml)
        lsec = _require(lb, "library", "library.yaml")
        src_raw = _require(lsec, "source", "library")
        prep_raw = lsec.get("prepare", {}) or {}
        loc_raw = _require(lsec, "local", "library")
        library = LibraryConfig(
            name=str(_require(lsec, "name", "library")),
            version=str(lsec.get("version", "v1")),
            source=LibrarySourceConfig(**{k: v for k, v in src_raw.items()
                                          if k in LibrarySourceConfig.__dataclass_fields__}),
            prepare=LibraryPrepConfig(**{k: v for k, v in prep_raw.items()
                                         if k in LibraryPrepConfig.__dataclass_fields__}),
            local=LibraryLocalConfig(**{k: expand(_require(loc_raw, k, "library.local"))
                                        for k in ("raw_json", "manifest", "prepared_dir", "index_csv")}),
        )
        if library.source.type not in {"chembl", "zinc", "local"}:
            raise ConfigError(f"library.source.type 无效: {library.source.type}（可选 chembl/zinc/local）")

    # 靶点配置（可选）
    target: TargetConfig | None = None
    gate: GateConfig | None = None
    docking: DockingParams | None = None
    if target_yaml:
        tg = _load_yaml(target_yaml)
        tsec = _require(tg, "target", "target.yaml")
        box_raw = tsec.get("box", {}) or {}
        box = BoxConfig(**{k: v for k, v in box_raw.items() if k in BoxConfig.__dataclass_fields__})
        if box.mode not in {"reference_ligand", "fpocket", "manual"}:
            raise ConfigError(f"target.box.mode 无效: {box.mode}")
        if box.mode == "manual" and not (box.center and box.size):
            raise ConfigError("target.box.mode=manual 时必须提供 center 与 size")
        target = TargetConfig(
            name=str(_require(tsec, "name", "target")),
            pdb_id=tsec.get("pdb_id"),
            structure_path=expand(tsec["structure_path"]) if tsec.get("structure_path") else None,
            chains=list(tsec.get("chains") or []),
            reference_ligand=tsec.get("reference_ligand"),
            box=box,
        )
        if not target.structure_path and not target.pdb_id:
            raise ConfigError("target 必须提供 pdb_id 或 structure_path")
        gate = GateConfig(**{k: v for k, v in (tg.get("gate", {}) or {}).items()
                             if k in GateConfig.__dataclass_fields__})
        docking = DockingParams(**{k: v for k, v in (tg.get("docking", {}) or {}).items()
                                   if k in DockingParams.__dataclass_fields__})

    return Config(
        framework=fw.get("framework", {}),
        paths=paths,
        tools=tools,
        runtime=runtime,
        library=library,  # type: ignore[arg-type]
        target=target,
        gate=gate,
        docking=docking,
    )
