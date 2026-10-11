"""配置体系测试：加载 / 校验 / 错误路径。"""

import pytest

from autovs.config import ConfigError, load_config

ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent


def test_load_full_config(cfg):
    assert cfg.framework["name"] == "AutoVS"
    assert cfg.library.source.type == "chembl"
    assert cfg.target.name == "EGFR_demo"
    assert cfg.docking.exhaustiveness == 8
    assert cfg.gate.score_window == 1.0


def test_tools_paths_exist(cfg):
    # 允许 fpocket 缺失（可选），其余工具必须存在
    missing = [m for m in cfg.tools.missing() if m != "fpocket"]
    assert not missing, f"工具缺失: {missing}"


def test_missing_config_file():
    with pytest.raises(ConfigError):
        load_config("/nonexistent/framework.yaml")


def test_invalid_box_mode(tmp_path):
    bad = tmp_path / "bad_target.yaml"
    bad.write_text(
        "target:\n  name: X\n  pdb_id: 1ABC\n  chains: [A]\n"
        "  reference_ligand: null\n  box:\n    mode: nonsense\n", encoding="utf-8"
    )
    with pytest.raises(ConfigError):
        load_config(ROOT / "config" / "framework.yaml", None, bad)
