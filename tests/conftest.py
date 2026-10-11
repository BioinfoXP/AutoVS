"""AutoVS 测试公共fixture：配置加载（全部本地，无需网络）。"""

from pathlib import Path

import pytest

from autovs.config import load_config

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def cfg():
    return load_config(ROOT / "config" / "framework.yaml", ROOT / "config" / "library.yaml",
                       ROOT / "config" / "targets" / "demo_egfr_1m17.yaml")
