"""FDA 库本地化测试：索引完整性 / 离线可用性（禁网）/ 抽样确定性。"""

import socket

from autovs.library import LocalLibrary


def test_index_loaded(cfg):
    lib = LocalLibrary(cfg.library.local.index_csv)
    stats = lib.stats()
    assert stats["total"] > 3000, "FDA 库应包含 3000+ 已上市药物"
    assert stats["dockable"] > stats["total"] * 0.95, "可对接比例应 ≥95%"


def test_entry_lookup(cfg):
    lib = LocalLibrary(cfg.library.local.index_csv)
    entry = next(iter(lib.iter_entries()))
    found = lib.get(entry.entry_id)
    assert found.entry_id == entry.entry_id
    assert found.smiles and found.pdbqt_path


def test_sample_deterministic(cfg):
    lib = LocalLibrary(cfg.library.local.index_csv)
    a = [e.entry_id for e in lib.sample(10, seed=42)]
    b = [e.entry_id for e in lib.sample(10, seed=42)]
    assert a == b, "同种子抽样必须确定性一致"


def test_offline_access(cfg, monkeypatch):
    """离线验证：禁网后本地库仍可读取且产物文件齐全。"""
    def _blocked(*args, **kwargs):
        raise RuntimeError("测试期间禁止网络访问")

    monkeypatch.setattr(socket, "socket", _blocked)
    monkeypatch.setattr(socket, "create_connection", _blocked)

    lib = LocalLibrary(cfg.library.local.index_csv)
    missing = 0
    for entry in lib.sample(20, seed=cfg.seed):
        if not entry.pdbqt_path or not entry.sdf_path:
            missing += 1
    assert missing == 0, "抽样条目应全部具备本地 SDF/PDBQT 产物"
