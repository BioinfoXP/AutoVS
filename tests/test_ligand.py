"""配体准备测试：质子化 / 3D 构象 / PDBQT 产物链。"""

from autovs.ligand import build_entry, enumerate_states


def test_enumerate_states_aspirin(cfg):
    states = enumerate_states("CC(=O)Oc1ccccc1C(=O)O", cfg.library.prepare)
    assert states, "应至少生成 1 个质子化状态"
    assert states[0]["smiles"]


def test_build_entry_outputs(cfg, tmp_path):
    res = build_entry("TEST00001_s1", "CC(=O)Oc1ccccc1C(=O)O", tmp_path,
                      cfg.library.prepare, cfg.tools, seed=cfg.seed)
    assert res["status"] in {"ok", "cached"}, res.get("error")
    assert (tmp_path / "TEST00001_s1.sdf").exists()
    assert (tmp_path / "TEST00001_s1.pdbqt").exists()


def test_build_entry_failure_isolated(cfg, tmp_path):
    """非法 SMILES 失败隔离：返回 failed 而不是抛异常。"""
    res = build_entry("BAD00001_s1", "not_a_smiles", tmp_path,
                      cfg.library.prepare, cfg.tools, seed=cfg.seed)
    assert res["status"] == "failed"
    assert res["error"]
