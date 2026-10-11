"""对接引擎测试：打分解析 / 质量门判据逻辑（纯本地）。"""

from pathlib import Path

from autovs.docking import vina_scores
from autovs.receptor import evaluate_redocking


def test_vina_scores_parse(tmp_path):
    pdbqt = tmp_path / "poses.pdbqt"
    pdbqt.write_text(
        "REMARK VINA RESULT:    -8.500      0.000      0.000\n"
        "REMARK VINA RESULT:    -7.200      1.500      3.100\n",
        encoding="utf-8",
    )
    assert vina_scores(pdbqt) == [-8.5, -7.2]


def test_evaluate_redocking_window_criterion(cfg, tmp_path):
    """打分窗口判据链路：RMSD 计算、窗口归组、PASS 判定。"""
    from rdkit import Chem

    ref = tmp_path / "ref.sdf"
    poses = tmp_path / "poses.sdf"
    # 用库内真实分子构造：ref 与 pose 相同（RMSD=0），验证 RMSD 计算链路
    lib_ref = Path(cfg.library.local.index_csv).parent.parent / "prepared"
    sample_sdf = next(lib_ref.rglob("*.sdf"))
    mol = next(m for m in Chem.SDMolSupplier(str(sample_sdf), removeHs=True) if m)
    for path in (ref, poses):
        w = Chem.SDWriter(str(path))
        w.write(mol)
        w.write(mol)  # 两个相同 pose
        w.close()

    ev = evaluate_redocking(ref, poses, [-9.0, -8.2], rmsd_threshold=2.0, score_window=1.0)
    assert ev["passed"] is True
    assert ev["min_rmsd_in_window"] == 0.0
    assert ev["top_pose_rmsd"] == 0.0
    assert len(ev["per_mode"]) == 2
