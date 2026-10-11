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


def test_run_batch_prefilter_outlier(cfg, tmp_path, monkeypatch):
    """离群大分子预过滤：超 MW/可旋转键阈值的条目提交前跳过并写 result.json。"""
    from autovs import docking as dock
    from autovs.library import LibraryEntry

    # 构造三个条目：正常 / MW 超限 / 可旋转键超限
    def make_entry(eid, mw, rot):
        return LibraryEntry(
            entry_id=eid, lib_id=eid, name=eid, smiles="", inchi_key="",
            status="ok",
            pdbqt_path=str(tmp_path / "fake.pdbqt"),  # 不存在 -> 正常条目也会因缺文件跳过
            sdf_path="", row={"mol_weight": str(mw), "rotatable_bonds": str(rot)},
        )

    entries = [make_entry("OK_ONE", 500, 8),
               make_entry("BIG_MW", 1620.69, 35),
               make_entry("BIG_ROT", 900, 50)]
    for e in entries:  # 全部带真实 pdbqt 文件：正常条目可进提交队列
        p = tmp_path / f"{e.entry_id}.pdbqt"
        p.write_text("REMARK\n")
        e.pdbqt_path = str(p)

    captured = []
    monkeypatch.setattr(dock, "_dock_entry_job", lambda job: captured.append(job) or
                        {"entry_id": job[0], "status": "ok", "config_hash": job[11],
                         "top_score": -9.0, "scores": [-9.0], "pose_pdbqt": "", "pose_sdf": ""})

    results = dock.run_batch(cfg, entries, tmp_path / "rec.pdbqt",
                             [0.0, 0.0, 0.0], [22.0, 22.0, 22.0], tmp_path / "out")

    statuses = {r["entry_id"]: r["status"] for r in results}
    assert statuses["BIG_MW"] == "skipped" and "MW" in \
        next(r["reason"] for r in results if r["entry_id"] == "BIG_MW")
    assert statuses["BIG_ROT"] == "skipped" and "rotatable_bonds" in \
        next(r["reason"] for r in results if r["entry_id"] == "BIG_ROT")
    assert statuses["OK_ONE"] != "skipped"
    # 被跳过的条目必须写 result.json 且不进入提交队列
    assert (tmp_path / "out" / "BIG_MW" / "result.json").exists()
    assert all(job[0] not in {"BIG_MW", "BIG_ROT"} for job in captured)


def test_dock_single_timeout_raises_docking_error(tmp_path, monkeypatch):
    """超时兜底：vina 超时被转换为 DockingError（失败隔离可捕获）。"""
    import subprocess as sp

    import pytest

    from autovs import docking as dock

    def fake_run_cmd(cmd, **kwargs):
        raise sp.TimeoutExpired(cmd, kwargs.get("timeout") or 1)

    monkeypatch.setattr(dock, "run_cmd", fake_run_cmd)

    class FakeTools:
        vina = "/bin/true"
        mk_export = "/bin/true"

    with pytest.raises(dock.DockingError, match="超时"):
        dock.dock_single(
            receptor_pdbqt="/dev/null", ligand_pdbqt="/dev/null",
            center=[0.0, 0.0, 0.0], size=[22.0, 22.0, 22.0],
            out_dir=tmp_path / "t", exhaustiveness=8, num_modes=9, seed=1,
            tools=FakeTools(), vina_threads=1, timeout_s=1800,
        )
