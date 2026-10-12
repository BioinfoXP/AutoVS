"""任务打包测试：任务目录 zip 交付包（纯本地）。"""

import zipfile

import pytest

from autovs.packaging import package_task


class _Paths:
    def __init__(self, root):
        self.output_figure_dir = root


class _Target:
    name = "DEMO"


class _Cfg:
    def __init__(self, root):
        self.paths = _Paths(root)
        self.target = _Target()


def test_package_task_zips_task_dir(tmp_path):
    """任务目录 -> zip：包内保留 <tag>/ 一级目录，文件计数正确，重复打包幂等。"""
    cfg = _Cfg(tmp_path)
    task = tmp_path / "DEMO" / "full"
    task.mkdir(parents=True)
    (task / "report.md").write_text("# demo", encoding="utf-8")
    (task / "full_ranking.csv").write_text("rank,name\n1,ASPIRIN\n", encoding="utf-8")
    sub = task / "fig"
    sub.mkdir()
    (sub / "pose.png").write_bytes(b"\x89PNG\r\n")

    r = package_task(cfg, "full")
    assert r["n_files"] == 3
    assert r["size_mb"] >= 0
    # 默认命名：<靶点>_<批次>_results_<YYYYMMDD>.zip
    zips = list((tmp_path / "DEMO").glob("DEMO_full_results_*.zip"))
    assert len(zips) == 1
    with zipfile.ZipFile(zips[0]) as zf:
        names = set(zf.namelist())
    assert {"full/report.md", "full/full_ranking.csv", "full/fig/pose.png"} <= names

    # 幂等：重复打包覆盖，不报错、不翻倍
    r2 = package_task(cfg, "full")
    assert r2["n_files"] == 3
    assert len(list((tmp_path / "DEMO").glob("DEMO_full_results_*.zip"))) == 1


def test_package_task_custom_out_and_missing_dir(tmp_path):
    """自定义 --out 路径生效；任务目录缺失时报 FileNotFoundError。"""
    cfg = _Cfg(tmp_path)
    with pytest.raises(FileNotFoundError):
        package_task(cfg, "nofull")

    task = tmp_path / "DEMO" / "screen"
    task.mkdir(parents=True)
    (task / "a.txt").write_text("x", encoding="utf-8")
    out = tmp_path / "elsewhere" / "pkg.zip"
    out.parent.mkdir()
    r = package_task(cfg, "screen", out=out)
    assert r["zip"] == str(out)
    assert out.exists()
