"""AutoVS 通用工具：日志 / 子进程 / JSON / 路径 / 缓存指纹。"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
from pathlib import Path
from typing import Any


# ---------- 日志 ----------
def setup_logging(level: str = "INFO") -> logging.Logger:
    """初始化全局日志（控制台输出，中文友好格式）。"""
    logging.basicConfig(
        level=getattr(logging, str(level).upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        force=True,
    )
    return logging.getLogger("autovs")


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"autovs.{name}")


# ---------- 路径 ----------
def ensure_dir(path: Path) -> Path:
    """创建目录并返回 Path。"""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def expand(path_str: str | Path) -> Path:
    """展开 ~ 并转为绝对路径。"""
    return Path(str(path_str)).expanduser().resolve()


# ---------- JSON ----------
def write_json(path: Path, data: Any) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


# ---------- 配置指纹（断点续跑判据） ----------
def config_hash(payload: Any) -> str:
    """对配置片段做稳定哈希，用于缓存有效性判断。"""
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=False)
    return hashlib.sha1(blob.encode("utf-8")).hexdigest()[:12]


# ---------- 子进程 ----------
class CommandError(RuntimeError):
    pass


def run_cmd(
    cmd: list[str],
    *,
    cwd: Path | None = None,
    log: Path | None = None,
    check: bool = True,
    timeout: int | None = None,
) -> subprocess.CompletedProcess:
    """执行外部命令，输出写入日志文件；失败时抛出含日志尾部的异常。"""
    log_path = Path(log) if log else None
    if log_path:
        ensure_dir(log_path.parent)
    proc = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if log_path:
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"$ {' '.join(str(c) for c in cmd)}\n")
            fh.write(proc.stdout or "")
            fh.write(proc.stderr or "")
            fh.write("\n")
    if check and proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-5:] + (proc.stderr or "").splitlines()[-5:])
        raise CommandError(f"命令失败（exit={proc.returncode}）: {' '.join(str(c) for c in cmd)}\n{tail}")
    return proc
