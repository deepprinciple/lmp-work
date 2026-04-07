from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any


def copy_if_needed(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    shutil.copy2(src, dst)


def prepare_branch_workdir(
    *,
    source_workdir: Path,
    branch_workdir: Path,
    source_restart_name: str = "equil_nvt.restart",
    source_data_name: str = "system.data",
) -> tuple[Path, Path]:
    branch_workdir.mkdir(parents=True, exist_ok=True)
    source_restart = source_workdir / source_restart_name
    source_data = source_workdir / source_data_name
    if not source_restart.exists():
        raise FileNotFoundError(f"Missing source restart: {source_restart}")
    if not source_data.exists():
        raise FileNotFoundError(f"Missing source system.data: {source_data}")

    copied_restart = branch_workdir / "equil_nvt.restart"
    copied_data = branch_workdir / "system.data"
    copy_if_needed(source_restart, copied_restart)
    copy_if_needed(source_data, copied_data)
    return copied_restart, copied_data


def write_branch_manifest(
    *,
    branch_workdir: Path,
    filename: str,
    payload: dict[str, Any],
) -> Path:
    out_file = branch_workdir / filename
    out_file.write_text(json.dumps(payload, indent=2))
    return out_file
