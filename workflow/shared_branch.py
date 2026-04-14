"""Branch manifest helpers used by the thermal mainline."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_branch_manifest(
    *,
    branch_workdir: Path,
    filename: str,
    payload: dict[str, Any],
) -> Path:
    out_file = branch_workdir / filename
    out_file.write_text(json.dumps(payload, indent=2))
    return out_file


__all__ = ["write_branch_manifest"]
