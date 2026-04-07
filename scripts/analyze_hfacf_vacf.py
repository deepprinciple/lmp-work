from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_main():
    source_path = Path(__file__).resolve().parents[1] / "analyze_hfacf_vacf.py"
    spec = importlib.util.spec_from_file_location("analyze_hfacf_vacf_legacy", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load legacy script: {source_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.main


if __name__ == "__main__":
    _load_main()()
