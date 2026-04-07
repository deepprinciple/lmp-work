from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_main():
    source_path = Path(__file__).resolve().parents[1] / "bench_runs" / "generate_5mol_baseline_configs.py"
    spec = importlib.util.spec_from_file_location("generate_5mol_baseline_configs_legacy", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load legacy script: {source_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.main


if __name__ == "__main__":
    _load_main()()
