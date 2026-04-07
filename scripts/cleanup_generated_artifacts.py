from __future__ import annotations

import argparse
import shutil
from pathlib import Path


DEFAULT_PATTERNS = [
    "**/run_*.log",
    "**/packmol.log",
    "**/packmol_try*.log",
    "**/*live_convergence*.json",
    "**/*live_latest.json",
    "**/*.dat",
    "**/*.pdb",
    "**/*.xyz",
    "**/.openff_cache",
]


def iter_matches(root: Path, patterns: list[str]) -> list[Path]:
    seen: set[Path] = set()
    matches: list[Path] = []
    for pattern in patterns:
        for path in root.glob(pattern):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            matches.append(path)
    return sorted(matches)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dry-run or delete generated bench-run artifacts.")
    parser.add_argument(
        "--root",
        default="/root/lmp-work/bench_runs",
        help="Root directory to scan for generated artifacts",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually delete matched files. Default is dry-run only.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    matches = iter_matches(root, DEFAULT_PATTERNS)
    print(f"scan_root: {root}")
    print(f"matches: {len(matches)}")
    for path in matches:
        print(path)
    if not args.apply:
        print("dry-run only; pass --apply to delete")
        return

    for path in matches:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    print("cleanup complete")


if __name__ == "__main__":
    main()
