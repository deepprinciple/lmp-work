"""CLI for HFACF thermal-conductivity post-processing."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.hfacf import analyze_hfacf_file, load_volume_from_gk_data
from scripts._cli import Argv, new_parser


def build_parser() -> argparse.ArgumentParser:
    parser = new_parser(
        "Analyze HFACF data and plot thermal-conductivity convergence."
    )
    parser.add_argument("--hfacf-file", required=True, help="LAMMPS HFACF file")
    parser.add_argument(
        "--gk-data-file",
        default=None,
        help="Optional gk_data.dat file for automatic average volume",
    )
    parser.add_argument(
        "--volume-a3",
        type=float,
        default=None,
        help="Average volume in A^3 when gk_data.dat is absent",
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        required=True,
        help="Simulation temperature in K",
    )
    parser.add_argument(
        "--t-max-ps",
        type=float,
        default=20.0,
        help="Maximum correlation time used for integration in ps",
    )
    parser.add_argument(
        "--plateau-start-ps",
        type=float,
        default=2.0,
        help="Start of the plateau analysis window in ps",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=11,
        help="Running-average window size",
    )
    parser.add_argument(
        "--json-output",
        default="thermal_hfacf_summary.json",
        help="JSON summary filename",
    )
    parser.add_argument(
        "--plot-output",
        default="thermal_hfacf_analysis.png",
        help="Plot filename",
    )
    return parser


def parse_args(argv: Argv = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main(argv: Argv = None) -> int:
    args = parse_args(argv)
    hfacf = Path(args.hfacf_file).resolve()

    if args.gk_data_file is not None:
        volume = load_volume_from_gk_data(Path(args.gk_data_file).resolve())
    else:
        if args.volume_a3 is None:
            raise ValueError("Either --gk-data-file or --volume-a3 must be provided.")
        volume = float(args.volume_a3)

    summary = analyze_hfacf_file(
        hfacf,
        temperature_k=float(args.temperature_k),
        volume_a3=volume,
        t_max_ps=float(args.t_max_ps),
        plateau_start_ps=float(args.plateau_start_ps),
        smooth_window=int(args.smooth_window),
        json_out=args.json_output,
        plot_out=args.plot_output,
    )
    print(json.dumps(summary, indent=2))
    return 0


__all__ = [
    "analyze_hfacf_file",
    "build_parser",
    "load_volume_from_gk_data",
    "main",
    "parse_args",
]


if __name__ == "__main__":
    raise SystemExit(main())
