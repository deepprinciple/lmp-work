"""CLI for generating LAMMPS input files."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from legacy.input_gk import build_replica_input_text
from scripts._cli import Argv, new_parser
from workflow.input_thermal import build_equilibration_input_text


def build_parser() -> argparse.ArgumentParser:
    parser = new_parser("Write LAMMPS input files for equilibration or GK replica runs.")
    parser.add_argument("--workdir", required=True, help="Case directory")
    parser.add_argument(
        "--mode",
        choices=["equilibration", "replica"],
        default="equilibration",
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=300.0,
        help="Simulation temperature in K",
    )
    parser.add_argument("--seed", type=int, default=20260316, help="Random seed")
    parser.add_argument(
        "--npt-pre-50-steps",
        type=int,
        default=200000,
        help="50 atm NPT steps",
    )
    parser.add_argument(
        "--npt-pre-10-steps",
        type=int,
        default=300000,
        help="10 atm NPT steps",
    )
    parser.add_argument(
        "--npt-1atm-steps",
        type=int,
        default=2000000,
        help="1 atm NPT steps",
    )
    parser.add_argument(
        "--nvt-steps",
        type=int,
        default=1000000,
        help="NVT equilibration steps",
    )
    parser.add_argument(
        "--decorrelation-steps",
        type=int,
        default=500000,
        help="Replica decorrelation steps",
    )
    parser.add_argument(
        "--gk-steps",
        type=int,
        default=10000000,
        help="NVE+GK production steps",
    )
    parser.add_argument(
        "--timestep-fs",
        type=float,
        default=1.0,
        help="Time step in fs",
    )
    parser.add_argument(
        "--sample-every-steps",
        type=int,
        default=5,
        help="ACF sampling interval in steps",
    )
    parser.add_argument(
        "--corr-points",
        type=int,
        default=8000,
        help="Number of correlation points",
    )
    parser.add_argument(
        "--equil-restart-relpath",
        default="../equil_nvt.restart",
        help="Relative path to equilibration restart for replica mode",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output input filename",
    )
    return parser


def parse_args(argv: Argv = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main(argv: Argv = None) -> int:
    args = parse_args(argv)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "equilibration":
        text = build_equilibration_input_text(
            temperature=args.temperature_k,
            seed=args.seed,
            npt_pre_50_steps=args.npt_pre_50_steps,
            npt_pre_10_steps=args.npt_pre_10_steps,
            npt_1atm_steps=args.npt_1atm_steps,
            nvt_steps=args.nvt_steps,
            timestep_fs=args.timestep_fs,
        )
        out_name = args.output or "in.equil.lammps"
    else:
        text = build_replica_input_text(
            temperature=args.temperature_k,
            replica_seed=args.seed,
            decorrelation_steps=args.decorrelation_steps,
            gk_steps=args.gk_steps,
            timestep_fs=args.timestep_fs,
            sample_every_steps=args.sample_every_steps,
            corr_points=args.corr_points,
            equil_restart_relpath=args.equil_restart_relpath,
        )
        out_name = args.output or "in.replica.lammps"

    out_file = workdir / out_name
    out_file.write_text(text)
    print(f"Wrote {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
