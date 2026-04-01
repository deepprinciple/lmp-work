"""
Minimal Green-Kubo analysis script for hfacf.dat + vacf.dat.

Example:
    python analyze_gk.py \
    --hfacf hfacf.dat \
    --vacf vacf.dat \
    --gk_data gk_data.dat \
    --T 300
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from gk_workflow.analysis import GreenKuboAnalyzer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze thermal conductivity and viscosity from hfacf.dat + vacf.dat"
    )
    parser.add_argument("--hfacf", required=True, help="LAMMPS heat-flux ACF file")
    parser.add_argument("--vacf", required=True, help="LAMMPS stress ACF file")
    parser.add_argument("--T", required=True, type=float, help="Temperature in K")
    parser.add_argument("--gk_data", default=None, help="Optional gk_data.dat for automatic average volume")
    parser.add_argument("--volume", type=float, default=None, help="Average volume in A^3 (used when --gk_data is absent)")
    parser.add_argument("--visc_t_max_ps", type=float, default=20.0, help="Viscosity integration cutoff (ps)")
    parser.add_argument("--kappa_t_max_ps", type=float, default=20.0, help="Thermal conductivity integration cutoff (ps)")
    parser.add_argument("--visc_plateau_start_ps", type=float, default=5.0, help="Viscosity plateau start time (ps)")
    parser.add_argument("--kappa_plateau_start_ps", type=float, default=2.0, help="Thermal conductivity plateau start time (ps)")
    parser.add_argument("--smooth_window", type=int, default=50, help="Running-integral smoothing window")
    parser.add_argument("--n_blocks", type=int, default=5, help="Block count for SEM estimation")
    parser.add_argument("--out", default="result.json", help="Output JSON filename")
    return parser.parse_args()


def load_volume_from_gk_data(gk_data_file: Path) -> float:
    rows = []
    with open(gk_data_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            rows.append(parts)

    if not rows:
        raise ValueError(f"No numeric rows found in gk_data file: {gk_data_file}")

    volumes = [float(row[2]) for row in rows]
    return sum(volumes) / len(volumes)


def main() -> None:
    args = parse_args()

    hfacf = Path(args.hfacf)
    vacf = Path(args.vacf)
    if not hfacf.exists():
        raise FileNotFoundError(f"hfacf file not found: {hfacf}")
    if not vacf.exists():
        raise FileNotFoundError(f"vacf file not found: {vacf}")

    if args.gk_data is not None:
        gk_data = Path(args.gk_data)
        if not gk_data.exists():
            raise FileNotFoundError(f"gk_data file not found: {gk_data}")
        volume = load_volume_from_gk_data(gk_data)
    else:
        if args.volume is None:
            raise ValueError("Either --gk_data or --volume must be provided")
        volume = args.volume
        gk_data = None

    workdir = hfacf.resolve().parent
    analyzer = GreenKuboAnalyzer(workdir)

    results = analyzer.analyze(
        corr_stress_file=vacf,
        corr_flux_file=hfacf,
        T=args.T,
        volume=volume,
        visc_t_max_ps=args.visc_t_max_ps,
        kappa_t_max_ps=args.kappa_t_max_ps,
        visc_plateau_start_ps=args.visc_plateau_start_ps,
        kappa_plateau_start_ps=args.kappa_plateau_start_ps,
        smooth_window=args.smooth_window,
        n_blocks=args.n_blocks,
        save_summary=False,
    )

    payload = {
        "inputs": {
            "hfacf": str(hfacf),
            "vacf": str(vacf),
            "gk_data": str(gk_data) if gk_data is not None else None,
            "temperature_K": args.T,
            "volume_A3": volume,
        },
        "transport": {
            "eta_cP": results["eta_cP"],
            "eta_sem_cP": results["eta_sem_cP"],
            "eta_Pa_s": results["eta_Pa_s"],
            "eta_sem_Pa_s": results["eta_sem_Pa_s"],
            "kappa_W_mK": results["kappa_W_mK"],
            "kappa_sem_W_mK": results["kappa_sem_W_mK"],
        },
    }

    out_file = workdir / args.out
    with open(out_file, "w") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

    print("=" * 70)
    print("HFACF / VACF analysis complete")
    print("=" * 70)
    print(f"  hfacf   : {hfacf}")
    print(f"  vacf    : {vacf}")
    if gk_data is not None:
        print(f"  gk_data : {gk_data}")
    print(f"  T       : {args.T:.2f} K")
    print(f"  V       : {volume:.2f} A^3")
    print(f"  eta     : {results['eta_cP']:.4f} +/- {results['eta_sem_cP']:.4f} cP")
    print(f"  kappa   : {results['kappa_W_mK']:.4f} +/- {results['kappa_sem_W_mK']:.4f} W/(m*K)")
    print(f"  json    : {out_file}")
    print("  note    : dielectric constant cannot be computed from hfacf.dat + vacf.dat alone")


if __name__ == "__main__":
    main()
