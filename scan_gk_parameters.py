"""Scan Green-Kubo analysis parameters on existing ACF outputs.

Example:
    conda run -n HTMD python scan_gk_parameters.py \
      --repdir /root/lmp-work/bench_runs/methanol_24core_accuracy/replica_01 \
      --temperature 298.15 \
      --gk-data gk_data_conv.dat \
      --vacf vacf_conv.dat \
      --hfacf hfacf_conv.dat \
      --kappa-t-max-ps 12 15 20 \
      --kappa-plateau-start-ps 2 4 6 \
      --smooth-window 11 21 \
      --kappa-rel-std-tol 0.20 0.25
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent))

from gk_workflow.analysis.greenkubo import (  # type: ignore
    KAPPA_CONV,
    VISC_CONV,
    _analyze_transport_running_integral,
)


def load_volume_from_gk_data(gk_data_file: Path) -> float:
    rows: list[list[str]] = []
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
        raise ValueError(f"No numeric rows found in {gk_data_file}")
    volumes = [float(row[2]) for row in rows]
    return sum(volumes) / len(volumes)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan GK analysis parameters on existing outputs")
    parser.add_argument("--repdir", required=True, help="Replica directory containing ACF files")
    parser.add_argument("--temperature", type=float, required=True, help="Simulation temperature in K")
    parser.add_argument("--gk-data", default="gk_data_conv.dat", help="Path to GK data file")
    parser.add_argument("--vacf", default="vacf_conv.dat", help="Path to VACF/stress correlation file")
    parser.add_argument("--hfacf", default="hfacf_conv.dat", help="Path to HFACF/heat-flux correlation file")
    parser.add_argument("--visc-t-max-ps", type=float, nargs="+", default=[20.0])
    parser.add_argument("--kappa-t-max-ps", type=float, nargs="+", default=[12.0, 15.0, 20.0])
    parser.add_argument("--visc-plateau-start-ps", type=float, nargs="+", default=[5.0])
    parser.add_argument("--kappa-plateau-start-ps", type=float, nargs="+", default=[2.0, 4.0, 6.0])
    parser.add_argument("--smooth-window", type=int, nargs="+", default=[11, 21])
    parser.add_argument("--visc-rel-std-tol", type=float, nargs="+", default=[0.12])
    parser.add_argument("--kappa-rel-std-tol", type=float, nargs="+", default=[0.20, 0.25])
    parser.add_argument("--visc-drift-tol", type=float, nargs="+", default=[0.12])
    parser.add_argument("--kappa-drift-tol", type=float, nargs="+", default=[0.20])
    parser.add_argument("--visc-plateau-end-ps", type=float, default=None)
    parser.add_argument("--kappa-plateau-end-ps", type=float, default=None)
    parser.add_argument("--n-blocks", type=int, default=4)
    parser.add_argument("--csv-out", default="gk_scan_results.csv", help="Output CSV path")
    parser.add_argument("--json-out", default="gk_scan_results.json", help="Output JSON path")
    return parser.parse_args()


def run_single_scan(
    *,
    vacf: Path,
    hfacf: Path,
    temperature: float,
    volume: float,
    visc_t_max_ps: float,
    kappa_t_max_ps: float,
    visc_plateau_start_ps: float,
    kappa_plateau_start_ps: float,
    smooth_window: int,
    visc_rel_std_tol: float,
    kappa_rel_std_tol: float,
    visc_drift_tol: float,
    kappa_drift_tol: float,
    visc_plateau_end_ps: float | None,
    kappa_plateau_end_ps: float | None,
) -> dict[str, Any]:
    visc = _analyze_transport_running_integral(
        corr_file=vacf,
        T=temperature,
        volume=volume,
        t_max_ps=visc_t_max_ps,
        plateau_start_ps=visc_plateau_start_ps,
        plateau_end_ps=visc_plateau_end_ps,
        smooth_window=smooth_window,
        rel_std_tol=visc_rel_std_tol,
        drift_tol=visc_drift_tol,
        min_ps_before_cutoff=max(0.2, visc_plateau_start_ps * 0.5),
        converter_scale=VISC_CONV,
        component_factor=volume / temperature,
        component_unit="Pa·s",
        result_key="eta_Pa_s",
    )
    kappa = _analyze_transport_running_integral(
        corr_file=hfacf,
        T=temperature,
        volume=volume,
        t_max_ps=kappa_t_max_ps,
        plateau_start_ps=kappa_plateau_start_ps,
        plateau_end_ps=kappa_plateau_end_ps,
        smooth_window=smooth_window,
        rel_std_tol=kappa_rel_std_tol,
        drift_tol=kappa_drift_tol,
        min_ps_before_cutoff=max(0.3, kappa_plateau_start_ps * 0.5),
        converter_scale=KAPPA_CONV,
        component_factor=volume / (temperature**2),
        component_unit="W/(m·K)",
        result_key="kappa_W_mK",
    )
    visc_diag = visc["diagnostics"]
    kappa_diag = kappa["diagnostics"]
    return {
        "visc_t_max_ps": visc_t_max_ps,
        "kappa_t_max_ps": kappa_t_max_ps,
        "visc_plateau_start_ps": visc_plateau_start_ps,
        "kappa_plateau_start_ps": kappa_plateau_start_ps,
        "smooth_window": smooth_window,
        "visc_rel_std_tol": visc_rel_std_tol,
        "kappa_rel_std_tol": kappa_rel_std_tol,
        "visc_drift_tol": visc_drift_tol,
        "kappa_drift_tol": kappa_drift_tol,
        "eta_cP": float(visc["eta_Pa_s"]) * 1e3,
        "eta_stable": bool(visc_diag["window_stable"]),
        "eta_reason": str(visc_diag["window_reason"]),
        "eta_window_start_ps": float(visc_diag["window_start_ps"]),
        "eta_window_end_ps": float(visc_diag["window_end_ps"]),
        "eta_window_rel_std": float(visc_diag["window_rel_std"]),
        "eta_window_drift_fraction": float(visc_diag["window_drift_fraction"]),
        "eta_window_slope_fraction": float(visc_diag["window_slope_fraction"]),
        "eta_component_rel_std": float(visc_diag["component_rel_std"]),
        "kappa_W_mK": float(kappa["kappa_W_mK"]),
        "kappa_stable": bool(kappa_diag["window_stable"]),
        "kappa_reason": str(kappa_diag["window_reason"]),
        "kappa_window_start_ps": float(kappa_diag["window_start_ps"]),
        "kappa_window_end_ps": float(kappa_diag["window_end_ps"]),
        "kappa_window_rel_std": float(kappa_diag["window_rel_std"]),
        "kappa_window_drift_fraction": float(kappa_diag["window_drift_fraction"]),
        "kappa_window_slope_fraction": float(kappa_diag["window_slope_fraction"]),
        "kappa_component_rel_std": float(kappa_diag["component_rel_std"]),
        "passes_both": bool(visc_diag["window_stable"]) and bool(kappa_diag["window_stable"]),
    }


def main() -> None:
    args = parse_args()
    repdir = Path(args.repdir).resolve()
    vacf = (repdir / args.vacf).resolve() if not Path(args.vacf).is_absolute() else Path(args.vacf).resolve()
    hfacf = (repdir / args.hfacf).resolve() if not Path(args.hfacf).is_absolute() else Path(args.hfacf).resolve()
    gk_data = (repdir / args.gk_data).resolve() if not Path(args.gk_data).is_absolute() else Path(args.gk_data).resolve()
    volume = load_volume_from_gk_data(gk_data)

    rows: list[dict[str, Any]] = []
    for values in itertools.product(
        args.visc_t_max_ps,
        args.kappa_t_max_ps,
        args.visc_plateau_start_ps,
        args.kappa_plateau_start_ps,
        args.smooth_window,
        args.visc_rel_std_tol,
        args.kappa_rel_std_tol,
        args.visc_drift_tol,
        args.kappa_drift_tol,
    ):
        row = run_single_scan(
            vacf=vacf,
            hfacf=hfacf,
            temperature=args.temperature,
            volume=volume,
            visc_t_max_ps=float(values[0]),
            kappa_t_max_ps=float(values[1]),
            visc_plateau_start_ps=float(values[2]),
            kappa_plateau_start_ps=float(values[3]),
            smooth_window=int(values[4]),
            visc_rel_std_tol=float(values[5]),
            kappa_rel_std_tol=float(values[6]),
            visc_drift_tol=float(values[7]),
            kappa_drift_tol=float(values[8]),
            visc_plateau_end_ps=args.visc_plateau_end_ps,
            kappa_plateau_end_ps=args.kappa_plateau_end_ps,
        )
        rows.append(row)

    rows.sort(
        key=lambda row: (
            not row["kappa_stable"],
            row["kappa_window_rel_std"],
            row["kappa_window_drift_fraction"],
            row["kappa_component_rel_std"],
        )
    )

    csv_out = (repdir / args.csv_out).resolve() if not Path(args.csv_out).is_absolute() else Path(args.csv_out).resolve()
    json_out = (repdir / args.json_out).resolve() if not Path(args.json_out).is_absolute() else Path(args.json_out).resolve()
    fieldnames = list(rows[0].keys()) if rows else []
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "repdir": str(repdir),
        "temperature_K": args.temperature,
        "volume_A3": volume,
        "rows": rows,
    }
    with open(json_out, "w") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)

    print(f"Scan rows: {len(rows)}")
    print(f"CSV: {csv_out}")
    print(f"JSON: {json_out}")
    if rows:
        best = rows[0]
        print(
            "Best kappa candidate:",
            f"kappa={best['kappa_W_mK']:.4f}",
            f"stable={best['kappa_stable']}",
            f"rel_std={best['kappa_window_rel_std']:.4f}",
            f"drift={best['kappa_window_drift_fraction']:.4f}",
            f"component_rel_std={best['kappa_component_rel_std']:.4f}",
        )


if __name__ == "__main__":
    main()
