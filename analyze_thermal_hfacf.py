"""Postprocess HFACF and visualize thermal conductivity convergence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gk_workflow.analysis.greenkubo import (
    compute_thermal_conductivity,
    parse_ave_correlate_detail,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze HFACF and plot thermal conductivity running integral."
    )
    parser.add_argument("--hfacf", required=True, help="LAMMPS HFACF file.")
    parser.add_argument(
        "--gk-data",
        default=None,
        help="Optional gk_data.dat file for automatic average volume.",
    )
    parser.add_argument(
        "--volume",
        type=float,
        default=None,
        help="Average volume in A^3 when gk_data.dat is absent.",
    )
    parser.add_argument("--temperature-k", type=float, required=True)
    parser.add_argument("--t-max-ps", type=float, default=20.0)
    parser.add_argument("--plateau-start-ps", type=float, default=2.0)
    parser.add_argument("--smooth-window", type=int, default=11)
    parser.add_argument("--json-out", default="thermal_hfacf_summary.json")
    parser.add_argument("--plot-out", default="thermal_hfacf_analysis.png")
    return parser.parse_args()


def load_volume_from_gk_data(gk_data_file: Path) -> float:
    rows = []
    with open(gk_data_file) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            rows.append(parts)

    if not rows:
        raise ValueError(f"No numeric rows found in gk_data file: {gk_data_file}")
    return sum(float(row[2]) for row in rows) / len(rows)


def main() -> None:
    args = parse_args()
    hfacf = Path(args.hfacf).resolve()
    if not hfacf.exists():
        raise FileNotFoundError(f"HFACF file not found: {hfacf}")

    if args.gk_data is not None:
        volume = load_volume_from_gk_data(Path(args.gk_data).resolve())
    else:
        if args.volume is None:
            raise ValueError("Either --gk-data or --volume must be provided.")
        volume = float(args.volume)

    detail = parse_ave_correlate_detail(hfacf)
    result = compute_thermal_conductivity(
        corr_flux_file=hfacf,
        T=float(args.temperature_k),
        volume=volume,
        t_max_ps=float(args.t_max_ps),
        plateau_start_ps=float(args.plateau_start_ps),
        smooth_window=int(args.smooth_window),
    )

    time_ps = np.asarray(result["time_ps"], dtype=float)
    running_kappa = np.asarray(result["running_integral_W_mK"], dtype=float)
    hfacf_time_ps = np.asarray(detail["time_fs"], dtype=float) / 1000.0
    hfacf_mean = np.mean(np.asarray(detail["corr"], dtype=float), axis=1)
    norm = np.max(np.abs(hfacf_mean))
    hfacf_norm = hfacf_mean / norm if norm > 0 else hfacf_mean

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0))
    axes[0].plot(hfacf_time_ps, hfacf_norm, lw=1.2)
    axes[0].axhline(0.0, color="gray", lw=0.8)
    axes[0].set_xlabel("Time (ps)")
    axes[0].set_ylabel("Normalized HFACF")
    axes[0].set_title("Heat-flux autocorrelation function")

    axes[1].plot(time_ps, running_kappa, lw=1.4, label="running kappa")
    axes[1].axvspan(
        float(result["plateau_start_ps"]),
        float(result["plateau_end_ps"]),
        color="tab:green",
        alpha=0.15,
        label="analysis window",
    )
    axes[1].axhline(
        float(result["kappa_W_mK"]),
        color="tab:red",
        ls="--",
        lw=1.2,
        label=f"kappa = {result['kappa_W_mK']:.4f} W/mK",
    )
    axes[1].set_xlabel("Time (ps)")
    axes[1].set_ylabel("Thermal conductivity (W/mK)")
    axes[1].set_title("Green-Kubo thermal conductivity convergence")
    axes[1].legend()

    fig.tight_layout()
    plot_out = hfacf.parent / args.plot_out
    fig.savefig(plot_out, dpi=180)
    plt.close(fig)

    summary = {
        "hfacf_file": str(hfacf),
        "temperature_K": float(args.temperature_k),
        "volume_A3": volume,
        "kappa_W_mK": float(result["kappa_W_mK"]),
        "kappa_components_W_mK": result["kappa_components_W_mK"],
        "plateau_start_ps": float(result["plateau_start_ps"]),
        "plateau_end_ps": float(result["plateau_end_ps"]),
        "effective_t_max_ps": float(result["effective_t_max_ps"]),
        "plot_file": str(plot_out),
    }
    json_out = hfacf.parent / args.json_out
    json_out.write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
