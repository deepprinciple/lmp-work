from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gk_workflow.analysis.rnemd import (  # noqa: E402
    _fit_profile_block,
    _load_exchange_file,
    _parse_ave_chunk_blocks,
    _select_stable_window,
)
from gk_workflow.utils.constants import (  # noqa: E402
    ANGSTROM_TO_METER,
    AVOGADRO,
    FEMTOSECOND_TO_SECOND,
    KCAL_MOL_TO_JOULE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write a live reverse-NEMD convergence JSON file.")
    parser.add_argument("--repdir", required=True, help="Replica directory containing temp_profile.dat and thermal_exchange.dat")
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--timestep-fs", type=float, default=1.0)
    parser.add_argument("--nbin", type=int, default=20)
    parser.add_argument("--edim", choices=("x", "y", "z"), default="z")
    parser.add_argument("--plateau-start-ps", type=float, default=20.0)
    parser.add_argument("--swap-buffer-bins", type=int, default=2)
    parser.add_argument("--kappa-rel-std-tol", type=float, default=0.15)
    parser.add_argument("--kappa-drift-tol", type=float, default=0.15)
    parser.add_argument("--min-window-points", type=int, default=5)
    parser.add_argument("--exchange-file", default="thermal_exchange.dat")
    parser.add_argument("--profile-file", default="temp_profile.dat")
    parser.add_argument("--interval-s", type=float, default=30.0)
    parser.add_argument("--max-cycles", type=int, default=0, help="0 means run until stopped")
    return parser.parse_args()


def _compute_live_summary(
    repdir: Path,
    *,
    timestep_fs: float,
    nbin: int,
    edim: str,
    plateau_start_ps: float,
    swap_buffer_bins: int,
    kappa_rel_std_tol: float,
    kappa_drift_tol: float,
    min_window_points: int,
    exchange_file: str,
    profile_file: str,
) -> dict:
    exchange = _load_exchange_file(repdir / exchange_file, timestep_fs=timestep_fs)
    profile_blocks = _parse_ave_chunk_blocks(repdir / profile_file, timestep_fs=timestep_fs)

    area_map = {
        "x": exchange["ly_A"] * exchange["lz_A"],
        "y": exchange["lx_A"] * exchange["lz_A"],
        "z": exchange["lx_A"] * exchange["ly_A"],
    }
    area_m2_series = area_map[edim] * (ANGSTROM_TO_METER ** 2)

    time_ps = []
    kappa_values = []
    left_r2 = []
    right_r2 = []
    gradients = []

    exch_steps = exchange["step"]
    for block in profile_blocks:
        step = float(block["step"][0])
        idx = int((exch_steps <= step).sum() - 1)
        if idx < 0:
            continue

        profile = _fit_profile_block(
            block,
            nbin=nbin,
            swap_buffer_bins=swap_buffer_bins,
        )
        elapsed_s = float(exchange["time_fs"][idx]) * FEMTOSECOND_TO_SECOND
        if elapsed_s <= 0.0:
            continue

        energy_j = float(exchange["exchanged_energy_kcal_mol"][idx]) * KCAL_MOL_TO_JOULE / AVOGADRO
        area_m2 = float(area_m2_series[idx])
        heat_flux_w_m2 = abs(energy_j) / (2.0 * area_m2 * elapsed_s)
        gradient_k_m = float(profile["abs_gradient_K_m"])
        if gradient_k_m <= 0.0:
            continue

        time_ps.append(float(block["time_ps"][0]))
        kappa_values.append(float(heat_flux_w_m2 / gradient_k_m))
        left_r2.append(float(profile["left_fit"]["r2"]))
        right_r2.append(float(profile["right_fit"]["r2"]))
        gradients.append(float(gradient_k_m))

    if not time_ps:
        raise ValueError("No overlapping profile/exchange samples available yet")

    time_ps_arr = __import__("numpy").array(time_ps, dtype=float)
    kappa_arr = __import__("numpy").array(kappa_values, dtype=float)
    window = _select_stable_window(
        time_ps_arr,
        kappa_arr,
        start_ps=plateau_start_ps,
        rel_std_tol=kappa_rel_std_tol,
        drift_tol=kappa_drift_tol,
        min_points=min_window_points,
    )

    latest_profile = _fit_profile_block(
        profile_blocks[-1],
        nbin=nbin,
        swap_buffer_bins=swap_buffer_bins,
    )
    latest_exchange_idx = len(exchange["step"]) - 1
    latest_heat_flux = (
        abs(float(exchange["exchanged_energy_kcal_mol"][latest_exchange_idx])) * KCAL_MOL_TO_JOULE / AVOGADRO
    ) / (
        2.0
        * float(area_m2_series[latest_exchange_idx])
        * float(exchange["time_fs"][latest_exchange_idx])
        * FEMTOSECOND_TO_SECOND
    )
    latest_gradient = float(latest_profile["abs_gradient_K_m"])
    latest_temperature = latest_profile["temperature_K"]

    return {
        "repdir": str(repdir),
        "updated_at_epoch_s": time.time(),
        "status": "running",
        "current_step": int(exchange["step"][latest_exchange_idx]),
        "current_time_ps": float(exchange["time_ps"][latest_exchange_idx]),
        "current_temperature_K": float(exchange["temperature_K"][latest_exchange_idx]),
        "current_kappa_W_mK": float(latest_heat_flux / latest_gradient) if latest_gradient > 0 else None,
        "window_kappa_W_mK": float(window["mean"]),
        "window_kappa_std_W_mK": float(window["std"]),
        "window_start_ps": float(window["start_ps"]),
        "window_end_ps": float(window["end_ps"]),
        "window_reason": str(window["reason"]),
        "window_stable": bool(window["stable"]),
        "latest_profile_deltaT_K": float(latest_temperature.max() - latest_temperature.min()),
        "latest_gradient_K_m": latest_gradient,
        "latest_left_r2": float(latest_profile["left_fit"]["r2"]),
        "latest_right_r2": float(latest_profile["right_fit"]["r2"]),
        "latest_min_count": float(latest_profile["min_count"]),
        "mean_left_r2": float(sum(left_r2) / len(left_r2)),
        "mean_right_r2": float(sum(right_r2) / len(right_r2)),
        "mean_gradient_K_m": float(sum(gradients) / len(gradients)),
        "running_time_ps": [float(x) for x in time_ps],
        "running_kappa_W_mK": [float(x) for x in kappa_values],
    }


def main() -> None:
    args = parse_args()
    repdir = Path(args.repdir).resolve()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    cycle = 0
    while True:
        cycle += 1
        try:
            result = _compute_live_summary(
                repdir,
                timestep_fs=args.timestep_fs,
                nbin=args.nbin,
                edim=args.edim,
                plateau_start_ps=args.plateau_start_ps,
                swap_buffer_bins=args.swap_buffer_bins,
                kappa_rel_std_tol=args.kappa_rel_std_tol,
                kappa_drift_tol=args.kappa_drift_tol,
                min_window_points=args.min_window_points,
                exchange_file=args.exchange_file,
                profile_file=args.profile_file,
            )
        except Exception as exc:
            result = {
                "repdir": str(repdir),
                "updated_at_epoch_s": time.time(),
                "status": "waiting_for_data",
                "error": str(exc),
            }

        tmp_path = output.with_suffix(output.suffix + ".tmp")
        tmp_path.write_text(json.dumps(result, indent=2))
        tmp_path.replace(output)

        if args.max_cycles > 0 and cycle >= args.max_cycles:
            break
        time.sleep(args.interval_s)


if __name__ == "__main__":
    main()
