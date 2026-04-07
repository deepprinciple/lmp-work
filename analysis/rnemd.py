"""
Reverse-NEMD thermal-conductivity analysis for `fix thermal/conductivity`.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - plotting is optional at runtime
    plt = None

from ..utils.constants import (
    ANGSTROM_TO_METER,
    AVOGADRO,
    FEMTOSECOND_TO_SECOND,
    KCAL_MOL_TO_JOULE,
)


def _load_exchange_file(path: Path, timestep_fs: float) -> Dict[str, np.ndarray]:
    rows: list[list[float]] = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            rows.append([float(x) for x in parts[:7]])

    if not rows:
        raise ValueError(f"No numeric rows found in exchange file: {path}")

    arr = np.asarray(rows, dtype=float)
    steps = arr[:, 0]
    step0 = float(steps[0])
    time_fs = (steps - step0) * timestep_fs
    return {
        "step": steps,
        "time_fs": time_fs,
        "time_ps": time_fs / 1000.0,
        "temperature_K": arr[:, 1],
        "volume_A3": arr[:, 2],
        "lx_A": arr[:, 3],
        "ly_A": arr[:, 4],
        "lz_A": arr[:, 5],
        "exchanged_energy_kcal_mol": arr[:, 6],
    }


def _parse_ave_chunk_blocks(path: Path, timestep_fs: float) -> list[dict[str, np.ndarray]]:
    blocks: list[dict[str, np.ndarray]] = []
    current_header: list[float] | None = None
    current_rows: list[list[float]] = []

    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                row = [float(x) for x in line.split()]
            except ValueError:
                continue

            if len(row) in {2, 3}:
                if current_header is not None and current_rows:
                    blocks.append(_build_chunk_block(current_header, current_rows, timestep_fs))
                current_header = row
                current_rows = []
                continue

            current_rows.append(row)

    if current_header is not None and current_rows:
        blocks.append(_build_chunk_block(current_header, current_rows, timestep_fs))

    if not blocks:
        raise ValueError(f"Failed to parse ave/chunk blocks from {path}")
    return blocks


def _build_chunk_block(header: list[float], rows: list[list[float]], timestep_fs: float) -> dict[str, np.ndarray]:
    header_arr = np.asarray(header, dtype=float)
    data = np.asarray(rows, dtype=float)
    timestep = float(header_arr[0])
    n_chunks = int(round(header_arr[1])) if len(header_arr) >= 2 else int(data.shape[0])
    return {
        "step": np.asarray([timestep], dtype=float),
        "time_fs": np.asarray([timestep * timestep_fs], dtype=float),
        "time_ps": np.asarray([timestep * timestep_fs / 1000.0], dtype=float),
        "n_chunks": np.asarray([n_chunks], dtype=float),
        "data": data,
    }


def _fit_line(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    slope, intercept = np.polyfit(x, y, deg=1)
    fit = slope * x + intercept
    ss_res = float(np.sum((y - fit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1.0e-30 and ss_res <= 1.0e-30 else (1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2),
    }


def _extract_profile_columns(data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if data.shape[1] < 4:
        raise ValueError("ave/chunk profile must contain at least 4 numeric columns")
    chunk_id = data[:, 0].astype(int)
    coord = data[:, 1].astype(float)
    count = data[:, 2].astype(float)
    temperature = data[:, -1].astype(float)
    return chunk_id, coord, count, temperature


def _fit_profile_block(
    block: dict[str, np.ndarray],
    *,
    nbin: int,
    swap_buffer_bins: int,
) -> dict[str, Any]:
    data = block["data"]
    chunk_id, coord, count, temperature = _extract_profile_columns(data)
    order = np.argsort(chunk_id)
    chunk_id = chunk_id[order]
    coord = coord[order]
    count = count[order]
    temperature = temperature[order]

    hot_idx = nbin // 2 + 1
    cold_idx = 1

    left_mask = (chunk_id >= cold_idx + swap_buffer_bins) & (chunk_id <= hot_idx - swap_buffer_bins)
    right_mask = (chunk_id >= hot_idx + swap_buffer_bins) & (chunk_id <= nbin - swap_buffer_bins + 1)

    if np.count_nonzero(left_mask) < 3 or np.count_nonzero(right_mask) < 3:
        raise ValueError("Not enough profile bins after excluding swap-adjacent layers")

    left_fit = _fit_line(coord[left_mask], temperature[left_mask])
    right_fit = _fit_line(coord[right_mask], temperature[right_mask])
    abs_grad_k_per_m = 0.5 * (
        abs(left_fit["slope"]) / ANGSTROM_TO_METER
        + abs(right_fit["slope"]) / ANGSTROM_TO_METER
    )

    return {
        "chunk_id": chunk_id,
        "coord_A": coord,
        "count": count,
        "temperature_K": temperature,
        "left_fit": left_fit,
        "right_fit": right_fit,
        "left_mask": left_mask,
        "right_mask": right_mask,
        "abs_gradient_K_m": float(abs_grad_k_per_m),
        "min_count": float(np.min(count)),
    }


def _select_stable_window(
    time_ps: np.ndarray,
    values: np.ndarray,
    *,
    start_ps: float,
    rel_std_tol: float,
    drift_tol: float,
    min_points: int,
) -> dict[str, Any]:
    start_idx = int(np.searchsorted(time_ps, start_ps, side="left"))
    best_mask = None
    best_score = None
    n = len(time_ps)

    for i in range(start_idx, max(start_idx + 1, n - min_points + 1)):
        for j in range(i + min_points, n + 1):
            mask = np.zeros(n, dtype=bool)
            mask[i:j] = True
            window = values[mask]
            if len(window) < min_points:
                continue
            mean_val = float(np.mean(window))
            if abs(mean_val) <= 1.0e-30:
                continue
            rel_std = float(np.std(window) / abs(mean_val))
            drift = float(abs(window[-1] - window[0]) / abs(mean_val))
            stable = rel_std <= rel_std_tol and drift <= drift_tol
            score = (stable, len(window), -rel_std, -drift)
            if best_score is None or score > best_score:
                best_score = score
                best_mask = mask

    if best_mask is None:
        best_mask = np.zeros(n, dtype=bool)
        tail_start = max(start_idx, n - max(min_points, n // 5))
        best_mask[tail_start:] = True
        reason = "fallback_tail_window"
    else:
        reason = "stable_window" if best_score[0] else "best_available_window"

    vals = values[best_mask]
    return {
        "mask": best_mask,
        "reason": reason,
        "stable": bool(best_score[0]) if best_score is not None else False,
        "mean": float(np.mean(vals)),
        "std": float(np.std(vals)),
        "start_ps": float(time_ps[best_mask][0]),
        "end_ps": float(time_ps[best_mask][-1]),
    }


def analyze_rnemd_replica(
    repdir: Path,
    *,
    timestep_fs: float,
    nbin: int,
    edim: str,
    plateau_start_ps: float = 20.0,
    swap_buffer_bins: int = 2,
    kappa_rel_std_tol: float = 0.15,
    kappa_drift_tol: float = 0.15,
    min_window_points: int = 5,
    profile_file: str = "temp_profile.dat",
    exchange_file: str = "thermal_exchange.dat",
    summary_name: str = "thermal_rnemd_summary.json",
    figure_profile_name: str = "thermal_profile_latest.png",
    figure_kappa_name: str = "thermal_kappa_running.png",
) -> Dict[str, Any]:
    repdir = Path(repdir)
    exchange = _load_exchange_file(repdir / exchange_file, timestep_fs=timestep_fs)
    profile_blocks = _parse_ave_chunk_blocks(repdir / profile_file, timestep_fs=timestep_fs)

    exch_steps = exchange["step"]
    time_ps = []
    kappa_values = []
    left_r2 = []
    right_r2 = []
    gradients = []
    profile_summaries = []

    area_map = {
        "x": exchange["ly_A"] * exchange["lz_A"],
        "y": exchange["lx_A"] * exchange["lz_A"],
        "z": exchange["lx_A"] * exchange["ly_A"],
    }
    area_m2_series = area_map[edim] * (ANGSTROM_TO_METER ** 2)

    for block in profile_blocks:
        step = float(block["step"][0])
        idx = int(np.searchsorted(exch_steps, step, side="right") - 1)
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
        grad = float(profile["abs_gradient_K_m"])
        if grad <= 0.0:
            continue
        kappa = heat_flux_w_m2 / grad

        time_ps.append(float(block["time_ps"][0]))
        kappa_values.append(float(kappa))
        left_r2.append(float(profile["left_fit"]["r2"]))
        right_r2.append(float(profile["right_fit"]["r2"]))
        gradients.append(float(grad))
        profile_summaries.append(
            {
                "time_ps": float(block["time_ps"][0]),
                "step": step,
                "left_r2": float(profile["left_fit"]["r2"]),
                "right_r2": float(profile["right_fit"]["r2"]),
                "abs_gradient_K_m": float(grad),
                "min_count": float(profile["min_count"]),
            }
        )

    if not time_ps:
        raise ValueError(f"No overlapping profile/exchange data could be analyzed in {repdir}")

    time_ps_arr = np.asarray(time_ps, dtype=float)
    kappa_arr = np.asarray(kappa_values, dtype=float)
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
    if plt is not None:
        profile_png = repdir / figure_profile_name
        plt.figure(figsize=(7.0, 4.8))
        coord = latest_profile["coord_A"]
        temp = latest_profile["temperature_K"]
        plt.plot(coord, temp, "o-", lw=1.2, label="temperature profile")
        for mask, fit, label in (
            (latest_profile["left_mask"], latest_profile["left_fit"], "left fit"),
            (latest_profile["right_mask"], latest_profile["right_fit"], "right fit"),
        ):
            x = coord[mask]
            y = fit["slope"] * x + fit["intercept"]
            plt.plot(x, y, "--", lw=1.6, label=f"{label} (R2={fit['r2']:.3f})")
        plt.xlabel(f"{edim} coordinate (A)")
        plt.ylabel("Temperature (K)")
        plt.title("Latest reverse-NEMD temperature profile")
        plt.legend()
        plt.tight_layout()
        plt.savefig(profile_png, dpi=180)
        plt.close()

        kappa_png = repdir / figure_kappa_name
        plt.figure(figsize=(7.0, 4.8))
        plt.plot(time_ps_arr, kappa_arr, lw=1.4, label="running kappa")
        mask = window["mask"]
        plt.axvspan(window["start_ps"], window["end_ps"], color="tab:green", alpha=0.15, label=window["reason"])
        plt.plot(time_ps_arr[mask], np.full(np.count_nonzero(mask), window["mean"]), "--", lw=1.6, label="window mean")
        plt.xlabel("Time (ps)")
        plt.ylabel("kappa (W/mK)")
        plt.title("reverse-NEMD thermal conductivity convergence")
        plt.legend()
        plt.tight_layout()
        plt.savefig(kappa_png, dpi=180)
        plt.close()

    result = {
        "repdir": str(repdir),
        "plateau_start_ps": plateau_start_ps,
        "swap_buffer_bins": swap_buffer_bins,
        "edim": edim,
        "nbin": nbin,
        "kappa_W_mK": float(window["mean"]),
        "kappa_window_std_W_mK": float(window["std"]),
        "kappa_window_start_ps": float(window["start_ps"]),
        "kappa_window_end_ps": float(window["end_ps"]),
        "kappa_window_stable": bool(window["stable"]),
        "kappa_window_reason": str(window["reason"]),
        "mean_left_r2": float(np.mean(left_r2)),
        "mean_right_r2": float(np.mean(right_r2)),
        "mean_abs_gradient_K_m": float(np.mean(gradients)),
        "running_time_ps": time_ps_arr.tolist(),
        "running_kappa_W_mK": kappa_arr.tolist(),
        "profile_samples": profile_summaries,
    }

    out_file = repdir / summary_name
    out_file.write_text(json.dumps(result, indent=2))
    return result
