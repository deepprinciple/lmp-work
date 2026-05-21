"""Generic correlation parsing and running-integral helpers."""
from __future__ import annotations

import warnings
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np


def _running_average(data: np.ndarray, window: int) -> np.ndarray:
    return np.convolve(data, np.ones(window) / window, mode="valid")


def _integrate_acf(
    time: np.ndarray,
    acf: np.ndarray,
    t_max: Optional[float] = None,
) -> tuple[np.ndarray, float]:
    if t_max is not None:
        mask = time <= t_max
        time = time[mask]
        acf = acf[mask]

    running = np.zeros(len(time))
    for i in range(1, len(time)):
        dt = time[i] - time[i - 1]
        running[i] = running[i - 1] + 0.5 * (acf[i - 1] + acf[i]) * dt

    return running, float(running[-1])


def _relative_metric(value: float, reference: float, fallback_scale: float) -> float:
    denom = max(abs(reference), fallback_scale, 1.0e-12)
    return abs(value) / denom


def parse_ave_correlate_detail(filepath: Path) -> Dict[str, np.ndarray | int | float | bool]:
    """Parse LAMMPS ``fix ave/correlate`` output and keep the last full block."""
    numeric_rows = []
    row_lengths = []
    with open(filepath) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                values = [float(x) for x in line.split()]
            except ValueError:
                continue
            numeric_rows.append(values)
            row_lengths.append(len(values))

    if not numeric_rows:
        raise ValueError(f"No numeric data found in {filepath}")

    count_by_len = Counter(row_lengths)
    data_len = max(count_by_len.items(), key=lambda kv: (kv[1], kv[0]))[0]
    if data_len < 5:
        raise ValueError(f"Unexpected ave/correlate column count: {data_len}")

    blocks = []
    current_block = []
    dropped_rows = 0

    for row in numeric_rows:
        if len(row) == data_len:
            current_block.append(row)
            continue

        if len(row) == 2:
            if current_block:
                blocks.append(np.asarray(current_block, dtype=float))
                current_block = []
            continue

        dropped_rows += 1

    if current_block:
        blocks.append(np.asarray(current_block, dtype=float))

    if not blocks:
        filtered = [row for row in numeric_rows if len(row) == data_len]
        if not filtered:
            raise ValueError(f"No valid correlation block found in {filepath}")
        blocks = [np.asarray(filtered, dtype=float)]

    if dropped_rows > 0:
        warnings.warn(
            f"{filepath} contains {dropped_rows} unrecognized rows; ignored",
            RuntimeWarning,
        )

    arr = blocks[-1]

    if arr.shape[1] >= 6:
        rounded = np.rint(arr[:, 0])
        diffs = np.diff(rounded)
        looks_index = (
            np.all(np.isclose(arr[:, 0], rounded))
            and len(diffs) > 0
            and np.all(diffs == 1)
            and rounded[0] in (0.0, 1.0)
        )
        if looks_index:
            time_delta = arr[:, 1]
            corr = arr[:, 3:]
        else:
            time_delta = arr[:, 0]
            corr = arr[:, 2:]
    else:
        looks_index = False
        time_delta = arr[:, 0]
        corr = arr[:, 2:]

    if corr.shape[1] < 3:
        raise ValueError(f"Need at least 3 correlation columns in {filepath}")

    ncount = arr[:, 2] if arr.shape[1] >= 6 and looks_index else arr[:, 1]

    valid_mask = np.isfinite(time_delta) & np.all(np.isfinite(corr), axis=1)
    positive_count_mask = valid_mask & (ncount > 0)
    data_mask = positive_count_mask if np.count_nonzero(positive_count_mask) >= 3 else valid_mask

    return {
        "time_fs": np.asarray(time_delta[data_mask], dtype=float),
        "corr": np.asarray(corr[data_mask], dtype=float),
        "ncount": np.asarray(ncount[data_mask], dtype=float),
        "raw_time_fs": np.asarray(time_delta, dtype=float),
        "raw_ncount": np.asarray(ncount, dtype=float),
        "selected_block_rows": int(arr.shape[0]),
        "selected_data_rows": int(np.count_nonzero(data_mask)),
        "n_blocks_detected": len(blocks),
        "used_positive_ncount_rows": bool(np.count_nonzero(positive_count_mask) >= 3),
    }


def _estimate_effective_cutoff_ps(
    time_ps: np.ndarray,
    corr_matrix: np.ndarray,
    requested_t_max_ps: Optional[float],
    *,
    min_ps_before_cutoff: float,
    zero_run_length: int = 3,
) -> float:
    if len(time_ps) == 0:
        return 0.0

    effective_end = time_ps[-1]
    if requested_t_max_ps is not None:
        effective_end = min(effective_end, requested_t_max_ps)

    mean_corr = np.mean(corr_matrix, axis=1)
    candidate_mask = time_ps >= min_ps_before_cutoff
    candidate_indices = np.where(candidate_mask)[0]
    if len(candidate_indices) >= zero_run_length:
        start_idx = int(candidate_indices[0])
        for idx in range(start_idx, len(mean_corr) - zero_run_length + 1):
            if np.all(mean_corr[idx : idx + zero_run_length] <= 0.0):
                effective_end = min(effective_end, float(time_ps[idx]))
                break

    return max(float(time_ps[0]), float(effective_end))


def _fallback_window_mask(n_points: int, min_points: int) -> np.ndarray:
    mask = np.zeros(n_points, dtype=bool)
    if n_points == 0:
        return mask
    start_idx = max(
        0,
        min(n_points - 1, n_points - max(min_points, int(n_points * 0.2))),
    )
    mask[start_idx:] = True
    if np.count_nonzero(mask) < min_points:
        mask[max(0, n_points - min_points) :] = True
    return mask


def _select_analysis_window(
    time_ps: np.ndarray,
    running_integral: np.ndarray,
    *,
    plateau_start_ps: float,
    plateau_end_ps: Optional[float],
    rel_std_tol: float,
    drift_tol: float,
    min_points: int = 8,
) -> Dict[str, object]:
    n_points = min(len(time_ps), len(running_integral))
    if n_points == 0:
        raise ValueError("Cannot select an analysis window from an empty series")
    if len(time_ps) != len(running_integral):
        time_ps = time_ps[:n_points]
        running_integral = running_integral[:n_points]

    analysis_end_ps = (
        float(time_ps[-1])
        if plateau_end_ps is None
        else min(float(plateau_end_ps), float(time_ps[-1]))
    )
    start_idx = int(np.searchsorted(time_ps, plateau_start_ps, side="left"))
    end_idx = int(np.searchsorted(time_ps, analysis_end_ps, side="right")) - 1
    if end_idx < start_idx:
        end_idx = n_points - 1

    min_points = min(min_points, n_points)
    scale = max(float(np.max(np.abs(running_integral))), 1.0e-12)
    best_start_idx: int | None = None
    best_stats: dict[str, float] | None = None

    for candidate_start in range(start_idx, max(start_idx + 1, end_idx - min_points + 2)):
        window = running_integral[candidate_start : end_idx + 1]
        if len(window) < min_points:
            continue

        time_window = time_ps[candidate_start : end_idx + 1]
        mean_val = float(window.mean())
        std_val = float(window.std())
        drift = float(window[-1] - window[0])
        rel_std = _relative_metric(std_val, mean_val, scale * 0.05)
        drift_fraction = _relative_metric(drift, mean_val, scale * 0.05)
        duration_ps = float(max(time_window[-1] - time_window[0], 0.0))

        if len(window) >= 2 and duration_ps > 0.0:
            slope = float(np.polyfit(time_window, window, 1)[0])
            slope_fraction = _relative_metric(slope * duration_ps, mean_val, scale * 0.05)
        else:
            slope = 0.0
            slope_fraction = 0.0

        stable = (
            rel_std <= rel_std_tol
            and drift_fraction <= drift_tol
            and slope_fraction <= drift_tol
        )
        if not stable:
            continue

        best_start_idx = candidate_start
        best_stats = {
            "mean": mean_val,
            "std": std_val,
            "rel_std": rel_std,
            "drift_fraction": drift_fraction,
            "slope_fraction": slope_fraction,
            "duration_ps": duration_ps,
            "slope": slope,
        }
        break

    if best_start_idx is None or best_stats is None:
        mask = _fallback_window_mask(n_points, min_points=min_points)
        window = running_integral[mask]
        time_window = time_ps[mask]
        mean_val = float(window.mean())
        std_val = float(window.std())
        drift = float(window[-1] - window[0]) if len(window) >= 2 else 0.0
        duration_ps = float(max(time_window[-1] - time_window[0], 0.0))
        rel_std = _relative_metric(std_val, mean_val, scale * 0.05)
        drift_fraction = _relative_metric(drift, mean_val, scale * 0.05)
        slope = (
            float(np.polyfit(time_window, window, 1)[0])
            if len(window) >= 2 and duration_ps > 0.0
            else 0.0
        )
        slope_fraction = (
            _relative_metric(slope * duration_ps, mean_val, scale * 0.05)
            if duration_ps > 0.0
            else 0.0
        )
        stable = (
            rel_std <= rel_std_tol
            and drift_fraction <= drift_tol
            and slope_fraction <= drift_tol
        )
        return {
            "mask": mask,
            "stable": stable,
            "reason": "fallback_tail_window_stable" if stable else "fallback_tail_window",
            "stats": {
                "mean": mean_val,
                "std": std_val,
                "rel_std": rel_std,
                "drift_fraction": drift_fraction,
                "slope_fraction": slope_fraction,
                "duration_ps": duration_ps,
                "slope": slope,
            },
        }

    mask = np.zeros(n_points, dtype=bool)
    mask[best_start_idx : end_idx + 1] = True
    return {
        "mask": mask,
        "stable": True,
        "reason": "auto_stable_window",
        "stats": best_stats,
    }


def _analyze_transport_running_integral(
    *,
    corr_file: Path,
    T: float,
    volume: float,
    t_max_ps: Optional[float],
    plateau_start_ps: float,
    plateau_end_ps: Optional[float],
    smooth_window: int,
    rel_std_tol: float,
    drift_tol: float,
    min_ps_before_cutoff: float,
    converter_scale: float,
    component_factor: float,
    component_unit: str,
    result_key: str,
) -> Dict[str, object]:
    detail = parse_ave_correlate_detail(corr_file)
    time_fs = detail["time_fs"]
    corr = detail["corr"]
    ncount = detail["ncount"]
    if corr.shape[1] < 3:
        raise ValueError(f"{corr_file} needs at least 3 correlation components")

    time_ps = time_fs / 1000.0
    effective_t_max_ps = _estimate_effective_cutoff_ps(
        time_ps,
        corr,
        t_max_ps,
        min_ps_before_cutoff=min_ps_before_cutoff,
    )
    t_max_fs = effective_t_max_ps * 1000.0

    running_components = []
    for col in range(3):
        run, _ = _integrate_acf(time_fs, corr[:, col], t_max=t_max_fs)
        running_components.append(run * converter_scale * component_factor)

    full_time_mask = time_ps <= effective_t_max_ps + 1.0e-12
    time_plot = time_ps[full_time_mask]
    running_components_arr = np.asarray(
        [component[: len(time_plot)] for component in running_components],
        dtype=float,
    )
    running_mean = np.mean(running_components_arr, axis=0)

    smooth_window = max(1, min(int(smooth_window), len(running_mean)))
    if smooth_window > 1:
        running_smooth = _running_average(running_mean, smooth_window)
        time_smooth = time_plot[smooth_window - 1 :]
        component_smooth = np.asarray(
            [_running_average(component, smooth_window) for component in running_components_arr],
            dtype=float,
        )
    else:
        running_smooth = running_mean
        time_smooth = time_plot
        component_smooth = running_components_arr

    n_smooth = min(len(time_smooth), len(running_smooth))
    if component_smooth.ndim == 2 and component_smooth.shape[1] > 0:
        n_smooth = min(n_smooth, component_smooth.shape[1])
    time_smooth = time_smooth[:n_smooth]
    running_smooth = running_smooth[:n_smooth]
    if component_smooth.ndim == 2:
        component_smooth = component_smooth[:, :n_smooth]

    window = _select_analysis_window(
        time_smooth,
        running_smooth,
        plateau_start_ps=plateau_start_ps,
        plateau_end_ps=plateau_end_ps,
        rel_std_tol=rel_std_tol,
        drift_tol=drift_tol,
    )
    window_mask = window["mask"]
    window_vals = running_smooth[window_mask]
    stats = window["stats"]

    component_values = [float(component[window_mask].mean()) for component in component_smooth]
    diagnostics = {
        "corr_file": str(corr_file),
        "max_lag_ps_available": float(time_ps[-1]) if len(time_ps) > 0 else 0.0,
        "effective_t_max_ps": float(effective_t_max_ps),
        "n_lags_used": int(len(time_plot)),
        "ncount_nonzero_fraction": (
            float(np.count_nonzero(ncount > 0) / len(ncount))
            if len(ncount) > 0
            else 0.0
        ),
        "min_ncount": float(np.min(ncount)) if len(ncount) > 0 else 0.0,
        "max_ncount": float(np.max(ncount)) if len(ncount) > 0 else 0.0,
        "window_reason": str(window["reason"]),
        "window_stable": bool(window["stable"]),
        "window_start_ps": float(time_smooth[window_mask][0]),
        "window_end_ps": float(time_smooth[window_mask][-1]),
        "window_points": int(np.count_nonzero(window_mask)),
        "window_duration_ps": float(stats["duration_ps"]),
        "window_rel_std": float(stats["rel_std"]),
        "window_drift_fraction": float(stats["drift_fraction"]),
        "window_slope_fraction": float(stats["slope_fraction"]),
        "component_values": component_values,
        "component_unit": component_unit,
        "component_mean": float(np.mean(component_values)),
        "component_std": float(np.std(component_values)),
        "component_rel_std": _relative_metric(
            float(np.std(component_values)),
            float(np.mean(component_values)),
            1.0e-12,
        ),
    }

    return {
        result_key: float(window_vals.mean()),
        "plateau_std": float(window_vals.std()),
        "component_values": component_values,
        "time_ps": time_smooth,
        "running_integral": running_smooth,
        "window_mask": window_mask,
        "effective_t_max_ps": float(effective_t_max_ps),
        "window_start_ps": float(time_smooth[window_mask][0]),
        "window_end_ps": float(time_smooth[window_mask][-1]),
        "diagnostics": diagnostics,
    }


__all__ = [
    "_analyze_transport_running_integral",
    "parse_ave_correlate_detail",
]
