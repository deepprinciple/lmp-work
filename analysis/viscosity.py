"""Green-Kubo shear viscosity analysis from LAMMPS stress outputs.

The shear viscosity is computed from the autocorrelation of the off-diagonal
stress tensor components (pxy, pxz, pyz) via the Green-Kubo formula:

    η = V / (k_B T) × (1/3) × Σ_{αβ} ∫₀^∞ ⟨P_{αβ}(0) P_{αβ}(t)⟩ dt

LAMMPS real-unit conversion
----------------------------
- Pressure  : atm   (1 atm = 101 325 Pa)
- Volume    : Å³    (1 Å³  = 10⁻³⁰ m³)
- Time      : fs    (1 fs  = 10⁻¹⁵ s)
- k_B       : 1.380 649 × 10⁻²³ J/K

The result is returned in mPa·s (= cP), which equals Pa·s × 10³.

Two post-processing modes are supported:

* ``acf`` reads LAMMPS ``fix ave/correlate`` output directly.
* ``pressure_blocks`` reads a raw ``step pxy pxz pyz`` pressure trace, builds
  per-block autocorrelation functions, and reports a block SEM.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

try:
    from ..analysis.correlation import (
        _analyze_transport_running_integral,
        parse_ave_correlate_detail,
    )
    from ..utils.constants import (
        ANGSTROM_TO_METER,
        ATM_TO_PASCAL,
        BOLTZMANN_J_K,
        FEMTOSECOND_TO_SECOND,
    )
except ImportError:
    from analysis.correlation import (
        _analyze_transport_running_integral,
        parse_ave_correlate_detail,
    )
    from utils.constants import (
        ANGSTROM_TO_METER,
        ATM_TO_PASCAL,
        BOLTZMANN_J_K,
        FEMTOSECOND_TO_SECOND,
    )


# ──────────────────────────────────────────────────────────────────────────
# Unit-conversion constant
# ──────────────────────────────────────────────────────────────────────────

# Derivation:
#   η [Pa·s] = V [m³] / (k_B [J/K] · T [K]) × ∫ ACF [Pa²] dt [s]
#
# In LAMMPS real units the ACF has units of atm², the integral is in fs, and
# V is in Å³.  Writing out the conversions:
#
#   η [Pa·s] = V[Å³]·(Å³→m³) / (k_B·T) × ACF[atm²]·(atm→Pa)² × dt[fs]·(fs→s)
#            = ∫ ACF[atm²·fs] dt  ×  (V[Å³] / T[K])  ×  ETA_CONV
#
# where ETA_CONV = (atm→Pa)² × (Å³→m³) × (fs→s) / k_B × 1e3  (×1e3 for mPa·s)

ETA_CONV: float = (
    ATM_TO_PASCAL ** 2         # (101 325)²  Pa²/atm²
    * ANGSTROM_TO_METER ** 3   # 10⁻³⁰       m³/Å³
    * FEMTOSECOND_TO_SECOND    # 10⁻¹⁵       s/fs
    / BOLTZMANN_J_K            # 1/k_B       K/J
    * 1e3                      # Pa·s → mPa·s
)


def _resolve_output_path(base_dir: Path, output: str | Path | None) -> Path | None:
    if output is None:
        return None
    path = Path(output)
    if path.is_absolute():
        return path
    return base_dir / path


def _mean_sem(values: np.ndarray, axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(values, axis=axis)
    n = values.shape[axis]
    if n < 2:
        sem = np.full_like(mean, np.nan, dtype=np.float64)
    else:
        sem = np.std(values, axis=axis, ddof=1) / math.sqrt(n)
    return mean, sem


def _check_viscosity_window(
    diagnostics: dict[str, Any],
    *,
    plateau_start_ps: float,
    min_window_ps: float,
    require_stable_window: bool,
) -> None:
    diagnostics["requested_plateau_start_ps"] = float(plateau_start_ps)
    diagnostics["min_window_ps"] = float(min_window_ps)
    diagnostics["window_duration_ok"] = (
        float(diagnostics.get("window_duration_ps", 0.0)) + 1.0e-12 >= float(min_window_ps)
    )
    diagnostics["plateau_start_beyond_effective_t_max"] = (
        float(plateau_start_ps) > float(diagnostics.get("effective_t_max_ps", 0.0)) + 1.0e-12
    )

    problems: list[str] = []
    if diagnostics["plateau_start_beyond_effective_t_max"]:
        problems.append(
            "requested plateau_start_ps is beyond the effective integration cutoff"
        )
    if not diagnostics["window_duration_ok"]:
        problems.append(
            f"selected window duration {float(diagnostics.get('window_duration_ps', 0.0)):.6g} ps "
            f"is shorter than min_window_ps={float(min_window_ps):.6g}"
        )
    if not bool(diagnostics.get("window_stable", False)):
        problems.append("selected viscosity window did not pass stability diagnostics")

    diagnostics["window_quality"] = "ok" if not problems else "warning"
    diagnostics["window_quality_messages"] = problems

    if require_stable_window and problems:
        raise ValueError("Viscosity analysis window rejected: " + "; ".join(problems))


# ──────────────────────────────────────────────────────────────────────────
# Core computation
# ──────────────────────────────────────────────────────────────────────────

def compute_viscosity(
    corr_stress_file: Path,
    T: float,
    volume: float,
    t_max_ps: Optional[float] = None,
    plateau_start_ps: float = 2.0,
    plateau_end_ps: Optional[float] = None,
    smooth_window: int = 1,
    min_window_ps: float = 0.0,
    require_stable_window: bool = False,
) -> Dict[str, Any]:
    """Compute shear viscosity from a stress autocorrelation file.

    The file must be the output of LAMMPS ``fix ave/correlate`` containing at
    least three correlation columns (pxy×pxy, pxz×pxz, pyz×pyz) in atm².

    Parameters
    ----------
    corr_stress_file:
        Path to ``stress_acf.dat`` (output of fix ave/correlate).
    T:
        Simulation temperature (K).
    volume:
        Simulation box volume (Å³) – use the production-averaged volume.
    t_max_ps:
        Maximum lag time to include in the integral (ps).
        ``None`` → auto-detect from first zero-crossing of the mean ACF.
    plateau_start_ps:
        Earliest integration time (ps) considered for plateau detection.
    plateau_end_ps:
        Latest integration time (ps) for plateau detection.  ``None`` → end.
    smooth_window:
        Running-average window applied to the running integral before
        plateau detection (steps).  1 = no smoothing.

    Returns
    -------
    dict with keys:

    * ``eta_mPas``               – shear viscosity in mPa·s (= cP)
    * ``eta_components_mPas``    – [η_xy, η_xz, η_yz] individual estimates
    * ``plateau_std_mPas``       – std deviation over the plateau window
    * ``time_ps``                – time axis of the running integral
    * ``running_integral_mPas``  – running η vs integration time
    * ``window_mask``            – boolean mask marking the plateau region
    * ``effective_t_max_ps``     – actual upper limit used in integration
    * ``window_start_ps``, ``window_end_ps``
    * ``diagnostics``            – detailed convergence metadata
    * ``temperature_K``, ``volume_A3``
    """
    analysis = _analyze_transport_running_integral(
        corr_file=Path(corr_stress_file),
        T=T,
        volume=volume,
        t_max_ps=t_max_ps,
        plateau_start_ps=plateau_start_ps,
        plateau_end_ps=plateau_end_ps,
        smooth_window=smooth_window,
        rel_std_tol=0.20,
        drift_tol=0.20,
        min_ps_before_cutoff=max(0.5, plateau_start_ps * 0.5),
        # η = ETA_CONV × (V/T) × ∫ ACF dt   (averaged over 3 components)
        converter_scale=ETA_CONV,
        component_factor=volume / T,
        component_unit="mPa·s",
        result_key="eta_mPas",
    )
    _check_viscosity_window(
        analysis["diagnostics"],
        plateau_start_ps=float(plateau_start_ps),
        min_window_ps=float(min_window_ps),
        require_stable_window=bool(require_stable_window),
    )

    return {
        "eta_mPas": float(analysis["eta_mPas"]),
        "eta_components_mPas": list(analysis["component_values"]),
        "plateau_std_mPas": float(analysis["plateau_std"]),
        "time_ps": analysis["time_ps"],
        "running_integral_mPas": analysis["running_integral"],
        "window_mask": analysis["window_mask"],
        "effective_t_max_ps": float(analysis["effective_t_max_ps"]),
        "window_start_ps": float(analysis["window_start_ps"]),
        "window_end_ps": float(analysis["window_end_ps"]),
        "diagnostics": analysis["diagnostics"],
        "temperature_K": T,
        "volume_A3": volume,
    }


def _read_pressure_tensor(
    pressure_file: Path,
    *,
    start_step: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Read ``step pxy pxz pyz`` samples written during production."""
    path = Path(pressure_file)
    raw = np.loadtxt(path, comments="#", dtype=np.float64)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    if raw.ndim != 2 or raw.shape[1] < 4:
        raise ValueError(f"Expected columns 'step pxy pxz pyz' in {path}")

    valid = np.all(np.isfinite(raw[:, :4]), axis=1)
    if start_step is not None:
        valid &= raw[:, 0] >= int(start_step)
    if not np.any(valid):
        raise ValueError(f"No pressure tensor samples found in {path}")

    steps = np.rint(raw[valid, 0]).astype(np.int64)
    pressure_atm = np.asarray(raw[valid, 1:4], dtype=np.float64)
    if steps.size < 2:
        raise ValueError(f"Need at least two pressure tensor samples in {path}")
    if np.any(np.diff(steps) <= 0):
        raise ValueError(f"Pressure tensor timesteps must be strictly increasing: {path}")
    return steps, pressure_atm


def _infer_sample_dt_fs(steps: np.ndarray, timestep_fs: float) -> tuple[float, int]:
    diffs = np.diff(steps)
    positive = diffs[diffs > 0]
    if positive.size == 0:
        raise ValueError("Cannot infer pressure tensor sample spacing")
    step_stride = int(round(float(np.median(positive))))
    if step_stride <= 0:
        raise ValueError("Pressure tensor sample spacing must be positive")
    return float(step_stride * timestep_fs), step_stride


def _autocorrelation_unbiased(values: np.ndarray, max_lag: int) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    values = values - np.mean(values)
    n_samples = values.size
    n_fft = 1 << (2 * n_samples - 1).bit_length()
    spectrum = np.fft.rfft(values, n=n_fft)
    corr = np.fft.irfft(spectrum * np.conjugate(spectrum), n=n_fft)[: max_lag + 1]
    counts = np.arange(n_samples, n_samples - max_lag - 1, -1, dtype=np.float64)
    return corr / counts


def _block_green_kubo_pressure(
    pressure_atm: np.ndarray,
    *,
    sample_dt_fs: float,
    volume_a3: float,
    temperature_k: float,
    block_size_ps: float,
    max_corr_ps: float,
    min_blocks: int,
) -> dict[str, Any]:
    sample_dt_ps = sample_dt_fs / 1000.0
    if sample_dt_ps <= 0.0:
        raise ValueError("sample_dt_fs must be positive")
    if max_corr_ps <= 0.0:
        raise ValueError("max_corr_ps must be positive")
    if block_size_ps <= 0.0:
        raise ValueError("block_size_ps must be positive")

    max_lag = int(round(max_corr_ps / sample_dt_ps))
    block_size = int(round(block_size_ps / sample_dt_ps))
    if max_lag < 1:
        raise ValueError("max_corr_ps is shorter than one pressure sample")
    if block_size <= max_lag:
        raise ValueError("analysis.block_size_ps must be larger than analysis.t_max_ps")

    n_blocks = pressure_atm.shape[0] // block_size
    if n_blocks < int(min_blocks):
        raise ValueError(
            f"Need at least {int(min_blocks)} pressure blocks; got {n_blocks}. "
            "Increase simulation.prod_steps or reduce analysis.block_size_ps."
        )

    trimmed = pressure_atm[: n_blocks * block_size]
    time_ps = np.arange(max_lag + 1, dtype=np.float64) * sample_dt_ps
    prefactor = ETA_CONV * volume_a3 / temperature_k

    block_eta_components = np.empty((n_blocks, 3, max_lag + 1), dtype=np.float64)
    block_acf_components = np.empty((n_blocks, 3, max_lag + 1), dtype=np.float64)

    for block_index in range(n_blocks):
        block = trimmed[block_index * block_size : (block_index + 1) * block_size]
        for component in range(3):
            acf = _autocorrelation_unbiased(block[:, component], max_lag)
            integral = np.zeros_like(acf)
            integral[1:] = np.cumsum(0.5 * (acf[1:] + acf[:-1]) * sample_dt_fs)
            block_acf_components[block_index, component] = acf
            block_eta_components[block_index, component] = integral * prefactor

    block_eta_mean = np.mean(block_eta_components, axis=1)
    block_acf_mean = np.mean(block_acf_components, axis=1)
    eta_mean, eta_sem = _mean_sem(block_eta_mean, axis=0)
    acf_mean, acf_sem = _mean_sem(block_acf_mean, axis=0)

    return {
        "time_ps": time_ps,
        "block_eta_mPas": block_eta_mean,
        "block_eta_components_mPas": block_eta_components,
        "eta_mean_mPas": eta_mean,
        "eta_sem_mPas": eta_sem,
        "acf_mean_atm2": acf_mean,
        "acf_sem_atm2": acf_sem,
        "component_eta_mean_mPas": np.mean(block_eta_components, axis=0),
        "block_size_samples": block_size,
        "max_lag_samples": max_lag,
        "n_blocks": n_blocks,
        "trimmed_samples": int(trimmed.shape[0]),
    }


def _pressure_window_stats(
    time_ps: np.ndarray,
    eta_mean: np.ndarray,
    eta_sem: np.ndarray,
    mask: np.ndarray,
) -> dict[str, float]:
    t_window = time_ps[mask]
    eta_window = eta_mean[mask]
    sem_window = eta_sem[mask]
    mean_val = float(np.mean(eta_window))
    std_val = float(np.std(eta_window))
    duration_ps = float(max(t_window[-1] - t_window[0], 0.0))
    drift = float(eta_window[-1] - eta_window[0]) if eta_window.size >= 2 else 0.0
    slope = (
        float(np.polyfit(t_window, eta_window, deg=1)[0])
        if eta_window.size >= 2 and duration_ps > 0.0
        else 0.0
    )
    finite_sem = sem_window[np.isfinite(sem_window)]
    mean_sem = float(np.mean(finite_sem)) if finite_sem.size else float("nan")
    eta_abs = max(abs(mean_val), 1.0e-12)
    return {
        "mean": mean_val,
        "std": std_val,
        "rel_std": std_val / eta_abs,
        "drift": drift,
        "drift_fraction": abs(drift) / eta_abs,
        "slope": slope,
        "slope_fraction": abs(slope) * duration_ps / eta_abs,
        "duration_ps": duration_ps,
        "mean_sem": mean_sem,
        "mean_sem_fraction": (mean_sem / eta_abs) if math.isfinite(mean_sem) else float("nan"),
    }


def _select_pressure_plateau_window(
    time_ps: np.ndarray,
    eta_mean: np.ndarray,
    eta_sem: np.ndarray,
    *,
    plateau_start_ps: float,
    plateau_end_ps: float | None,
    auto_plateau: bool,
    auto_window_ps: float,
    auto_step_ps: float,
    min_window_ps: float,
    rel_std_tol: float,
    drift_tol: float,
    min_points: int = 3,
) -> dict[str, Any]:
    if time_ps.size == 0:
        raise ValueError("Cannot select a viscosity window from an empty curve")
    min_window_ps = max(0.0, float(min_window_ps))
    min_points = max(1, int(min_points))
    max_time_ps = float(time_ps[-1])

    def _mask_for(start_ps: float, end_ps: float) -> np.ndarray:
        return (time_ps >= start_ps - 1.0e-12) & (time_ps <= end_ps + 1.0e-12)

    if auto_plateau:
        window_ps = float(auto_window_ps)
        step_ps = float(auto_step_ps)
        if window_ps <= 0.0:
            raise ValueError("analysis.auto_window_ps must be positive")
        if step_ps <= 0.0:
            raise ValueError("analysis.auto_step_ps must be positive")
        if window_ps + 1.0e-12 < min_window_ps:
            raise ValueError("analysis.auto_window_ps must be >= analysis.min_window_ps")

        max_start_ps = max_time_ps - window_ps
        if plateau_end_ps is not None:
            max_start_ps = min(max_start_ps, float(plateau_end_ps) - window_ps)
        if float(plateau_start_ps) > max_start_ps + 1.0e-12:
            raise ValueError("No room for an automatic viscosity plateau window")

        best: dict[str, Any] | None = None
        for start_ps in np.arange(float(plateau_start_ps), max_start_ps + 0.5 * step_ps, step_ps):
            end_ps = float(start_ps + window_ps)
            mask = _mask_for(float(start_ps), end_ps)
            if np.count_nonzero(mask) < min_points:
                continue
            stats = _pressure_window_stats(time_ps, eta_mean, eta_sem, mask)
            if stats["duration_ps"] + 1.0e-12 < min_window_ps:
                continue
            score = (
                3.0 * stats["slope_fraction"]
                + 2.0 * stats["drift_fraction"]
                + stats["rel_std"]
                + (stats["mean_sem_fraction"] if math.isfinite(stats["mean_sem_fraction"]) else 0.0)
            )
            candidate = {
                "mask": mask,
                "reason": "auto_block_window",
                "stats": stats,
                "score": float(score),
            }
            if best is None or candidate["score"] < best["score"]:
                best = candidate

        if best is None:
            raise ValueError("Automatic viscosity plateau selection found no valid window")
        stats = best["stats"]
        stable = (
            stats["rel_std"] <= rel_std_tol
            and stats["drift_fraction"] <= drift_tol
            and stats["slope_fraction"] <= drift_tol
        )
        best["stable"] = stable
        return best

    end_ps = max_time_ps if plateau_end_ps is None else min(float(plateau_end_ps), max_time_ps)
    start_ps = float(plateau_start_ps)
    if start_ps > end_ps + 1.0e-12:
        raise ValueError("analysis.plateau_start_ps is beyond the available integration range")
    mask = _mask_for(start_ps, end_ps)
    if np.count_nonzero(mask) < min_points:
        raise ValueError("Manual viscosity plateau window has too few points")
    stats = _pressure_window_stats(time_ps, eta_mean, eta_sem, mask)
    if stats["duration_ps"] + 1.0e-12 < min_window_ps:
        raise ValueError(
            f"Manual viscosity plateau window is {stats['duration_ps']:.6g} ps, "
            f"shorter than min_window_ps={min_window_ps:.6g}"
        )
    stable = (
        stats["rel_std"] <= rel_std_tol
        and stats["drift_fraction"] <= drift_tol
        and stats["slope_fraction"] <= drift_tol
    )
    return {
        "mask": mask,
        "reason": "manual_block_window",
        "stats": stats,
        "score": float("nan"),
        "stable": stable,
    }


def compute_viscosity_from_pressure_file(
    pressure_file: Path,
    *,
    T: float,
    volume: float,
    timestep_fs: float,
    t_max_ps: float = 20.0,
    block_size_ps: float = 200.0,
    plateau_start_ps: float = 5.0,
    plateau_end_ps: float | None = None,
    auto_plateau: bool = True,
    auto_window_ps: float = 10.0,
    auto_step_ps: float = 1.0,
    min_window_ps: float = 10.0,
    min_blocks: int = 3,
    rel_std_tol: float = 0.25,
    drift_tol: float = 0.25,
    start_step: int | None = None,
    require_stable_window: bool = False,
) -> Dict[str, Any]:
    """Compute shear viscosity from a raw ``pxy/pxz/pyz`` pressure trace.

    The pressure trace is split into independent blocks.  Each block gets its
    own FFT autocorrelation and Green-Kubo running integral, making the final
    plateau estimate report a block SEM instead of only within-curve roughness.
    """
    pressure_path = Path(pressure_file).resolve()
    if not pressure_path.exists():
        raise FileNotFoundError(f"Pressure tensor file not found: {pressure_path}")

    steps, pressure_atm = _read_pressure_tensor(pressure_path, start_step=start_step)
    sample_dt_fs, sample_step_stride = _infer_sample_dt_fs(steps, timestep_fs)
    block = _block_green_kubo_pressure(
        pressure_atm,
        sample_dt_fs=sample_dt_fs,
        volume_a3=float(volume),
        temperature_k=float(T),
        block_size_ps=float(block_size_ps),
        max_corr_ps=float(t_max_ps),
        min_blocks=int(min_blocks),
    )

    window = _select_pressure_plateau_window(
        block["time_ps"],
        block["eta_mean_mPas"],
        block["eta_sem_mPas"],
        plateau_start_ps=float(plateau_start_ps),
        plateau_end_ps=plateau_end_ps,
        auto_plateau=bool(auto_plateau),
        auto_window_ps=float(auto_window_ps),
        auto_step_ps=float(auto_step_ps),
        min_window_ps=float(min_window_ps),
        rel_std_tol=float(rel_std_tol),
        drift_tol=float(drift_tol),
    )
    window_mask = np.asarray(window["mask"], dtype=bool)
    block_plateau = np.mean(block["block_eta_mPas"][:, window_mask], axis=1)
    eta_mpas = float(np.mean(block_plateau))
    eta_sem = (
        float(np.std(block_plateau, ddof=1) / math.sqrt(block_plateau.size))
        if block_plateau.size > 1
        else float("nan")
    )
    component_values = [
        float(np.mean(block["component_eta_mean_mPas"][component, window_mask]))
        for component in range(3)
    ]
    stats = dict(window["stats"])
    window_time = block["time_ps"][window_mask]

    diagnostics: dict[str, Any] = {
        "analysis_method": "pressure_blocks",
        "pressure_file": str(pressure_path),
        "selected_start_step": int(steps[0]),
        "selected_end_step": int(steps[-1]),
        "n_samples": int(pressure_atm.shape[0]),
        "sample_step_stride": int(sample_step_stride),
        "sample_dt_fs": float(sample_dt_fs),
        "max_corr_ps": float(t_max_ps),
        "block_size_ps": float(block_size_ps),
        "block_size_samples": int(block["block_size_samples"]),
        "max_lag_samples": int(block["max_lag_samples"]),
        "n_blocks": int(block["n_blocks"]),
        "min_blocks": int(min_blocks),
        "trimmed_samples": int(block["trimmed_samples"]),
        "window_reason": str(window["reason"]),
        "window_stable": bool(window["stable"]),
        "window_start_ps": float(window_time[0]),
        "window_end_ps": float(window_time[-1]),
        "window_points": int(np.count_nonzero(window_mask)),
        "window_duration_ps": float(stats["duration_ps"]),
        "window_rel_std": float(stats["rel_std"]),
        "window_drift_fraction": float(stats["drift_fraction"]),
        "window_slope_fraction": float(stats["slope_fraction"]),
        "window_mean_sem_mPas": float(stats["mean_sem"]),
        "window_mean_sem_fraction": float(stats["mean_sem_fraction"]),
        "window_score": float(window["score"]),
        "component_values": component_values,
        "component_unit": "mPa·s",
        "component_mean": float(np.mean(component_values)),
        "component_std": float(np.std(component_values)),
        "component_rel_std": float(np.std(component_values) / max(abs(np.mean(component_values)), 1.0e-12)),
    }
    diagnostics["window_quality"] = "ok" if diagnostics["window_stable"] else "warning"
    diagnostics["window_quality_messages"] = (
        [] if diagnostics["window_stable"] else ["selected viscosity window did not pass stability diagnostics"]
    )
    if require_stable_window and not diagnostics["window_stable"]:
        raise ValueError("Viscosity analysis window rejected: selected window did not pass stability diagnostics")

    return {
        "eta_mPas": eta_mpas,
        "eta_sem_mPas": eta_sem,
        "eta_components_mPas": component_values,
        "plateau_std_mPas": float(np.std(block["eta_mean_mPas"][window_mask])),
        "block_eta_plateau_mPas": block_plateau,
        "time_ps": block["time_ps"],
        "running_integral_mPas": block["eta_mean_mPas"],
        "running_integral_sem_mPas": block["eta_sem_mPas"],
        "acf_mean_atm2": block["acf_mean_atm2"],
        "acf_sem_atm2": block["acf_sem_atm2"],
        "window_mask": window_mask,
        "effective_t_max_ps": float(block["time_ps"][-1]),
        "window_start_ps": float(window_time[0]),
        "window_end_ps": float(window_time[-1]),
        "diagnostics": diagnostics,
        "temperature_K": T,
        "volume_A3": volume,
    }


# ──────────────────────────────────────────────────────────────────────────
# Volume helper
# ──────────────────────────────────────────────────────────────────────────

def load_volume_from_thermo(thermo_file: Path) -> float:
    """Return the mean production volume (Å³) from *gk_thermo.dat*.

    Expected format (written by ``fix print`` in ``build_gk_input``)::

        # step temp press vol
        1000  300.1  -5.3  28341.7
        ...

    Parameters
    ----------
    thermo_file:
        Path to the thermo log written during the production run.

    Returns
    -------
    float
        Mean volume in Å³.
    """
    volumes: list[float] = []
    with open(thermo_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    volumes.append(float(parts[3]))
                except ValueError:
                    continue
    if not volumes:
        raise ValueError(f"No volume data found in {thermo_file}")
    return float(np.mean(volumes))


# ──────────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────────

def _import_matplotlib_pyplot():
    """Import matplotlib.pyplot with a clear error for NumPy 2.x mismatches."""
    import matplotlib

    mpl_parts = matplotlib.__version__.split(".")[:2]
    np_parts = np.__version__.split(".")[:2]
    mpl_major_minor = tuple(int(part) for part in mpl_parts)
    np_major_minor = tuple(int(part) for part in np_parts)
    if np_major_minor >= (2, 0) and mpl_major_minor < (3, 9):
        raise ImportError(
            "matplotlib>=3.9 is required with NumPy 2.x "
            f"(found matplotlib {matplotlib.__version__}, "
            f"numpy {np.__version__}). "
            "Install with: python -m pip install 'matplotlib>=3.9'"
        )
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _mark_plateau_window(
    ax: Any,
    *,
    start_ps: float,
    end_ps: float,
    label: str = "plateau window",
) -> None:
    """Mark an analysis window without relying on matplotlib patch spans."""
    start = float(start_ps)
    end = float(end_ps)
    if not (math.isfinite(start) and math.isfinite(end)):
        return
    ax.axvline(start, color="tab:green", ls="--", lw=1.0, alpha=0.85, label=label)
    if abs(end - start) > 1.0e-12:
        ax.axvline(end, color="tab:green", ls=":", lw=1.0, alpha=0.85)


def _write_viscosity_plot(
    detail: Dict[str, Any],
    result: Dict[str, Any],
    *,
    plot_out: Path,
) -> None:
    plt = _import_matplotlib_pyplot()

    time_ps = np.asarray(result["time_ps"], dtype=float)
    running_eta = np.asarray(result["running_integral_mPas"], dtype=float)

    acf_time_ps = np.asarray(detail["time_fs"], dtype=float) / 1000.0
    acf_mat = np.asarray(detail["corr"], dtype=float)
    acf_mean = np.mean(acf_mat, axis=1)
    norm = np.max(np.abs(acf_mean))
    acf_norm = acf_mean / norm if norm > 0 else acf_mean

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0))

    # Panel 1: stress ACF
    labels = ["pxy", "pxz", "pyz"]
    for col, lbl in enumerate(labels):
        col_norm = acf_mat[:, col] / norm if norm > 0 else acf_mat[:, col]
        axes[0].plot(acf_time_ps, col_norm, lw=0.8, alpha=0.55, label=lbl)
    axes[0].plot(acf_time_ps, acf_norm, lw=1.6, color="k", label="mean")
    axes[0].axhline(0.0, color="gray", lw=0.8, ls="--")
    axes[0].set_xlabel("Lag time (ps)")
    axes[0].set_ylabel("Normalised stress ACF")
    axes[0].set_title("Stress autocorrelation function")
    axes[0].legend(fontsize=8)

    # Panel 2: running viscosity
    eta_components = result["eta_components_mPas"]
    comp_labels = ["η_xy", "η_xz", "η_yz"]
    # (running integrals per component are not stored separately; show mean only)
    axes[1].plot(time_ps, running_eta, lw=1.6, color="tab:blue", label="running η (mean)")
    _mark_plateau_window(
        axes[1],
        start_ps=float(result["window_start_ps"]),
        end_ps=float(result["window_end_ps"]),
    )
    axes[1].axhline(
        float(result["eta_mPas"]),
        color="tab:red",
        ls="--",
        lw=1.4,
        label=f"η = {result['eta_mPas']:.3f} mPa·s",
    )
    axes[1].set_xlabel("Integration time (ps)")
    axes[1].set_ylabel("Shear viscosity (mPa·s)")
    axes[1].set_title("Green-Kubo viscosity convergence")
    axes[1].legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(plot_out, dpi=180)
    plt.close(fig)


def _write_pressure_block_plot(
    result: Dict[str, Any],
    *,
    plot_out: Path,
) -> None:
    plt = _import_matplotlib_pyplot()

    time_ps = np.asarray(result["time_ps"], dtype=float)
    running_eta = np.asarray(result["running_integral_mPas"], dtype=float)
    running_sem = np.asarray(result["running_integral_sem_mPas"], dtype=float)
    acf_mean = np.asarray(result["acf_mean_atm2"], dtype=float)
    acf_sem = np.asarray(result["acf_sem_atm2"], dtype=float)
    norm = np.max(np.abs(acf_mean))
    acf_norm = acf_mean / norm if norm > 0 else acf_mean
    acf_sem_norm = acf_sem / norm if norm > 0 else acf_sem

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0))

    axes[0].plot(time_ps, acf_norm, lw=1.5, color="k", label="mean block ACF")
    if np.any(np.isfinite(acf_sem_norm)):
        axes[0].fill_between(
            time_ps,
            acf_norm - acf_sem_norm,
            acf_norm + acf_sem_norm,
            alpha=0.2,
            linewidth=0,
            label="SEM",
        )
    axes[0].axhline(0.0, color="gray", lw=0.8, ls="--")
    axes[0].set_xlabel("Lag time (ps)")
    axes[0].set_ylabel("Normalised stress ACF")
    axes[0].set_title("Block-averaged stress autocorrelation")
    axes[0].legend(fontsize=8)

    axes[1].plot(time_ps, running_eta, lw=1.6, color="tab:blue", label="running η")
    if np.any(np.isfinite(running_sem)):
        axes[1].fill_between(
            time_ps,
            running_eta - running_sem,
            running_eta + running_sem,
            alpha=0.22,
            linewidth=0,
            label="block SEM",
        )
    _mark_plateau_window(
        axes[1],
        start_ps=float(result["window_start_ps"]),
        end_ps=float(result["window_end_ps"]),
    )
    axes[1].axhline(
        float(result["eta_mPas"]),
        color="tab:red",
        ls="--",
        lw=1.4,
        label=f"η = {result['eta_mPas']:.3f} ± {result['eta_sem_mPas']:.3f} mPa·s",
    )
    axes[1].set_xlabel("Integration time (ps)")
    axes[1].set_ylabel("Shear viscosity (mPa·s)")
    axes[1].set_title("Green-Kubo viscosity convergence")
    axes[1].legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(plot_out, dpi=180)
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────
# High-level entry point
# ──────────────────────────────────────────────────────────────────────────

def analyze_viscosity_file(
    stress_acf_file: Path,
    *,
    temperature_k: float,
    volume_a3: float,
    t_max_ps: float = 20.0,
    plateau_start_ps: float = 2.0,
    plateau_end_ps: float | None = None,
    smooth_window: int = 11,
    min_window_ps: float = 0.0,
    require_stable_window: bool = False,
    json_out: str | Path | None = "viscosity_summary.json",
    plot_out: str | Path | None = "viscosity_analysis.png",
) -> Dict[str, Any]:
    """Parse a stress ACF file, compute η, and write plots + JSON summary.

    Parameters
    ----------
    stress_acf_file:
        Path to the ``stress_acf.dat`` file from LAMMPS.
    temperature_k:
        Simulation temperature in Kelvin.
    volume_a3:
        Production-averaged volume in Å³ (from ``load_volume_from_thermo``).
    t_max_ps:
        Maximum ACF lag time to integrate (ps).
    plateau_start_ps:
        Earliest time (ps) to search for a stable plateau.
    smooth_window:
        Smoothing window for the running integral plot.
    json_out:
        Where to write the JSON summary.  Relative paths are anchored to
        the directory of *stress_acf_file*.  ``None`` to skip.
    plot_out:
        Where to write the PNG plot.  Same anchoring rules.  ``None`` to skip.

    Returns
    -------
    dict
        Summary dictionary (same content as the JSON file).
    """
    acf_path = Path(stress_acf_file).resolve()
    if not acf_path.exists():
        raise FileNotFoundError(f"Stress ACF file not found: {acf_path}")

    detail = parse_ave_correlate_detail(acf_path)
    result = compute_viscosity(
        corr_stress_file=acf_path,
        T=float(temperature_k),
        volume=float(volume_a3),
        t_max_ps=float(t_max_ps),
        plateau_start_ps=float(plateau_start_ps),
        plateau_end_ps=plateau_end_ps,
        smooth_window=int(smooth_window),
        min_window_ps=float(min_window_ps),
        require_stable_window=bool(require_stable_window),
    )

    plot_path = _resolve_output_path(acf_path.parent, plot_out)
    plot_error: str | None = None
    summary: Dict[str, Any] = {
        "analysis_method": "acf",
        "stress_acf_file": str(acf_path),
        "temperature_K": float(temperature_k),
        "volume_A3": float(volume_a3),
        "eta_mPas": float(result["eta_mPas"]),
        "eta_cP": float(result["eta_mPas"]),   # 1 mPa·s ≡ 1 cP
        "eta_components_mPas": result["eta_components_mPas"],
        "plateau_start_ps": float(result["window_start_ps"]),
        "plateau_end_ps": float(result["window_end_ps"]),
        "effective_t_max_ps": float(result["effective_t_max_ps"]),
        "plot_file": str(plot_path) if plot_path is not None else None,
        "diagnostics": result["diagnostics"],
    }

    if plot_path is not None:
        try:
            _write_viscosity_plot(detail, result, plot_out=plot_path)
        except Exception as exc:  # noqa: BLE001 - plotting is optional output.
            plot_error = f"{type(exc).__name__}: {exc}"
            summary["plot_file"] = None
            summary["plot_error"] = plot_error
            print(f"  [warn] viscosity plot failed: {plot_error}")

    json_path = _resolve_output_path(acf_path.parent, json_out)
    if json_path is not None:
        json_path.write_text(json.dumps(summary, indent=2))

    return summary


def analyze_viscosity_pressure_file(
    pressure_file: Path,
    *,
    temperature_k: float,
    volume_a3: float,
    timestep_fs: float,
    t_max_ps: float = 20.0,
    block_size_ps: float = 200.0,
    plateau_start_ps: float = 5.0,
    plateau_end_ps: float | None = None,
    auto_plateau: bool = True,
    auto_window_ps: float = 10.0,
    auto_step_ps: float = 1.0,
    min_window_ps: float = 10.0,
    min_blocks: int = 3,
    rel_std_tol: float = 0.25,
    drift_tol: float = 0.25,
    start_step: int | None = None,
    require_stable_window: bool = False,
    json_out: str | Path | None = "viscosity_summary.json",
    plot_out: str | Path | None = "viscosity_analysis.png",
    running_out: str | Path | None = "viscosity_running.csv",
    blocks_out: str | Path | None = "viscosity_blocks.csv",
) -> Dict[str, Any]:
    """Analyze a raw pressure tensor trace using block Green-Kubo estimates."""
    pressure_path = Path(pressure_file).resolve()
    result = compute_viscosity_from_pressure_file(
        pressure_path,
        T=float(temperature_k),
        volume=float(volume_a3),
        timestep_fs=float(timestep_fs),
        t_max_ps=float(t_max_ps),
        block_size_ps=float(block_size_ps),
        plateau_start_ps=float(plateau_start_ps),
        plateau_end_ps=plateau_end_ps,
        auto_plateau=bool(auto_plateau),
        auto_window_ps=float(auto_window_ps),
        auto_step_ps=float(auto_step_ps),
        min_window_ps=float(min_window_ps),
        min_blocks=int(min_blocks),
        rel_std_tol=float(rel_std_tol),
        drift_tol=float(drift_tol),
        start_step=start_step,
        require_stable_window=bool(require_stable_window),
    )

    running_path = _resolve_output_path(pressure_path.parent, running_out)
    if running_path is not None:
        curve = np.column_stack(
            (
                np.asarray(result["time_ps"], dtype=float),
                np.asarray(result["running_integral_mPas"], dtype=float),
                np.asarray(result["running_integral_sem_mPas"], dtype=float),
                np.asarray(result["acf_mean_atm2"], dtype=float),
                np.asarray(result["acf_sem_atm2"], dtype=float),
            )
        )
        np.savetxt(
            running_path,
            curve,
            delimiter=",",
            header="time_ps,eta_mPas_mean,eta_mPas_sem,acf_atm2_mean,acf_atm2_sem",
            comments="",
        )

    block_path = _resolve_output_path(pressure_path.parent, blocks_out)
    if block_path is not None:
        block_eta = np.asarray(result["block_eta_plateau_mPas"], dtype=float)
        np.savetxt(
            block_path,
            np.column_stack((np.arange(1, block_eta.size + 1), block_eta)),
            delimiter=",",
            header="block,eta_mPas_plateau",
            comments="",
        )

    plot_path = _resolve_output_path(pressure_path.parent, plot_out)
    plot_error: str | None = None
    summary: Dict[str, Any] = {
        "analysis_method": "pressure_blocks",
        "pressure_file": str(pressure_path),
        "temperature_K": float(temperature_k),
        "volume_A3": float(volume_a3),
        "eta_mPas": float(result["eta_mPas"]),
        "eta_cP": float(result["eta_mPas"]),
        "eta_sem_mPas": float(result["eta_sem_mPas"]),
        "eta_components_mPas": result["eta_components_mPas"],
        "plateau_start_ps": float(result["window_start_ps"]),
        "plateau_end_ps": float(result["window_end_ps"]),
        "effective_t_max_ps": float(result["effective_t_max_ps"]),
        "running_file": str(running_path) if running_path is not None else None,
        "blocks_file": str(block_path) if block_path is not None else None,
        "plot_file": str(plot_path) if plot_path is not None else None,
        "diagnostics": result["diagnostics"],
    }

    if plot_path is not None:
        try:
            _write_pressure_block_plot(result, plot_out=plot_path)
        except Exception as exc:  # noqa: BLE001 - plotting is optional output.
            plot_error = f"{type(exc).__name__}: {exc}"
            summary["plot_file"] = None
            summary["plot_error"] = plot_error
            print(f"  [warn] viscosity plot failed: {plot_error}")

    json_path = _resolve_output_path(pressure_path.parent, json_out)
    if json_path is not None:
        json_path.write_text(json.dumps(summary, indent=2))

    return summary


__all__ = [
    "analyze_viscosity_pressure_file",
    "analyze_viscosity_file",
    "compute_viscosity_from_pressure_file",
    "compute_viscosity",
    "load_volume_from_thermo",
    "ETA_CONV",
]
