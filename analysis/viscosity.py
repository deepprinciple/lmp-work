"""Green-Kubo shear viscosity analysis from LAMMPS ``fix ave/correlate`` output.

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

This module re-uses the generic ``_analyze_transport_running_integral``
engine from ``analysis.hfacf``, overriding only the unit-conversion
constants and physical prefactors appropriate for viscosity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

try:
    from ..analysis.hfacf import (
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
    from analysis.hfacf import (
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
        Simulation box volume (Å³) – use the NVE-averaged volume.
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


# ──────────────────────────────────────────────────────────────────────────
# Volume helper
# ──────────────────────────────────────────────────────────────────────────

def load_volume_from_thermo(thermo_file: Path) -> float:
    """Return the mean NVE volume (Å³) from *gk_thermo.dat*.

    Expected format (written by ``fix print`` in ``build_gk_input``)::

        # step temp press vol
        1000  300.1  -5.3  28341.7
        ...

    Parameters
    ----------
    thermo_file:
        Path to the thermo log written during the NVE production run.

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

def _write_viscosity_plot(
    detail: Dict[str, Any],
    result: Dict[str, Any],
    *,
    plot_out: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
    axes[1].axvspan(
        float(result["window_start_ps"]),
        float(result["window_end_ps"]),
        color="tab:green",
        alpha=0.15,
        label="plateau window",
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
    smooth_window: int = 11,
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
        NVE-averaged volume in Å³ (from ``load_volume_from_thermo``).
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

    def _resolve(p: str | Path | None) -> Path | None:
        if p is None:
            return None
        pp = Path(p)
        return pp if pp.is_absolute() else acf_path.parent / pp

    detail = parse_ave_correlate_detail(acf_path)
    result = compute_viscosity(
        corr_stress_file=acf_path,
        T=float(temperature_k),
        volume=float(volume_a3),
        t_max_ps=float(t_max_ps),
        plateau_start_ps=float(plateau_start_ps),
        smooth_window=int(smooth_window),
    )

    plot_path = _resolve(plot_out)
    if plot_path is not None:
        _write_viscosity_plot(detail, result, plot_out=plot_path)

    summary: Dict[str, Any] = {
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

    json_path = _resolve(json_out)
    if json_path is not None:
        json_path.write_text(json.dumps(summary, indent=2))

    return summary


__all__ = [
    "analyze_viscosity_file",
    "compute_viscosity",
    "load_volume_from_thermo",
    "ETA_CONV",
]
