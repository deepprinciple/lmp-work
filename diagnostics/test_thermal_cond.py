#!/usr/bin/env python3
"""
Standalone thermal conductivity diagnostic script for EMD trajectories.

Purpose:
- Re-process existing EMD heat-flux outputs without modifying the main workflow
- Test whether alternative post-processing can reduce the large kappa bias
- Keep all logic local to this file

Methods implemented:
1. Einstein / Helfand relation on the integrated heat-flux density time series
2. Double-exponential fit to the HCACF with analytic integration

Expected inputs:
- gk_data.dat : columns "# step T V pxy pxz pyz Jx Jy Jz Mx My Mz"
- hfacf.dat   : LAMMPS fix ave/correlate output ("ave running" compatible)

Notes on units:
- This workflow writes heat-flux density j = J / V in LAMMPS real units:
    kcal / (mol * A^2 * fs)
- Therefore both diagnostics use the same prefactor as the existing GK analysis:
    kappa = KAPPA_CONV * V / T^2 * integral_term
  where integral_term has units of:
    kcal^2 / (mol^2 * A^4) * fs
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# SI constants
BOLTZMANN_J_K = 1.380649e-23
AVOGADRO = 6.02214076e23
ANGSTROM_TO_METER = 1.0e-10
FEMTOSECOND_TO_SECOND = 1.0e-15
KCAL_MOL_TO_JOULE = 4184.0


# Heat-flux density conversion:
# 1 kcal / (mol * A^2 * fs) -> W / m^2
HF_TO_SI = (
    KCAL_MOL_TO_JOULE
    / AVOGADRO
    / (ANGSTROM_TO_METER ** 2)
    / FEMTOSECOND_TO_SECOND
)

# Integral conversion:
# [kcal^2 / (mol^2 * A^4) * fs] -> [W^2 / m^4 * s]
# then multiply by V[A^3] / (kB * T^2)
KAPPA_CONV = (
    (ANGSTROM_TO_METER ** 3)
    * (HF_TO_SI ** 2)
    * FEMTOSECOND_TO_SECOND
    / BOLTZMANN_J_K
)


GK_COLUMNS = [
    "step",
    "T",
    "V",
    "pxy",
    "pxz",
    "pyz",
    "Jx",
    "Jy",
    "Jz",
    "Mx",
    "My",
    "Mz",
]


def positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("value must be > 0")
    return ivalue


def nonnegative_float(value: str) -> float:
    fvalue = float(value)
    if fvalue < 0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return fvalue


def moving_average(y: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return y.copy()
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(y, kernel, mode="same")


def fit_r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if ss_tot <= 0.0:
        return 1.0 if ss_res <= 1.0e-30 else float("nan")
    return 1.0 - ss_res / ss_tot


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    slope, intercept = np.polyfit(x, y, deg=1)
    y_fit = slope * x + intercept
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": fit_r2(y, y_fit),
    }


def double_exponential(t_ps: np.ndarray, a: float, tau1_ps: float, b: float, tau2_ps: float) -> np.ndarray:
    return a * np.exp(-t_ps / tau1_ps) + b * np.exp(-t_ps / tau2_ps)


def parse_gk_data(gk_data_path: Path, timestep_fs: float) -> pd.DataFrame:
    df = pd.read_csv(
        gk_data_path,
        comment="#",
        delim_whitespace=True,
        header=None,
        names=GK_COLUMNS,
    )
    if df.empty:
        raise ValueError(f"No numeric rows in {gk_data_path}")
    df["step"] = df["step"].astype(float)
    step0 = float(df["step"].iloc[0])
    df["time_fs"] = (df["step"] - step0) * timestep_fs
    df["time_ps"] = df["time_fs"] / 1000.0
    return df


def parse_ave_correlate_last_block(filepath: Path, timestep_fs: float) -> dict[str, np.ndarray]:
    numeric_rows: list[list[float]] = []
    row_lengths: list[int] = []

    with open(filepath) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                row = [float(x) for x in parts]
            except ValueError:
                continue
            numeric_rows.append(row)
            row_lengths.append(len(row))

    if not numeric_rows:
        raise ValueError(f"No numeric rows found in {filepath}")

    count_by_len = Counter(row_lengths)
    data_len = max(count_by_len.items(), key=lambda item: (item[1], item[0]))[0]
    if data_len < 5:
        raise ValueError(f"Unexpected ave/correlate data width {data_len} in {filepath}")

    blocks: list[np.ndarray] = []
    current_block: list[list[float]] = []

    for row in numeric_rows:
        if len(row) == data_len:
            current_block.append(row)
            continue
        if len(row) == 2:
            if current_block:
                blocks.append(np.asarray(current_block, dtype=float))
                current_block = []
            continue

    if current_block:
        blocks.append(np.asarray(current_block, dtype=float))

    if not blocks:
        filtered = [row for row in numeric_rows if len(row) == data_len]
        if not filtered:
            raise ValueError(f"No valid correlation rows found in {filepath}")
        blocks = [np.asarray(filtered, dtype=float)]

    arr = blocks[-1]

    # Compatible with:
    # [Index, TimeDelta, Ncount, corr1, corr2, corr3]
    # [TimeDelta, Ncount, corr1, corr2, corr3]
    if arr.shape[1] >= 6:
        col0 = arr[:, 0]
        rounded = np.rint(col0)
        diffs = np.diff(rounded)
        looks_like_index = (
            np.all(np.isclose(col0, rounded))
            and len(diffs) > 0
            and np.all(diffs == 1)
        )
        if looks_like_index:
            lag_steps = arr[:, 1]
            ncount = arr[:, 2]
            corr = arr[:, 3:6]
        else:
            lag_steps = arr[:, 0]
            ncount = arr[:, 1]
            corr = arr[:, 2:5]
    else:
        lag_steps = arr[:, 0]
        ncount = arr[:, 1]
        corr = arr[:, 2:5]

    return {
        "lag_steps": lag_steps.astype(float),
        "lag_fs": lag_steps.astype(float) * timestep_fs,
        "lag_ps": lag_steps.astype(float) * timestep_fs / 1000.0,
        "ncount": ncount.astype(float),
        "corr": corr.astype(float),
        "n_blocks": len(blocks),
    }


def autodetect_repdirs(case_dir: Path) -> list[Path]:
    repdirs = []
    for gk_path in sorted(case_dir.glob("**/gk_data.dat")):
        repdir = gk_path.parent
        if (repdir / "hfacf.dat").exists():
            repdirs.append(repdir)
    if not repdirs:
        raise FileNotFoundError(f"No replica directories with gk_data.dat + hfacf.dat found under {case_dir}")
    return repdirs


def rep_label(repdir: Path, case_dir: Path | None = None) -> str:
    if case_dir is not None:
        try:
            rel = repdir.relative_to(case_dir)
            parts = rel.parts
            if len(parts) >= 2 and parts[-1] == "replica_01" and parts[-2].startswith("replica_"):
                return parts[-2]
            return "__".join(parts)
        except ValueError:
            pass
    return repdir.name


def compute_einstein_diagnostic(
    gk_df: pd.DataFrame,
    temperature_k: float,
    volume_a3: float,
    max_lag_ps: float,
    fit_start_ps: float,
    fit_end_ps: float,
    subtract_mean_flux: bool,
    output_png: Path,
) -> dict[str, float | list[float]]:
    flux = gk_df[["Jx", "Jy", "Jz"]].to_numpy(dtype=float)
    if subtract_mean_flux:
        flux = flux - np.mean(flux, axis=0, keepdims=True)

    time_fs = gk_df["time_fs"].to_numpy(dtype=float)
    if len(time_fs) < 3:
        raise ValueError("Not enough gk_data points for Einstein diagnostic")

    dt_fs = float(np.median(np.diff(time_fs)))
    if not np.allclose(np.diff(time_fs), dt_fs, rtol=1.0e-6, atol=1.0e-6):
        raise ValueError("gk_data time spacing is not uniform enough for this diagnostic")

    w = np.cumsum(flux, axis=0) * dt_fs
    max_lag_idx = int(round(max_lag_ps * 1000.0 / dt_fs))
    max_lag_idx = max(2, min(max_lag_idx, len(w) - 1))

    lag_fs = np.arange(max_lag_idx + 1, dtype=float) * dt_fs
    lag_ps = lag_fs / 1000.0
    msd_vec = np.zeros(max_lag_idx + 1, dtype=float)

    for lag in range(1, max_lag_idx + 1):
        dw = w[lag:] - w[:-lag]
        msd_vec[lag] = np.mean(np.sum(dw * dw, axis=1))

    fit_mask = (lag_ps >= fit_start_ps) & (lag_ps <= fit_end_ps)
    if np.count_nonzero(fit_mask) < 3:
        raise ValueError("Einstein fit window contains fewer than 3 points")

    fit = linear_fit(lag_fs[fit_mask], msd_vec[fit_mask])
    slope_vec_fs = float(fit["slope"])
    kappa_w_mk = KAPPA_CONV * (volume_a3 / (temperature_k ** 2)) * (slope_vec_fs / 6.0)

    line = fit["slope"] * lag_fs + fit["intercept"]
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(lag_ps, msd_vec, label="Einstein MSD", lw=1.6)
    plt.plot(lag_ps[fit_mask], line[fit_mask], "--", label="Linear fit", lw=1.8)
    plt.xlabel("Lag time (ps)")
    plt.ylabel(r"MSD of integrated heat flux density")
    plt.title("Einstein/Helfand thermal-conductivity diagnostic")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png, dpi=180)
    plt.close()

    return {
        "dt_fs": dt_fs,
        "max_lag_ps": float(lag_ps[-1]),
        "fit_start_ps": fit_start_ps,
        "fit_end_ps": fit_end_ps,
        "slope_vector_per_fs": slope_vec_fs,
        "intercept": float(fit["intercept"]),
        "r2": float(fit["r2"]),
        "kappa_W_mK": float(kappa_w_mk),
        "mean_flux_components": np.mean(gk_df[["Jx", "Jy", "Jz"]].to_numpy(dtype=float), axis=0).tolist(),
        "subtract_mean_flux": bool(subtract_mean_flux),
    }


def compute_biexponential_diagnostic(
    hfacf_detail: dict[str, np.ndarray],
    temperature_k: float,
    volume_a3: float,
    fit_start_ps: float,
    fit_end_ps: float,
    smooth_window: int,
    output_png: Path,
) -> dict[str, float | int]:
    lag_ps = hfacf_detail["lag_ps"]
    lag_fs = hfacf_detail["lag_fs"]
    corr = hfacf_detail["corr"]
    cavg = np.mean(corr, axis=1)
    cavg_smooth = moving_average(cavg, smooth_window)

    fit_mask = (lag_ps >= fit_start_ps) & (lag_ps <= fit_end_ps)
    if np.count_nonzero(fit_mask) < 5:
        raise ValueError("Biexponential fit window contains fewer than 5 points")

    x = lag_ps[fit_mask]
    y = cavg_smooth[fit_mask]

    y0 = float(max(np.max(y), 1.0e-30))
    p0 = [0.7 * y0, 0.02, 0.3 * y0, 0.5]
    lower = [0.0, 1.0e-6, 0.0, 1.0e-6]
    upper = [np.inf, max(10.0, fit_end_ps * 10.0), np.inf, max(50.0, fit_end_ps * 20.0)]

    params, _ = curve_fit(
        double_exponential,
        x,
        y,
        p0=p0,
        bounds=(lower, upper),
        maxfev=50000,
    )
    a, tau1_ps, b, tau2_ps = [float(v) for v in params]
    y_fit = double_exponential(x, *params)
    r2 = fit_r2(y, y_fit)

    integral_ps = a * tau1_ps + b * tau2_ps
    integral_fs = integral_ps * 1000.0
    kappa_w_mk = KAPPA_CONV * (volume_a3 / (temperature_k ** 2)) * integral_fs

    full_fit = double_exponential(lag_ps, *params)
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(lag_ps, cavg, label="Raw HCACF avg", lw=1.0, alpha=0.55)
    if smooth_window > 1:
        plt.plot(lag_ps, cavg_smooth, label=f"Smoothed HCACF (window={smooth_window})", lw=1.2)
    plt.plot(lag_ps[fit_mask], full_fit[fit_mask], "--", label="Double-exp fit", lw=1.8)
    plt.xlim(0.0, max(fit_end_ps * 1.2, fit_end_ps + 0.2))
    plt.xlabel("Lag time (ps)")
    plt.ylabel("HCACF")
    plt.title("HCACF double-exponential diagnostic")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png, dpi=180)
    plt.close()

    return {
        "n_blocks_in_source": int(hfacf_detail["n_blocks"]),
        "fit_start_ps": fit_start_ps,
        "fit_end_ps": fit_end_ps,
        "smooth_window": int(smooth_window),
        "A": a,
        "tau1_ps": tau1_ps,
        "B": b,
        "tau2_ps": tau2_ps,
        "integral_ps": float(integral_ps),
        "integral_fs": float(integral_fs),
        "r2": float(r2),
        "kappa_W_mK": float(kappa_w_mk),
        "lag_dt_fs": float(lag_fs[1] - lag_fs[0]) if len(lag_fs) > 1 else float("nan"),
    }


def process_replica(
    repdir: Path,
    outdir: Path,
    timestep_fs: float,
    temperature_k_arg: float | None,
    einstein_max_lag_ps: float,
    einstein_fit_start_ps: float,
    einstein_fit_end_ps: float,
    subtract_mean_flux: bool,
    biexp_fit_start_ps: float,
    biexp_fit_end_ps: float,
    biexp_smooth_window: int,
) -> dict[str, object]:
    gk_data_path = repdir / "gk_data.dat"
    hfacf_path = repdir / "hfacf.dat"
    if not gk_data_path.exists() or not hfacf_path.exists():
        raise FileNotFoundError(f"Missing gk_data.dat or hfacf.dat in {repdir}")

    gk_df = parse_gk_data(gk_data_path, timestep_fs=timestep_fs)
    hfacf_detail = parse_ave_correlate_last_block(hfacf_path, timestep_fs=timestep_fs)

    temperature_k = float(temperature_k_arg if temperature_k_arg is not None else gk_df["T"].mean())
    volume_a3 = float(gk_df["V"].mean())

    outdir.mkdir(parents=True, exist_ok=True)
    einstein_png = outdir / "einstein_msd.png"
    biexp_png = outdir / "hfacf_biexp_fit.png"

    einstein = compute_einstein_diagnostic(
        gk_df=gk_df,
        temperature_k=temperature_k,
        volume_a3=volume_a3,
        max_lag_ps=einstein_max_lag_ps,
        fit_start_ps=einstein_fit_start_ps,
        fit_end_ps=einstein_fit_end_ps,
        subtract_mean_flux=subtract_mean_flux,
        output_png=einstein_png,
    )
    biexp = compute_biexponential_diagnostic(
        hfacf_detail=hfacf_detail,
        temperature_k=temperature_k,
        volume_a3=volume_a3,
        fit_start_ps=biexp_fit_start_ps,
        fit_end_ps=biexp_fit_end_ps,
        smooth_window=biexp_smooth_window,
        output_png=biexp_png,
    )

    payload = {
        "repdir": str(repdir),
        "temperature_K": temperature_k,
        "volume_A3": volume_a3,
        "gk_rows": int(len(gk_df)),
        "hfacf_last_block_rows": int(len(hfacf_detail["lag_ps"])),
        "einstein": einstein,
        "biexp": biexp,
        "outputs": {
            "einstein_plot": str(einstein_png),
            "biexp_plot": str(biexp_png),
        },
    }

    with open(outdir / "summary.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    return payload


def summarize_batch(results: Iterable[dict[str, object]], outpath: Path) -> None:
    rows = []
    for item in results:
        rows.append(
            {
                "repdir": item["repdir"],
                "temperature_K": item["temperature_K"],
                "volume_A3": item["volume_A3"],
                "einstein_kappa_W_mK": item["einstein"]["kappa_W_mK"],
                "einstein_r2": item["einstein"]["r2"],
                "biexp_kappa_W_mK": item["biexp"]["kappa_W_mK"],
                "biexp_r2": item["biexp"]["r2"],
            }
        )

    df = pd.DataFrame(rows)
    summary = {
        "replicas": rows,
        "aggregate": {
            "einstein_mean_W_mK": float(df["einstein_kappa_W_mK"].mean()),
            "einstein_std_W_mK": float(df["einstein_kappa_W_mK"].std(ddof=1)) if len(df) > 1 else 0.0,
            "biexp_mean_W_mK": float(df["biexp_kappa_W_mK"].mean()),
            "biexp_std_W_mK": float(df["biexp_kappa_W_mK"].std(ddof=1)) if len(df) > 1 else 0.0,
            "mean_einstein_r2": float(df["einstein_r2"].mean()),
            "mean_biexp_r2": float(df["biexp_r2"].mean()),
        },
    }
    with open(outpath, "w") as handle:
        json.dump(summary, handle, indent=2)


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standalone EMD thermal-conductivity diagnostics")
    parser.add_argument("--case-dir", type=Path, default=None, help="Case directory containing replica subdirectories")
    parser.add_argument("--repdir", type=Path, action="append", default=None, help="Replica directory with gk_data.dat and hfacf.dat")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for diagnostic outputs")
    parser.add_argument("--timestep-fs", type=float, required=True, help="LAMMPS MD timestep in fs")
    parser.add_argument("--temperature-k", type=float, default=None, help="Override average temperature from gk_data.dat")

    parser.add_argument("--einstein-max-lag-ps", type=nonnegative_float, default=20.0)
    parser.add_argument("--einstein-fit-start-ps", type=nonnegative_float, default=2.0)
    parser.add_argument("--einstein-fit-end-ps", type=nonnegative_float, default=10.0)
    parser.add_argument(
        "--no-subtract-mean-flux",
        action="store_true",
        help="Disable subtraction of the average heat-flux density before integration",
    )

    parser.add_argument("--biexp-fit-start-ps", type=nonnegative_float, default=0.0)
    parser.add_argument("--biexp-fit-end-ps", type=nonnegative_float, default=5.0)
    parser.add_argument("--biexp-smooth-window", type=positive_int, default=1)

    return parser


def main() -> None:
    parser = build_argparser()
    args = parser.parse_args()

    if args.case_dir is None and not args.repdir:
        parser.error("Provide either --case-dir or at least one --repdir")

    if args.einstein_fit_end_ps <= args.einstein_fit_start_ps:
        parser.error("Einstein fit end must be greater than fit start")
    if args.biexp_fit_end_ps <= args.biexp_fit_start_ps:
        parser.error("Biexp fit end must be greater than fit start")

    if args.case_dir is not None:
        repdirs = autodetect_repdirs(args.case_dir.resolve())
        case_dir = args.case_dir.resolve()
    else:
        repdirs = [path.resolve() for path in args.repdir]
        case_dir = None

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for repdir in repdirs:
        label = rep_label(repdir, case_dir=case_dir)
        rep_outdir = output_dir / label
        result = process_replica(
            repdir=repdir,
            outdir=rep_outdir,
            timestep_fs=float(args.timestep_fs),
            temperature_k_arg=args.temperature_k,
            einstein_max_lag_ps=float(args.einstein_max_lag_ps),
            einstein_fit_start_ps=float(args.einstein_fit_start_ps),
            einstein_fit_end_ps=float(args.einstein_fit_end_ps),
            subtract_mean_flux=not args.no_subtract_mean_flux,
            biexp_fit_start_ps=float(args.biexp_fit_start_ps),
            biexp_fit_end_ps=float(args.biexp_fit_end_ps),
            biexp_smooth_window=int(args.biexp_smooth_window),
        )
        results.append(result)

        print(f"[{label}] Einstein kappa = {result['einstein']['kappa_W_mK']:.6f} W/mK | R2 = {result['einstein']['r2']:.4f}")
        print(f"[{label}] Biexp    kappa = {result['biexp']['kappa_W_mK']:.6f} W/mK | R2 = {result['biexp']['r2']:.4f}")

    summarize_batch(results, output_dir / "batch_summary.json")
    print(f"Wrote outputs to {output_dir}")


if __name__ == "__main__":
    main()
