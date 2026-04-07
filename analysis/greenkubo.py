"""
Green-Kubo 分析模块

从 LAMMPS fix ave/correlate 输出文件计算：
  - 剪切粘度 η（via 应力张量自相关函数）
  - 热导率 λ（via 热流自相关函数）

单位制：LAMMPS real units
  - 压力: atm
  - 热流: kcal/(mol·Å²·fs)  (已除以体积)
  - 时间: fs
  - 体积: Å³
  - 温度: K
"""
import json
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional
from collections import Counter
import warnings

from ..utils.constants import (
    BOLTZMANN_J_K,
    ATM_TO_PASCAL,
    FEMTOSECOND_TO_SECOND,
    ANGSTROM_TO_METER,
    KCAL_MOL_TO_JOULE,
    AVOGADRO,
)
from .statistics import integrate_acf, block_average, running_average


# ──────────────────────────────────────────────────────────────────────────────
# 单位转换系数
# ──────────────────────────────────────────────────────────────────────────────

# 粘度: η [Pa·s] = V[Å³] / (kB[J/K] × T[K]) × ∫C_P[atm²·fs]
#
# ∫C_P [atm²·fs] → SI [Pa²·s]:  × ATM_TO_PASCAL² × FEMTOSECOND_TO_SECOND
# V [Å³] → m³:                  × ANGSTROM_TO_METER³
# ──────────────────────────────────────────────────────────────────
# η_CONV = ANGSTROM_TO_METER³ × ATM_TO_PASCAL² × FEMTOSECOND_TO_SECOND / BOLTZMANN_J_K
VISC_CONV = (
    ANGSTROM_TO_METER ** 3
    * ATM_TO_PASCAL ** 2
    * FEMTOSECOND_TO_SECOND
    / BOLTZMANN_J_K
)

# 热导率: λ [W/(m·K)] = V[Å³] / (kB[J/K] × T²[K²]) × ∫C_J[kcal²/(mol²·Å⁴·fs)]
#
# heat flux unit: kcal/(mol·Å²·fs) → SI: W/m²
#   1 kcal/(mol·Å²·fs) = KCAL_MOL_TO_JOULE / AVOGADRO / ANGSTROM_TO_METER² / FEMTOSECOND_TO_SECOND
HF_TO_SI = KCAL_MOL_TO_JOULE / AVOGADRO / ANGSTROM_TO_METER ** 2 / FEMTOSECOND_TO_SECOND

# ∫C_J [kcal²/(mol²·Å⁴·fs)] → SI [W²/m⁴·s]:
#   × HF_TO_SI² × FEMTOSECOND_TO_SECOND
# V [Å³] → m³: × ANGSTROM_TO_METER³
# λ_CONV = ANGSTROM_TO_METER³ × HF_TO_SI² × FEMTOSECOND_TO_SECOND / BOLTZMANN_J_K
KAPPA_CONV = (
    ANGSTROM_TO_METER ** 3
    * HF_TO_SI ** 2
    * FEMTOSECOND_TO_SECOND
    / BOLTZMANN_J_K
)


def _relative_metric(value: float, reference: float, fallback_scale: float) -> float:
    denom = max(abs(reference), fallback_scale, 1.0e-12)
    return abs(value) / denom


# ──────────────────────────────────────────────────────────────────────────────
# 文件解析
# ──────────────────────────────────────────────────────────────────────────────

def parse_ave_correlate(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    detail = parse_ave_correlate_detail(filepath)
    return detail["time_fs"], detail["corr"]


def parse_ave_correlate_detail(filepath: Path) -> Dict[str, np.ndarray | int | float | bool]:
    """
    解析 LAMMPS fix ave/correlate 输出文件。

    兼容两类常见输出：
    1. 单块输出:
       [Index, TimeDelta, Ncount, corr1, corr2, ...]
       或
       [TimeDelta, Ncount, corr1, corr2, ...]
    2. ``ave running`` 多块输出:
       多个 block 依次写到同一文件中，每个 block 之前会有一行
       ``<timestep> <nrepeat>``，随后是 nrepeat 行相关函数数据。

    对于多块输出，分析时应使用“最后一个完整 block”，否则会把不同时间点的
    running average 拼接在一起，甚至把 ``Ncount`` 误读为相关函数，导致
    热导率/粘度积分出现天文数字。

    Args:
        filepath: ave/correlate 输出文件路径

    Returns:
        包含 time_fs / corr / ncount 等解析细节的字典
    """
    numeric_rows = []
    row_lengths = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                vals = [float(x) for x in line.split()]
            except ValueError:
                continue
            numeric_rows.append(vals)
            row_lengths.append(len(vals))

    if not numeric_rows:
        raise ValueError(f"文件无数据: {filepath}")

    count_by_len = Counter(row_lengths)
    data_len = max(count_by_len.items(), key=lambda kv: (kv[1], kv[0]))[0]
    if data_len < 5:
        raise ValueError(
            f"ave/correlate 数据列数异常: {data_len}（期望至少5列）"
        )

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
        filtered = [r for r in numeric_rows if len(r) == data_len]
        if not filtered:
            raise ValueError(f"未找到有效相关函数数据: {filepath}")
        blocks = [np.asarray(filtered, dtype=float)]

    if dropped_rows > 0:
        warnings.warn(
            f"{filepath} 中有 {dropped_rows} 行无法识别，已忽略",
            RuntimeWarning,
        )

    arr = blocks[-1]

    # 兼容:
    #   [Index, TimeDelta, Ncount, corr1, corr2, ...]
    #   [TimeDelta, Ncount, corr1, corr2, ...]
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
        time_delta = arr[:, 0]
        corr = arr[:, 2:]

    if corr.shape[1] < 3:
        raise ValueError(f"相关函数列不足: {corr.shape[1]}（文件: {filepath}）")

    if arr.shape[1] >= 6 and looks_index:
        ncount = arr[:, 2]
    else:
        ncount = arr[:, 1]

    valid_mask = np.isfinite(time_delta) & np.all(np.isfinite(corr), axis=1)
    positive_count_mask = valid_mask & (ncount > 0)
    if np.count_nonzero(positive_count_mask) >= 3:
        data_mask = positive_count_mask
    else:
        data_mask = valid_mask

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
    start_idx = max(0, min(n_points - 1, n_points - max(min_points, int(n_points * 0.2))))
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
        raise ValueError("无法从空时间序列中选择分析窗口")
    if len(time_ps) != len(running_integral):
        time_ps = time_ps[:n_points]
        running_integral = running_integral[:n_points]

    analysis_end_ps = float(time_ps[-1]) if plateau_end_ps is None else min(float(plateau_end_ps), float(time_ps[-1]))
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

        stable = rel_std <= rel_std_tol and drift_fraction <= drift_tol and slope_fraction <= drift_tol
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
        slope = float(np.polyfit(time_window, window, 1)[0]) if len(window) >= 2 and duration_ps > 0.0 else 0.0
        slope_fraction = _relative_metric(slope * duration_ps, mean_val, scale * 0.05) if duration_ps > 0.0 else 0.0
        stable = rel_std <= rel_std_tol and drift_fraction <= drift_tol and slope_fraction <= drift_tol
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
        raise ValueError(f"{corr_file} 需要至少 3 个相关函数分量，实际 {corr.shape[1]}")

    time_ps = time_fs / 1000.0
    effective_t_max_ps = _estimate_effective_cutoff_ps(
        time_ps,
        corr,
        t_max_ps,
        min_ps_before_cutoff=min_ps_before_cutoff,
    )
    t_max_fs = effective_t_max_ps * 1000.0 if effective_t_max_ps is not None else None

    running_components = []
    for col in range(3):
        run, _ = integrate_acf(time_fs, corr[:, col], t_max=t_max_fs)
        running_components.append(run * converter_scale * component_factor)

    full_time_mask = time_ps <= effective_t_max_ps + 1.0e-12
    time_plot = time_ps[full_time_mask]
    running_components_arr = np.asarray([comp[: len(time_plot)] for comp in running_components], dtype=float)
    running_mean = np.mean(running_components_arr, axis=0)

    smooth_window = max(1, min(int(smooth_window), len(running_mean)))
    if smooth_window > 1:
        running_smooth = running_average(running_mean, smooth_window)
        time_smooth = time_plot[smooth_window - 1 :]
        component_smooth = np.asarray([running_average(comp, smooth_window) for comp in running_components_arr], dtype=float)
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

    component_values = [float(comp[window_mask].mean()) for comp in component_smooth]
    diagnostics = {
        "corr_file": str(corr_file),
        "max_lag_ps_available": float(time_ps[-1]) if len(time_ps) > 0 else 0.0,
        "effective_t_max_ps": float(effective_t_max_ps),
        "n_lags_used": int(len(time_plot)),
        "ncount_nonzero_fraction": float(np.count_nonzero(ncount > 0) / len(ncount)) if len(ncount) > 0 else 0.0,
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
        "component_rel_std": _relative_metric(float(np.std(component_values)), float(np.mean(component_values)), 1.0e-12),
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


# ──────────────────────────────────────────────────────────────────────────────
# 粘度
# ──────────────────────────────────────────────────────────────────────────────

def compute_viscosity(
    corr_stress_file: Path,
    T: float,
    volume: float,
    t_max_ps: Optional[float] = None,
    plateau_start_ps: float = 5.0,
    plateau_end_ps: Optional[float] = None,
    smooth_window: int = 1,
) -> Dict:
    """
    从应力自相关函数计算剪切粘度

    Green-Kubo 公式：
        η = V / (kB·T) × ∫₀^∞ <Pij(0) Pij(t)> dt

    取三个非对角分量（pxy, pxz, pyz）的平均值。

    Args:
        corr_stress_file: fix ave/correlate 应力输出文件
                          （包含 pxy*pxy, pxz*pxz, pyz*pyz 三列）
        T: 温度 (K)
        volume: 盒子体积 (Å³)
        t_max_ps: ACF 积分截止时间 (ps)，None = 全程
        plateau_start_ps: 积分平台区起始时间 (ps)，用于取平均
        plateau_end_ps: 积分平台区终止时间 (ps)，None = 全程末尾 20%
        smooth_window: 对 running integral 的滑动平均窗口（奇数步）

    Returns:
        结果字典：
          eta_Pa_s   — 粘度 [Pa·s]
          eta_cP     — 粘度 [cP = mPa·s]
          eta_components — 三分量各自的粘度
          running_integral — 每个时间点的累积积分（已转换到 Pa·s）
          time_ps    — 对应时间轴 (ps)
          plateau_mean, plateau_std — 平台区均值和标准差
    """
    analysis = _analyze_transport_running_integral(
        corr_file=corr_stress_file,
        T=T,
        volume=volume,
        t_max_ps=t_max_ps,
        plateau_start_ps=plateau_start_ps,
        plateau_end_ps=plateau_end_ps,
        smooth_window=smooth_window,
        rel_std_tol=0.12,
        drift_tol=0.12,
        min_ps_before_cutoff=max(0.2, plateau_start_ps * 0.5),
        converter_scale=VISC_CONV,
        component_factor=volume / T,
        component_unit="Pa·s",
        result_key="eta_Pa_s",
    )
    eta = float(analysis["eta_Pa_s"])
    eta_std = float(analysis["plateau_std"])
    eta_components = list(analysis["component_values"])

    return {
        "eta_Pa_s": eta,
        "eta_cP": eta * 1e3,          # 1 Pa·s = 1000 cP
        "eta_components_Pa_s": eta_components,
        "eta_components_cP": [x * 1e3 for x in eta_components],
        "plateau_mean_Pa_s": eta,
        "plateau_std_Pa_s": eta_std,
        "plateau_start_ps": float(analysis["window_start_ps"]),
        "plateau_end_ps": float(analysis["window_end_ps"]),
        "running_integral_Pa_s": analysis["running_integral"],
        "time_ps": analysis["time_ps"],
        "analysis_window_mask": analysis["window_mask"],
        "effective_t_max_ps": float(analysis["effective_t_max_ps"]),
        "diagnostics": analysis["diagnostics"],
        "temperature": T,
        "volume_A3": volume,
    }


# ──────────────────────────────────────────────────────────────────────────────
# 热导率
# ──────────────────────────────────────────────────────────────────────────────

def compute_thermal_conductivity(
    corr_flux_file: Path,
    T: float,
    volume: float,
    t_max_ps: Optional[float] = None,
    plateau_start_ps: float = 2.0,
    plateau_end_ps: Optional[float] = None,
    smooth_window: int = 1,
) -> Dict:
    """
    从热流自相关函数计算热导率

    Green-Kubo 公式：
        λ = V / (kB·T²) × ∫₀^∞ <Jα(0) Jα(t)> dt

    取三分量（Jx, Jy, Jz）的平均值。

    注意：生产运行中 Jα = c_myFlux[α] / vol，已除以体积，
    因此积分公式前置因子为 V（而非 1/V²）。

    Args:
        corr_flux_file: fix ave/correlate 热流输出文件
                        （包含 Jx*Jx, Jy*Jy, Jz*Jz 三列）
        T: 温度 (K)
        volume: 盒子体积 (Å³)
        t_max_ps: ACF 积分截止时间 (ps)，None = 全程
        plateau_start_ps: 平台区起始时间 (ps)
        plateau_end_ps: 平台区终止时间 (ps)
        smooth_window: 滑动平均窗口

    Returns:
        结果字典：
          kappa_W_mK   — 热导率 [W/(m·K)]
          running_integral — 每个时间点的累积积分
          time_ps
          plateau_mean, plateau_std
    """
    analysis = _analyze_transport_running_integral(
        corr_file=corr_flux_file,
        T=T,
        volume=volume,
        t_max_ps=t_max_ps,
        plateau_start_ps=plateau_start_ps,
        plateau_end_ps=plateau_end_ps,
        smooth_window=smooth_window,
        rel_std_tol=0.20,
        drift_tol=0.20,
        min_ps_before_cutoff=max(0.3, plateau_start_ps * 0.5),
        converter_scale=KAPPA_CONV,
        component_factor=volume / (T ** 2),
        component_unit="W/(m·K)",
        result_key="kappa_W_mK",
    )
    kappa = float(analysis["kappa_W_mK"])
    kappa_std = float(analysis["plateau_std"])
    kappa_components = list(analysis["component_values"])

    return {
        "kappa_W_mK": kappa,
        "kappa_components_W_mK": kappa_components,
        "plateau_mean_W_mK": kappa,
        "plateau_std_W_mK": kappa_std,
        "plateau_start_ps": float(analysis["window_start_ps"]),
        "plateau_end_ps": float(analysis["window_end_ps"]),
        "running_integral_W_mK": analysis["running_integral"],
        "time_ps": analysis["time_ps"],
        "analysis_window_mask": analysis["window_mask"],
        "effective_t_max_ps": float(analysis["effective_t_max_ps"]),
        "diagnostics": analysis["diagnostics"],
        "temperature": T,
        "volume_A3": volume,
    }


# ──────────────────────────────────────────────────────────────────────────────
# 高层接口：同时计算两者
# ──────────────────────────────────────────────────────────────────────────────

class GreenKuboAnalyzer:
    """
    Green-Kubo 分析器

    从 GreenKuboRunner 产生的输出文件一次性计算
    粘度和热导率，附带误差估计。

    Example::

        analyzer = GreenKuboAnalyzer(workdir=Path("./run1"))
        results = analyzer.analyze(
            corr_stress_file=Path("greenkubo_corr_stress.dat"),
            corr_flux_file=Path("greenkubo_corr_flux.dat"),
            T=298.0,
            volume=18500.0,
        )
        print(f"η = {results['eta_cP']:.3f} cP")
        print(f"λ = {results['kappa_W_mK']:.4f} W/(m·K)")
    """

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)

    def analyze(
        self,
        corr_stress_file: Path,
        corr_flux_file: Path,
        T: float,
        volume: float,
        visc_t_max_ps: Optional[float] = 20.0,
        kappa_t_max_ps: Optional[float] = 20.0,
        visc_plateau_start_ps: float = 5.0,
        kappa_plateau_start_ps: float = 2.0,
        visc_plateau_end_ps: Optional[float] = None,
        kappa_plateau_end_ps: Optional[float] = None,
        smooth_window: int = 50,
        n_blocks: int = 5,
        save_summary: bool = True,
        save_diagnostics: bool = True,
    ) -> Dict:
        """
        完整 Green-Kubo 分析

        Args:
            corr_stress_file: 应力 ACF 文件
            corr_flux_file: 热流 ACF 文件
            T: 温度 (K)
            volume: 盒子体积 (Å³)
            visc_t_max_ps: 粘度积分截止 (ps)
            kappa_t_max_ps: 热导率积分截止 (ps)
            visc_plateau_start_ps: 粘度平台起始 (ps)
            kappa_plateau_start_ps: 热导率平台起始 (ps)
            visc_plateau_end_ps: 粘度平台终止 (ps)
            kappa_plateau_end_ps: 热导率平台终止 (ps)
            smooth_window: running integral 滑动平均窗口
            n_blocks: 平台区分块数，用于统计误差
            save_summary: 是否保存文本摘要
            save_diagnostics: 是否保存结构化诊断文件

        Returns:
            包含粘度和热导率结果的字典
        """
        print("=" * 70)
        print("Green-Kubo 分析")
        print("=" * 70)
        print(f"  温度: {T} K")
        print(f"  体积: {volume:.2f} Å³")

        # ── 粘度 ──────────────────────────────────────────────────────────────
        print("\n  计算剪切粘度...")
        visc = compute_viscosity(
            corr_stress_file=corr_stress_file,
            T=T,
            volume=volume,
            t_max_ps=visc_t_max_ps,
            plateau_start_ps=visc_plateau_start_ps,
            plateau_end_ps=visc_plateau_end_ps,
            smooth_window=smooth_window,
        )

        # 用 block averaging 估计平台误差
        plateau_mask = visc["analysis_window_mask"]
        plateau_data = visc["running_integral_Pa_s"][plateau_mask]
        try:
            _, visc_sem, visc_rel = block_average(plateau_data, n_blocks=n_blocks)
        except ValueError:
            visc_sem, visc_rel = visc["plateau_std_Pa_s"], 0.0

        print(f"  η = {visc['eta_cP']:.4f} ± {visc_sem*1e3:.4f} cP")
        print(f"  η = {visc['eta_Pa_s']:.4e} ± {visc_sem:.4e} Pa·s")

        # ── 热导率 ────────────────────────────────────────────────────────────
        print("\n  计算热导率...")
        kappa = compute_thermal_conductivity(
            corr_flux_file=corr_flux_file,
            T=T,
            volume=volume,
            t_max_ps=kappa_t_max_ps,
            plateau_start_ps=kappa_plateau_start_ps,
            plateau_end_ps=kappa_plateau_end_ps,
            smooth_window=smooth_window,
        )

        plateau_mask_k = kappa["analysis_window_mask"]
        plateau_data_k = kappa["running_integral_W_mK"][plateau_mask_k]
        try:
            _, kappa_sem, kappa_rel = block_average(plateau_data_k, n_blocks=n_blocks)
        except ValueError:
            kappa_sem, kappa_rel = kappa["plateau_std_W_mK"], 0.0

        print(f"  λ = {kappa['kappa_W_mK']:.4f} ± {kappa_sem:.4f} W/(m·K)")

        results = {
            # 粘度
            "eta_Pa_s": visc["eta_Pa_s"],
            "eta_cP": visc["eta_cP"],
            "eta_sem_Pa_s": visc_sem,
            "eta_sem_cP": visc_sem * 1e3,
            "eta_rel_error": visc_rel,
            "eta_components_cP": visc["eta_components_cP"],
            "visc_time_ps": visc["time_ps"],
            "visc_running_Pa_s": visc["running_integral_Pa_s"],
            # 热导率
            "kappa_W_mK": kappa["kappa_W_mK"],
            "kappa_sem_W_mK": kappa_sem,
            "kappa_rel_error": kappa_rel,
            "kappa_components_W_mK": kappa["kappa_components_W_mK"],
            "kappa_time_ps": kappa["time_ps"],
            "kappa_running_W_mK": kappa["running_integral_W_mK"],
            # 参数
            "temperature": T,
            "volume_A3": volume,
            "diagnostics": {
                "viscosity": visc["diagnostics"],
                "thermal_conductivity": kappa["diagnostics"],
            },
        }

        if save_summary:
            self._save_summary(results)
        if save_diagnostics:
            self._save_diagnostics(results)

        print("\n  完成！")
        return results

    def _save_summary(self, results: Dict):
        """保存文本摘要到工作目录"""
        summary_file = self.workdir / "greenkubo_results.txt"
        with open(summary_file, "w") as f:
            f.write("=" * 70 + "\n")
            f.write("Green-Kubo 分析结果\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"温度: {results['temperature']:.2f} K\n")
            f.write(f"体积: {results['volume_A3']:.2f} Å³\n\n")

            f.write("剪切粘度:\n")
            f.write(f"  η = {results['eta_cP']:.4f} ± {results['eta_sem_cP']:.4f} cP\n")
            f.write(f"  η = {results['eta_Pa_s']:.4e} ± {results['eta_sem_Pa_s']:.4e} Pa·s\n")
            f.write(f"  相对误差: {results['eta_rel_error']*100:.1f}%\n")
            comps = results["eta_components_cP"]
            f.write(f"  分量(pxy,pxz,pyz): {comps[0]:.4f}, {comps[1]:.4f}, {comps[2]:.4f} cP\n\n")

            f.write("热导率:\n")
            f.write(f"  λ = {results['kappa_W_mK']:.4f} ± {results['kappa_sem_W_mK']:.4f} W/(m·K)\n")
            f.write(f"  相对误差: {results['kappa_rel_error']*100:.1f}%\n")
            comps_k = results["kappa_components_W_mK"]
            f.write(f"  分量(Jx,Jy,Jz): {comps_k[0]:.4f}, {comps_k[1]:.4f}, {comps_k[2]:.4f} W/(m·K)\n")

        print(f"  结果摘要: {summary_file}")

    def _save_diagnostics(self, results: Dict):
        diagnostics_file = self.workdir / "greenkubo_diagnostics.json"
        with open(diagnostics_file, "w") as f:
            json.dump(results["diagnostics"], f, indent=2, ensure_ascii=False)
