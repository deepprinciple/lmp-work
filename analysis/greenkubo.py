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


# ──────────────────────────────────────────────────────────────────────────────
# 文件解析
# ──────────────────────────────────────────────────────────────────────────────

def parse_ave_correlate(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
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
        (time_delta_fs, corr_matrix)
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

    return time_delta, corr


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
    time_fs, corr = parse_ave_correlate(corr_stress_file)
    if corr.shape[1] < 3:
        raise ValueError(
            f"应力 ACF 文件应包含 3 列（pxy, pxz, pyz），实际 {corr.shape[1]} 列"
        )

    time_ps = time_fs / 1000.0
    t_max_fs = t_max_ps * 1000.0 if t_max_ps is not None else None

    # 对三分量分别积分，取平均
    integrals = []
    running_list = []
    for col in range(3):
        run, _ = integrate_acf(time_fs, corr[:, col], t_max=t_max_fs)
        # 转换为 Pa·s
        run_si = run * VISC_CONV * volume / T
        running_list.append(run_si)
        integrals.append(run_si)

    running_mean = np.mean(running_list, axis=0)  # shape: (nrepeat,)

    # 截断到 t_max
    time_plot = time_ps.copy()
    if t_max_ps is not None:
        mask = time_plot <= t_max_ps
        time_plot = time_plot[mask]
        running_mean = running_mean[: len(time_plot)]

    # 可选平滑
    smooth_window = max(1, min(int(smooth_window), len(running_mean)))
    if smooth_window > 1:
        running_smooth = running_average(running_mean, smooth_window)
        time_smooth = time_plot[smooth_window - 1 :]
    else:
        running_smooth = running_mean
        time_smooth = time_plot

    # 平台区取均值
    p_end = plateau_end_ps if plateau_end_ps is not None else time_smooth[-1] * 0.8
    plateau_mask = (time_smooth >= plateau_start_ps) & (time_smooth <= p_end)
    plateau_vals = running_smooth[plateau_mask]

    if len(plateau_vals) == 0:
        # fallback: 后 20% 数据
        n = len(running_smooth)
        plateau_vals = running_smooth[int(n * 0.8) :]

    eta = float(plateau_vals.mean())
    eta_std = float(plateau_vals.std())

    # 各分量
    eta_components = []
    for run_si in running_list:
        run_arr = run_si[: len(time_plot)]
        if smooth_window > 1:
            run_arr = running_average(run_arr, smooth_window)
        plateau_c = run_arr[plateau_mask] if plateau_mask.sum() > 0 else run_arr[int(len(run_arr) * 0.8):]
        eta_components.append(float(plateau_c.mean()))

    return {
        "eta_Pa_s": eta,
        "eta_cP": eta * 1e3,          # 1 Pa·s = 1000 cP
        "eta_components_Pa_s": eta_components,
        "eta_components_cP": [x * 1e3 for x in eta_components],
        "plateau_mean_Pa_s": eta,
        "plateau_std_Pa_s": eta_std,
        "plateau_start_ps": plateau_start_ps,
        "plateau_end_ps": float(p_end),
        "running_integral_Pa_s": running_smooth,
        "time_ps": time_smooth,
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
    time_fs, corr = parse_ave_correlate(corr_flux_file)
    if corr.shape[1] < 3:
        raise ValueError(
            f"热流 ACF 文件应包含 3 列（Jx, Jy, Jz），实际 {corr.shape[1]} 列"
        )

    time_ps = time_fs / 1000.0
    t_max_fs = t_max_ps * 1000.0 if t_max_ps is not None else None

    running_list = []
    for col in range(3):
        run, _ = integrate_acf(time_fs, corr[:, col], t_max=t_max_fs)
        # 转换为 W/(m·K)
        run_si = run * KAPPA_CONV * volume / (T ** 2)
        running_list.append(run_si)

    running_mean = np.mean(running_list, axis=0)

    time_plot = time_ps.copy()
    if t_max_ps is not None:
        mask = time_plot <= t_max_ps
        time_plot = time_plot[mask]
        running_mean = running_mean[: len(time_plot)]

    smooth_window = max(1, min(int(smooth_window), len(running_mean)))
    if smooth_window > 1:
        running_smooth = running_average(running_mean, smooth_window)
        time_smooth = time_plot[smooth_window - 1 :]
    else:
        running_smooth = running_mean
        time_smooth = time_plot

    p_end = plateau_end_ps if plateau_end_ps is not None else time_smooth[-1] * 0.8
    plateau_mask = (time_smooth >= plateau_start_ps) & (time_smooth <= p_end)
    plateau_vals = running_smooth[plateau_mask]

    if len(plateau_vals) == 0:
        n = len(running_smooth)
        plateau_vals = running_smooth[int(n * 0.8):]

    kappa = float(plateau_vals.mean())
    kappa_std = float(plateau_vals.std())

    kappa_components = []
    for run_si in running_list:
        run_arr = run_si[: len(time_plot)]
        if smooth_window > 1:
            run_arr = running_average(run_arr, smooth_window)
        plateau_c = run_arr[plateau_mask] if plateau_mask.sum() > 0 else run_arr[int(len(run_arr) * 0.8):]
        kappa_components.append(float(plateau_c.mean()))

    return {
        "kappa_W_mK": kappa,
        "kappa_components_W_mK": kappa_components,
        "plateau_mean_W_mK": kappa,
        "plateau_std_W_mK": kappa_std,
        "plateau_start_ps": plateau_start_ps,
        "plateau_end_ps": float(p_end),
        "running_integral_W_mK": running_smooth,
        "time_ps": time_smooth,
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
        smooth_window: int = 50,
        n_blocks: int = 5,
        save_summary: bool = True,
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
            smooth_window: running integral 滑动平均窗口
            n_blocks: 平台区分块数，用于统计误差
            save_summary: 是否保存文本摘要

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
            smooth_window=smooth_window,
        )

        # 用 block averaging 估计平台误差
        plateau_mask = (visc["time_ps"] >= visc_plateau_start_ps)
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
            smooth_window=smooth_window,
        )

        plateau_mask_k = (kappa["time_ps"] >= kappa_plateau_start_ps)
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
        }

        if save_summary:
            self._save_summary(results)

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
