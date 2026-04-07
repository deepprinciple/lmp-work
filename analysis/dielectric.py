"""
介电常数计算模块

从 GreenKuboRunner 产生的时间序列数据（gk_data.dat）计算静态介电常数。

方法：偶极矩涨落公式（Kirkwood-Fröhlich）
    ε = 1 + <ΔM²> / (3 ε₀ V kB T)

其中：
  <ΔM²> = <M²> - <M>²  （总偶极矩均方涨落，单位 e²·Å²）
  V     — 盒子平均体积（Å³）
  T     — 温度（K）

LAMMPS real 单位：
  - 偶极矩: e·Å（elementary charge × Angstrom）
  - 体积: Å³
"""
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple

from ..utils.constants import (
    BOLTZMANN_J_K,
    ELEMENTARY_CHARGE,
    VACUUM_PERMITTIVITY,
    ANGSTROM_TO_METER,
)
from .statistics import block_average, running_average, equilibration_check


# ──────────────────────────────────────────────────────────────────────────────
# 单位转换系数
# ──────────────────────────────────────────────────────────────────────────────

# ε - 1 = <ΔM²>[e²·Å²] × DIEL_CONV / (V[Å³] × T[K])
#
# (e·Å)² → C²·m²: × (ELEMENTARY_CHARGE × ANGSTROM_TO_METER)²
# ε₀·kB  → F/m × J/K = C²/(N·m) × J/K
# V·T 前置因子: Å³ → m³: × ANGSTROM_TO_METER³
#
# DIEL_CONV = (e·Å in SI)² / (3 × ε₀ × kB × ANGSTROM_TO_METER³)
#           = (ELEMENTARY_CHARGE × ANGSTROM_TO_METER)²
#             / (3 × VACUUM_PERMITTIVITY × BOLTZMANN_J_K × ANGSTROM_TO_METER³)

_e_ang_SI = ELEMENTARY_CHARGE * ANGSTROM_TO_METER      # C·m per e·Å
DIEL_CONV = _e_ang_SI ** 2 / (
    3.0 * VACUUM_PERMITTIVITY * BOLTZMANN_J_K * ANGSTROM_TO_METER ** 3
)
# DIEL_CONV ≈ 7.0e5  Å³·K / (e²·Å²)
# 即: ε - 1 = <ΔM²> × DIEL_CONV / (V × T)


# ──────────────────────────────────────────────────────────────────────────────
# 数据解析
# ──────────────────────────────────────────────────────────────────────────────

def parse_gk_data(
    gk_data_file: Path,
) -> Dict[str, np.ndarray]:
    """
    解析 GreenKuboRunner.run_greenkubo_sampling() 产生的时间序列文件

    文件列格式（以 # 开头的行为注释）：
        step  T  V  pxy  pxz  pyz  Jx  Jy  Jz  Mx  My  Mz

    Args:
        gk_data_file: gk_data.dat 文件路径

    Returns:
        字典，键为列名，值为 np.ndarray
    """
    data = np.loadtxt(gk_data_file, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 12:
        raise ValueError(
            f"期望 12 列（step T V pxy pxz pyz Jx Jy Jz Mx My Mz），"
            f"实际 {data.shape[1]} 列"
        )

    return {
        "step": data[:, 0],
        "T":    data[:, 1],
        "V":    data[:, 2],
        "pxy":  data[:, 3],
        "pxz":  data[:, 4],
        "pyz":  data[:, 5],
        "Jx":   data[:, 6],
        "Jy":   data[:, 7],
        "Jz":   data[:, 8],
        "Mx":   data[:, 9],
        "My":   data[:, 10],
        "Mz":   data[:, 11],
    }


# ──────────────────────────────────────────────────────────────────────────────
# 介电常数计算
# ──────────────────────────────────────────────────────────────────────────────

def compute_dielectric_constant(
    gk_data_file: Path,
    T: Optional[float] = None,
    discard_fraction: float = 0.0,
    n_blocks: int = 10,
) -> Dict:
    """
    从时间序列计算静态介电常数

    ε = 1 + <ΔM²> / (3 ε₀ V kB T)

    使用偶极矩涨落：<ΔM²> = <Mx²+My²+Mz²> - (<Mx>²+<My>²+<Mz>²)

    Args:
        gk_data_file: gk_data.dat 文件
        T: 温度 (K)，None 则从数据中取时间平均
        discard_fraction: 丢弃开头的比例（预平衡数据），0.0–0.5
        n_blocks: 误差估计的块数

    Returns:
        结果字典：
          epsilon        — 静态介电常数
          epsilon_sem    — 块标准误差
          M_rms          — 均方根偶极矩 (e·Å)
          delta_M2       — <ΔM²> (e²·Å²)
          volume         — 平均体积 (Å³)
          temperature    — 使用的温度 (K)
    """
    series = parse_gk_data(gk_data_file)
    n = len(series["step"])

    # 丢弃预平衡部分
    start = int(n * discard_fraction)
    if start > 0:
        for key in series:
            series[key] = series[key][start:]

    Mx = series["Mx"]
    My = series["My"]
    Mz = series["Mz"]
    V_arr = series["V"]
    T_arr = series["T"]

    # 使用数据中的平均温度和体积（除非用户指定）
    V_mean = float(V_arr.mean())
    T_use = T if T is not None else float(T_arr.mean())

    # 偶极矩涨落
    M2 = Mx ** 2 + My ** 2 + Mz ** 2          # M²(t)
    M2_mean = float(M2.mean())
    Mx_mean = float(Mx.mean())
    My_mean = float(My.mean())
    Mz_mean = float(Mz.mean())
    delta_M2 = M2_mean - (Mx_mean ** 2 + My_mean ** 2 + Mz_mean ** 2)

    epsilon = 1.0 + delta_M2 * DIEL_CONV / (V_mean * T_use)

    # 误差估计：对 M²(t) 做 block average
    try:
        _, M2_sem, _ = block_average(M2, n_blocks=n_blocks)
    except ValueError:
        M2_sem = float(M2.std()) / np.sqrt(len(M2))

    # ε 相对于 <M²> 的不确定度传播
    # ε ≈ 1 + <M²> × DIEL_CONV / (V·T)  （忽略 <M>² 的不确定度）
    epsilon_sem = M2_sem * DIEL_CONV / (V_mean * T_use)

    M_rms = float(np.sqrt(M2_mean))

    return {
        "epsilon": epsilon,
        "epsilon_sem": epsilon_sem,
        "delta_M2_e2A2": delta_M2,
        "M2_mean_e2A2": M2_mean,
        "M_rms_eA": M_rms,
        "Mx_mean_eA": Mx_mean,
        "My_mean_eA": My_mean,
        "Mz_mean_eA": Mz_mean,
        "volume_A3": V_mean,
        "temperature": T_use,
        "n_samples": len(Mx),
        "discard_fraction": discard_fraction,
    }


# ──────────────────────────────────────────────────────────────────────────────
# 运行平均（收敛性检验）
# ──────────────────────────────────────────────────────────────────────────────

def dielectric_convergence(
    gk_data_file: Path,
    T: Optional[float] = None,
    window_ps: Optional[float] = None,
    dt_sample_fs: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算介电常数关于采样时间的收敛曲线

    逐步累积数据，每一步用所有已采集数据重新估计 ε，
    从而判断模拟是否已收敛。

    Args:
        gk_data_file: gk_data.dat 文件
        T: 温度 (K)，None 则从数据中取平均
        window_ps: 每个收敛点之间的步长（ps），None 则取总点数 / 100
        dt_sample_fs: 采样间隔（fs），用于计算时间轴

    Returns:
        (time_ps, epsilon_vs_time)
    """
    series = parse_gk_data(gk_data_file)
    Mx = series["Mx"]
    My = series["My"]
    Mz = series["Mz"]
    V_arr = series["V"]
    T_arr = series["T"]
    n = len(Mx)

    T_use = T if T is not None else float(T_arr.mean())

    step_size = max(1, int(window_ps * 1000.0 / dt_sample_fs)) if window_ps else max(1, n // 100)
    checkpoints = list(range(step_size, n + 1, step_size))
    if not checkpoints or checkpoints[-1] < n:
        checkpoints.append(n)

    eps_list = []
    t_list = []

    for end in checkpoints:
        mx = Mx[:end]
        my = My[:end]
        mz = Mz[:end]
        V_m = float(V_arr[:end].mean())

        M2 = (mx ** 2 + my ** 2 + mz ** 2).mean()
        dM2 = float(M2) - (float(mx.mean()) ** 2 + float(my.mean()) ** 2 + float(mz.mean()) ** 2)
        eps = 1.0 + dM2 * DIEL_CONV / (V_m * T_use)
        eps_list.append(eps)
        t_list.append(end * dt_sample_fs / 1000.0)  # ps

    return np.array(t_list), np.array(eps_list)


# ──────────────────────────────────────────────────────────────────────────────
# 高层接口
# ──────────────────────────────────────────────────────────────────────────────

class DielectricAnalyzer:
    """
    介电常数分析器

    Example::

        analyzer = DielectricAnalyzer(workdir=Path("./run1"))
        results = analyzer.analyze(
            gk_data_file=Path("greenkubo_data.dat"),
            T=298.0,
        )
        print(f"ε = {results['epsilon']:.2f} ± {results['epsilon_sem']:.2f}")
    """

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)

    def analyze(
        self,
        gk_data_file: Path,
        T: Optional[float] = None,
        discard_fraction: float = 0.1,
        n_blocks: int = 10,
        dt_sample_fs: float = 10.0,
        save_summary: bool = True,
    ) -> Dict:
        """
        完整介电常数分析

        Args:
            gk_data_file: gk_data.dat 文件路径
            T: 温度 (K)，None 则从数据自动获取
            discard_fraction: 丢弃开头比例（去除预平衡）
            n_blocks: 误差块数
            dt_sample_fs: 采样时间间隔（fs），用于收敛曲线
            save_summary: 是否保存摘要文件

        Returns:
            结果字典
        """
        print("=" * 70)
        print("介电常数分析")
        print("=" * 70)

        results = compute_dielectric_constant(
            gk_data_file=gk_data_file,
            T=T,
            discard_fraction=discard_fraction,
            n_blocks=n_blocks,
        )

        print(f"  温度: {results['temperature']:.2f} K")
        print(f"  体积: {results['volume_A3']:.2f} Å³")
        print(f"  样本数: {results['n_samples']}")
        print(f"\n  ε = {results['epsilon']:.3f} ± {results['epsilon_sem']:.3f}")
        print(f"  <ΔM²> = {results['delta_M2_e2A2']:.4f} e²·Å²")
        print(f"  M_rms = {results['M_rms_eA']:.4f} e·Å")

        # 收敛曲线
        t_conv, eps_conv = dielectric_convergence(
            gk_data_file, T=results["temperature"], dt_sample_fs=dt_sample_fs
        )
        results["convergence_time_ps"] = t_conv
        results["convergence_epsilon"] = eps_conv

        if save_summary:
            self._save_summary(results)

        print("\n  完成！")
        return results

    def _save_summary(self, results: Dict):
        summary_file = self.workdir / "dielectric_results.txt"
        with open(summary_file, "w") as f:
            f.write("=" * 70 + "\n")
            f.write("介电常数分析结果\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"温度: {results['temperature']:.2f} K\n")
            f.write(f"体积: {results['volume_A3']:.2f} Å³\n")
            f.write(f"样本数: {results['n_samples']}\n")
            f.write(f"丢弃比例: {results['discard_fraction']*100:.0f}%\n\n")

            f.write("结果:\n")
            f.write(f"  静态介电常数 ε = {results['epsilon']:.3f} ± {results['epsilon_sem']:.3f}\n")
            f.write(f"  <ΔM²>         = {results['delta_M2_e2A2']:.6f} e²·Å²\n")
            f.write(f"  均方根偶极矩  = {results['M_rms_eA']:.4f} e·Å\n")
            f.write(f"  <Mx>          = {results['Mx_mean_eA']:.4f} e·Å\n")
            f.write(f"  <My>          = {results['My_mean_eA']:.4f} e·Å\n")
            f.write(f"  <Mz>          = {results['Mz_mean_eA']:.4f} e·Å\n")

        print(f"  结果摘要: {summary_file}")
