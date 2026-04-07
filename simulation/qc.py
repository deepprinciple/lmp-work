"""
模拟质量控制模块（QC）

对各阶段模拟输出文件进行自动化质量检验：
  - NVT 平衡：温度收敛
  - NPT 平衡：温度收敛、密度收敛（可与实验值对比）
  - NVE 生产：能量守恒（漂移率检验）
  - NVT/NPT 生产：温度/密度稳定性
  - GK 采样：温度和体积稳定性

用法::

    from gk_workflow.simulation.qc import SimulationQC

    qc = SimulationQC()

    # 单阶段检查
    checks = qc.check_npt_equilibration(
        log_file=Path("work/npt_equilibration.log"),
        T_target=298.0,
        rho_exp=0.867,
    )
    SimulationQC.print_report("NPT平衡", checks)

    # 全流程一键检查
    qc.full_workflow_check(workdir=Path("work"), T_target=298.0, rho_exp=0.867)
"""
import re
import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ..analysis.statistics import block_average


# ──────────────────────────────────────────────────────────────────────────────
# 数据结构
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class QCCheck:
    """单项检查结果"""
    name: str           # 检查名称
    passed: bool        # 是否通过
    value: float        # 测量值（如相对偏差 %）
    threshold: float    # 判断阈值（同单位）
    unit: str           # 单位
    message: str        # 说明文字
    details: dict = field(default_factory=dict)  # 原始数值等


# ──────────────────────────────────────────────────────────────────────────────
# 文件解析
# ──────────────────────────────────────────────────────────────────────────────

# LAMMPS log 文件中的列名 → 统一名称映射
_LOG_COL_MAP = {
    "step": "step",
    "time": "time",
    "temp": "temp",
    "press": "press",
    "pe": "pe",
    "ke": "ke",
    "toteng": "etotal",   # LAMMPS log 把 etotal 显示为 TotEng
    "vol": "vol",
    "density": "density",
    "v_rho": "density",   # NPT 平衡中的密度变量
    "cpu": "cpu",
}


def parse_thermo_dat(filepath: Path) -> Dict[str, np.ndarray]:
    """
    解析 fix print 产生的热力学数据文件

    文件首行为注释格式列名，如::

        # step time(ps) temp(K) press(atm) pe(kcal/mol) ke(kcal/mol) etotal(kcal/mol) vol(A^3)

    列名中的单位括号会被自动去除。

    Args:
        filepath: thermo dat 文件路径

    Returns:
        {列名: np.ndarray}
    """
    with open(filepath) as f:
        header_line = f.readline()

    # 去掉 # 和单位括号，得到纯列名
    header_line = header_line.lstrip("# ").strip()
    raw_names = header_line.split()
    col_names = [re.sub(r"\(.*?\)", "", n).lower() for n in raw_names]

    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    return {
        name: data[:, i]
        for i, name in enumerate(col_names)
        if i < data.shape[1]
    }


def parse_lammps_log(log_file: Path) -> Dict[str, np.ndarray]:
    """
    解析 LAMMPS log 文件，提取所有热力学数据（拼接所有 run 段）

    LAMMPS log 中 thermo 块结构::

        Step Temp Press PE KE TotEng Vol Density CPU
        0    298.0  1.0  -12345  456  -11889  18500  0.867  0.0
        1000 298.5  0.9  ...
        Loop time of ...

    列名经过统一化映射（如 TotEng → etotal）。

    Args:
        log_file: LAMMPS log 文件路径

    Returns:
        {列名: np.ndarray}（多段 run 自动拼接）
    """
    all_data: Dict[str, List[float]] = {}
    col_names: List[str] = []
    in_thermo = False

    with open(log_file) as f:
        for line in f:
            line_stripped = line.strip()

            # 检测列名行：以 "Step" 开头
            if re.match(r"^Step\b", line_stripped, re.IGNORECASE):
                raw_cols = line_stripped.split()
                col_names = [_LOG_COL_MAP.get(c.lower(), c.lower()) for c in raw_cols]
                in_thermo = True
                continue

            # 检测 thermo 段结束
            if in_thermo and line_stripped.startswith("Loop time"):
                in_thermo = False
                continue

            # 读取数据行
            if in_thermo and line_stripped:
                parts = line_stripped.split()
                # 跳过非数字行（如 WARNING 等）
                try:
                    values = [float(v) for v in parts]
                except ValueError:
                    in_thermo = False
                    continue

                if len(values) == len(col_names):
                    for name, val in zip(col_names, values):
                        all_data.setdefault(name, []).append(val)

    return {k: np.array(v) for k, v in all_data.items()}


def parse_ave_time_density(filepath: Path) -> np.ndarray:
    """
    解析 fix ave/time 产生的密度文件（density_npt.dat）

    跳过所有 # 注释行，返回最后一列（密度值）的时间序列。

    Args:
        filepath: density_npt.dat 文件路径

    Returns:
        密度时间序列 (g/cm³)
    """
    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # 当前 workflow 写出的是:
    # TimeStep temp press density vol
    # 因而倒数第二列是密度，最后一列是体积。
    if data.shape[1] >= 4:
        return data[:, -2]
    return data[:, -1]


def parse_gk_data(filepath: Path) -> Dict[str, np.ndarray]:
    """
    解析 GK 采样时间序列文件

    列格式：step T V pxy pxz pyz Jx Jy Jz Mx My Mz

    Args:
        filepath: greenkubo_data.dat 文件路径

    Returns:
        {列名: np.ndarray}
    """
    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    col_names = ["step", "temp", "vol", "pxy", "pxz", "pyz",
                 "jx", "jy", "jz", "mx", "my", "mz"]
    return {
        name: data[:, i]
        for i, name in enumerate(col_names)
        if i < data.shape[1]
    }


# ──────────────────────────────────────────────────────────────────────────────
# 底层检验函数
# ──────────────────────────────────────────────────────────────────────────────

def _tail(arr: np.ndarray, frac: float = 0.2) -> np.ndarray:
    """取时间序列末尾 frac 比例的数据"""
    n = len(arr)
    if n <= 2:
        return arr
    start = max(0, int(n * (1 - frac)))
    if n - start < 2:
        start = n - 2
    return arr[start:]


def _check_mean_bias(
    series: np.ndarray,
    target: float,
    name: str,
    unit: str,
    tol: float,
    tail_frac: float = 0.2,
    n_blocks: int = 5,
) -> QCCheck:
    """均值偏差检验：|mean(tail) - target| / |target| < tol"""
    tail = _tail(series, tail_frac)
    tail = tail[np.isfinite(tail)]
    if len(tail) == 0:
        return QCCheck(
            name=f"{name}均值偏差",
            passed=False,
            value=float("inf"),
            threshold=tol * 100,
            unit="%",
            message="有效数据不足，无法计算均值偏差",
        )
    mean_val = float(tail.mean())
    rel_bias = abs(mean_val - target) / abs(target) if target != 0 else abs(mean_val)

    # 用 block average 估计均值误差
    try:
        _, sem, _ = block_average(tail, n_blocks=n_blocks)
    except ValueError:
        sem = float(tail.std()) / np.sqrt(len(tail))

    return QCCheck(
        name=f"{name}均值偏差",
        passed=rel_bias <= tol,
        value=rel_bias * 100,
        threshold=tol * 100,
        unit="%",
        message=(
            f"均值 {mean_val:.5g} {unit}，目标 {target:.5g} {unit}，"
            f"偏差 {rel_bias*100:.2f}%，SEM ±{sem:.4g}"
        ),
        details={"mean": mean_val, "target": target, "sem": sem},
    )


def _check_stability(
    series: np.ndarray,
    name: str,
    unit: str,
    tol: float,
    tail_frac: float = 0.2,
) -> QCCheck:
    """稳定性检验：std(tail) / |mean(tail)| < tol（相对标准差）"""
    tail = _tail(series, tail_frac)
    tail = tail[np.isfinite(tail)]
    if len(tail) == 0:
        return QCCheck(
            name=f"{name}稳定性（RSD）",
            passed=False,
            value=float("inf"),
            threshold=tol * 100,
            unit="%",
            message="有效数据不足，无法计算稳定性",
        )
    mean_val = float(tail.mean())
    std_val = float(tail.std())
    rel_std = std_val / abs(mean_val) if mean_val != 0 else std_val

    return QCCheck(
        name=f"{name}稳定性（RSD）",
        passed=rel_std <= tol,
        value=rel_std * 100,
        threshold=tol * 100,
        unit="%",
        message=(
            f"末尾 {tail_frac*100:.0f}% 相对标准差 {rel_std*100:.2f}%"
            f"（mean={mean_val:.5g}, std={std_val:.4g} {unit}）"
        ),
        details={"mean": mean_val, "std": std_val, "rel_std": rel_std},
    )


def _check_drift(
    series: np.ndarray,
    name: str,
    unit: str,
    tol: float,
    tail_frac: float = 0.2,
) -> QCCheck:
    """
    线性漂移检验：用线性拟合尾部数据，检测是否仍在趋势性变化中

    tol: |slope * duration / mean| 的阈值（相对量纲）
    """
    tail = _tail(series, tail_frac)
    tail = tail[np.isfinite(tail)]
    if len(tail) < 2:
        return QCCheck(
            name=f"{name}残余漂移",
            passed=False,
            value=float("inf"),
            threshold=tol * 100,
            unit="%/段",
            message="有效数据不足，无法计算漂移",
        )
    n = len(tail)
    x = np.arange(n, dtype=float)
    mean_val = float(tail.mean())

    if n == 2:
        slope = float(tail[1] - tail[0])
    else:
        try:
            coeffs = np.polyfit(x, tail, 1)
            slope = float(coeffs[0])  # 单位/点
        except np.linalg.LinAlgError:
            slope = float((tail[-1] - tail[0]) / max(1, n - 1))
    drift = abs(slope) * n / abs(mean_val) if mean_val != 0 else abs(slope) * n

    return QCCheck(
        name=f"{name}残余漂移",
        passed=drift <= tol,
        value=drift * 100,
        threshold=tol * 100,
        unit="%/段",
        message=(
            f"尾部线性漂移量 {drift*100:.2f}%"
            f"（slope={slope:.4g} {unit}/step）"
        ),
        details={"slope": slope, "drift_fraction": drift, "mean": mean_val},
    )


def _check_energy_conservation(
    etotal: np.ndarray,
    time_ps: np.ndarray,
    tol: float = 1e-4,
) -> QCCheck:
    """
    NVE 能量守恒检验：线性拟合总能量，计算漂移率

    drift_rate = |slope [kcal/mol/ns]| / |mean(E) [kcal/mol]|  (per ns)
    """
    t_ns = (time_ps - time_ps[0]) / 1000.0
    mean_e = float(etotal.mean())

    coeffs = np.polyfit(t_ns, etotal, 1)
    slope = float(coeffs[0])  # kcal/mol per ns
    drift_rate = abs(slope) / abs(mean_e) if mean_e != 0 else 0.0

    # 计算 R²，检测是否有明显线性趋势
    e_pred = np.polyval(coeffs, t_ns)
    ss_res = float(np.sum((etotal - e_pred) ** 2))
    ss_tot = float(np.sum((etotal - mean_e) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return QCCheck(
        name="能量守恒（漂移率）",
        passed=drift_rate <= tol,
        value=drift_rate,
        threshold=tol,
        unit="/ns",
        message=(
            f"漂移率 {drift_rate:.2e}/ns"
            f"（slope={slope:.4f} kcal/mol/ns，"
            f"mean={mean_e:.2f} kcal/mol，R²={r2:.4f}）"
        ),
        details={
            "slope_kcal_ns": slope,
            "drift_rate_per_ns": drift_rate,
            "mean_energy": mean_e,
            "r2": r2,
        },
    )


def summarize_checks(checks: List[QCCheck]) -> dict:
    """Return a machine-readable summary for a list of QC checks."""
    return {
        "passed": bool(checks) and all(check.passed for check in checks),
        "checks": checks,
    }


def _check_window_mean(
    series: np.ndarray,
    *,
    lower: float,
    upper: float,
    name: str,
    unit: str,
    tail_frac: float = 0.5,
) -> QCCheck:
    tail = _tail(series, tail_frac)
    tail = tail[np.isfinite(tail)]
    if len(tail) == 0:
        return QCCheck(
            name=f"{name}窗口均值",
            passed=False,
            value=float("nan"),
            threshold=upper,
            unit=unit,
            message="有效数据不足，无法计算窗口均值",
        )
    mean_val = float(tail.mean())
    passed = lower <= mean_val <= upper
    return QCCheck(
        name=f"{name}窗口均值",
        passed=passed,
        value=mean_val,
        threshold=upper,
        unit=unit,
        message=f"末尾均值 {mean_val:.5g} {unit}，目标窗口 [{lower:.5g}, {upper:.5g}] {unit}",
        details={"mean": mean_val, "lower": lower, "upper": upper},
    )


def evaluate_npt_density_window(
    *,
    log_file: Path,
    density_target: float,
    density_file: Optional[Path] = None,
    density_window_lo: float = 0.95,
    density_window_hi: float = 1.05,
    tail_frac: float = 0.5,
) -> dict:
    data = parse_lammps_log(log_file)
    rho_series: Optional[np.ndarray] = None
    if "density" in data and len(data["density"]) > 3:
        rho_series = data["density"]
    elif density_file is not None and density_file.exists():
        rho_series = parse_ave_time_density(density_file)

    checks: List[QCCheck] = []
    if rho_series is None:
        checks.append(
            QCCheck(
                name="密度窗口均值",
                passed=False,
                value=float("nan"),
                threshold=density_target * density_window_hi,
                unit="g/cm³",
                message="缺少密度时间序列，无法判断 NPT 是否达到目标密度窗口",
            )
        )
        return summarize_checks(checks)

    checks.append(
        _check_window_mean(
            rho_series,
            lower=density_target * density_window_lo,
            upper=density_target * density_window_hi,
            name="密度",
            unit="g/cm³",
            tail_frac=tail_frac,
        )
    )
    return summarize_checks(checks)


def evaluate_nvt_thermo_stage(
    *,
    log_file: Path,
    temperature_target: float,
    tail_frac: float = 0.5,
    temperature_mean_tol: float = 0.10,
    temperature_std_tol: float = 0.10,
    etotal_std_tol: float = 0.10,
    pe_std_tol: float = 0.10,
    ke_std_tol: float = 0.01,
) -> dict:
    data = parse_lammps_log(log_file)
    checks: List[QCCheck] = []

    if "temp" in data:
        temp = data["temp"]
        checks.append(_check_mean_bias(temp, temperature_target, "温度", "K", temperature_mean_tol, tail_frac))
        checks.append(_check_stability(temp, "温度", "K", temperature_std_tol, tail_frac))
    if "etotal" in data:
        checks.append(_check_stability(data["etotal"], "总能量", "kcal/mol", etotal_std_tol, tail_frac))
    if "pe" in data:
        checks.append(_check_stability(data["pe"], "总势能", "kcal/mol", pe_std_tol, tail_frac))
    if "ke" in data:
        checks.append(_check_stability(data["ke"], "总动能", "kcal/mol", ke_std_tol, tail_frac))

    return summarize_checks(checks)


def evaluate_npt_segment(
    *,
    log_file: Path,
    temperature_target: float,
    density_target: Optional[float] = None,
    density_file: Optional[Path] = None,
    tail_frac: float = 0.5,
    temperature_mean_tol: float = 0.02,
    temperature_std_tol: float = 0.05,
    density_mean_tol: float = 0.05,
    density_std_tol: float = 0.01,
    density_drift_tol: float = 0.01,
) -> dict:
    qc = SimulationQC()
    checks = qc.check_npt_equilibration(
        log_file=log_file,
        T_target=temperature_target,
        rho_exp=density_target,
        density_file=density_file,
        T_mean_tol=temperature_mean_tol,
        rho_mean_tol=density_mean_tol,
        rho_std_tol=density_std_tol,
        rho_drift_tol=density_drift_tol,
        tail_frac=tail_frac,
    )
    if "temp" in parse_lammps_log(log_file):
        for idx, check in enumerate(checks):
            if check.name == "温度稳定性（RSD）":
                checks[idx] = _check_stability(
                    parse_lammps_log(log_file)["temp"],
                    "温度",
                    "K",
                    temperature_std_tol,
                    tail_frac,
                )
                break
    return summarize_checks(checks)


def evaluate_nvt_segment(
    *,
    log_file: Path,
    temperature_target: float,
    tail_frac: float = 0.5,
    temperature_mean_tol: float = 0.02,
    temperature_std_tol: float = 0.05,
    temperature_drift_tol: float = 0.02,
    energy_key: str = "etotal",
    energy_fluctuation_tol: float = 0.01,
) -> dict:
    data = parse_lammps_log(log_file)
    checks: List[QCCheck] = []

    if "temp" in data:
        temp = data["temp"]
        checks.append(_check_mean_bias(temp, temperature_target, "温度", "K", temperature_mean_tol, tail_frac))
        checks.append(_check_stability(temp, "温度", "K", temperature_std_tol, tail_frac))
        checks.append(_check_drift(temp, "温度", "K", temperature_drift_tol, tail_frac))

    if energy_key in data:
        checks.append(
            _check_stability(
                data[energy_key],
                f"{energy_key}能量",
                "kcal/mol",
                energy_fluctuation_tol,
                tail_frac,
            )
        )

    return summarize_checks(checks)


# ──────────────────────────────────────────────────────────────────────────────
# 高层 QC 类
# ──────────────────────────────────────────────────────────────────────────────

class SimulationQC:
    """
    模拟质量控制器

    对各阶段的 LAMMPS 输出文件进行结构化质量检验。

    Example::

        from gk_workflow.simulation.qc import SimulationQC
        from pathlib import Path

        qc = SimulationQC()

        # 一键全流程检查
        qc.full_workflow_check(
            workdir=Path("./work"),
            T_target=298.0,
            rho_exp=0.867,      # 实验密度（可选）
        )
    """

    # ── NVT 平衡 ──────────────────────────────────────────────────────────────

    def check_nvt_equilibration(
        self,
        log_file: Path,
        T_target: float,
        T_mean_tol: float = 0.02,    # 均值偏差容差 2%
        T_std_tol: float = 0.05,     # 稳定性容差 5%（小体系波动大）
        T_drift_tol: float = 0.02,   # 尾部漂移容差 2%
        tail_frac: float = 0.3,      # 取末尾 30% 做统计
    ) -> List[QCCheck]:
        """
        检查 NVT 平衡质量（从 nvt_equilibration.log 读取）

        Args:
            log_file: nvt_equilibration.log 路径
            T_target: 目标温度 (K)
            T_mean_tol: 温度均值偏差容差（相对）
            T_std_tol: 温度相对标准差容差（NVT 系综波动正常）
            T_drift_tol: 尾部线性漂移容差（相对）
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_lammps_log(log_file)
        checks: List[QCCheck] = []

        if "temp" not in data:
            return checks

        temp = data["temp"]
        checks.append(_check_mean_bias(temp, T_target, "温度", "K", T_mean_tol, tail_frac))
        checks.append(_check_stability(temp, "温度", "K", T_std_tol, tail_frac))
        checks.append(_check_drift(temp, "温度", "K", T_drift_tol, tail_frac))

        return checks

    # ── NPT 平衡 ──────────────────────────────────────────────────────────────

    def check_npt_equilibration(
        self,
        log_file: Path,
        T_target: float,
        rho_exp: Optional[float] = None,
        density_file: Optional[Path] = None,
        T_mean_tol: float = 0.02,
        rho_mean_tol: float = 0.05,   # 密度与实验值偏差容差 5%
        rho_std_tol: float = 0.01,    # 密度稳定性容差 1%
        rho_drift_tol: float = 0.01,  # 密度漂移容差 1%
        tail_frac: float = 0.3,
    ) -> List[QCCheck]:
        """
        检查 NPT 平衡质量（从 npt_equilibration.log 和 density_npt.dat 读取）

        Args:
            log_file: npt_equilibration.log 路径
            T_target: 目标温度 (K)
            rho_exp: 实验密度 (g/cm³)，None 则跳过与实验值对比
            density_file: density_npt.dat 路径（可选，密度时间序列更细）
            T_mean_tol: 温度均值偏差容差
            rho_mean_tol: 密度与实验值偏差容差
            rho_std_tol: 密度稳定性容差
            rho_drift_tol: 密度漂移容差
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_lammps_log(log_file)
        checks: List[QCCheck] = []

        # ── 温度
        if "temp" in data:
            temp = data["temp"]
            checks.append(_check_mean_bias(temp, T_target, "温度", "K", T_mean_tol, tail_frac))
            checks.append(_check_stability(temp, "温度", "K", 0.05, tail_frac))

        # ── 密度
        rho_series: Optional[np.ndarray] = None

        if "density" in data and len(data["density"]) > 5:
            rho_series = data["density"]
        elif density_file is not None and density_file.exists():
            try:
                rho_series = parse_ave_time_density(density_file)
            except Exception:
                pass

        if rho_series is not None:
            if rho_exp is not None:
                checks.append(
                    _check_mean_bias(rho_series, rho_exp, "密度", "g/cm³", rho_mean_tol, tail_frac)
                )
            checks.append(_check_stability(rho_series, "密度", "g/cm³", rho_std_tol, tail_frac))
            checks.append(_check_drift(rho_series, "密度", "g/cm³", rho_drift_tol, tail_frac))

        return checks

    # ── NVE 生产 ──────────────────────────────────────────────────────────────

    def check_nve_production(
        self,
        thermo_file: Path,
        T_target: Optional[float] = None,
        drift_tol: float = 1e-4,    # 能量漂移率容差（/ns）
        T_std_tol: float = 0.03,    # 温度稳定性（NVE 波动较大）
        tail_frac: float = 0.3,
    ) -> List[QCCheck]:
        """
        检查 NVE 生产运行质量（从 nve_prod_thermo.dat 读取）

        重点检查能量守恒。NVE 系综中温度会有较大统计波动，属正常现象。

        Args:
            thermo_file: nve_prod_thermo.dat 路径
            T_target: 目标温度 (K)，None 则只检查稳定性
            drift_tol: 能量漂移率容差（per ns），推荐 1e-4
            T_std_tol: 温度相对标准差容差
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_thermo_dat(thermo_file)
        checks: List[QCCheck] = []

        # ── 能量守恒（核心检查）
        if "etotal" in data and "time" in data:
            checks.append(
                _check_energy_conservation(data["etotal"], data["time"], tol=drift_tol)
            )

        # ── 温度稳定性（NVE 波动允许较大）
        if "temp" in data:
            temp = data["temp"]
            if T_target is not None:
                checks.append(
                    _check_mean_bias(temp, T_target, "温度", "K", 0.03, tail_frac)
                )
            checks.append(_check_stability(temp, "温度", "K", T_std_tol, tail_frac))

        return checks

    # ── NVT 生产 ──────────────────────────────────────────────────────────────

    def check_nvt_production(
        self,
        thermo_file: Path,
        T_target: float,
        T_mean_tol: float = 0.02,
        T_std_tol: float = 0.05,
        tail_frac: float = 0.3,
    ) -> List[QCCheck]:
        """
        检查 NVT 生产运行质量（从 nvt_prod_thermo.dat 读取）

        Args:
            thermo_file: nvt_prod_thermo.dat 路径
            T_target: 目标温度 (K)
            T_mean_tol: 温度均值偏差容差
            T_std_tol: 温度稳定性容差
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_thermo_dat(thermo_file)
        checks: List[QCCheck] = []

        if "temp" in data:
            temp = data["temp"]
            checks.append(_check_mean_bias(temp, T_target, "温度", "K", T_mean_tol, tail_frac))
            checks.append(_check_stability(temp, "温度", "K", T_std_tol, tail_frac))
            checks.append(_check_drift(temp, "温度", "K", 0.01, tail_frac))

        return checks

    # ── NPT 生产 ──────────────────────────────────────────────────────────────

    def check_npt_production(
        self,
        thermo_file: Path,
        T_target: float,
        rho_exp: Optional[float] = None,
        T_mean_tol: float = 0.02,
        rho_mean_tol: float = 0.05,
        rho_std_tol: float = 0.01,
        tail_frac: float = 0.3,
    ) -> List[QCCheck]:
        """
        检查 NPT 生产运行质量（从 npt_prod_thermo.dat 读取）

        Args:
            thermo_file: npt_prod_thermo.dat 路径
            T_target: 目标温度 (K)
            rho_exp: 实验密度 (g/cm³)（可选）
            T_mean_tol: 温度均值偏差容差
            rho_mean_tol: 密度偏差容差（与实验值）
            rho_std_tol: 密度稳定性容差
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_thermo_dat(thermo_file)
        checks: List[QCCheck] = []

        if "temp" in data:
            checks.append(
                _check_mean_bias(data["temp"], T_target, "温度", "K", T_mean_tol, tail_frac)
            )
            checks.append(_check_stability(data["temp"], "温度", "K", 0.05, tail_frac))

        if "density" in data:
            rho = data["density"]
            if rho_exp is not None:
                checks.append(
                    _check_mean_bias(rho, rho_exp, "密度", "g/cm³", rho_mean_tol, tail_frac)
                )
            checks.append(_check_stability(rho, "密度", "g/cm³", rho_std_tol, tail_frac))

        return checks

    # ── GK 采样 ───────────────────────────────────────────────────────────────

    def check_gk_sampling(
        self,
        gk_data_file: Path,
        T_target: float,
        T_mean_tol: float = 0.03,
        T_std_tol: float = 0.04,
        vol_std_tol: float = 0.002,  # NVE 体积应极其稳定
        tail_frac: float = 0.3,
    ) -> List[QCCheck]:
        """
        检查 Green-Kubo 采样运行质量（从 greenkubo_data.dat 读取）

        NVE 系综：温度围绕目标值波动，体积保持恒定。

        Args:
            gk_data_file: greenkubo_data.dat 路径
            T_target: 目标温度 (K)
            T_mean_tol: 温度均值偏差容差
            T_std_tol: 温度稳定性容差
            vol_std_tol: 体积稳定性容差（NVE 应 < 0.2%）
            tail_frac: 统计尾部比例

        Returns:
            QCCheck 列表
        """
        data = parse_gk_data(gk_data_file)
        checks: List[QCCheck] = []

        if "temp" in data:
            temp = data["temp"]
            checks.append(_check_mean_bias(temp, T_target, "温度", "K", T_mean_tol, tail_frac))
            checks.append(_check_stability(temp, "温度", "K", T_std_tol, tail_frac))

        if "vol" in data:
            checks.append(_check_stability(data["vol"], "体积", "Å³", vol_std_tol, tail_frac))

        return checks

    # ── 报告打印 ──────────────────────────────────────────────────────────────

    @staticmethod
    def print_report(stage_name: str, checks: List[QCCheck]) -> None:
        """打印单阶段 QC 报告"""
        print("=" * 70)
        print(f"QC 报告 — {stage_name}")
        print("=" * 70)

        if not checks:
            print("  （无可用数据，跳过）")
            print()
            return

        n_pass = sum(1 for c in checks if c.passed)
        n_total = len(checks)

        for c in checks:
            icon = "✓" if c.passed else "✗"
            status = "PASS" if c.passed else "FAIL"
            print(f"  [{icon} {status}]  {c.name}")
            print(f"           {c.message}")
            print(f"           测量值 {c.value:.3g} {c.unit}，阈值 {c.threshold:.3g} {c.unit}")

        overall = "全部通过 ✓" if n_pass == n_total else f"{n_pass}/{n_total} 项通过"
        print()
        print(f"  总结：{overall}")
        print("=" * 70 + "\n")

    # ── 全流程一键检查 ────────────────────────────────────────────────────────

    def full_workflow_check(
        self,
        workdir: Path,
        T_target: float,
        P_target: float = 1.0,
        rho_exp: Optional[float] = None,
        nve_prefix: str = "nve_prod",
        nvt_prefix: str = "nvt_prod",
        npt_prefix: str = "npt_prod",
        gk_prefix: str = "greenkubo",
    ) -> Dict[str, List[QCCheck]]:
        """
        对工作目录中所有标准输出文件执行完整 QC 检查

        自动查找并检查（文件存在时）：
          - nvt_equilibration.log
          - npt_equilibration.log / density_npt.dat
          - {nve_prefix}_thermo.dat
          - {nvt_prefix}_thermo.dat
          - {npt_prefix}_thermo.dat
          - {gk_prefix}_data.dat

        Args:
            workdir: 工作目录
            T_target: 目标温度 (K)
            P_target: 目标压力 (atm)
            rho_exp: 实验密度 (g/cm³)（可选）
            nve_prefix / nvt_prefix / npt_prefix / gk_prefix: 文件前缀

        Returns:
            {阶段名: [QCCheck, ...]}
        """
        workdir = Path(workdir)
        all_results: Dict[str, List[QCCheck]] = {}

        _stages = [
            ("NVT平衡", workdir / "nvt_equilibration.log",
             lambda f: self.check_nvt_equilibration(f, T_target)),
            ("NPT平衡", workdir / "npt_equilibration.log",
             lambda f: self.check_npt_equilibration(
                 f, T_target, rho_exp=rho_exp,
                 density_file=workdir / "density_npt.dat"
             )),
            ("NVE生产", workdir / f"{nve_prefix}_thermo.dat",
             lambda f: self.check_nve_production(f, T_target)),
            ("NVT生产", workdir / f"{nvt_prefix}_thermo.dat",
             lambda f: self.check_nvt_production(f, T_target)),
            ("NPT生产", workdir / f"{npt_prefix}_thermo.dat",
             lambda f: self.check_npt_production(f, T_target, rho_exp=rho_exp)),
            ("GK采样", workdir / f"{gk_prefix}_data.dat",
             lambda f: self.check_gk_sampling(f, T_target)),
        ]

        for stage_name, filepath, checker in _stages:
            if filepath.exists():
                try:
                    checks = checker(filepath)
                    all_results[stage_name] = checks
                    self.print_report(stage_name, checks)
                except Exception as e:
                    print(f"  [警告] {stage_name} QC 失败: {e}\n")

        self._print_summary(all_results)
        return all_results

    @staticmethod
    def _print_summary(all_results: Dict[str, List[QCCheck]]) -> None:
        """打印总体汇总"""
        if not all_results:
            return

        print("=" * 70)
        print("总体 QC 汇总")
        print("=" * 70)

        any_fail = False
        for stage, checks in all_results.items():
            if not checks:
                print(f"  -  {stage}：无数据")
                continue
            n_pass = sum(1 for c in checks if c.passed)
            n_total = len(checks)
            icon = "✓" if n_pass == n_total else "✗"
            if n_pass < n_total:
                any_fail = True
            print(f"  {icon}  {stage}：{n_pass}/{n_total} 项通过")

            # 列出未通过的项
            failed = [c for c in checks if not c.passed]
            for c in failed:
                print(f"       → FAIL: {c.name}（{c.message}）")

        print()
        if any_fail:
            print("  结论：存在质量问题，请检查上方 FAIL 项目后再进行分析。")
        else:
            print("  结论：所有检查通过，模拟质量良好，可继续进行分析。")
        print("=" * 70 + "\n")
