"""
统计分析模块：Block平均、误差估计、自相关时间

用于评估MD模拟结果的统计不确定度。
"""
import numpy as np
from typing import Dict, Tuple, Optional


def block_average(
    data: np.ndarray,
    n_blocks: int = 10,
) -> Tuple[float, float, float]:
    """
    Block平均法估计统计误差

    将数据切分为n_blocks个连续块，计算块均值的标准误差。
    当块大小远大于相关时间时，各块独立，SEM估计有效。

    Args:
        data: 1D时间序列数据
        n_blocks: 块数（通常10–20）

    Returns:
        (均值, 块标准误差, 相对误差)
    """
    n = len(data)
    block_size = n // n_blocks
    if block_size < 5:
        raise ValueError(
            f"块大小过小 ({block_size})，请减少 n_blocks 或增加数据量"
        )

    # 只使用能被整除的部分
    data_use = data[: n_blocks * block_size]
    block_means = data_use.reshape(n_blocks, block_size).mean(axis=1)

    mean = float(block_means.mean())
    sem = float(block_means.std(ddof=1) / np.sqrt(n_blocks))
    rel_error = abs(sem / mean) if mean != 0.0 else 0.0

    return mean, sem, rel_error


def block_average_scan(
    data: np.ndarray,
    min_blocks: int = 5,
    max_blocks: int = 50,
) -> Dict[str, list]:
    """
    扫描不同块大小，寻找误差收敛的平台区

    当块大小超过数据相关时间后，SEM趋于稳定（平台）。
    查看 errors vs block_sizes 图像，平台处的误差即为真实统计误差。

    Args:
        data: 1D时间序列数据
        min_blocks: 最小块数
        max_blocks: 最大块数

    Returns:
        字典: {block_sizes, means, errors, rel_errors}
    """
    n = len(data)
    result: Dict[str, list] = {
        "block_sizes": [],
        "means": [],
        "errors": [],
        "rel_errors": [],
    }

    for nb in range(min_blocks, min(max_blocks + 1, n // 5 + 1)):
        try:
            mean, sem, rel_err = block_average(data, n_blocks=nb)
            result["block_sizes"].append(n // nb)
            result["means"].append(mean)
            result["errors"].append(sem)
            result["rel_errors"].append(rel_err)
        except ValueError:
            break

    return result


def statistical_inefficiency(
    data: np.ndarray,
    max_lag: Optional[int] = None,
) -> float:
    """
    计算统计低效性 g（statistical inefficiency）

    g = 1 + 2 * Σ_{t=1}^{τ} C(t)/C(0)

    其中 C(t) 为归一化自相关函数，截断条件为 C(t) ≤ 0。
    有效独立样本数 N_eff = N / g。

    Args:
        data: 1D时间序列数据
        max_lag: 最大滞后步数，None 则取 N//4

    Returns:
        g ≥ 1
    """
    n = len(data)
    if max_lag is None:
        max_lag = n // 4

    data_c = data - data.mean()
    c0 = float(np.dot(data_c, data_c)) / n
    if c0 == 0.0:
        return 1.0

    g = 1.0
    for lag in range(1, min(max_lag + 1, n)):
        c_lag = float(np.dot(data_c[: n - lag], data_c[lag:])) / (n - lag)
        rho = c_lag / c0
        if rho <= 0.0:
            break
        g += 2.0 * rho

    return max(1.0, g)


def autocorrelation_time(
    data: np.ndarray,
    max_lag: Optional[int] = None,
) -> float:
    """
    计算积分自相关时间 τ = (g − 1) / 2

    Args:
        data: 1D时间序列数据
        max_lag: 最大滞后步数

    Returns:
        自相关时间（数据点单位）
    """
    return (statistical_inefficiency(data, max_lag=max_lag) - 1.0) / 2.0


def running_average(data: np.ndarray, window: int) -> np.ndarray:
    """
    滑动平均（等权重）

    Args:
        data: 1D数据
        window: 窗口大小

    Returns:
        长度为 len(data) - window + 1 的数组
    """
    return np.convolve(data, np.ones(window) / window, mode="valid")


def equilibration_check(
    data: np.ndarray,
    n_fracs: int = 5,
) -> Tuple[int, float]:
    """
    简单平衡性检验：将数据按比例分段，找前后两半均值差 < 5% 的最早起点

    Args:
        data: 1D时间序列数据
        n_fracs: 检验的分割比例数量（0.1, 0.2, ..., 0.5）

    Returns:
        (建议丢弃的点数 discard_idx, 平衡后的均值)
    """
    n = len(data)
    best_start = 0

    for i in range(1, n_fracs + 1):
        frac = i / (2 * n_fracs)  # 0.1, 0.2, ..., 0.5
        start = int(n * frac)
        tail = data[start:]
        if len(tail) < 20:
            break

        half = len(tail) // 2
        m1 = tail[:half].mean()
        m2 = tail[half:].mean()
        m_total = tail.mean()

        diff = abs(m1 - m2) / abs(m_total) if m_total != 0 else abs(m1 - m2)
        if diff < 0.05:
            best_start = start
            break

    return best_start, float(data[best_start:].mean())


def integrate_acf(
    time: np.ndarray,
    acf: np.ndarray,
    t_max: Optional[float] = None,
) -> Tuple[np.ndarray, float]:
    """
    数值积分自相关函数（梯形法）

    Args:
        time: 时间数组（同单位）
        acf: 自相关函数值
        t_max: 积分截止时间，None 则积分到全程

    Returns:
        (running_integral, total_integral)
        running_integral 为每个时间点的累积积分值
    """
    if t_max is not None:
        mask = time <= t_max
        time = time[mask]
        acf = acf[mask]

    running = np.zeros(len(time))
    for i in range(1, len(time)):
        dt = time[i] - time[i - 1]
        running[i] = running[i - 1] + 0.5 * (acf[i - 1] + acf[i]) * dt

    return running, float(running[-1])
