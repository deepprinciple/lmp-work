"""
分析模块：Green-Kubo 分析、介电常数与统计工具
"""
from .greenkubo import GreenKuboAnalyzer, compute_viscosity, compute_thermal_conductivity
from .dielectric import DielectricAnalyzer, compute_dielectric_constant
from .statistics import (
    block_average,
    block_average_scan,
    statistical_inefficiency,
    autocorrelation_time,
    running_average,
    equilibration_check,
    integrate_acf,
)

__all__ = [
    # Green-Kubo
    "GreenKuboAnalyzer",
    "compute_viscosity",
    "compute_thermal_conductivity",
    # 介电常数
    "DielectricAnalyzer",
    "compute_dielectric_constant",
    # 统计工具
    "block_average",
    "block_average_scan",
    "statistical_inefficiency",
    "autocorrelation_time",
    "running_average",
    "equilibration_check",
    "integrate_acf",
]
