"""Thermal mainline analysis exports."""

from __future__ import annotations

from importlib import import_module
from typing import Any

__all__ = [
    "analyze_hfacf_file",
    "analyze_rnemd_replica",
    "compute_thermal_conductivity",
    "load_volume_from_gk_data",
    "parse_ave_correlate_detail",
]


def __getattr__(name: str) -> Any:
    if name in {
        "analyze_hfacf_file",
        "compute_thermal_conductivity",
        "load_volume_from_gk_data",
        "parse_ave_correlate_detail",
    }:
        module = import_module(".hfacf", __name__)
        return getattr(module, name)
    if name == "analyze_rnemd_replica":
        module = import_module(".rnemd", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
