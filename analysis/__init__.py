"""Analysis exports for the BAMBOO viscosity demo."""

from __future__ import annotations

from importlib import import_module
from typing import Any

__all__ = [
    "analyze_viscosity_file",
    "analyze_viscosity_pressure_file",
    "compute_viscosity",
    "compute_viscosity_from_pressure_file",
    "load_volume_from_thermo",
    "parse_ave_correlate_detail",
]


def __getattr__(name: str) -> Any:
    if name in {
        "analyze_viscosity_file",
        "analyze_viscosity_pressure_file",
        "compute_viscosity",
        "compute_viscosity_from_pressure_file",
        "load_volume_from_thermo",
    }:
        module = import_module(".viscosity", __name__)
        return getattr(module, name)
    if name == "parse_ave_correlate_detail":
        module = import_module(".correlation", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
