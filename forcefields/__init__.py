"""Force-field implementations."""
from .ligpargen import LigParGen, get_ligpargen_lammps_settings
from .openff import OpenFF, get_openff_lammps_settings


def get_lammps_settings_for_engine(engine: str) -> dict:
    """Return default LAMMPS settings for a supported force-field engine."""
    token = str(engine).strip().lower()
    if token == "openff":
        return get_openff_lammps_settings()
    if token == "ligpargen":
        return get_ligpargen_lammps_settings()
    raise ValueError(f"Unsupported forcefield engine: {engine}")


__all__ = [
    "LigParGen",
    "OpenFF",
    "get_lammps_settings_for_engine",
    "get_ligpargen_lammps_settings",
    "get_openff_lammps_settings",
]
