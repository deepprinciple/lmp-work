"""Built-in geometries and partial charges for common electrolyte ions.

OpenFF AM1-BCC is unreliable (or undefined) for monatomic and some polyatomic
ions, so the workflow ships hand-built geometries plus known per-atom charges
that match the BAMBOO reference data file:

- Li+   single atom, charge +1.0
- PF6-  regular octahedron (P-F = 1.60 A), P:+1.34, F:-0.39 each
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class IonPreset:
    """Static description of an ion as geometry + per-atom charges."""

    smiles: str
    geometry: tuple[tuple[str, float, float, float], ...]
    charges: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.geometry) != len(self.charges):
            raise ValueError(
                f"geometry and charges length mismatch: "
                f"{len(self.geometry)} vs {len(self.charges)}"
            )


# Regular octahedron for PF6-; P-F bond ~1.60 A (experimental ~1.59 A).
_PF6_BOND = 1.60
PF6_MINUS = IonPreset(
    smiles="F[P-](F)(F)(F)(F)F",
    geometry=(
        ("P", 0.000, 0.000, 0.000),
        ("F",  _PF6_BOND, 0.0, 0.0),
        ("F", -_PF6_BOND, 0.0, 0.0),
        ("F", 0.0,  _PF6_BOND, 0.0),
        ("F", 0.0, -_PF6_BOND, 0.0),
        ("F", 0.0, 0.0,  _PF6_BOND),
        ("F", 0.0, 0.0, -_PF6_BOND),
    ),
    charges=(1.34, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39),
)


ION_PRESETS: dict[str, IonPreset] = {
    "Li+": IonPreset(smiles="[Li+]", geometry=(("Li", 0.0, 0.0, 0.0),), charges=(1.0,)),
    "Na+": IonPreset(smiles="[Na+]", geometry=(("Na", 0.0, 0.0, 0.0),), charges=(1.0,)),
    "K+":  IonPreset(smiles="[K+]",  geometry=(("K",  0.0, 0.0, 0.0),), charges=(1.0,)),
    "PF6-": PF6_MINUS,
}


__all__ = ["IonPreset", "ION_PRESETS", "PF6_MINUS"]
