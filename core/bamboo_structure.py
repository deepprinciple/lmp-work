"""3D single-molecule structure generation.

Two paths:

1. ``generate_molecule(smiles, name)`` for organic solvents — RDKit ETKDG
   embedding plus MMFF (fallback UFF) optimisation. Returns the RDKit Mol
   together with the per-atom ``(element, x, y, z)`` tuples in RDKit atom
   order (this same order is what the OpenFF charge generator and the
   data-file writer will assume downstream).

2. ``load_ion_preset(name)`` for ions covered by ``forcefields/ions.py`` —
   skips RDKit / MMFF entirely and returns the built-in template geometry.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from forcefields.ions import ION_PRESETS, IonPreset
from utils.constants import ATOMIC_MASS

XyzAtoms = list[tuple[str, float, float, float]]


def generate_molecule(
    smiles: str,
    name: str,
    *,
    optimize: bool = True,
    max_iters: int = 500,
    seed: int = 42,
) -> tuple[Chem.Mol, XyzAtoms, float]:
    """SMILES → ``(rdkit_mol, xyz_atoms, molecular_weight_g_per_mol)``.

    The RDKit Mol carries explicit hydrogens; ``xyz_atoms`` is a list of
    ``(element_symbol, x, y, z)`` in the same RDKit atom order.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {name!r}: {smiles}")
    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol, randomSeed=seed)
    if status != 0:
        # ETKDG occasionally fails; one retry with a different seed is usually
        # enough for the small organics we care about.
        status = AllChem.EmbedMolecule(mol, randomSeed=seed + 7919, useRandomCoords=True)
    if status != 0:
        raise RuntimeError(f"RDKit failed to embed a 3D conformer for {name} ({smiles})")

    if optimize:
        optimized = False
        if AllChem.MMFFHasAllMoleculeParams(mol):
            try:
                rc = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
                optimized = rc >= 0
            except Exception:
                optimized = False
        if not optimized:
            # UFF is more permissive — use it for charged or unusual species.
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
            except Exception:
                pass  # accept the embedded geometry as-is

    conf = mol.GetConformer()
    xyz_atoms: XyzAtoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        xyz_atoms.append((atom.GetSymbol(), float(pos.x), float(pos.y), float(pos.z)))

    mw = float(Descriptors.MolWt(mol))
    return mol, xyz_atoms, mw


def load_ion_preset(preset_name: str) -> tuple[XyzAtoms, list[float], float]:
    """Return ``(xyz_atoms, charges, mw)`` for a built-in ion preset."""
    if preset_name not in ION_PRESETS:
        raise KeyError(
            f"Unknown ion preset: {preset_name!r}. "
            f"Available: {sorted(ION_PRESETS)}"
        )
    preset: IonPreset = ION_PRESETS[preset_name]
    xyz_atoms: XyzAtoms = [tuple(a) for a in preset.geometry]  # type: ignore[misc]
    charges = list(preset.charges)
    mw = sum(ATOMIC_MASS[elem] for elem, *_ in xyz_atoms)
    return xyz_atoms, charges, mw


def write_xyz(
    path: Path,
    xyz_atoms: Iterable[tuple[str, float, float, float]],
    comment: str = "",
) -> Path:
    """Write a plain XYZ file (one molecule per file, used as Packmol input)."""
    atoms = list(xyz_atoms)
    lines = [str(len(atoms)), comment]
    for elem, x, y, z in atoms:
        lines.append(f"{elem:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    path.write_text("\n".join(lines) + "\n")
    return path


__all__ = ["generate_molecule", "load_ion_preset", "write_xyz", "XyzAtoms"]
