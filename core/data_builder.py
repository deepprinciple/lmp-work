"""Write a BAMBOO-compatible LAMMPS data file (``in.data``).

BAMBOO expects ``atom_style full`` with per-atom partial charges. The data file
does *not* contain Bonds / Angles / Dihedrals / Impropers — BAMBOO learns the
short-range physics from coordinates + element types, and PPPM handles
long-range electrostatics via the charges we write here.

Atom-type ordering: we use ascending atomic number for whichever elements are
present in the system. The element list is returned alongside the data path so
that ``workflow.input_viscosity`` can render the matching ``pair_coeff`` line.

The element symbols used in Masses comments, Atoms comments, and the
``pair_coeff`` line are all upper-cased (``Li → LI``) to match the BAMBOO
sample data file convention.
"""
from __future__ import annotations

from pathlib import Path

from utils.constants import ATOMIC_MASS, ATOMIC_NUMBER

AtomRecord = tuple[str, float, float, float, float, int]
# (element_symbol, partial_charge, x, y, z, mol_id)


def build_bamboo_data(
    atoms: list[AtomRecord],
    box_lengths: tuple[float, float, float],
    output: Path,
    *,
    title: str = "LAMMPS data file - BAMBOO system",
    neutrality_tolerance: float = 1e-3,
) -> tuple[Path, list[str]]:
    """Write ``in.data`` and return ``(output_path, elements_for_pair_coeff)``.

    Parameters
    ----------
    atoms:
        Per-atom records ``(element, charge, x, y, z, mol_id)`` in the order
        they appear in the packed box. ``mol_id`` should be sequentially
        assigned, one unique value per molecule.
    box_lengths:
        Orthorhombic box dimensions ``(Lx, Ly, Lz)`` in A. Origin at ``(0,0,0)``.
    output:
        Destination path.
    title:
        First-line title comment.
    neutrality_tolerance:
        Abort if the box's net charge exceeds this — PPPM cannot converge on
        a non-neutral system.

    Returns
    -------
    output_path:
        Absolute path to the data file.
    elements:
        Element symbols (UPPER-cased) ordered so that the index corresponds to
        the atom_type number used in the data file. Pass this through to the
        ``pair_coeff`` writer.
    """
    if not atoms:
        raise ValueError("Empty atom list — nothing to write.")

    Lx, Ly, Lz = box_lengths
    if min(Lx, Ly, Lz) <= 0:
        raise ValueError(f"Box must have positive lengths, got {box_lengths}")

    # Neutrality check (whole box) before any I/O.
    total_charge = sum(charge for _, charge, *_ in atoms)
    if abs(total_charge) > neutrality_tolerance:
        raise RuntimeError(
            f"System is not charge-neutral: net charge = {total_charge:+.6f} e "
            f"(tolerance {neutrality_tolerance}). PPPM cannot run on a charged "
            "cell — check composition / charge assignment."
        )

    # Determine present elements and assign deterministic type IDs.
    present = {elem for elem, *_ in atoms}
    unknown = [e for e in present if e not in ATOMIC_NUMBER]
    if unknown:
        raise RuntimeError(
            f"Unknown element(s) {unknown}. Extend ATOMIC_NUMBER / ATOMIC_MASS "
            "in utils/constants.py to include them."
        )
    elements_internal = sorted(present, key=lambda e: ATOMIC_NUMBER[e])
    elements_display = [e.upper() for e in elements_internal]
    type_map = {e: i + 1 for i, e in enumerate(elements_internal)}

    lines: list[str] = [
        f"# {title}",
        "",
        f"{len(atoms)} atoms",
        f"{len(elements_internal)} atom types",
        "",
        f"0.0 {Lx:.6f} xlo xhi",
        f"0.0 {Ly:.6f} ylo yhi",
        f"0.0 {Lz:.6f} zlo zhi",
        "",
        "Masses",
        "",
    ]
    for e in elements_internal:
        lines.append(f"{type_map[e]} {ATOMIC_MASS[e]:.4f} # {e.upper()}")

    lines += ["", "Atoms # full", ""]
    for idx, (elem, charge, x, y, z, mol_id) in enumerate(atoms, start=1):
        lines.append(
            f"{idx} {mol_id} {type_map[elem]} {charge:.6f} "
            f"{x:.6f} {y:.6f} {z:.6f} # {elem.upper()}"
        )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")

    print(
        f"  wrote {output}  "
        f"({len(atoms)} atoms, {len(elements_internal)} types: "
        f"{' '.join(elements_display)})"
    )
    return output.resolve(), elements_display


__all__ = ["AtomRecord", "build_bamboo_data"]
