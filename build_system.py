"""
Build system.data from SMILES (minimal workflow).

Steps:
1) RDKit 3D structure
2) Packmol liquid box
3) OpenFF single-molecule parameterization
4) Build multi-molecule system.data
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Import through the package root so the script works from the project directory.
sys.path.insert(0, str(Path(__file__).parent.parent))

from gk_workflow.core import MoleculeStructure, PackmolBuilder, LAMMPSBuilder
from gk_workflow.forcefields import OpenFF


def build_system_data_only(
    smiles: str,
    name: str,
    workdir: Path,
    n_molecules: int,
    density: float,
    *,
    packmol_seed: int = 192911,
    packmol_max_attempts: int = 3,
    packmol_seed_step: int = 97,
    packmol_nloop: int | None = None,
    packmol_full_box: bool = False,
    packmol_margin: float | None = None,
    packmol_strict: bool = False,
) -> dict:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Build System Workflow (up to system.data)")
    print("=" * 70)
    print(f"  molecule : {name} ({smiles})")
    print(f"  workdir  : {workdir}")
    print(f"  n_mol    : {n_molecules}")
    print(f"  density  : {density:.6f} g/cm^3")
    print("=" * 70 + "\n")

    # Step 1: 3D structure
    structure = MoleculeStructure(smiles, name)
    mol = structure.generate(optimize=True)
    mol_xyz = structure.save_xyz(workdir / f"{name}.xyz")

    # Step 2: Packmol
    packer = PackmolBuilder(workdir)
    system_xyz, box_length = packer.build_box(
        mol_xyz=mol_xyz,
        n_molecules=n_molecules,
        density=density,
        mol_weight=structure.molecular_weight,
        seed=packmol_seed,
        max_attempts=packmol_max_attempts,
        seed_step=packmol_seed_step,
        nloop=packmol_nloop,
        full_box=packmol_full_box,
        margin=packmol_margin,
        allow_imperfect=not packmol_strict,
    )

    # Step 3: OpenFF parameterization
    forcefield = OpenFF(workdir)
    single_data = forcefield.generate_parameters(mol, smiles, name)

    # Step 4: Build system.data
    builder = LAMMPSBuilder(workdir)
    system_data = builder.build_system(
        single_data=single_data,
        system_xyz=system_xyz,
        n_molecules=n_molecules,
        n_atoms_per_mol=structure.n_atoms,
        box_length=box_length,
        forcefield_name=forcefield.name,
    )

    result = {
        "smiles": smiles,
        "name": name,
        "workdir": str(workdir),
        "n_molecules": n_molecules,
        "n_atoms_per_molecule": structure.n_atoms,
        "molecular_weight_g_mol": structure.molecular_weight,
        "density_init_g_cm3": density,
        "packmol_seed": packmol_seed,
        "packmol_max_attempts": packmol_max_attempts,
        "packmol_seed_step": packmol_seed_step,
        "packmol_nloop": packmol_nloop,
        "packmol_full_box": packmol_full_box,
        "packmol_margin_A": packmol_margin,
        "packmol_strict": packmol_strict,
        "box_length_A": box_length,
        "single_molecule_data": str(single_data),
        "system_xyz": str(system_xyz),
        "system_data": str(system_data),
    }

    report_file = workdir / "build_system_report.json"
    with open(report_file, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print("=" * 70)
    print("Build finished")
    print("=" * 70)
    print(f"  single molecule data : {single_data}")
    print(f"  system xyz           : {system_xyz}")
    print(f"  system data          : {system_data}")
    print(f"  report               : {report_file}")

    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build system.data from SMILES")
    parser.add_argument("--smiles", required=True, help="SMILES string")
    parser.add_argument("--name", required=True, help="Molecule name")
    parser.add_argument("--workdir", required=True, help="Working directory")
    parser.add_argument("--n_mol", required=True, type=int, help="Number of molecules")
    parser.add_argument("--density", required=True, type=float, help="Initial packing density (g/cm^3)")
    parser.add_argument("--packmol_seed", type=int, default=192911, help="Packmol random seed")
    parser.add_argument("--packmol_max_attempts", type=int, default=3, help="Packmol maximum attempts")
    parser.add_argument("--packmol_seed_step", type=int, default=97, help="Packmol seed increment per retry")
    parser.add_argument("--packmol_nloop", type=int, default=None, help="Optional Packmol nloop")
    parser.add_argument("--packmol_full_box", action="store_true", help="Use inside box instead of inside cube")
    parser.add_argument("--packmol_margin", type=float, default=None, help="Margin for inside cube mode (A)")
    parser.add_argument("--packmol_strict", action="store_true", help="Require perfect packing only")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_system_data_only(
        smiles=args.smiles,
        name=args.name,
        workdir=Path(args.workdir),
        n_molecules=args.n_mol,
        density=args.density,
        packmol_seed=args.packmol_seed,
        packmol_max_attempts=args.packmol_max_attempts,
        packmol_seed_step=args.packmol_seed_step,
        packmol_nloop=args.packmol_nloop,
        packmol_full_box=args.packmol_full_box,
        packmol_margin=args.packmol_margin,
        packmol_strict=args.packmol_strict,
    )


if __name__ == "__main__":
    main()
