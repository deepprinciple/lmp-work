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

# Import through the project root so the script works from this directory.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gk_workflow.core import MoleculeStructure, PackmolBuilder, LAMMPSBuilder
from gk_workflow.forcefields import LigParGen, OpenFF


def build_system_data_only(
    smiles: str,
    name: str,
    workdir: Path,
    n_molecules: int,
    density: float,
    *,
    packmol_density_scale: float = 0.85,
    packmol_tolerance: float = 2.0,
    packmol_seed: int = 192911,
    packmol_max_attempts: int = 3,
    packmol_seed_step: int = 97,
    packmol_nloop: int | None = None,
    packmol_full_box: bool = False,
    packmol_margin: float | None = None,
    packmol_strict: bool = False,
    box_aspect_ratio: tuple[float, float, float] | None = None,
    box_lengths: tuple[float, float, float] | None = None,
    structure_output_formats: tuple[str, ...] = ("xyz",),
    forcefield_engine: str = "openff",
    forcefield_version: str = "openff-2.0.0",
    forcefield_strict_stereo: bool = True,
    forcefield_charge_method: str = "am1bcc",
    forcefield_charge_fallback: str | None = "gasteiger",
    forcefield_charge_file: str | None = None,
    forcefield_use_cache: bool = True,
    forcefield_input_mode: str = "generated_pdb",
    forcefield_residue_name: str = "MOL",
    forcefield_n_optimizations: int = 0,
    forcefield_wrapper_script: str | None = None,
    forcefield_debug: bool = False,
) -> dict:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Build System Workflow (up to system.data)")
    print("=" * 70)
    print(f"  molecule : {name} ({smiles})")
    print(f"  workdir  : {workdir}")
    print(f"  n_mol    : {n_molecules}")
    print(f"  target density : {density:.6f} g/cm^3")
    print(f"  pack density   : {density * packmol_density_scale:.6f} g/cm^3")
    print("=" * 70 + "\n")

    if packmol_density_scale <= 0:
        raise ValueError("packmol_density_scale 必须 > 0")
    if packmol_tolerance <= 0:
        raise ValueError("packmol_tolerance 必须 > 0")

    # Step 1: 3D structure
    structure = MoleculeStructure(smiles, name)
    mol = structure.generate(optimize=True)
    mol_xyz = structure.save_xyz(workdir / f"{name}.xyz")
    saved_structures: dict[str, str] = {"xyz": str(mol_xyz)}

    for fmt in structure_output_formats:
        fmt_norm = str(fmt).strip().lower()
        if fmt_norm == "xyz":
            continue
        if fmt_norm == "pdb":
            saved_structures["pdb"] = str(structure.save_pdb(workdir / f"{name}.pdb"))
            continue
        if fmt_norm == "mol":
            saved_structures["mol"] = str(structure.save_mol(workdir / f"{name}.mol"))
            continue
        raise ValueError(f"Unsupported structure output format: {fmt}")

    # Step 2: Packmol
    packer = PackmolBuilder(workdir)
    system_xyz, box_lengths_out = packer.build_box(
        mol_xyz=mol_xyz,
        n_molecules=n_molecules,
        density=density * packmol_density_scale,
        mol_weight=structure.molecular_weight,
        tolerance=packmol_tolerance,
        seed=packmol_seed,
        max_attempts=packmol_max_attempts,
        seed_step=packmol_seed_step,
        nloop=packmol_nloop,
        full_box=packmol_full_box,
        margin=packmol_margin,
        allow_imperfect=not packmol_strict,
        box_aspect_ratio=box_aspect_ratio,
        box_lengths=box_lengths,
    )

    # Step 3: single-molecule parameterization
    engine = str(forcefield_engine).strip().lower()
    if engine == "openff":
        forcefield = OpenFF(
            workdir,
            version=forcefield_version,
            strict_stereo=forcefield_strict_stereo,
            charge_method=forcefield_charge_method,
            charge_fallback=forcefield_charge_fallback,
            charge_file=forcefield_charge_file,
            use_cache=forcefield_use_cache,
        )
    elif engine == "ligpargen":
        forcefield = LigParGen(
            workdir,
            charge_method=forcefield_charge_method,
            input_mode=forcefield_input_mode,
            residue_name=forcefield_residue_name,
            n_optimizations=forcefield_n_optimizations,
            wrapper_script=forcefield_wrapper_script,
            debug=forcefield_debug,
        )
    else:
        raise ValueError(f"Unsupported forcefield engine: {forcefield_engine}")
    single_data = forcefield.generate_parameters(mol, smiles, name)

    # Step 4: Build system.data
    builder = LAMMPSBuilder(workdir)
    system_data = builder.build_system(
        single_data=single_data,
        system_xyz=system_xyz,
        n_molecules=n_molecules,
        n_atoms_per_mol=structure.n_atoms,
        box_lengths=box_lengths_out,
        forcefield_name=forcefield.name,
    )

    result = {
        "smiles": smiles,
        "name": name,
        "workdir": str(workdir),
        "n_molecules": n_molecules,
        "n_atoms_per_molecule": structure.n_atoms,
        "molecular_weight_g_mol": structure.molecular_weight,
        "density_target_g_cm3": density,
        "density_init_g_cm3": density * packmol_density_scale,
        "packmol_density_scale": packmol_density_scale,
        "packmol_tolerance_A": packmol_tolerance,
        "packmol_seed": packmol_seed,
        "packmol_max_attempts": packmol_max_attempts,
        "packmol_seed_step": packmol_seed_step,
        "packmol_nloop": packmol_nloop,
        "packmol_full_box": packmol_full_box,
        "packmol_margin_A": packmol_margin,
        "packmol_strict": packmol_strict,
        "box_aspect_ratio": list(box_aspect_ratio) if box_aspect_ratio is not None else None,
        "box_lengths_A_requested": list(box_lengths) if box_lengths is not None else None,
        "structure_output_formats": list(structure_output_formats),
        "structure_files": saved_structures,
        "forcefield_engine": engine,
        "forcefield_version": forcefield_version,
        "forcefield_strict_stereo": forcefield_strict_stereo,
        "forcefield_charge_method": forcefield_charge_method,
        "forcefield_charge_fallback": forcefield_charge_fallback,
        "forcefield_charge_file": forcefield_charge_file,
        "forcefield_use_cache": forcefield_use_cache,
        "forcefield_input_mode": forcefield_input_mode,
        "forcefield_residue_name": forcefield_residue_name,
        "forcefield_n_optimizations": forcefield_n_optimizations,
        "forcefield_wrapper_script": forcefield_wrapper_script,
        "forcefield_debug": forcefield_debug,
        "box_length_A": box_lengths_out[0] if abs(box_lengths_out[0] - box_lengths_out[1]) < 1.0e-12 and abs(box_lengths_out[1] - box_lengths_out[2]) < 1.0e-12 else None,
        "box_lengths_A": list(box_lengths_out),
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
    parser.add_argument(
        "--packmol_density_scale",
        type=float,
        default=0.85,
        help="Scale factor applied to target density when building the initial Packmol box",
    )
    parser.add_argument(
        "--packmol_tolerance",
        type=float,
        default=2.0,
        help="Packmol tolerance in angstrom",
    )
    parser.add_argument("--packmol_seed", type=int, default=192911, help="Packmol random seed")
    parser.add_argument("--packmol_max_attempts", type=int, default=3, help="Packmol maximum attempts")
    parser.add_argument("--packmol_seed_step", type=int, default=97, help="Packmol seed increment per retry")
    parser.add_argument("--packmol_nloop", type=int, default=None, help="Optional Packmol nloop")
    parser.add_argument("--packmol_full_box", action="store_true", help="Use inside box instead of inside cube")
    parser.add_argument("--packmol_margin", type=float, default=None, help="Margin for inside cube mode (A)")
    parser.add_argument("--packmol_strict", action="store_true", help="Require perfect packing only")
    parser.add_argument(
        "--box_aspect_ratio",
        nargs=3,
        type=float,
        default=None,
        metavar=("RX", "RY", "RZ"),
        help="Orthorhombic box aspect ratio used to derive (Lx, Ly, Lz) from total target volume",
    )
    parser.add_argument(
        "--box_lengths_A",
        nargs=3,
        type=float,
        default=None,
        metavar=("LX", "LY", "LZ"),
        help="Explicit orthorhombic box lengths in angstrom",
    )
    parser.add_argument(
        "--structure_formats",
        nargs="*",
        default=("xyz",),
        help="Structure files to export alongside the internal XYZ (xyz/mol/pdb)",
    )
    parser.add_argument(
        "--forcefield_engine",
        default="openff",
        choices=("openff", "ligpargen"),
        help="Force-field engine used for single-molecule parameterization",
    )
    parser.add_argument("--forcefield_version", default="openff-2.0.0", help="OpenFF force-field version")
    parser.add_argument("--forcefield_charge_method", default="am1bcc", help="OpenFF partial charge method")
    parser.add_argument("--forcefield_charge_fallback", default="gasteiger", help="Fallback partial charge method")
    parser.add_argument("--forcefield_charge_file", default=None, help="External partial charges file (JSON or plain text)")
    parser.add_argument("--forcefield_no_cache", action="store_true", help="Disable OpenFF local cache")
    parser.add_argument("--forcefield_relaxed_stereo", action="store_true", help="Allow undefined stereochemistry")
    parser.add_argument(
        "--forcefield_input_mode",
        default="generated_pdb",
        help="LigParGen input mode (smiles or generated_pdb)",
    )
    parser.add_argument(
        "--forcefield_residue_name",
        default="MOL",
        help="LigParGen residue name for generated outputs",
    )
    parser.add_argument(
        "--forcefield_n_optimizations",
        type=int,
        default=0,
        help="LigParGen geometry optimization count",
    )
    parser.add_argument(
        "--forcefield_wrapper_script",
        default=None,
        help="LigParGen wrapper script path",
    )
    parser.add_argument(
        "--forcefield_debug",
        action="store_true",
        help="Keep LigParGen/BOSS intermediate files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_system_data_only(
        smiles=args.smiles,
        name=args.name,
        workdir=Path(args.workdir),
        n_molecules=args.n_mol,
        density=args.density,
        packmol_density_scale=args.packmol_density_scale,
        packmol_tolerance=args.packmol_tolerance,
        packmol_seed=args.packmol_seed,
        packmol_max_attempts=args.packmol_max_attempts,
        packmol_seed_step=args.packmol_seed_step,
        packmol_nloop=args.packmol_nloop,
        packmol_full_box=args.packmol_full_box,
        packmol_margin=args.packmol_margin,
        packmol_strict=args.packmol_strict,
        box_aspect_ratio=tuple(args.box_aspect_ratio) if args.box_aspect_ratio is not None else None,
        box_lengths=tuple(args.box_lengths_A) if args.box_lengths_A is not None else None,
        structure_output_formats=tuple(args.structure_formats),
        forcefield_engine=args.forcefield_engine,
        forcefield_version=args.forcefield_version,
        forcefield_strict_stereo=not args.forcefield_relaxed_stereo,
        forcefield_charge_method=args.forcefield_charge_method,
        forcefield_charge_fallback=args.forcefield_charge_fallback,
        forcefield_charge_file=args.forcefield_charge_file,
        forcefield_use_cache=not args.forcefield_no_cache,
        forcefield_input_mode=args.forcefield_input_mode,
        forcefield_residue_name=args.forcefield_residue_name,
        forcefield_n_optimizations=args.forcefield_n_optimizations,
        forcefield_wrapper_script=args.forcefield_wrapper_script,
        forcefield_debug=args.forcefield_debug,
    )


if __name__ == "__main__":
    main()
