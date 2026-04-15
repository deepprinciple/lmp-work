"""ANI Green-Kubo viscosity workflow entry point.

Usage
-----
    python md_viscosity.py --config configs/viscosity.yaml

Workflow stages
---------------
1. ``build_system``  – SMILES → 3-D geometry → Packmol box → LAMMPS data file
2. ``write_input``   – generate equilibration and Green-Kubo LAMMPS input files
3. ``run_lammps``    – run equilibration then GK production with ANI pair style
4. ``analyze``       – integrate stress ACF → shear viscosity + convergence plot

Any stage can be skipped by setting its flag to ``false`` in the config.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

# ── project-local imports (all relative to lmp-work root) ─────────────────
try:
    from core.data_builder import AniDataBuilder
    from core.packing import PackmolBuilder
    from core.structure import MoleculeStructure
    from workflow.config import get_required, get_section, load_yaml
    from workflow.input_ani_viscosity import build_equil_input, build_gk_input
    from analysis.viscosity import analyze_viscosity_file, load_volume_from_thermo
except ImportError as exc:  # pragma: no cover
    sys.exit(
        f"Import error: {exc}\n"
        "Run this script from the lmp-work project root:\n"
        "  python md_viscosity.py --config configs/viscosity.yaml"
    )


# ──────────────────────────────────────────────────────────────────────────
# LAMMPS execution helper
# ──────────────────────────────────────────────────────────────────────────

def _run_lammps(
    workdir: Path,
    lammps_cfg: dict[str, Any],
    input_file: Path,
    log_file: str,
) -> None:
    """Execute LAMMPS for an ANI simulation.

    The ANI pair style manages its own GPU via PyTorch, so we do NOT pass
    LAMMPS GPU-package flags (-sf gpu).  The command is simply::

        mpirun [--allow-run-as-root] -np N lmp_mpi -in input.in -log log.out
    """
    executable = lammps_cfg.get("executable", "lmp_mpi")
    mpi_command = lammps_cfg.get("mpi_command", "mpirun")
    mpi_ranks = int(lammps_cfg.get("mpi_ranks", 1))
    allow_root = bool(lammps_cfg.get("allow_run_as_root", True))

    cmd: list[str] = [mpi_command]
    if allow_root:
        cmd.append("--allow-run-as-root")
    cmd += ["-np", str(mpi_ranks), executable, "-in", input_file.name, "-log", log_file]

    print("=" * 70)
    print(f"Run LAMMPS  [{input_file.name}]")
    print("=" * 70)
    print(f"  cwd    : {workdir}")
    print(f"  cmd    : {' '.join(cmd)}")

    env = os.environ.copy()
    subprocess.run(cmd, cwd=workdir, env=env, check=True)


# ──────────────────────────────────────────────────────────────────────────
# Stage helpers
# ──────────────────────────────────────────────────────────────────────────

def stage_build_system(
    cfg: dict[str, Any],
    workdir: Path,
) -> tuple[Path, tuple[float, float, float]]:
    """SMILES → geometry → Packmol box → LAMMPS data file.

    Returns
    -------
    (data_file_path, (Lx, Ly, Lz))
    """
    case_cfg = get_section(cfg, "case")
    model_cfg = get_section(cfg, "model")
    struct_cfg = get_section(cfg, "structure")

    smiles = get_required(cfg, "case", "smiles")
    name = get_required(cfg, "case", "name")

    print("\n── Stage 1: build system ──────────────────────────────────────")

    # 1a. Generate 3-D single-molecule geometry.
    mol = MoleculeStructure(smiles=smiles, name=name)
    mol.generate(optimize=True)
    for fmt in struct_cfg.get("output_formats", ["xyz"]):
        getattr(mol, f"save_{fmt}")(workdir / f"{name}.{fmt}")

    # 1b. Pack molecules into a box.
    n_mol = int(model_cfg.get("n_molecules", 500))
    density = float(model_cfg.get("density_g_cm3", 1.0))
    packer = PackmolBuilder(workdir=workdir)
    system_xyz, box_lengths = packer.build_box(
        mol_xyz=str(workdir / f"{name}.xyz"),
        n_molecules=n_mol,
        density=density,
        mol_weight=mol.molecular_weight,
        tolerance=float(model_cfg.get("packmol_tolerance_A", 2.0)),
        seed=int(model_cfg.get("packmol_seed", 12345)),
        nloop=int(model_cfg.get("packmol_nloop", 200)),
        max_attempts=int(model_cfg.get("packmol_max_attempts", 5)),
        seed_step=int(model_cfg.get("packmol_seed_step", 17)),
        full_box=bool(model_cfg.get("packmol_full_box", False)),
        margin=model_cfg.get("packmol_margin_A"),
        allow_imperfect=not bool(model_cfg.get("packmol_strict", False)),
        density_scale=float(model_cfg.get("packmol_density_scale", 0.85)),
    )

    # 1c. Convert XYZ → LAMMPS atomic data file.
    data_file = workdir / "system.data"
    builder = AniDataBuilder()
    builder.build(system_xyz, box_lengths, data_file)

    Lx, Ly, Lz = box_lengths
    V = Lx * Ly * Lz
    print(f"  system.data : {n_mol} molecules  |  box {Lx:.1f}×{Ly:.1f}×{Lz:.1f} Å  |  V={V:.0f} Å³")
    return data_file, box_lengths


def stage_write_input(
    cfg: dict[str, Any],
    workdir: Path,
    data_file: Path,
) -> tuple[Path, Path]:
    """Write equilibration and GK LAMMPS input files.

    Returns
    -------
    (equil_input_path, gk_input_path)
    """
    sim_cfg = get_section(cfg, "simulation")
    ani_cfg = get_section(cfg, "ani")

    print("\n── Stage 2: write LAMMPS inputs ───────────────────────────────")

    common = dict(
        model_file=get_required(cfg, "ani", "model_file"),
        temperature=float(sim_cfg.get("temperature_K", 300.0)),
        timestep_fs=float(sim_cfg.get("timestep_fs", 0.5)),
        device=str(ani_cfg.get("device", "cuda")),
        num_models=int(ani_cfg.get("num_models", 1)),
    )

    restart_file = str(sim_cfg.get("equil_restart_file", "equil_nvt.restart"))

    equil_text = build_equil_input(
        data_file=data_file.name,
        nvt_steps=int(sim_cfg.get("nvt_steps", 40000)),
        seed=int(sim_cfg.get("seed", 12345)),
        restart_file=restart_file,
        **common,
    )
    equil_path = workdir / str(sim_cfg.get("equil_input_filename", "in.equil.lammps"))
    equil_path.write_text(equil_text)
    print(f"  wrote: {equil_path.name}")

    gk_text = build_gk_input(
        restart_file=restart_file,
        prod_steps=int(sim_cfg.get("prod_steps", 2000000)),
        corr_length=int(sim_cfg.get("corr_length", 40000)),
        sample_every=int(sim_cfg.get("sample_every", 1)),
        acf_file=str(sim_cfg.get("acf_file", "stress_acf.dat")),
        thermo_file=str(sim_cfg.get("thermo_file", "gk_thermo.dat")),
        **common,
    )
    gk_path = workdir / str(sim_cfg.get("gk_input_filename", "in.gk.lammps"))
    gk_path.write_text(gk_text)
    print(f"  wrote: {gk_path.name}")

    return equil_path, gk_path


def stage_run_lammps(
    cfg: dict[str, Any],
    workdir: Path,
    equil_path: Path,
    gk_path: Path,
) -> None:
    """Run equilibration then GK production."""
    lammps_cfg = get_section(cfg, "lammps")
    sim_cfg = get_section(cfg, "simulation")

    print("\n── Stage 3: run LAMMPS ─────────────────────────────────────────")

    _run_lammps(workdir, lammps_cfg, equil_path,
                str(sim_cfg.get("equil_log_file", "run_equil.log")))
    _run_lammps(workdir, lammps_cfg, gk_path,
                str(sim_cfg.get("gk_log_file", "run_gk.log")))


def stage_analyze(
    cfg: dict[str, Any],
    workdir: Path,
) -> dict[str, Any]:
    """Integrate stress ACF → shear viscosity + convergence plot."""
    sim_cfg = get_section(cfg, "simulation")
    analysis_cfg = get_section(cfg, "analysis")

    print("\n── Stage 4: analyze ────────────────────────────────────────────")

    acf_file = workdir / str(sim_cfg.get("acf_file", "stress_acf.dat"))
    thermo_file = workdir / str(sim_cfg.get("thermo_file", "gk_thermo.dat"))

    volume_a3 = load_volume_from_thermo(thermo_file)
    print(f"  NVE mean volume : {volume_a3:.1f} Å³")

    summary = analyze_viscosity_file(
        acf_file,
        temperature_k=float(get_section(cfg, "simulation").get("temperature_K", 300.0)),
        volume_a3=volume_a3,
        t_max_ps=float(analysis_cfg.get("t_max_ps", 20.0)),
        plateau_start_ps=float(analysis_cfg.get("plateau_start_ps", 5.0)),
        smooth_window=int(analysis_cfg.get("smooth_window", 11)),
        json_out=str(analysis_cfg.get("summary_file", "viscosity_summary.json")),
        plot_out=str(analysis_cfg.get("plot_file", "viscosity_analysis.png")),
    )

    eta = summary["eta_mPas"]
    comps = summary["eta_components_mPas"]
    print(f"\n  ┌─ Result ─────────────────────────────────────────")
    print(f"  │  η = {eta:.3f} mPa·s (cP)")
    print(f"  │  components: xy={comps[0]:.3f}  xz={comps[1]:.3f}  yz={comps[2]:.3f} mPa·s")
    print(f"  │  plateau   : {summary['plateau_start_ps']:.1f} – {summary['plateau_end_ps']:.1f} ps")
    print(f"  └──────────────────────────────────────────────────")

    return summary


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="ANI Green-Kubo viscosity workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args(argv)

    cfg = load_yaml(args.config)
    run_flags = get_section(cfg, "run")

    workdir = Path(get_required(cfg, "case", "workdir")).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    print(f"\nWorkdir: {workdir}")

    # Track paths across stages so skipped stages can still read outputs.
    data_file: Path | None = None
    equil_path: Path | None = None
    gk_path: Path | None = None

    if run_flags.get("build_system", True):
        data_file, _ = stage_build_system(cfg, workdir)
    else:
        data_file = workdir / "system.data"
        print(f"\n[skip] build_system – expecting {data_file}")

    if run_flags.get("write_input", True):
        assert data_file is not None
        equil_path, gk_path = stage_write_input(cfg, workdir, data_file)
    else:
        sim_cfg = get_section(cfg, "simulation")
        equil_path = workdir / str(sim_cfg.get("equil_input_filename", "in.equil.lammps"))
        gk_path = workdir / str(sim_cfg.get("gk_input_filename", "in.gk.lammps"))
        print(f"\n[skip] write_input – expecting {equil_path.name} / {gk_path.name}")

    if run_flags.get("run_lammps", True):
        assert equil_path is not None and gk_path is not None
        stage_run_lammps(cfg, workdir, equil_path, gk_path)
    else:
        print("\n[skip] run_lammps")

    if run_flags.get("analyze", True):
        stage_analyze(cfg, workdir)
    else:
        print("\n[skip] analyze")

    print("\nDone.")


if __name__ == "__main__":
    main()
