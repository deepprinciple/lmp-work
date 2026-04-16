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
import datetime as dt
import hashlib
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
    from workflow.config import (
        REPO_ROOT,
        apply_default_lammps_ani_environment,
        apply_md_viscosity_path_resolution,
        get_required,
        get_section,
        load_yaml,
    )
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

def _build_lammps_env(lammps_cfg: dict[str, Any]) -> dict[str, str]:
    """Build the runtime environment for ``lmp_mpi``.

    ``LAMMPS_PLUGIN_PATH`` priority:
    1. Current process environment
    2. ``lammps.env.LAMMPS_PLUGIN_PATH`` from YAML
    3. ``$LAMMPS_ANI_ROOT/build`` when ``ani_plugin.so`` is present
    """
    env = os.environ.copy()
    extra_env = lammps_cfg.get("env") or {}
    if not isinstance(extra_env, dict):
        raise TypeError("lammps.env must be a mapping when provided")

    extra = {str(key): str(value) for key, value in extra_env.items()}
    plugin_yaml = os.path.expandvars(extra.pop("LAMMPS_PLUGIN_PATH", "").strip())

    for key, value in extra.items():
        env[key] = value

    plugin_os = os.environ.get("LAMMPS_PLUGIN_PATH", "").strip()
    if plugin_os:
        env["LAMMPS_PLUGIN_PATH"] = plugin_os
        return env

    if plugin_yaml and "${" not in plugin_yaml:
        env["LAMMPS_PLUGIN_PATH"] = plugin_yaml
        return env

    ani_root = env.get("LAMMPS_ANI_ROOT", "").strip()
    if ani_root:
        plugin_dir = Path(ani_root) / "build"
        if (plugin_dir / "ani_plugin.so").is_file():
            env["LAMMPS_PLUGIN_PATH"] = str(plugin_dir.resolve())

    return env


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
    use_hwthread_cpus = bool(lammps_cfg.get("use_hwthread_cpus", False))
    oversubscribe = bool(lammps_cfg.get("oversubscribe", False))

    cmd: list[str] = [mpi_command]
    if allow_root:
        cmd.append("--allow-run-as-root")
    if use_hwthread_cpus:
        cmd.append("--use-hwthread-cpus")
    if oversubscribe:
        cmd.append("--oversubscribe")
    cmd += ["-np", str(mpi_ranks), executable, "-in", input_file.name, "-log", log_file]

    print("=" * 70)
    print(f"Run LAMMPS  [{input_file.name}]")
    print("=" * 70)
    print(f"  cwd    : {workdir}")
    print(f"  cmd    : {' '.join(cmd)}")

    env = _build_lammps_env(lammps_cfg)
    subprocess.run(cmd, cwd=workdir, env=env, check=True)


def _validate_npt_density_trace(
    workdir: Path,
    sim_cfg: dict[str, Any],
    model_cfg: dict[str, Any],
) -> None:
    """Check that the tail of the NPT density trace is close to target."""
    if bool(sim_cfg.get("equil_skip_density_check", False)):
        print("  [skip] NPT density check (equil_skip_density_check: true)")
        return

    trace_name = str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat"))
    trace_path = workdir / trace_name
    if not trace_path.exists():
        raise FileNotFoundError(f"NPT density trace missing: {trace_path}")

    target_density = float(model_cfg.get("density_g_cm3", 0.0))
    if target_density <= 0:
        raise ValueError("model.density_g_cm3 must be > 0 for NPT density check")

    npt_steps = int(sim_cfg.get("npt_steps", 1))
    last_steps = int(sim_cfg.get("equil_density_check_last_steps", 50000))
    rel_tol = float(sim_cfg.get("equil_density_rel_tol", 0.05))
    min_step = max(1, npt_steps - last_steps + 1)

    samples: list[tuple[int, float]] = []
    for line in trace_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            step = int(float(parts[0]))
            rho = float(parts[1])
        except ValueError:
            continue
        if step >= min_step:
            samples.append((step, rho))

    if len(samples) < 3:
        raise ValueError(
            f"Not enough NPT density samples in {trace_path} for steps >= {min_step} "
            f"(got {len(samples)} rows)"
        )

    densities = [rho for _, rho in samples]
    mean_density = sum(densities) / len(densities)
    rel_err = abs(mean_density - target_density) / target_density
    if rel_err > rel_tol:
        raise RuntimeError(
            f"NPT density check failed: mean rho={mean_density:.6f} g/cm^3 vs target "
            f"{target_density:.6f} g/cm^3 (relative error {rel_err:.4f} > {rel_tol}) "
            f"over steps >= {min_step} (n={len(densities)} samples)"
        )

    print(
        f"  NPT density check OK: mean rho={mean_density:.4f} g/cm^3 vs target "
        f"{target_density:.4f} g/cm^3 (|Δ|/ρ <= {rel_tol}; {len(densities)} samples)"
    )


def _state_paths(workdir: Path) -> tuple[Path, Path]:
    """Return paths for stage marker and machine-readable state."""
    return workdir / "stage.done", workdir / "state.json"


def _utc_now_iso() -> str:
    """Return a compact UTC timestamp for state tracking."""
    return dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z")


def _load_state(workdir: Path) -> dict[str, Any]:
    """Load existing ``state.json`` or return a fresh state mapping."""
    _, state_file = _state_paths(workdir)
    if state_file.exists():
        try:
            data = json.loads(state_file.read_text())
            if isinstance(data, dict):
                return data
        except Exception:
            pass
    return {
        "status": "running",
        "stages": {},
    }


def _save_state(workdir: Path, state: dict[str, Any]) -> None:
    """Persist stage state to ``state.json``."""
    _, state_file = _state_paths(workdir)
    state_file.write_text(json.dumps(state, indent=2, sort_keys=True))


def _mark_stage(
    workdir: Path,
    state: dict[str, Any],
    stage: str,
    *,
    status: str,
    extra: dict[str, Any] | None = None,
) -> None:
    """Update ``stage.done`` + ``state.json`` for a stage transition."""
    stage_file, _ = _state_paths(workdir)
    entry: dict[str, Any] = {
        "status": status,
        "time_utc": _utc_now_iso(),
    }
    if extra:
        entry.update(extra)
    state.setdefault("stages", {})[stage] = entry
    state["last_stage"] = stage
    _save_state(workdir, state)
    with stage_file.open("a", encoding="utf-8") as fh:
        fh.write(f"{entry['time_utc']} {stage} {status}\n")


def _workflow_paths(cfg: dict[str, Any], workdir: Path) -> dict[str, Path]:
    """Collect key workflow file paths under *workdir*."""
    sim_cfg = get_section(cfg, "simulation")
    analysis_cfg = get_section(cfg, "analysis")
    return {
        "data_file": workdir / "system.data",
        "equil_input": workdir / str(sim_cfg.get("equil_input_filename", "in.equil.lammps")),
        "gk_input": workdir / str(sim_cfg.get("gk_input_filename", "in.gk.lammps")),
        "equil_restart": workdir / str(sim_cfg.get("equil_restart_file", "equil_nvt.restart")),
        "npt_thermo": workdir / str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat")),
        "acf_file": workdir / str(sim_cfg.get("acf_file", "stress_acf.dat")),
        "thermo_file": workdir / str(sim_cfg.get("thermo_file", "gk_thermo.dat")),
        "summary_file": workdir / str(analysis_cfg.get("summary_file", "viscosity_summary.json")),
        "plot_file": workdir / str(analysis_cfg.get("plot_file", "viscosity_analysis.png")),
    }


def _relative_paths(paths: list[Path], workdir: Path) -> list[str]:
    """Return paths relative to *workdir* when possible."""
    out: list[str] = []
    for path in paths:
        try:
            out.append(str(path.relative_to(workdir)))
        except ValueError:
            out.append(str(path))
    return out


def _paths_exist(paths: list[Path]) -> bool:
    """Return True when every path in *paths* exists."""
    return all(path.exists() for path in paths)


def _outputs_up_to_date(outputs: list[Path], inputs: list[Path]) -> bool:
    """Return True when outputs exist and are newer than all inputs."""
    if not outputs or not inputs:
        return False
    if not _paths_exist(outputs) or not _paths_exist(inputs):
        return False
    oldest_output = min(path.stat().st_mtime for path in outputs)
    newest_input = max(path.stat().st_mtime for path in inputs)
    return oldest_output >= newest_input


def _stage_signature(cfg: dict[str, Any], stage: str) -> str | None:
    """Build a compact hash for the config subset relevant to *stage*."""
    if stage == "build_system":
        payload: dict[str, Any] = {
            "case": get_section(cfg, "case"),
            "structure": get_section(cfg, "structure"),
            "model": get_section(cfg, "model"),
        }
    elif stage == "write_input":
        payload = {
            "ani": get_section(cfg, "ani"),
            "simulation": get_section(cfg, "simulation"),
        }
    elif stage == "analyze":
        sim_cfg = get_section(cfg, "simulation")
        payload = {
            "analysis": get_section(cfg, "analysis"),
            "simulation": {
                "temperature_K": sim_cfg.get("temperature_K"),
                "acf_file": sim_cfg.get("acf_file"),
                "thermo_file": sim_cfg.get("thermo_file"),
            },
        }
    else:
        return None

    raw = json.dumps(payload, sort_keys=True, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]


def _can_reuse_stage(
    state: dict[str, Any],
    stage: str,
    *,
    outputs: list[Path],
    inputs: list[Path] | None = None,
    signature: str | None = None,
    allow_artifact_fallback: bool = False,
) -> bool:
    """Return True when an existing stage output set is safe to reuse."""
    if not _paths_exist(outputs):
        return False
    if inputs is not None and not _outputs_up_to_date(outputs, inputs):
        return False
    if signature is None:
        return True

    stages = state.get("stages")
    entry = stages.get(stage) if isinstance(stages, dict) else None
    if not isinstance(entry, dict):
        return allow_artifact_fallback

    recorded = entry.get("signature")
    if recorded is None:
        return allow_artifact_fallback
    return recorded == signature


def _print_reuse(label: str, workdir: Path, outputs: list[Path]) -> None:
    """Log that a stage is reusing up-to-date artifacts."""
    files = ", ".join(_relative_paths(outputs, workdir))
    print(f"  [reuse] {label} – using {files}")


def _rollup_status(*statuses: str) -> str:
    """Collapse sub-step statuses into one stage status."""
    active = [status for status in statuses if status != "skipped"]
    if not active:
        return "skipped"
    if all(status == "reused" for status in active):
        return "reused"
    return "done"


def _load_summary_file(summary_file: Path) -> dict[str, Any]:
    """Load a JSON summary produced by the analysis stage."""
    data = json.loads(summary_file.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Summary file must contain a JSON object: {summary_file}")
    return data


def _print_viscosity_result(summary: dict[str, Any]) -> None:
    """Pretty-print the viscosity summary."""
    eta = float(summary["eta_mPas"])
    comps = [float(x) for x in summary["eta_components_mPas"]]
    print(f"\n  ┌─ Result ─────────────────────────────────────────")
    print(f"  │  η = {eta:.3f} mPa·s (cP)")
    print(f"  │  components: xy={comps[0]:.3f}  xz={comps[1]:.3f}  yz={comps[2]:.3f} mPa·s")
    print(f"  │  plateau   : {float(summary['plateau_start_ps']):.1f} – {float(summary['plateau_end_ps']):.1f} ps")
    print(f"  └──────────────────────────────────────────────────")


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
    model_cfg = get_section(cfg, "model")
    struct_cfg = get_section(cfg, "structure")

    smiles = get_required(cfg, "case", "smiles")
    name = get_required(cfg, "case", "name")

    print("\n── Stage 1: build system ──────────────────────────────────────")

    # 1a. Generate 3-D single-molecule geometry.
    mol = MoleculeStructure(smiles=smiles, name=name)
    mol.generate(optimize=True)
    mol_xyz = workdir / f"{name}.xyz"
    mol.save_xyz(mol_xyz)
    for fmt in struct_cfg.get("output_formats", []):
        fmt_name = str(fmt).strip().lower()
        if fmt_name == "xyz":
            continue
        save_fn = getattr(mol, f"save_{fmt_name}", None)
        if save_fn is None:
            raise ValueError(f"Unsupported structure output format: {fmt}")
        save_fn(workdir / f"{name}.{fmt_name}")

    # 1b. Pack molecules into a box.
    n_mol = int(model_cfg.get("n_molecules", 500))
    density = float(model_cfg.get("density_g_cm3", 1.0))
    packmol_density_scale = float(model_cfg.get("packmol_density_scale", 0.85))
    if packmol_density_scale <= 0:
        raise ValueError("model.packmol_density_scale must be > 0")
    packer = PackmolBuilder(workdir=workdir)
    system_xyz, box_lengths = packer.build_box(
        mol_xyz=mol_xyz,
        n_molecules=n_mol,
        density=density * packmol_density_scale,
        mol_weight=mol.molecular_weight,
        tolerance=float(model_cfg.get("packmol_tolerance_A", 2.0)),
        seed=int(model_cfg.get("packmol_seed", 12345)),
        nloop=int(model_cfg.get("packmol_nloop", 200)),
        max_attempts=int(model_cfg.get("packmol_max_attempts", 5)),
        seed_step=int(model_cfg.get("packmol_seed_step", 17)),
        full_box=bool(model_cfg.get("packmol_full_box", False)),
        margin=model_cfg.get("packmol_margin_A"),
        allow_imperfect=not bool(model_cfg.get("packmol_strict", False)),
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
        npt_steps=int(sim_cfg.get("npt_steps", 200000)),
        npt_thermo_file=str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat")),
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
    *,
    run_equil: bool = True,
    run_gk: bool = True,
) -> dict[str, Any]:
    """Run equilibration and/or GK production."""
    lammps_cfg = get_section(cfg, "lammps")
    model_cfg = get_section(cfg, "model")
    sim_cfg = get_section(cfg, "simulation")
    paths = _workflow_paths(cfg, workdir)
    result: dict[str, Any] = {
        "run_equil": run_equil,
        "run_gk": run_gk,
        "equil_outputs": _relative_paths([paths["equil_restart"], paths["npt_thermo"]], workdir),
        "gk_outputs": _relative_paths([paths["acf_file"], paths["thermo_file"]], workdir),
    }

    print("\n── Stage 3: run LAMMPS ─────────────────────────────────────────")

    if run_equil:
        equil_outputs = [paths["equil_restart"], paths["npt_thermo"]]
        if _outputs_up_to_date(equil_outputs, [equil_path]):
            _print_reuse("equilibration run", workdir, equil_outputs)
            result["equil_status"] = "reused"
        else:
            _run_lammps(
                workdir,
                lammps_cfg,
                equil_path,
                str(sim_cfg.get("equil_log_file", "run_equil.log")),
            )
            result["equil_status"] = "done"
        _validate_npt_density_trace(workdir, sim_cfg, model_cfg)
    else:
        print("  [skip] equilibration run")
        result["equil_status"] = "skipped"

    if run_gk:
        gk_outputs = [paths["acf_file"], paths["thermo_file"]]
        if _outputs_up_to_date(gk_outputs, [gk_path, paths["equil_restart"]]):
            _print_reuse("Green-Kubo run", workdir, gk_outputs)
            result["gk_status"] = "reused"
        else:
            _run_lammps(
                workdir,
                lammps_cfg,
                gk_path,
                str(sim_cfg.get("gk_log_file", "run_gk.log")),
            )
            result["gk_status"] = "done"
    else:
        print("  [skip] Green-Kubo run")
        result["gk_status"] = "skipped"

    return result


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

    _print_viscosity_result(summary)
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

    config_path = Path(args.config).resolve()
    cfg = load_yaml(config_path)
    lammps_cfg = cfg.get("lammps")
    if not isinstance(lammps_cfg, dict):
        lammps_cfg = {}
    if lammps_cfg.get("auto_environment", True):
        apply_default_lammps_ani_environment()
    apply_md_viscosity_path_resolution(cfg, REPO_ROOT)
    run_flags = get_section(cfg, "run")

    workdir = Path(get_required(cfg, "case", "workdir")).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    paths = _workflow_paths(cfg, workdir)
    print(f"\nWorkdir: {workdir}")
    state = _load_state(workdir)
    state.update(
        {
            "status": "running",
            "config_path": str(config_path),
            "started_at_utc": _utc_now_iso(),
        }
    )
    state.pop("error", None)
    state.pop("ended_at_utc", None)
    state.pop("failed_stage", None)
    _save_state(workdir, state)

    # Track paths across stages so skipped stages can still read outputs.
    data_file: Path | None = None
    equil_path: Path | None = None
    gk_path: Path | None = None
    current_stage: str | None = None

    try:
        if run_flags.get("build_system", True):
            current_stage = "build_system"
            build_signature = _stage_signature(cfg, "build_system")
            if _can_reuse_stage(
                state,
                "build_system",
                outputs=[paths["data_file"]],
                signature=build_signature,
                allow_artifact_fallback=True,
            ):
                data_file = paths["data_file"]
                _print_reuse("build_system", workdir, [data_file])
                _mark_stage(
                    workdir,
                    state,
                    "build_system",
                    status="reused",
                    extra={
                        "signature": build_signature,
                        "outputs": _relative_paths([data_file], workdir),
                    },
                )
            else:
                data_file, _ = stage_build_system(cfg, workdir)
                _mark_stage(
                    workdir,
                    state,
                    "build_system",
                    status="done",
                    extra={
                        "signature": build_signature,
                        "outputs": _relative_paths([data_file], workdir),
                    },
                )
        else:
            data_file = paths["data_file"]
            print(f"\n[skip] build_system – expecting {data_file}")
            _mark_stage(workdir, state, "build_system", status="skipped")

        if run_flags.get("write_input", True):
            current_stage = "write_input"
            assert data_file is not None
            write_signature = _stage_signature(cfg, "write_input")
            write_outputs = [paths["equil_input"], paths["gk_input"]]
            if _can_reuse_stage(
                state,
                "write_input",
                outputs=write_outputs,
                inputs=[data_file],
                signature=write_signature,
                allow_artifact_fallback=True,
            ):
                equil_path = paths["equil_input"]
                gk_path = paths["gk_input"]
                _print_reuse("write_input", workdir, write_outputs)
                _mark_stage(
                    workdir,
                    state,
                    "write_input",
                    status="reused",
                    extra={
                        "signature": write_signature,
                        "outputs": _relative_paths(write_outputs, workdir),
                    },
                )
            else:
                equil_path, gk_path = stage_write_input(cfg, workdir, data_file)
                _mark_stage(
                    workdir,
                    state,
                    "write_input",
                    status="done",
                    extra={
                        "signature": write_signature,
                        "outputs": _relative_paths([equil_path, gk_path], workdir),
                    },
                )
        else:
            equil_path = paths["equil_input"]
            gk_path = paths["gk_input"]
            print(f"\n[skip] write_input – expecting {equil_path.name} / {gk_path.name}")
            _mark_stage(workdir, state, "write_input", status="skipped")

        if run_flags.get("run_lammps", True):
            current_stage = "run_lammps"
            assert equil_path is not None and gk_path is not None
            run_equil = bool(run_flags.get("run_equil", True))
            run_gk = bool(run_flags.get("run_gk", True))
            run_info = stage_run_lammps(
                cfg,
                workdir,
                equil_path,
                gk_path,
                run_equil=run_equil,
                run_gk=run_gk,
            )
            _mark_stage(
                workdir,
                state,
                "run_lammps",
                status=_rollup_status(
                    str(run_info.get("equil_status", "skipped")),
                    str(run_info.get("gk_status", "skipped")),
                ),
                extra=run_info,
            )
        else:
            print("\n[skip] run_lammps")
            _mark_stage(workdir, state, "run_lammps", status="skipped")

        if run_flags.get("analyze", True):
            current_stage = "analyze"
            analyze_signature = _stage_signature(cfg, "analyze")
            analyze_outputs = [paths["summary_file"], paths["plot_file"]]
            if _can_reuse_stage(
                state,
                "analyze",
                outputs=analyze_outputs,
                inputs=[paths["acf_file"], paths["thermo_file"]],
                signature=analyze_signature,
            ):
                try:
                    summary = _load_summary_file(paths["summary_file"])
                except Exception:
                    summary = stage_analyze(cfg, workdir)
                    analyze_status = "done"
                else:
                    _print_reuse("analyze", workdir, analyze_outputs)
                    _print_viscosity_result(summary)
                    analyze_status = "reused"
            else:
                summary = stage_analyze(cfg, workdir)
                analyze_status = "done"
            _mark_stage(
                workdir,
                state,
                "analyze",
                status=analyze_status,
                extra={
                    "signature": analyze_signature,
                    "eta_mPas": summary.get("eta_mPas"),
                    "outputs": _relative_paths(analyze_outputs, workdir),
                },
            )
        else:
            print("\n[skip] analyze")
            _mark_stage(workdir, state, "analyze", status="skipped")

        state["status"] = "succeeded"
        state["ended_at_utc"] = _utc_now_iso()
        _save_state(workdir, state)
        print("\nDone.")
    except Exception as exc:
        state["status"] = "failed"
        if current_stage is not None:
            state["failed_stage"] = current_stage
        state["error"] = {
            "type": exc.__class__.__name__,
            "message": str(exc),
        }
        state["ended_at_utc"] = _utc_now_iso()
        _save_state(workdir, state)
        raise


if __name__ == "__main__":
    main()
