"""BAMBOO Green-Kubo viscosity workflow entry point.

Usage
-----
    python md_viscosity.py --config configs/viscosity.yaml

Workflow stages
---------------
1. ``build_system``  – composition → 3-D geometries → Packmol box → BAMBOO data
2. ``write_input``   – generate equilibration and Green-Kubo LAMMPS input files
3. ``run_lammps``    – run equilibration then GK production with BAMBOO pair style
4. ``analyze``       – post-process stress traces → shear viscosity + convergence plot

Any stage can be skipped by setting its flag to ``false`` in the config.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

# ── project-local imports (all relative to lmp-work root) ─────────────────
try:
    from workflow.config import (
        REPO_ROOT,
        apply_md_viscosity_path_resolution,
        get_required,
        get_section,
        load_yaml,
    )
    from workflow.bamboo_system import build_bamboo_system, load_bamboo_build_report
    from workflow.input_viscosity import build_equil_input, build_gk_input
    from analysis.viscosity import (
        analyze_viscosity_file,
        analyze_viscosity_pressure_file,
        load_volume_from_thermo,
    )
except ImportError as exc:  # pragma: no cover
    sys.exit(
        f"Import error: {exc}\n"
        "Run this script from the lmp-work project root:\n"
        "  python md_viscosity.py --config configs/viscosity.yaml"
    )


# ──────────────────────────────────────────────────────────────────────────
# LAMMPS execution helper
# ──────────────────────────────────────────────────────────────────────────

_ENV_VAR_RE = re.compile(r"\$(?:\{([^}]+)\}|([A-Za-z_][A-Za-z0-9_]*))")


def _expand_runtime_env_value(value: str, env: dict[str, str]) -> str:
    """Expand ``$VAR`` / ``${VAR}`` against the environment being built."""

    def repl(match: re.Match[str]) -> str:
        name = match.group(1) or match.group(2) or ""
        return env.get(name, "")

    return _ENV_VAR_RE.sub(repl, os.path.expanduser(value))


def _build_lammps_env(lammps_cfg: dict[str, Any]) -> dict[str, str]:
    """Build the runtime environment for BAMBOO LAMMPS."""
    env = os.environ.copy()
    extra_env = lammps_cfg.get("env") or {}
    if not isinstance(extra_env, dict):
        raise TypeError("lammps.env must be a mapping when provided")

    for key, value in extra_env.items():
        env[str(key)] = _expand_runtime_env_value(str(value), env)
    return env


def _run_lammps(
    workdir: Path,
    lammps_cfg: dict[str, Any],
    input_file: Path,
    log_file: str,
) -> None:
    """Execute BAMBOO-enabled LAMMPS."""
    executable = lammps_cfg.get("executable", "lmp_mpi")
    mpi_command = lammps_cfg.get("mpi_command", "mpirun")
    mpi_ranks = int(lammps_cfg.get("mpi_ranks", 1))
    allow_root = bool(lammps_cfg.get("allow_run_as_root", True))
    use_hwthread_cpus = bool(lammps_cfg.get("use_hwthread_cpus", False))
    oversubscribe = bool(lammps_cfg.get("oversubscribe", False))
    kokkos = bool(lammps_cfg.get("kokkos", True))
    kokkos_gpus = int(lammps_cfg.get("kokkos_gpus", 1))
    suffix = str(lammps_cfg.get("suffix", "kk")).strip()
    extra_args = lammps_cfg.get("extra_args", [])
    if extra_args is None:
        extra_args = []
    if not isinstance(extra_args, list):
        raise TypeError("lammps.extra_args must be a list when provided")

    cmd: list[str] = [mpi_command]
    if allow_root:
        cmd.append("--allow-run-as-root")
    if use_hwthread_cpus:
        cmd.append("--use-hwthread-cpus")
    if oversubscribe:
        cmd.append("--oversubscribe")
    cmd += ["-np", str(mpi_ranks), executable]
    if kokkos:
        cmd += ["-k", "on", "g", str(kokkos_gpus)]
    if suffix:
        cmd += ["-sf", suffix]
    cmd += [str(arg) for arg in extra_args]
    cmd += ["-in", input_file.name, "-log", log_file]

    print("=" * 70)
    print(f"Run LAMMPS  [{input_file.name}]")
    print("=" * 70)
    print(f"  cwd    : {workdir}")
    print(f"  cmd    : {' '.join(cmd)}")

    env = _build_lammps_env(lammps_cfg)
    subprocess.run(cmd, cwd=workdir, env=env, check=True)


def _validate_equilibration_trace(
    workdir: Path,
    sim_cfg: dict[str, Any],
    model_cfg: dict[str, Any],
) -> None:
    """Check that the equilibration tail is close to target density/temperature."""
    trace_name = str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat"))
    trace_path = workdir / trace_name
    if not trace_path.exists():
        raise FileNotFoundError(f"Equilibration trace missing: {trace_path}")

    samples: list[tuple[int, float | None, float]] = []
    for line in trace_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            step = int(float(parts[0]))
        except ValueError:
            continue
        try:
            if len(parts) >= 3:
                temp = float(parts[1])
                rho = float(parts[2])
            else:
                temp = None
                rho = float(parts[1])
        except ValueError:
            continue
        samples.append((step, temp, rho))

    if len(samples) < 3:
        raise ValueError(f"Not enough equilibration samples in {trace_path} (got {len(samples)} rows)")

    nvt_steps = int(sim_cfg.get("nvt_steps", 0))
    npt_steps = int(sim_cfg.get("npt_steps", 1))
    final_nvt_steps = int(sim_cfg.get("final_nvt_steps", 0))

    if not bool(sim_cfg.get("equil_skip_density_check", False)):
        target_density = float(model_cfg.get("density_g_cm3", 0.0))
        if target_density <= 0:
            raise ValueError("model.density_g_cm3 must be > 0 for NPT density check")

        last_steps = int(sim_cfg.get("equil_density_check_last_steps", 50000))
        rel_tol = float(sim_cfg.get("equil_density_rel_tol", 0.05))
        npt_end = nvt_steps + npt_steps
        min_step = max(nvt_steps + 1, npt_end - last_steps + 1)
        densities = [rho for step, _, rho in samples if min_step <= step <= npt_end]
        if len(densities) < 3:
            raise ValueError(
                f"Not enough NPT density samples in {trace_path} for steps {min_step}..{npt_end} "
                f"(got {len(densities)} rows)"
            )

        mean_density = sum(densities) / len(densities)
        rel_err = abs(mean_density - target_density) / target_density
        if rel_err > rel_tol:
            raise RuntimeError(
                f"NPT density check failed: mean rho={mean_density:.6f} g/cm^3 vs target "
                f"{target_density:.6f} g/cm^3 (relative error {rel_err:.4f} > {rel_tol}) "
                f"over steps {min_step}..{npt_end} (n={len(densities)})"
            )

        print(
            f"  NPT density check OK: mean rho={mean_density:.4f} g/cm^3 vs target "
            f"{target_density:.4f} g/cm^3 (|Δ|/ρ <= {rel_tol}; {len(densities)} samples)"
        )
    else:
        print("  [skip] NPT density check (equil_skip_density_check: true)")

    if bool(sim_cfg.get("equil_skip_temperature_check", False)):
        print("  [skip] final NVT temperature check (equil_skip_temperature_check: true)")
        return
    if final_nvt_steps <= 0:
        print("  [skip] final NVT temperature check (final_nvt_steps <= 0)")
        return

    target_temp = float(sim_cfg.get("temperature_K", 0.0))
    if target_temp <= 0:
        raise ValueError("simulation.temperature_K must be > 0 for equilibration temperature check")

    last_steps = int(sim_cfg.get("equil_temp_check_last_steps", 10000))
    abs_tol = float(sim_cfg.get("equil_temp_abs_tol_K", 15.0))
    final_nvt_end = nvt_steps + npt_steps + final_nvt_steps
    min_step = max(nvt_steps + npt_steps + 1, final_nvt_end - last_steps + 1)
    temps = [temp for step, temp, _ in samples if temp is not None and min_step <= step <= final_nvt_end]
    if len(temps) < 3:
        raise ValueError(
            f"Not enough final NVT temperature samples in {trace_path} for steps {min_step}..{final_nvt_end} "
            f"(got {len(temps)} rows)"
        )

    mean_temp = sum(temps) / len(temps)
    abs_err = abs(mean_temp - target_temp)
    if abs_err > abs_tol:
        raise RuntimeError(
            f"Final NVT temperature check failed: mean T={mean_temp:.3f} K vs target "
            f"{target_temp:.3f} K (|ΔT|={abs_err:.3f} K > {abs_tol}) "
            f"over steps {min_step}..{final_nvt_end} (n={len(temps)})"
        )

    print(
        f"  Final NVT temperature check OK: mean T={mean_temp:.2f} K vs target "
        f"{target_temp:.2f} K (|ΔT| <= {abs_tol}; {len(temps)} samples)"
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
    bamboo_cfg = get_section(cfg, "bamboo")
    return {
        "data_file": workdir / str(bamboo_cfg.get("data_file", "in.data")),
        "build_report": workdir / str(bamboo_cfg.get("build_report_file", "build_report.json")),
        "equil_input": workdir / str(sim_cfg.get("equil_input_filename", "in.equil.lammps")),
        "gk_input": workdir / str(sim_cfg.get("gk_input_filename", "in.gk.lammps")),
        "equil_restart": workdir / str(sim_cfg.get("equil_restart_file", "equil_nvt.restart")),
        "npt_thermo": workdir / str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat")),
        "acf_file": workdir / str(sim_cfg.get("acf_file", "stress_acf.dat")),
        "pressure_file": workdir / str(sim_cfg.get("pressure_file", "pressure_tensor.dat")),
        "thermo_file": workdir / str(sim_cfg.get("thermo_file", "gk_thermo.dat")),
        "summary_file": workdir / str(analysis_cfg.get("summary_file", "viscosity_summary.json")),
        "plot_file": workdir / str(analysis_cfg.get("plot_file", "viscosity_analysis.png")),
        "running_file": workdir / str(analysis_cfg.get("running_file", "viscosity_running.csv")),
        "blocks_file": workdir / str(analysis_cfg.get("blocks_file", "viscosity_blocks.csv")),
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
            "components": cfg.get("components"),
            "composition": cfg.get("composition"),
            "structure": get_section(cfg, "structure"),
            "model": get_section(cfg, "model"),
            "bamboo": {
                "data_file": get_section(cfg, "bamboo").get("data_file"),
                "build_report_file": get_section(cfg, "bamboo").get("build_report_file"),
                "neutrality_tolerance": get_section(cfg, "bamboo").get("neutrality_tolerance"),
            },
        }
    elif stage == "write_input":
        payload = {
            "bamboo": get_section(cfg, "bamboo"),
            "simulation": get_section(cfg, "simulation"),
        }
    elif stage == "analyze":
        sim_cfg = get_section(cfg, "simulation")
        payload = {
            "analysis": get_section(cfg, "analysis"),
            "simulation": {
                "temperature_K": sim_cfg.get("temperature_K"),
                "timestep_fs": sim_cfg.get("timestep_fs"),
                "acf_file": sim_cfg.get("acf_file"),
                "pressure_file": sim_cfg.get("pressure_file"),
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
    method = str(summary.get("analysis_method", "acf"))
    sem = summary.get("eta_sem_mPas")
    sem_text = f" ± {float(sem):.3f}" if sem is not None else ""
    print(f"\n  ┌─ Result ─────────────────────────────────────────")
    print(f"  │  η = {eta:.3f}{sem_text} mPa·s (cP)")
    print(f"  │  method     : {method}")
    print(f"  │  components: xy={comps[0]:.3f}  xz={comps[1]:.3f}  yz={comps[2]:.3f} mPa·s")
    print(f"  │  plateau   : {float(summary['plateau_start_ps']):.2f} – {float(summary['plateau_end_ps']):.2f} ps")
    print(f"  └──────────────────────────────────────────────────")


def _resolve_workdir(cfg: dict[str, Any]) -> Path:
    """Return the absolute case workdir from the config."""
    return Path(get_required(cfg, "case", "workdir")).resolve()


# ──────────────────────────────────────────────────────────────────────────
# Stage helpers
# ──────────────────────────────────────────────────────────────────────────

def stage_build_system(
    cfg: dict[str, Any],
    workdir: Path,
) -> tuple[Path, tuple[float, float, float]]:
    """Composition → geometries → Packmol box → BAMBOO LAMMPS data file.

    Returns
    -------
    (data_file_path, (Lx, Ly, Lz))
    """
    result = build_bamboo_system(cfg, workdir)
    box_lengths = result.box_lengths
    Lx, Ly, Lz = box_lengths
    V = Lx * Ly * Lz
    print(
        f"  in.data : {result.n_molecules} molecules  |  "
        f"box {Lx:.1f} x {Ly:.1f} x {Lz:.1f} A  |  V={V:.0f} A^3"
    )
    return result.data_file, box_lengths


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
    bamboo_cfg = get_section(cfg, "bamboo")
    build_report = load_bamboo_build_report(workdir, cfg)
    elements = build_report.get("elements_in_pair_coeff")
    if not isinstance(elements, list) or not elements:
        raise ValueError("BAMBOO build report is missing elements_in_pair_coeff")

    print("\n-- Stage 2: write BAMBOO Green-Kubo inputs --------------------")
    pppm_mesh = bamboo_cfg.get("pppm_mesh")
    if pppm_mesh is not None:
        pppm_mesh = tuple(int(value) for value in pppm_mesh)

    common = dict(
        model_file=get_required(cfg, "bamboo", "model_file"),
        elements=[str(element) for element in elements],
        pair_style_args=tuple(bamboo_cfg.get("pair_style_args", (5.0, 5.0, 10.0, 1))),
        pppm_accuracy=float(bamboo_cfg.get("pppm_accuracy", 1.0e-6)),
        pppm_mesh=pppm_mesh,
        temperature=float(sim_cfg.get("temperature_K", 300.0)),
        pressure_atm=float(sim_cfg.get("pressure_atm", 1.0)),
        timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
        thermo_every=int(sim_cfg.get("thermo_every", 1000)),
    )

    restart_file = str(sim_cfg.get("equil_restart_file", "equil_nvt.restart"))

    equil_text = build_equil_input(
        data_file=data_file.name,
        nvt_steps=int(sim_cfg.get("nvt_steps", 40000)),
        npt_steps=int(sim_cfg.get("npt_steps", 200000)),
        final_nvt_steps=int(sim_cfg.get("final_nvt_steps", 20000)),
        npt_thermo_file=str(sim_cfg.get("npt_thermo_file", "npt_thermo.dat")),
        seed=int(sim_cfg.get("seed", 12345)),
        restart_file=restart_file,
        minimize=bool(sim_cfg.get("minimize", True)),
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
        pressure_file=str(sim_cfg.get("pressure_file", "pressure_tensor.dat")),
        pressure_every=int(sim_cfg.get("pressure_every", sim_cfg.get("sample_every", 1))),
        thermo_file=str(sim_cfg.get("thermo_file", "gk_thermo.dat")),
        model_file=common["model_file"],
        elements=common["elements"],
        pair_style_args=common["pair_style_args"],
        pppm_accuracy=common["pppm_accuracy"],
        pppm_mesh=common["pppm_mesh"],
        temperature=common["temperature"],
        timestep_fs=common["timestep_fs"],
        thermo_every=common["thermo_every"],
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
    gk_outputs = [paths["acf_file"], paths["thermo_file"]]
    if sim_cfg.get("pressure_file", "pressure_tensor.dat"):
        gk_outputs.append(paths["pressure_file"])
    result: dict[str, Any] = {
        "run_equil": run_equil,
        "run_gk": run_gk,
        "equil_outputs": _relative_paths([paths["equil_restart"], paths["npt_thermo"]], workdir),
        "gk_outputs": _relative_paths(gk_outputs, workdir),
    }

    print("\n── Stage 3: run LAMMPS ─────────────────────────────────────────")

    if run_equil:
        equil_outputs = [paths["equil_restart"], paths["npt_thermo"]]
        reuse_equil = _outputs_up_to_date(equil_outputs, [equil_path])
        if reuse_equil:
            try:
                _validate_equilibration_trace(workdir, sim_cfg, model_cfg)
            except Exception as exc:
                print(f"  [rerun] equilibration outputs failed validation: {exc}")
                reuse_equil = False

        if reuse_equil:
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
            _validate_equilibration_trace(workdir, sim_cfg, model_cfg)
    else:
        print("  [skip] equilibration run")
        result["equil_status"] = "skipped"

    if run_gk:
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
    """Analyze the Green-Kubo production traces."""
    sim_cfg = get_section(cfg, "simulation")
    analysis_cfg = get_section(cfg, "analysis")

    print("\n── Stage 4: analyze ────────────────────────────────────────────")

    acf_file = workdir / str(sim_cfg.get("acf_file", "stress_acf.dat"))
    pressure_file = workdir / str(sim_cfg.get("pressure_file", "pressure_tensor.dat"))
    thermo_file = workdir / str(sim_cfg.get("thermo_file", "gk_thermo.dat"))

    volume_a3 = load_volume_from_thermo(thermo_file)
    print(f"  NVE mean volume : {volume_a3:.1f} Å³")

    method = str(analysis_cfg.get("method", "acf")).strip().lower()
    if method not in {"acf", "pressure_blocks", "auto"}:
        raise ValueError("analysis.method must be one of: acf, pressure_blocks, auto")

    if method == "auto":
        method = "pressure_blocks" if pressure_file.exists() else "acf"

    print(f"  analysis method : {method}")
    common = dict(
        temperature_k=float(sim_cfg.get("temperature_K", 300.0)),
        volume_a3=volume_a3,
        t_max_ps=float(analysis_cfg.get("t_max_ps", 20.0)),
        plateau_start_ps=float(analysis_cfg.get("plateau_start_ps", 5.0)),
        plateau_end_ps=(
            None
            if analysis_cfg.get("plateau_end_ps") is None
            else float(analysis_cfg.get("plateau_end_ps"))
        ),
        min_window_ps=float(analysis_cfg.get("min_window_ps", 0.0)),
        require_stable_window=bool(analysis_cfg.get("require_stable_window", False)),
        json_out=str(analysis_cfg.get("summary_file", "viscosity_summary.json")),
        plot_out=str(analysis_cfg.get("plot_file", "viscosity_analysis.png")),
    )

    if method == "pressure_blocks":
        summary = analyze_viscosity_pressure_file(
            pressure_file,
            timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
            block_size_ps=float(analysis_cfg.get("block_size_ps", 200.0)),
            auto_plateau=bool(analysis_cfg.get("auto_plateau", True)),
            auto_window_ps=float(analysis_cfg.get("auto_window_ps", 10.0)),
            auto_step_ps=float(analysis_cfg.get("auto_step_ps", 1.0)),
            min_blocks=int(analysis_cfg.get("min_blocks", 3)),
            rel_std_tol=float(analysis_cfg.get("rel_std_tol", 0.25)),
            drift_tol=float(analysis_cfg.get("drift_tol", 0.25)),
            start_step=(
                None
                if analysis_cfg.get("start_step") is None
                else int(analysis_cfg.get("start_step"))
            ),
            running_out=str(analysis_cfg.get("running_file", "viscosity_running.csv")),
            blocks_out=str(analysis_cfg.get("blocks_file", "viscosity_blocks.csv")),
            **common,
        )
    else:
        summary = analyze_viscosity_file(
            acf_file,
            smooth_window=int(analysis_cfg.get("smooth_window", 11)),
            **common,
        )

    _print_viscosity_result(summary)
    return summary


# ──────────────────────────────────────────────────────────────────────────
# Single-case runner
# ──────────────────────────────────────────────────────────────────────────

def _run_single_case(cfg: dict[str, Any], config_path: Path) -> dict[str, Any] | None:
    """Run one viscosity workflow case and return its analysis summary when available."""
    run_flags = get_section(cfg, "run")
    workdir = _resolve_workdir(cfg)
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
    summary: dict[str, Any] | None = None

    try:
        if run_flags.get("build_system", True):
            current_stage = "build_system"
            build_signature = _stage_signature(cfg, "build_system")
            if _can_reuse_stage(
                state,
                "build_system",
                outputs=[paths["data_file"], paths["build_report"]],
                signature=build_signature,
                allow_artifact_fallback=True,
            ):
                data_file = paths["data_file"]
                build_outputs = [paths["data_file"], paths["build_report"]]
                _print_reuse("build_system", workdir, build_outputs)
                _mark_stage(
                    workdir,
                    state,
                    "build_system",
                    status="reused",
                    extra={
                        "signature": build_signature,
                        "outputs": _relative_paths(build_outputs, workdir),
                    },
                )
            else:
                data_file, _ = stage_build_system(cfg, workdir)
                build_outputs = [data_file, paths["build_report"]]
                _mark_stage(
                    workdir,
                    state,
                    "build_system",
                    status="done",
                    extra={
                        "signature": build_signature,
                        "outputs": _relative_paths(build_outputs, workdir),
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
            analysis_cfg = get_section(cfg, "analysis")
            analyze_method = str(analysis_cfg.get("method", "acf")).strip().lower()
            if analyze_method == "auto":
                analyze_method = "pressure_blocks" if paths["pressure_file"].exists() else "acf"
            if analyze_method == "pressure_blocks":
                analyze_outputs.extend([paths["running_file"], paths["blocks_file"]])
                analyze_inputs = [paths["pressure_file"], paths["thermo_file"]]
            else:
                analyze_inputs = [paths["acf_file"], paths["thermo_file"]]
            if _can_reuse_stage(
                state,
                "analyze",
                outputs=analyze_outputs,
                inputs=analyze_inputs,
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
        return summary
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


def _run_workflow(
    cfg: dict[str, Any],
    config_path: Path,
) -> dict[str, Any] | None:
    """Run one viscosity workflow case."""
    return _run_single_case(cfg, config_path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="BAMBOO Green-Kubo viscosity workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args(argv)

    config_path = Path(args.config).resolve()
    cfg = load_yaml(config_path)
    apply_md_viscosity_path_resolution(cfg, REPO_ROOT)
    _run_workflow(cfg, config_path)


if __name__ == "__main__":
    main()
