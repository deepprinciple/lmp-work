from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Any

from gk_workflow.forcefields import get_lammps_settings_for_engine


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for thermal workflow configs."
        ) from exc

    with open(path) as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def get_required(cfg: dict[str, Any], section: str, key: str) -> Any:
    if section not in cfg or not isinstance(cfg[section], dict):
        raise KeyError(f"Missing section '{section}' in config")
    if key not in cfg[section]:
        raise KeyError(f"Missing key '{section}.{key}' in config")
    return cfg[section][key]


def get_section(cfg: dict[str, Any], section: str) -> dict[str, Any]:
    value = cfg.get(section, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise TypeError(f"Section '{section}' must be a mapping")
    return value


def build_forcefield_settings(forcefield_cfg: dict[str, Any]) -> dict[str, Any]:
    engine = str(forcefield_cfg.get("engine", "openff")).lower()
    return get_lammps_settings_for_engine(engine)


def build_lammps_command(
    lammps_cfg: dict[str, Any],
    input_file: Path,
) -> list[str]:
    executable = lammps_cfg.get("executable", "/opt/lammps/bin/lmp")
    mpi_command = lammps_cfg.get("mpi_command", "mpirun")
    mpi_ranks = int(lammps_cfg.get("mpi_ranks", 4))
    allow_run_as_root = bool(lammps_cfg.get("allow_run_as_root", True))
    use_hwthread_cpus = bool(lammps_cfg.get("use_hwthread_cpus", False))
    use_gpu = bool(lammps_cfg.get("use_gpu", False))
    gpu_count = int(lammps_cfg.get("gpu_count", 0))

    cmd = [mpi_command]
    if allow_run_as_root:
        cmd.append("--allow-run-as-root")
    if use_hwthread_cpus:
        cmd.append("--use-hwthread-cpus")
    cmd += ["-np", str(mpi_ranks), executable]
    if use_gpu:
        cmd += ["-sf", "gpu", "-pk", "gpu", str(gpu_count)]
    return cmd + ["-in", input_file.name]


def run_lammps(
    workdir: Path,
    lammps_cfg: dict[str, Any],
    input_file: Path,
    log_file: str,
) -> None:
    cmd = build_lammps_command(lammps_cfg, input_file) + ["-log", log_file]
    env = os.environ.copy()
    if "omp_threads" in lammps_cfg:
        env["OMP_NUM_THREADS"] = str(lammps_cfg["omp_threads"])

    print("=" * 70)
    print("Run LAMMPS")
    print("=" * 70)
    print("  cwd     :", workdir)
    print("  command :", " ".join(cmd))
    if "OMP_NUM_THREADS" in env:
        print("  OMP     :", env["OMP_NUM_THREADS"])

    subprocess.run(cmd, cwd=workdir, env=env, check=True)
