"""One-shot entrypoint for the LAMMPS Green-Kubo workflow.

Usage:
    python runlmp.py --input config.yaml
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent))

from gk_workflow.analysis import GreenKuboAnalyzer
from gk_workflow.build_system import build_system_data_only
from gk_workflow.write_lammps_input import (
    build_equilibration_input_text,
    build_replica_input_text,
)


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required for runlmp.py. Install it with: pip install pyyaml"
        ) from exc

    with open(path) as f:
        data = yaml.safe_load(f)
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


def load_volume_from_gk_data(gk_data_file: Path) -> float:
    rows = []
    with open(gk_data_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            rows.append(parts)

    if not rows:
        raise ValueError(f"No numeric rows found in gk_data file: {gk_data_file}")

    volumes = [float(row[2]) for row in rows]
    return sum(volumes) / len(volumes)


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def sem(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    m = mean(values)
    var = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return math.sqrt(var) / math.sqrt(len(values))


def write_equilibration_input(workdir: Path, sim_cfg: dict[str, Any]) -> Path:
    out_name = sim_cfg.get("equil_input_filename", "in.equil.lammps")
    out_file = workdir / out_name
    out_file.write_text(
        build_equilibration_input_text(
            temperature=float(sim_cfg.get("temperature_K", 300.0)),
            seed=int(sim_cfg.get("seed", 20260316)),
            npt_pre_50_steps=int(sim_cfg.get("npt_pre_50_steps", 200000)),
            npt_pre_10_steps=int(sim_cfg.get("npt_pre_10_steps", 300000)),
            npt_1atm_steps=int(sim_cfg.get("npt_1atm_steps", 2000000)),
            nvt_steps=int(sim_cfg.get("nvt_steps", 1000000)),
            timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
        )
    )
    return out_file


def write_replica_inputs(workdir: Path, sim_cfg: dict[str, Any]) -> list[tuple[int, Path, Path]]:
    n_replicas = int(sim_cfg.get("n_replicas", 3))
    replica_seed_base = int(sim_cfg.get("replica_seed_base", 310000))
    replica_seed_step = int(sim_cfg.get("replica_seed_step", 137))
    replica_dir_prefix = str(sim_cfg.get("replica_dir_prefix", "replica_"))
    replica_input_name = str(sim_cfg.get("replica_input_filename", "in.replica.lammps"))
    replica_decorrelation_steps = int(sim_cfg.get("replica_decorrelation_steps", 500000))

    outputs: list[tuple[int, Path, Path]] = []
    for idx in range(1, n_replicas + 1):
        repdir = workdir / f"{replica_dir_prefix}{idx:02d}"
        repdir.mkdir(parents=True, exist_ok=True)
        input_file = repdir / replica_input_name
        seed = replica_seed_base + (idx - 1) * replica_seed_step
        input_file.write_text(
            build_replica_input_text(
                temperature=float(sim_cfg.get("temperature_K", 300.0)),
                replica_seed=seed,
                decorrelation_steps=replica_decorrelation_steps,
                gk_steps=int(sim_cfg.get("gk_steps", 10000000)),
                timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
                sample_every_steps=int(sim_cfg.get("acf_sample_every_steps", 5)),
                corr_points=int(sim_cfg.get("acf_corr_points", 8000)),
                equil_restart_relpath="../equil_nvt.restart",
            )
        )
        outputs.append((idx, repdir, input_file))
    return outputs


def build_lammps_command(lammps_cfg: dict[str, Any], input_file: Path) -> list[str]:
    executable = lammps_cfg.get("executable", "/opt/lammps/bin/lmp")
    mpi_command = lammps_cfg.get("mpi_command", "mpirun")
    mpi_ranks = int(lammps_cfg.get("mpi_ranks", 8))
    allow_run_as_root = bool(lammps_cfg.get("allow_run_as_root", True))
    use_gpu = bool(lammps_cfg.get("use_gpu", True))
    gpu_count = int(lammps_cfg.get("gpu_count", 1))

    cmd = [mpi_command]
    if allow_run_as_root:
        cmd.append("--allow-run-as-root")
    cmd += ["-np", str(mpi_ranks), executable]
    if use_gpu:
        cmd += ["-sf", "gpu", "-pk", "gpu", str(gpu_count)]
    return cmd + ["-in", input_file.name]


def run_lammps(workdir: Path, lammps_cfg: dict[str, Any], input_file: Path, log_file: str) -> None:
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


def analyze_replica(repdir: Path, sim_cfg: dict[str, Any], analysis_cfg: dict[str, Any]) -> dict[str, Any]:
    hfacf = repdir / analysis_cfg.get("hfacf_file", "hfacf.dat")
    vacf = repdir / analysis_cfg.get("vacf_file", "vacf.dat")
    gk_data = repdir / analysis_cfg.get("gk_data_file", "gk_data.dat")

    if not hfacf.exists():
        raise FileNotFoundError(f"Missing HFACF file: {hfacf}")
    if not vacf.exists():
        raise FileNotFoundError(f"Missing VACF file: {vacf}")

    if gk_data.exists():
        volume = load_volume_from_gk_data(gk_data)
    else:
        manual_volume = analysis_cfg.get("volume_A3")
        if manual_volume is None:
            raise ValueError(f"gk_data.dat not found in {repdir} and analysis.volume_A3 is not set")
        volume = float(manual_volume)

    analyzer = GreenKuboAnalyzer(repdir)
    results = analyzer.analyze(
        corr_stress_file=vacf,
        corr_flux_file=hfacf,
        T=float(sim_cfg.get("temperature_K", 300.0)),
        volume=volume,
        visc_t_max_ps=float(analysis_cfg.get("visc_t_max_ps", 20.0)),
        kappa_t_max_ps=float(analysis_cfg.get("kappa_t_max_ps", 20.0)),
        visc_plateau_start_ps=float(analysis_cfg.get("visc_plateau_start_ps", 5.0)),
        kappa_plateau_start_ps=float(analysis_cfg.get("kappa_plateau_start_ps", 2.0)),
        smooth_window=int(analysis_cfg.get("smooth_window", 50)),
        n_blocks=int(analysis_cfg.get("n_blocks", 5)),
        save_summary=False,
    )
    return {
        "replica_dir": str(repdir),
        "volume_A3": volume,
        "eta_cP": results["eta_cP"],
        "eta_sem_cP_internal": results["eta_sem_cP"],
        "eta_Pa_s": results["eta_Pa_s"],
        "kappa_W_mK": results["kappa_W_mK"],
        "kappa_sem_W_mK_internal": results["kappa_sem_W_mK"],
    }


def write_result_csv(workdir: Path, replica_results: list[dict[str, Any]], output_name: str) -> Path:
    eta_vals = [row["eta_cP"] for row in replica_results]
    kappa_vals = [row["kappa_W_mK"] for row in replica_results]
    volume_vals = [row["volume_A3"] for row in replica_results]
    eta_internal = [row["eta_sem_cP_internal"] for row in replica_results]
    kappa_internal = [row["kappa_sem_W_mK_internal"] for row in replica_results]

    summary = {
        "row_type": "summary",
        "replica": "all",
        "volume_A3": mean(volume_vals),
        "eta_cP": mean(eta_vals),
        "eta_sem_cP": sem(eta_vals),
        "eta_sem_cP_internal_mean": mean(eta_internal),
        "kappa_W_mK": mean(kappa_vals),
        "kappa_sem_W_mK": sem(kappa_vals),
        "kappa_sem_W_mK_internal_mean": mean(kappa_internal),
        "replica_dir": str(workdir),
    }

    out_file = workdir / output_name
    fieldnames = [
        "row_type",
        "replica",
        "replica_dir",
        "volume_A3",
        "eta_cP",
        "eta_sem_cP",
        "eta_sem_cP_internal_mean",
        "kappa_W_mK",
        "kappa_sem_W_mK",
        "kappa_sem_W_mK_internal_mean",
    ]

    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)
        for idx, row in enumerate(replica_results, start=1):
            writer.writerow(
                {
                    "row_type": "replica",
                    "replica": idx,
                    "replica_dir": row["replica_dir"],
                    "volume_A3": row["volume_A3"],
                    "eta_cP": row["eta_cP"],
                    "eta_sem_cP": row["eta_sem_cP_internal"],
                    "eta_sem_cP_internal_mean": row["eta_sem_cP_internal"],
                    "kappa_W_mK": row["kappa_W_mK"],
                    "kappa_sem_W_mK": row["kappa_sem_W_mK_internal"],
                    "kappa_sem_W_mK_internal_mean": row["kappa_sem_W_mK_internal"],
                }
            )
    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the full LAMMPS Green-Kubo workflow from a YAML config")
    parser.add_argument("--input", required=True, help="YAML config file")
    args = parser.parse_args()

    cfg = load_yaml(Path(args.input))
    run_cfg = get_section(cfg, "run")
    case_cfg = get_section(cfg, "case")
    model_cfg = get_section(cfg, "model")
    sim_cfg = get_section(cfg, "simulation")
    lammps_cfg = get_section(cfg, "lammps")
    analysis_cfg = get_section(cfg, "analysis")

    workdir = Path(get_required(cfg, "case", "workdir"))
    workdir.mkdir(parents=True, exist_ok=True)

    do_build = bool(run_cfg.get("build_system", True))
    do_write_input = bool(run_cfg.get("write_input", True))
    do_run_lammps = bool(run_cfg.get("run_lammps", True))
    do_analyze = bool(run_cfg.get("analyze", True))

    if do_build:
        build_system_data_only(
            smiles=str(get_required(cfg, "case", "smiles")),
            name=str(get_required(cfg, "case", "name")),
            workdir=workdir,
            n_molecules=int(get_required(cfg, "model", "n_molecules")),
            density=float(get_required(cfg, "model", "density_g_cm3")),
            packmol_seed=int(model_cfg.get("packmol_seed", 192911)),
            packmol_max_attempts=int(model_cfg.get("packmol_max_attempts", 3)),
            packmol_seed_step=int(model_cfg.get("packmol_seed_step", 97)),
            packmol_nloop=model_cfg.get("packmol_nloop"),
            packmol_full_box=bool(model_cfg.get("packmol_full_box", False)),
            packmol_margin=model_cfg.get("packmol_margin_A"),
            packmol_strict=bool(model_cfg.get("packmol_strict", False)),
        )
    else:
        system_data = workdir / "system.data"
        if not system_data.exists():
            raise FileNotFoundError(f"build_system is disabled but system.data is missing: {system_data}")

    equil_input = workdir / sim_cfg.get("equil_input_filename", "in.equil.lammps")
    replica_specs: list[tuple[int, Path, Path]] = []
    if do_write_input:
        equil_input = write_equilibration_input(workdir, sim_cfg)
        print(f"Wrote {equil_input}")
        replica_specs = write_replica_inputs(workdir, sim_cfg)
        for idx, _, input_file in replica_specs:
            print(f"Wrote replica {idx} input: {input_file}")
    else:
        if not equil_input.exists():
            raise FileNotFoundError(f"write_input is disabled but equilibration input is missing: {equil_input}")
        n_replicas = int(sim_cfg.get("n_replicas", 3))
        replica_dir_prefix = str(sim_cfg.get("replica_dir_prefix", "replica_"))
        replica_input_name = str(sim_cfg.get("replica_input_filename", "in.replica.lammps"))
        replica_specs = [
            (idx, workdir / f"{replica_dir_prefix}{idx:02d}", workdir / f"{replica_dir_prefix}{idx:02d}" / replica_input_name)
            for idx in range(1, n_replicas + 1)
        ]
        for _, repdir, input_file in replica_specs:
            if not input_file.exists():
                raise FileNotFoundError(f"Replica input is missing: {input_file}")
            repdir.mkdir(parents=True, exist_ok=True)

    if do_run_lammps:
        equil_log = str(lammps_cfg.get("equil_log_file", "run_equil.log"))
        run_lammps(workdir, lammps_cfg, equil_input, equil_log)
        equil_restart = workdir / "equil_nvt.restart"
        if not equil_restart.exists():
            raise FileNotFoundError(f"Equilibration finished but restart not found: {equil_restart}")

        replica_log_prefix = str(lammps_cfg.get("replica_log_prefix", "run_replica"))
        for idx, repdir, input_file in replica_specs:
            log_name = f"{replica_log_prefix}{idx}.log"
            run_lammps(repdir, lammps_cfg, input_file, log_name)

    if do_analyze:
        replica_results = []
        for idx, repdir, _ in replica_specs:
            result = analyze_replica(repdir, sim_cfg, analysis_cfg)
            result["replica_index"] = idx
            replica_results.append(result)

        out_name = str(analysis_cfg.get("output", "result.csv"))
        result_file = write_result_csv(workdir, replica_results, out_name)
        print("=" * 70)
        print("Workflow finished")
        print("=" * 70)
        print(f"  replicas : {len(replica_results)}")
        print(f"  result   : {result_file}")


if __name__ == "__main__":
    main()
