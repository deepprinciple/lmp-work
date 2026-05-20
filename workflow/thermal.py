"""Internal reverse-NEMD thermal-conductivity workflow."""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from scripts.build_system import build_system_data_only
from scripts._cli import Argv, new_parser
from workflow.config import get_required, get_section, load_yaml
from workflow.input_thermal import (
    build_equilibration_input_text,
    build_rnemd_replica_input_text,
)
from workflow.lammps_runner import build_forcefield_settings, run_lammps
from workflow.shared_branch import write_branch_manifest

DEFAULT_THERMAL_REPLICA_DIR_PREFIX = "thermal_conductivity/method_rnemd/replica_"
DEFAULT_THERMAL_REPLICA_INPUT_FILENAME = "in.thermal_rnemd.lammps"
DEFAULT_THERMAL_REPLICA_LOG_PREFIX = "../../../run_thermal_case_replica"
DEFAULT_THERMAL_SUMMARY_FILE = "thermal_rnemd_summary.json"
DEFAULT_THERMAL_BATCH_SUMMARY_FILE = "thermal_rnemd_results.json"


def _to_triplet(
    value: Any,
    *,
    default: tuple[float, float, float],
) -> tuple[float, float, float]:
    if value is None:
        return default
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError("Expected a list/tuple with 3 numbers")
    return tuple(float(v) for v in value)


def _resolve_thermal_layout(thermal_cfg: dict[str, Any]) -> dict[str, str]:
    return {
        "replica_dir_prefix": str(
            thermal_cfg.get(
                "replica_dir_prefix",
                DEFAULT_THERMAL_REPLICA_DIR_PREFIX,
            )
        ),
        "replica_input_filename": str(
            thermal_cfg.get(
                "replica_input_filename",
                DEFAULT_THERMAL_REPLICA_INPUT_FILENAME,
            )
        ),
        "replica_log_prefix": str(
            thermal_cfg.get(
                "replica_log_prefix",
                DEFAULT_THERMAL_REPLICA_LOG_PREFIX,
            )
        ),
        "summary_file": str(
            thermal_cfg.get(
                "summary_file",
                DEFAULT_THERMAL_SUMMARY_FILE,
            )
        ),
        "batch_summary_file": str(
            thermal_cfg.get(
                "batch_summary_file",
                DEFAULT_THERMAL_BATCH_SUMMARY_FILE,
            )
        ),
    }


def write_shared_equil_input(
    workdir: Path,
    sim_cfg: dict[str, Any],
    visualization_cfg: dict[str, Any],
    lammps_settings: dict[str, Any],
) -> Path:
    out_name = str(sim_cfg.get("equil_input_filename", "in.equil.lammps"))
    out_file = workdir / out_name
    out_file.write_text(
        build_equilibration_input_text(
            temperature=float(sim_cfg.get("temperature_K", 300.0)),
            seed=int(sim_cfg.get("seed", 20260316)),
            npt_pre_50_steps=int(sim_cfg.get("npt_pre_50_steps", 0)),
            npt_pre_10_steps=int(sim_cfg.get("npt_pre_10_steps", 0)),
            npt_1atm_steps=int(sim_cfg.get("npt_1atm_steps", 2000000)),
            nvt_steps=int(sim_cfg.get("nvt_steps", 1000000)),
            timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
            render_cfg=visualization_cfg,
            lammps_settings=lammps_settings,
        )
    )
    return out_file


def write_rnemd_replica_inputs(
    workdir: Path,
    sim_cfg: dict[str, Any],
    thermal_cfg: dict[str, Any],
    visualization_cfg: dict[str, Any],
    lammps_settings: dict[str, Any],
) -> list[tuple[int, Path, Path]]:
    layout = _resolve_thermal_layout(thermal_cfg)
    n_replicas = int(sim_cfg.get("n_replicas", 1))
    replica_seed_base = int(sim_cfg.get("replica_seed_base", 310000))
    replica_seed_step = int(sim_cfg.get("replica_seed_step", 137))

    outputs: list[tuple[int, Path, Path]] = []
    equil_restart = workdir / "equil_nvt.restart"
    for idx in range(1, n_replicas + 1):
        repdir = workdir / f"{layout['replica_dir_prefix']}{idx:02d}"
        repdir.mkdir(parents=True, exist_ok=True)
        input_file = repdir / layout["replica_input_filename"]
        restart_rel = os.path.relpath(equil_restart, repdir)
        replica_seed = replica_seed_base + (idx - 1) * replica_seed_step
        input_file.write_text(
            build_rnemd_replica_input_text(
                temperature=float(sim_cfg.get("temperature_K", 300.0)),
                replica_seed=replica_seed,
                decorrelation_steps=int(sim_cfg.get("replica_decorrelation_steps", 0)),
                run_steps=int(
                    thermal_cfg.get(
                        "run_steps",
                        sim_cfg.get("fit_sample_steps", 200000),
                    )
                ),
                timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
                swap_every_steps=int(thermal_cfg.get("swap_every_steps", 100)),
                nbin=int(thermal_cfg.get("nbin", 20)),
                edim=str(thermal_cfg.get("edim", "z")),
                nswap=int(thermal_cfg.get("nswap", 1)),
                profile_every_steps=int(thermal_cfg.get("profile_every_steps", 10)),
                profile_repeat=int(thermal_cfg.get("profile_repeat", 100)),
                profile_freq=int(thermal_cfg.get("profile_freq", 1000)),
                profile_file=str(thermal_cfg.get("profile_file", "temp_profile.dat")),
                exchange_file=str(
                    thermal_cfg.get("exchange_file", "thermal_exchange.dat")
                ),
                thermo_file=str(
                    thermal_cfg.get("thermo_file", "thermal_rnemd_thermo.dat")
                ),
                equil_restart_relpath=restart_rel,
                render_cfg=visualization_cfg,
                lammps_settings=lammps_settings,
            )
        )
        outputs.append((idx, repdir, input_file))
    return outputs


def build_branch_manifest_payload(
    replica_specs: list[tuple[int, Path, Path]],
    thermal_cfg: dict[str, Any],
) -> dict[str, Any]:
    return {
        "workflow": "thermal_rnemd_only",
        "replicas": [
            {
                "index": idx,
                "replica_dir": str(repdir),
                "input_file": str(input_file),
            }
            for idx, repdir, input_file in replica_specs
        ],
        "thermal_rnemd": thermal_cfg,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = new_parser("Run the thermal-only reverse-NEMD conductivity workflow.")
    parser.add_argument("--config", required=True, help="YAML workflow config file")
    return parser


def parse_args(argv: Argv = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _build_system(
    cfg: dict[str, Any],
    workdir: Path,
    structure_cfg: dict[str, Any],
    model_cfg: dict[str, Any],
    forcefield_cfg: dict[str, Any],
    thermal_cfg: dict[str, Any],
) -> None:
    build_system_data_only(
        smiles=str(get_required(cfg, "case", "smiles")),
        name=str(get_required(cfg, "case", "name")),
        workdir=workdir,
        n_molecules=int(get_required(cfg, "model", "n_molecules")),
        density=float(get_required(cfg, "model", "density_g_cm3")),
        packmol_density_scale=float(model_cfg.get("packmol_density_scale", 0.85)),
        packmol_tolerance=float(model_cfg.get("packmol_tolerance_A", 2.0)),
        packmol_seed=int(model_cfg.get("packmol_seed", 192911)),
        packmol_max_attempts=int(model_cfg.get("packmol_max_attempts", 3)),
        packmol_seed_step=int(model_cfg.get("packmol_seed_step", 97)),
        packmol_nloop=model_cfg.get("packmol_nloop"),
        packmol_full_box=bool(model_cfg.get("packmol_full_box", False)),
        packmol_margin=model_cfg.get("packmol_margin_A"),
        packmol_strict=bool(model_cfg.get("packmol_strict", False)),
        box_aspect_ratio=_to_triplet(
            thermal_cfg.get("box_aspect_ratio", model_cfg.get("box_aspect_ratio")),
            default=(1.0, 1.0, 3.0),
        ),
        box_lengths=_to_triplet(
            model_cfg.get("box_lengths_A"),
            default=(0.0, 0.0, 0.0),
        )
        if model_cfg.get("box_lengths_A") is not None
        else None,
        structure_output_formats=tuple(
            structure_cfg.get("output_formats", ["xyz", "pdb"])
        ),
        forcefield_engine=str(forcefield_cfg.get("engine", "openff")),
        forcefield_version=str(forcefield_cfg.get("version", "openff-2.0.0")),
        forcefield_strict_stereo=bool(forcefield_cfg.get("strict_stereo", True)),
        forcefield_charge_method=str(forcefield_cfg.get("charge_method", "am1bcc")),
        forcefield_charge_fallback=forcefield_cfg.get("charge_fallback", "gasteiger"),
        forcefield_charge_file=forcefield_cfg.get("charge_file"),
        forcefield_use_cache=bool(forcefield_cfg.get("use_cache", True)),
        forcefield_input_mode=str(forcefield_cfg.get("input_mode", "generated_pdb")),
        forcefield_residue_name=str(forcefield_cfg.get("residue_name", "MOL")),
        forcefield_n_optimizations=int(forcefield_cfg.get("n_optimizations", 0)),
        forcefield_wrapper_script=forcefield_cfg.get("wrapper_script"),
        forcefield_debug=bool(forcefield_cfg.get("debug", False)),
    )


def _analyze_replica(repdir: Path, *, sim_cfg: dict[str, Any], thermal_cfg: dict[str, Any]) -> dict[str, Any]:
    from analysis.rnemd import analyze_rnemd_replica

    return analyze_rnemd_replica(
        repdir,
        timestep_fs=float(sim_cfg.get("timestep_fs", 1.0)),
        nbin=int(thermal_cfg.get("nbin", 20)),
        edim=str(thermal_cfg.get("edim", "z")),
        plateau_start_ps=float(thermal_cfg.get("plateau_start_ps", 20.0)),
        swap_buffer_bins=int(thermal_cfg.get("swap_buffer_bins", 2)),
        kappa_rel_std_tol=float(thermal_cfg.get("kappa_rel_std_tol", 0.15)),
        kappa_drift_tol=float(thermal_cfg.get("kappa_drift_tol", 0.15)),
        min_window_points=int(thermal_cfg.get("min_window_points", 5)),
        profile_file=str(thermal_cfg.get("profile_file", "temp_profile.dat")),
        exchange_file=str(thermal_cfg.get("exchange_file", "thermal_exchange.dat")),
        summary_name=_resolve_thermal_layout(thermal_cfg)["summary_file"],
        figure_profile_name=str(
            thermal_cfg.get("profile_plot_file", "thermal_profile_latest.png")
        ),
        figure_kappa_name=str(
            thermal_cfg.get("kappa_plot_file", "thermal_kappa_running.png")
        ),
    )


def run_config(config_path: str | Path) -> int:
    """Run the thermal-conductivity workflow for one YAML config."""
    cfg = load_yaml(Path(config_path).resolve())
    run_cfg = get_section(cfg, "run")
    structure_cfg = get_section(cfg, "structure")
    model_cfg = get_section(cfg, "model")
    forcefield_cfg = get_section(cfg, "forcefield")
    sim_cfg = get_section(cfg, "simulation")
    lammps_cfg = get_section(cfg, "lammps")
    visualization_cfg = get_section(cfg, "visualization")
    thermal_cfg = get_section(cfg, "thermal_rnemd")

    lammps_settings = build_forcefield_settings(forcefield_cfg)
    workdir = Path(get_required(cfg, "case", "workdir")).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    if bool(run_cfg.get("build_system", True)):
        _build_system(
            cfg=cfg,
            workdir=workdir,
            structure_cfg=structure_cfg,
            model_cfg=model_cfg,
            forcefield_cfg=forcefield_cfg,
            thermal_cfg=thermal_cfg,
        )
    elif not (workdir / "system.data").exists():
        raise FileNotFoundError(
            f"build_system is disabled but system.data is missing: {workdir / 'system.data'}"
        )

    if bool(run_cfg.get("write_input", True)):
        equil_input = write_shared_equil_input(
            workdir,
            sim_cfg,
            visualization_cfg,
            lammps_settings,
        )
        replica_specs = write_rnemd_replica_inputs(
            workdir,
            sim_cfg,
            thermal_cfg,
            visualization_cfg,
            lammps_settings,
        )
        manifest = write_branch_manifest(
            branch_workdir=workdir,
            filename="thermal_branch_manifest.json",
            payload={
                "workdir": str(workdir),
                **build_branch_manifest_payload(replica_specs, thermal_cfg),
            },
        )
        print(f"Wrote shared equilibration input: {equil_input}")
        print(f"Wrote branch manifest: {manifest}")
    else:
        layout = _resolve_thermal_layout(thermal_cfg)
        replica_specs = [
            (
                idx,
                workdir / f"{layout['replica_dir_prefix']}{idx:02d}",
                workdir / f"{layout['replica_dir_prefix']}{idx:02d}"
                / layout["replica_input_filename"],
            )
            for idx in range(1, int(sim_cfg.get("n_replicas", 1)) + 1)
        ]

    if bool(run_cfg.get("run_lammps", False)):
        equil_input = workdir / str(sim_cfg.get("equil_input_filename", "in.equil.lammps"))
        equil_log = str(lammps_cfg.get("equil_log_file", "run_equil.log"))
        run_lammps(workdir, lammps_cfg, equil_input, equil_log)

        layout = _resolve_thermal_layout(thermal_cfg)
        replica_log_prefix = layout["replica_log_prefix"]
        parallel_replicas = bool(sim_cfg.get("parallel_replicas", False))
        if parallel_replicas:
            max_workers = int(
                sim_cfg.get("replica_parallel_workers", len(replica_specs))
            )
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_map = {
                    executor.submit(
                        run_lammps,
                        repdir,
                        lammps_cfg,
                        input_file,
                        f"{replica_log_prefix}{idx}.log",
                    ): idx
                    for idx, repdir, input_file in replica_specs
                }
                for future in as_completed(future_map):
                    future.result()
        else:
            for idx, repdir, input_file in replica_specs:
                run_lammps(
                    repdir,
                    lammps_cfg,
                    input_file,
                    f"{replica_log_prefix}{idx}.log",
                )

    if bool(run_cfg.get("analyze", False)):
        layout = _resolve_thermal_layout(thermal_cfg)
        results = []
        for idx, repdir, _ in replica_specs:
            result = _analyze_replica(
                repdir,
                sim_cfg=sim_cfg,
                thermal_cfg={
                    **thermal_cfg,
                    "summary_file": layout["summary_file"],
                },
            )
            result["replica_index"] = idx
            results.append(result)

        summary = {
            "workflow": "thermal_rnemd_only",
            "workdir": str(workdir),
            "replicas": results,
        }
        if results:
            kappas = [float(row["kappa_W_mK"]) for row in results]
            summary["aggregate"] = {
                "kappa_mean_W_mK": float(sum(kappas) / len(kappas)),
                "kappa_min_W_mK": float(min(kappas)),
                "kappa_max_W_mK": float(max(kappas)),
            }

        summary_file = workdir / layout["batch_summary_file"]
        summary_file.write_text(json.dumps(summary, indent=2))
        print(f"Wrote thermal branch summary: {summary_file}")
    return 0


def main(argv: Argv = None) -> int:
    args = parse_args(argv)
    return run_config(args.config)


if __name__ == "__main__":
    raise SystemExit(main())
