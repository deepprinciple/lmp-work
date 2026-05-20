"""Unified transport-property workflow dispatcher.

Usage
-----
    python md_transport.py --config configs/viscosity.yaml
    python md_transport.py --config configs/config.yaml

The actual physics workflows remain in their focused entry points:

* ``md_viscosity.py`` for BAMBOO Green-Kubo viscosity
* ``md_run.py`` for reverse-NEMD thermal conductivity

This dispatcher reads ``job.type`` from the YAML config and forwards the run to
the matching workflow.  If ``job.type`` is missing, it infers the route from the
sections present in the config so existing YAML files remain usable.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from workflow.config import get_section, load_yaml


THERMAL_JOB_TYPES = {
    "thermal",
    "thermal_conductivity",
    "conductivity",
    "rnemd",
}
VISCOSITY_JOB_TYPES = {
    "viscosity",
    "eta",
}


def _normalize_job_type(value: Any) -> str:
    return str(value).strip().lower().replace("-", "_")


def _infer_job_type(cfg: dict[str, Any]) -> str:
    job_cfg = get_section(cfg, "job")
    explicit = (
        job_cfg.get("type")
        or job_cfg.get("jobtype")
        or cfg.get("jobtype")
        or cfg.get("job_type")
    )
    if explicit:
        return _normalize_job_type(explicit)

    if "thermal_rnemd" in cfg:
        return "thermal_conductivity"
    if "bamboo" in cfg and "components" in cfg:
        return "viscosity"

    raise ValueError(
        "Could not infer job.type. Set job.type to 'thermal_conductivity' or 'viscosity'."
    )


def _job_method(cfg: dict[str, Any]) -> str:
    job_cfg = get_section(cfg, "job")
    method = job_cfg.get("method") or ""
    return str(method).strip().lower().replace("-", "_")


def _analysis_method(cfg: dict[str, Any]) -> str:
    analysis_cfg = get_section(cfg, "analysis")
    method = analysis_cfg.get("method") or ""
    return str(method).strip().lower().replace("-", "_")


def _resolve_workflow(cfg: dict[str, Any]) -> tuple[str, str, str]:
    job_type = _infer_job_type(cfg)
    job_method = _job_method(cfg)
    analysis_method = _analysis_method(cfg)

    if job_type in THERMAL_JOB_TYPES:
        if job_method and job_method not in {"rnemd", "reverse_nemd"}:
            raise ValueError(
                f"Unsupported thermal_conductivity job.method '{job_method}'. "
                "The unified dispatcher currently supports rnemd."
            )
        if analysis_method and analysis_method not in {"rnemd", "reverse_nemd"}:
            raise ValueError(
                f"Unsupported thermal_conductivity analysis.method '{analysis_method}'. "
                "The unified dispatcher currently supports rNEMD post-processing."
            )
        return "thermal_conductivity", job_method or "rnemd", analysis_method or "rnemd"

    if job_type in VISCOSITY_JOB_TYPES:
        if job_method and job_method not in {"green_kubo", "gk"}:
            raise ValueError(
                f"Unsupported viscosity job.method '{job_method}'. Use green_kubo."
            )
        if analysis_method and analysis_method not in {"acf", "pressure_blocks", "auto"}:
            raise ValueError(
                f"Unsupported viscosity analysis.method '{analysis_method}'. "
                "Use pressure_blocks, acf, or auto."
            )
        return "viscosity", job_method or "green_kubo", analysis_method or "pressure_blocks"

    raise ValueError(
        f"Unsupported job.type '{job_type}'. "
        "Use 'thermal_conductivity' or 'viscosity'."
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Unified dispatcher for lmp-work transport workflows",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", required=True, help="YAML workflow config file")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print the selected workflow without running it",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config_path = Path(args.config).resolve()
    cfg = load_yaml(config_path)
    workflow, method, analysis_method = _resolve_workflow(cfg)

    print(f"Transport workflow : {workflow}")
    print(f"Method             : {method}")
    print(f"Analysis method    : {analysis_method}")
    print(f"Config             : {config_path}")

    if args.dry_run:
        return 0

    if workflow == "thermal_conductivity":
        from md_run import main as run_thermal

        return int(run_thermal(["--config", str(config_path)]) or 0)

    if workflow == "viscosity":
        from md_viscosity import main as run_viscosity

        run_viscosity(["--config", str(config_path)])
        return 0

    raise AssertionError(f"Unhandled workflow: {workflow}")


if __name__ == "__main__":
    raise SystemExit(main())
