"""YAML configuration helpers for workflow entrypoints."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any

# Repository root (this file lives in <repo>/workflow/config.py).
REPO_ROOT: Path = Path(__file__).resolve().parent.parent


def resolve_config_file_path(value: Any, repo_root: Path | None = None) -> str:
    """Expand env vars and resolve file paths against ``repo_root``."""
    root = repo_root if repo_root is not None else REPO_ROOT
    text = os.path.expandvars(os.path.expanduser(str(value).strip()))
    if not text:
        return ""
    path = Path(text)
    if path.is_absolute():
        return str(path.resolve())
    return str((root / path).resolve())


def resolve_lammps_executable(value: Any, repo_root: Path | None = None) -> str:
    """Expand env vars while keeping bare command names for PATH lookup."""
    root = repo_root if repo_root is not None else REPO_ROOT
    text = os.path.expandvars(os.path.expanduser(str(value).strip()))
    if not text:
        return "lmp_mpi"
    if "/" not in text and not text.startswith((".", "~")):
        return text
    path = Path(text)
    if path.is_absolute():
        return str(path.resolve())
    return str((root / path).resolve())


def _expand_lammps_env_value(value: Any, repo_root: Path) -> str:
    """Expand YAML ``lammps.env`` values without mangling PATH-like strings."""
    text = os.path.expanduser(str(value))
    if not text.strip():
        return text
    if "$" in text:
        return text
    if ":" in text:
        return text
    if "/" in text or text.startswith("."):
        path = Path(text)
        if not path.is_absolute():
            return str((repo_root / path).resolve())
        return str(path.resolve())
    return text


def apply_md_viscosity_path_resolution(cfg: dict[str, Any], repo_root: Path | None = None) -> None:
    """Mutate viscosity config in place to resolve file paths and env values."""
    root = repo_root if repo_root is not None else REPO_ROOT

    bamboo_cfg = cfg.get("bamboo")
    if isinstance(bamboo_cfg, dict) and bamboo_cfg.get("model_file") is not None:
        bamboo_cfg["model_file"] = resolve_config_file_path(bamboo_cfg["model_file"], root)

    components_cfg = cfg.get("components")
    if isinstance(components_cfg, list):
        for component in components_cfg:
            if isinstance(component, dict) and component.get("charge_file") is not None:
                component["charge_file"] = resolve_config_file_path(component["charge_file"], root)

    lammps_cfg = cfg.get("lammps")
    if not isinstance(lammps_cfg, dict):
        return
    if lammps_cfg.get("executable") is not None:
        lammps_cfg["executable"] = resolve_lammps_executable(lammps_cfg["executable"], root)
    env_cfg = lammps_cfg.get("env")
    if isinstance(env_cfg, dict):
        for key in list(env_cfg.keys()):
            env_cfg[str(key)] = _expand_lammps_env_value(env_cfg[key], root)


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required to load workflow configs.") from exc

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


__all__ = [
    "REPO_ROOT",
    "apply_md_viscosity_path_resolution",
    "get_required",
    "get_section",
    "load_yaml",
    "resolve_config_file_path",
    "resolve_lammps_executable",
]
