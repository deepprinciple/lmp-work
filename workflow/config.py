"""YAML configuration helpers for workflow entrypoints."""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

# Repository root (this file lives in <repo>/workflow/config.py).
REPO_ROOT: Path = Path(__file__).resolve().parent.parent


def _default_lammps_ani_root(home: Path) -> str:
    """Pick a tree that actually contains ``build/ani_plugin.so``."""
    plugin = Path("build") / "ani_plugin.so"
    candidates = [
        home / "lammps-ani-src",
        home / "src" / "lammps-ani",
        Path("/root/src/lammps-ani"),
        REPO_ROOT.parent / "lammps-ani-src",
    ]
    seen: set[str] = set()
    for candidate in candidates:
        try:
            resolved = candidate.expanduser().resolve()
        except (OSError, RuntimeError):
            continue
        key = str(resolved)
        if key in seen:
            continue
        seen.add(key)
        if (resolved / plugin).is_file():
            return key
    return str((home / "lammps-ani-src").resolve())


def apply_default_lammps_ani_environment() -> None:
    """Populate conservative default env vars for the ANI workflow."""
    if os.environ.get("LAMMPS_SKIP_AUTO_ENV", "").strip() == "1":
        return

    home = Path.home()

    if not os.environ.get("LAMMPS_PREFIX", "").strip():
        os.environ["LAMMPS_PREFIX"] = str(home / ".local-lammps-ani")

    if not os.environ.get("LAMMPS_ANI_ROOT", "").strip():
        os.environ["LAMMPS_ANI_ROOT"] = _default_lammps_ani_root(home)

    prefix = Path(os.environ["LAMMPS_PREFIX"])
    bin_dir = prefix / "bin"
    if bin_dir.is_dir():
        path = os.environ.get("PATH", "")
        resolved_bin = str(bin_dir.resolve())
        parts = [part for part in path.split(":") if part]
        if resolved_bin not in parts:
            os.environ["PATH"] = f"{resolved_bin}:{path}" if path else resolved_bin

    prepend_ld: list[str] = []
    lib_dir = prefix / "lib"
    if lib_dir.is_dir():
        prepend_ld.append(str(lib_dir.resolve()))

    try:
        import torch

        torch_lib = Path(torch.__file__).resolve().parent / "lib"
        if torch_lib.is_dir():
            prepend_ld.append(str(torch_lib))
        nccl_lib = Path(torch.__file__).resolve().parent.parent / "nvidia" / "nccl" / "lib"
        if nccl_lib.is_dir():
            prepend_ld.append(str(nccl_lib.resolve()))
    except ImportError:
        pass

    cuda_lib = Path("/usr/local/cuda/lib64")
    if cuda_lib.is_dir():
        prepend_ld.append(str(cuda_lib.resolve()))

    conda_or_venv = os.environ.get("CONDA_PREFIX", "").strip() or sys.prefix
    if conda_or_venv:
        env_lib = Path(conda_or_venv) / "lib"
        if env_lib.is_dir():
            prepend_ld.append(str(env_lib.resolve()))

    if not prepend_ld:
        return

    existing = os.environ.get("LD_LIBRARY_PATH", "").strip()
    existing_parts = [part for part in existing.split(":") if part] if existing else []
    for entry in prepend_ld:
        if entry not in existing_parts:
            existing_parts.insert(0, entry)
    os.environ["LD_LIBRARY_PATH"] = ":".join(existing_parts)


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
    text = os.path.expandvars(os.path.expanduser(str(value)))
    if not text.strip():
        return text
    if "${" in text:
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

    ani_cfg = cfg.get("ani")
    if isinstance(ani_cfg, dict) and ani_cfg.get("model_file") is not None:
        ani_cfg["model_file"] = resolve_config_file_path(ani_cfg["model_file"], root)

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
    "apply_default_lammps_ani_environment",
    "apply_md_viscosity_path_resolution",
    "get_required",
    "get_section",
    "load_yaml",
    "resolve_config_file_path",
    "resolve_lammps_executable",
]
