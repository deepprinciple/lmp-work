"""Configuration helpers for thermal workflow entrypoints."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required for thermal workflow configs.") from exc

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


__all__ = ["get_required", "get_section", "load_yaml"]
