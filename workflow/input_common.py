"""Shared helpers for generating LAMMPS input files."""
from __future__ import annotations

from typing import Any


def _default_lammps_settings() -> dict[str, Any]:
    return {
        "pair_style": "lj/cut/coul/long 12.0",
        "pair_modify": "mix arithmetic tail yes",
        "bond_style": "harmonic",
        "angle_style": "harmonic",
        "dihedral_style": "fourier",
        "improper_style": "cvff",
        "kspace_style": "pppm 1.0e-5",
        "special_bonds": {
            "lj": [0.0, 0.0, 0.5],
            "coul": [0.0, 0.0, 0.8333],
        },
    }


def _resolve_lammps_settings(
    lammps_settings: dict[str, Any] | None,
) -> dict[str, Any]:
    merged = _default_lammps_settings()
    if not lammps_settings:
        return merged
    merged.update({k: v for k, v in lammps_settings.items() if k != "special_bonds"})
    if "special_bonds" in lammps_settings:
        merged["special_bonds"] = lammps_settings["special_bonds"]
    return merged


def _build_forcefield_block(lammps_settings: dict[str, Any] | None) -> str:
    settings = _resolve_lammps_settings(lammps_settings)
    lj = settings["special_bonds"]["lj"]
    coul = settings["special_bonds"]["coul"]
    return "\n".join(
        [
            f"pair_style      {settings['pair_style']}",
            f"bond_style      {settings['bond_style']}",
            f"angle_style     {settings['angle_style']}",
            f"dihedral_style  {settings['dihedral_style']}",
            f"improper_style  {settings['improper_style']}",
            "",
            f"pair_modify     {settings['pair_modify']}",
            f"kspace_style    {settings['kspace_style']}",
            f"special_bonds   lj {lj[0]} {lj[1]} {lj[2]} coul {coul[0]} {coul[1]} {coul[2]}",
        ]
    )


def _build_restart_forcefield_block(lammps_settings: dict[str, Any] | None) -> str:
    settings = _resolve_lammps_settings(lammps_settings)
    lj = settings["special_bonds"]["lj"]
    coul = settings["special_bonds"]["coul"]
    return "\n".join(
        [
            f"pair_modify     {settings['pair_modify']}",
            f"kspace_style    {settings['kspace_style']}",
            f"special_bonds   lj {lj[0]} {lj[1]} {lj[2]} coul {coul[0]} {coul[1]} {coul[2]}",
        ]
    )


def _render_line(prefix: str, parts: list[str]) -> str:
    lines: list[str] = []
    current = prefix
    for part in parts:
        candidate = f"{current} {part}".rstrip()
        if len(candidate) <= 88:
            current = candidate
            continue
        lines.append(f"{current} &")
        current = "                " + part
    lines.append(current)
    return "\n".join(lines)


def build_render_last_frame_text(*, render_cfg: dict[str, Any], filename: str) -> str:
    if not bool(render_cfg.get("render_last_frame", False)):
        return ""

    width = int(render_cfg.get("image_width", 1400))
    height = int(render_cfg.get("image_height", 1000))
    view_theta = float(render_cfg.get("view_theta_deg", 70.0))
    view_phi = float(render_cfg.get("view_phi_deg", 30.0))
    center_x = float(render_cfg.get("center_x_frac", 0.5))
    center_y = float(render_cfg.get("center_y_frac", 0.5))
    center_z = float(render_cfg.get("center_z_frac", 0.5))
    zoom = float(render_cfg.get("zoom", 1.6))
    draw_box = "yes" if bool(render_cfg.get("draw_box", True)) else "no"
    box_line_diam = float(render_cfg.get("box_line_diam_frac", 0.02))
    atom_diameter = float(render_cfg.get("atom_diameter_A", 1.2))
    background_color = str(render_cfg.get("background_color", "white"))
    box_color = str(render_cfg.get("box_color", "gray"))

    parts = [
        "type",
        "type",
        "size",
        str(width),
        str(height),
        "view",
        f"{view_theta:.3f}",
        f"{view_phi:.3f}",
        "center",
        "s",
        f"{center_x:.3f}",
        f"{center_y:.3f}",
        f"{center_z:.3f}",
        "zoom",
        f"{zoom:.3f}",
        "box",
        draw_box,
        f"{box_line_diam:.4f}",
        "modify",
        "backcolor",
        background_color,
        "boxcolor",
        box_color,
        "adiam",
        "*",
        f"{atom_diameter:.3f}",
    ]

    atom_colors = render_cfg.get("atom_colors", {})
    if isinstance(atom_colors, dict):
        for atom_type, color in atom_colors.items():
            parts.extend(["acolor", str(atom_type), str(color)])

    atom_diameters = render_cfg.get("atom_diameters", {})
    if isinstance(atom_diameters, dict):
        for atom_type, diameter in atom_diameters.items():
            parts.extend(["adiam", str(atom_type), str(float(diameter))])

    command = _render_line(f"write_dump      all image {filename}", parts)
    return (
        "# ================================================================\n"
        "# Optional final-frame render\n"
        "# Color-by-type uses visualization.atom_colors from config.yaml\n"
        "# ================================================================\n"
        f"{command}\n"
    )


__all__ = [
    "_build_forcefield_block",
    "_build_restart_forcefield_block",
    "build_render_last_frame_text",
]
