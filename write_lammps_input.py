"""Generate LAMMPS input files for the Green-Kubo workflow."""
from __future__ import annotations

import argparse
from pathlib import Path
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


def _resolve_lammps_settings(lammps_settings: dict[str, Any] | None) -> dict[str, Any]:
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


def build_equilibration_input_text(
    *,
    temperature: float,
    seed: int,
    npt_pre_50_steps: int,
    npt_pre_10_steps: int,
    npt_1atm_steps: int,
    nvt_steps: int,
    timestep_fs: float,
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("equil_image_file", "equil_final.png")),
    )
    return f"""# ================================================================
# LAMMPS Input: common equilibration
# Generated automatically by write_lammps_input.py
# Flow: system.data -> minimization -> NPT(1 atm) -> NVT
# ================================================================

units           real
atom_style      full
boundary        p p p

{forcefield_block}

read_data       system.data

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {seed}

thermo          500
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

# ================================================================
# Stage 1: Energy minimization
# ================================================================
min_style       cg
minimize        1.0e-6 1.0e-8 50000 500000
write_restart   stage1_minimized.restart
write_data      stage1_minimized.data
reset_timestep  0

velocity        all create ${{T}} ${{seed}} rot yes dist gaussian
fix             mom all momentum 100 linear 1 1 1

variable        t_inst   equal temp
variable        p_inst   equal press
variable        rho_inst equal density
variable        vol_inst equal vol
variable        pe_inst  equal pe
variable        ke_inst  equal ke
variable        et_inst  equal etotal

fix             nptlog all ave/time 100 10 1000 v_t_inst v_p_inst v_rho_inst v_vol_inst file density_npt.dat

# ================================================================
# Stage 2: NPT equilibration at 1 atm
# ================================================================
fix             f_npt_1 all npt temp ${{T}} ${{T}} 100.0 iso 1.0 1.0 1000.0 drag 2.0
run             {npt_1atm_steps}
unfix           f_npt_1
write_restart   equil_npt.restart
write_data      equil_npt.data

# ================================================================
# Stage 3: NVT equilibration at fixed volume
# ================================================================
fix             nvtlog all ave/time 100 10 1000 v_t_inst v_p_inst v_rho_inst v_vol_inst v_pe_inst v_ke_inst v_et_inst file thermo_nvt.dat
fix             f_nvt_eq all nvt temp ${{T}} ${{T}} 100.0
run             {nvt_steps}
unfix           f_nvt_eq
unfix           nvtlog
write_restart   equil_nvt.restart
write_data      equil_nvt.data
{render_block}unfix           nptlog
"""


def build_replica_input_text(
    *,
    temperature: float,
    replica_seed: int,
    decorrelation_steps: int,
    gk_steps: int,
    timestep_fs: float,
    sample_every_steps: int,
    corr_points: int,
    equil_restart_relpath: str = "../equil_nvt.restart",
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )
    return f"""# ================================================================
# LAMMPS Input: replica Green-Kubo production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> NVE + GK
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}

thermo          2000
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

fix             mom all momentum 100 linear 1 1 1
velocity        all create ${{T}} ${{seed}} rot yes dist gaussian

# ================================================================
# Stage R1: short NVT decorrelation
# ================================================================
fix             f_nvt_decor all nvt temp ${{T}} ${{T}} 100.0
run             {decorrelation_steps}
unfix           f_nvt_decor
write_restart   replica_init.restart
write_data      replica_init.data
reset_timestep  0

# ================================================================
# Stage R2: NVE + Green-Kubo production
# ================================================================
unfix           mom
fix             f_nve all nve

compute         myKE     all ke/atom
compute         myPE     all pe/atom
compute         myStress all stress/atom NULL virial
compute         myFlux   all heat/flux myKE myPE myStress

variable        Jx equal c_myFlux[1]/vol
variable        Jy equal c_myFlux[2]/vol
variable        Jz equal c_myFlux[3]/vol

variable        Nev_k  equal {sample_every_steps}
variable        Nrep_k equal {corr_points}
variable        Nfq_k  equal ${{Nev_k}}*${{Nrep_k}}

fix             hfacf all ave/correlate ${{Nev_k}} ${{Nrep_k}} ${{Nfq_k}} &
                v_Jx v_Jy v_Jz &
                type auto file hfacf.dat ave running

variable        pxy_v equal pxy
variable        pxz_v equal pxz
variable        pyz_v equal pyz

variable        Nev_v  equal {sample_every_steps}
variable        Nrep_v equal {corr_points}
variable        Nfq_v  equal ${{Nev_v}}*${{Nrep_v}}

fix             vacf all ave/correlate ${{Nev_v}} ${{Nrep_v}} ${{Nfq_v}} &
                v_pxy_v v_pxz_v v_pyz_v &
                type auto file vacf.dat ave running

compute         dip all dipole
variable        Mx equal c_dip[1]
variable        My equal c_dip[2]
variable        Mz equal c_dip[3]

variable        sample_every equal 20
fix             gkdata all print ${{sample_every}} &
                "$(step) $(temp) $(vol) $(v_pxy_v) $(v_pxz_v) $(v_pyz_v) $(v_Jx) $(v_Jy) $(v_Jz) $(v_Mx) $(v_My) $(v_Mz)" &
                file gk_data.dat screen no title "# step T V pxy pxz pyz Jx Jy Jz Mx My Mz"

run             {gk_steps}

write_restart   prod_final.restart
write_data      prod_final.data
{render_block}"""


def build_initial_equilibration_input_text(
    *,
    temperature: float,
    seed: int,
    npt_pre_50_steps: int,
    npt_pre_10_steps: int,
    timestep_fs: float,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_forcefield_block(lammps_settings)
    return f"""# ================================================================
# LAMMPS Input: adaptive equilibration pre-stage
# Flow: system.data -> minimization
# ================================================================

units           real
atom_style      full
boundary        p p p

{forcefield_block}

read_data       system.data

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {seed}

thermo          500
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

min_style       cg
minimize        1.0e-6 1.0e-8 50000 500000
write_restart   stage1_minimized.restart
write_data      stage1_minimized.data
reset_timestep  0

velocity        all create ${{T}} ${{seed}} rot yes dist gaussian
fix             mom all momentum 100 linear 1 1 1
write_restart   stage1_minimized.restart
"""


def build_npt_segment_input_text(
    *,
    temperature: float,
    timestep_fs: float,
    run_steps: int,
    input_restart_relpath: str,
    output_restart_name: str,
    density_file: str = "density_npt_segment.dat",
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    return f"""# ================================================================
# LAMMPS Input: adaptive NPT segment
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {input_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}

thermo          250
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

variable        t_inst   equal temp
variable        p_inst   equal press
variable        rho_inst equal density
variable        vol_inst equal vol
fix             nptlog all ave/time 50 10 500 v_t_inst v_p_inst v_rho_inst v_vol_inst file {density_file}
fix             f_npt_1 all npt temp ${{T}} ${{T}} 100.0 iso 1.0 1.0 1000.0 drag 2.0
run             {run_steps}
unfix           f_npt_1
unfix           nptlog
write_restart   {output_restart_name}
write_data      equil_npt.data
"""


def build_nvt_segment_input_text(
    *,
    temperature: float,
    timestep_fs: float,
    run_steps: int,
    input_restart_relpath: str,
    output_restart_name: str,
    thermo_file: str = "thermo_nvt_segment.dat",
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    return f"""# ================================================================
# LAMMPS Input: adaptive NVT segment
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {input_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}

thermo          250
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

variable        t_inst   equal temp
variable        pe_inst  equal pe
variable        ke_inst  equal ke
variable        et_inst  equal etotal
fix             nvtlog all ave/time 50 10 500 v_t_inst v_pe_inst v_ke_inst v_et_inst file {thermo_file}
fix             f_nvt_eq all nvt temp ${{T}} ${{T}} 100.0
run             {run_steps}
unfix           f_nvt_eq
unfix           nvtlog
write_restart   {output_restart_name}
write_data      equil_nvt.data
"""


def build_nve_gk_segment_input_text(
    *,
    temperature: float,
    replica_seed: int,
    timestep_fs: float,
    run_steps: int,
    sample_every_steps: int,
    corr_points: int,
    input_restart_relpath: str,
    output_restart_name: str,
    hfacf_file: str,
    vacf_file: str,
    gk_data_file: str,
    thermo_file: str,
    initialize_velocities: bool = False,
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )
    init_velocity_block = ""
    if initialize_velocities:
        init_velocity_block = (
            "velocity        all create ${T} ${seed} rot yes dist gaussian\n"
            "fix             mom all momentum 100 linear 1 1 1\n"
            "run             0\n"
            "unfix           mom\n"
        )
    return f"""# ================================================================
# LAMMPS Input: replica NVE Green-Kubo segment
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {input_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}

thermo          1000
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

{init_velocity_block}fix             f_nve all nve

compute         myKE     all ke/atom
compute         myPE     all pe/atom
compute         myStress all stress/atom NULL virial
compute         myFlux   all heat/flux myKE myPE myStress

variable        Jx equal c_myFlux[1]/vol
variable        Jy equal c_myFlux[2]/vol
variable        Jz equal c_myFlux[3]/vol

variable        Nev_k  equal {sample_every_steps}
variable        Nrep_k equal {corr_points}
variable        Nfq_k  equal ${{Nev_k}}*${{Nrep_k}}

fix             hfacf all ave/correlate ${{Nev_k}} ${{Nrep_k}} ${{Nfq_k}} &
                v_Jx v_Jy v_Jz &
                type auto file {hfacf_file} ave running

variable        pxy_v equal pxy
variable        pxz_v equal pxz
variable        pyz_v equal pyz

variable        Nev_v  equal {sample_every_steps}
variable        Nrep_v equal {corr_points}
variable        Nfq_v  equal ${{Nev_v}}*${{Nrep_v}}

fix             vacf all ave/correlate ${{Nev_v}} ${{Nrep_v}} ${{Nfq_v}} &
                v_pxy_v v_pxz_v v_pyz_v &
                type auto file {vacf_file} ave running

compute         dip all dipole
variable        Mx equal c_dip[1]
variable        My equal c_dip[2]
variable        Mz equal c_dip[3]

variable        sample_every equal 20
fix             gkdata all print ${{sample_every}} &
                "$(step) $(temp) $(vol) $(v_pxy_v) $(v_pxz_v) $(v_pyz_v) $(v_Jx) $(v_Jy) $(v_Jz) $(v_Mx) $(v_My) $(v_Mz)" &
                file {gk_data_file} screen no title "# step T V pxy pxz pyz Jx Jy Jz Mx My Mz"

fix             thermodat all print 200 &
                "$(step) $(temp) $(pe) $(ke) $(etotal) $(vol)" &
                file {thermo_file} screen no title "# step temp pe ke etotal vol"

run             {run_steps}

unfix           f_nve
unfix           hfacf
unfix           vacf
unfix           gkdata
unfix           thermodat
write_restart   {output_restart_name}
write_data      prod_final.data
{render_block}"""


def build_rnemd_replica_input_text(
    *,
    temperature: float,
    replica_seed: int,
    decorrelation_steps: int,
    run_steps: int,
    timestep_fs: float,
    swap_every_steps: int,
    nbin: int,
    edim: str = "z",
    nswap: int = 1,
    profile_every_steps: int = 10,
    profile_repeat: int = 100,
    profile_freq: int = 1000,
    profile_file: str = "temp_profile.dat",
    exchange_file: str = "thermal_exchange.dat",
    thermo_file: str = "thermal_rnemd_thermo.dat",
    equil_restart_relpath: str = "../equil_nvt.restart",
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    if edim not in {"x", "y", "z"}:
        raise ValueError("edim 必须是 x / y / z")
    if nbin % 2 != 0:
        raise ValueError("fix thermal/conductivity 要求 nbin 为偶数")

    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )
    return f"""# ================================================================
# LAMMPS Input: replica reverse-NEMD thermal conductivity production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> NVE + fix thermal/conductivity
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}
variable        Nbin equal {nbin}
variable        kB   equal 0.0019872041

thermo          1000
thermo_style    custom step temp press density vol lx ly lz pe ke etotal
thermo_modify   flush yes

fix             mom all momentum 100 linear 1 1 1
velocity        all create ${{T}} ${{seed}} rot yes dist gaussian

# ================================================================
# Stage R1: short NVT decorrelation
# ================================================================
fix             f_nvt_decor all nvt temp ${{T}} ${{T}} 100.0
run             {decorrelation_steps}
unfix           f_nvt_decor
write_restart   replica_init.restart
write_data      replica_init.data
reset_timestep  0

# ================================================================
# Stage R2: NVE + reverse NEMD thermal conductivity production
# ================================================================
unfix           mom
fix             f_nve all nve
fix             tc all thermal/conductivity {swap_every_steps} {edim} {nbin} swap {nswap}

compute         ke_atom all ke/atom
variable        temp_atom atom c_ke_atom/(1.5*${{kB}})
variable        binw equal l{edim}/${{Nbin}}
compute         layers all chunk/atom bin/1d {edim} lower ${{binw}} units box

fix             tprof all ave/chunk {profile_every_steps} {profile_repeat} {profile_freq} layers v_temp_atom file {profile_file}

variable        e_xfer equal f_tc
variable        sample_every equal {swap_every_steps}
fix             tcprint all print ${{sample_every}} &
                "$(step) $(temp) $(vol) $(lx) $(ly) $(lz) $(v_e_xfer)" &
                file {exchange_file} screen no title "# step T V lx ly lz exchanged_energy"

fix             thermodat all print 200 &
                "$(step) $(temp) $(pe) $(ke) $(etotal) $(vol) $(lx) $(ly) $(lz)" &
                file {thermo_file} screen no title "# step temp pe ke etotal vol lx ly lz"

run             {run_steps}

unfix           f_nve
unfix           tc
unfix           tprof
unfix           tcprint
unfix           thermodat
write_restart   prod_final.restart
write_data      prod_final.data
{render_block}"""


def build_viscosity_rnemd_input_text(
    *,
    temperature: float,
    replica_seed: int,
    decorrelation_steps: int,
    run_steps: int,
    timestep_fs: float,
    swap_every_steps: int,
    nbin: int,
    vdim: str = "x",
    pdim: str = "z",
    nswap: int = 1,
    vtarget: float | None = None,
    profile_every_steps: int = 10,
    profile_repeat: int = 100,
    profile_freq: int = 1000,
    profile_file: str = "velocity_profile.dat",
    exchange_file: str = "momentum_exchange.dat",
    thermo_file: str = "viscosity_rnemd_thermo.dat",
    equil_restart_relpath: str = "../equil_nvt.restart",
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    if vdim not in {"x", "y", "z"}:
        raise ValueError("vdim 必须是 x / y / z")
    if pdim not in {"x", "y", "z"}:
        raise ValueError("pdim 必须是 x / y / z")
    if vdim == pdim:
        raise ValueError("fix viscosity 要求 vdim 和 pdim 不同")
    if nbin % 2 != 0:
        raise ValueError("fix viscosity 要求 nbin 为偶数")

    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )
    velocity_component = {"x": "vx", "y": "vy", "z": "vz"}[vdim]
    fix_viscosity = f"fix             mp all viscosity {swap_every_steps} {vdim} {pdim} {nbin} swap {nswap}"
    if vtarget is not None:
        fix_viscosity += f" vtarget {float(vtarget):.6f}"

    return f"""# ================================================================
# LAMMPS Input: replica reverse-NEMD viscosity production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> NVE + fix viscosity
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}
variable        Nbin equal {nbin}

thermo          1000
thermo_style    custom step temp press density vol lx ly lz pe ke etotal
thermo_modify   flush yes

fix             mom all momentum 100 linear 1 1 1
velocity        all create ${{T}} ${{seed}} rot yes dist gaussian

# ================================================================
# Stage R1: short NVT decorrelation
# ================================================================
fix             f_nvt_decor all nvt temp ${{T}} ${{T}} 100.0
run             {decorrelation_steps}
unfix           f_nvt_decor
write_restart   replica_init.restart
write_data      replica_init.data
reset_timestep  0

# ================================================================
# Stage R2: NVE + reverse NEMD viscosity production
# ================================================================
unfix           mom
fix             f_nve all nve
{fix_viscosity}

variable        binw equal l{pdim}/${{Nbin}}
compute         layers all chunk/atom bin/1d {pdim} lower ${{binw}} units box
fix             vprof all ave/chunk {profile_every_steps} {profile_repeat} {profile_freq} layers {velocity_component} file {profile_file}

variable        p_xfer equal f_mp
variable        sample_every equal {swap_every_steps}
fix             mpprint all print ${{sample_every}} &
                "$(step) $(temp) $(vol) $(lx) $(ly) $(lz) $(v_p_xfer)" &
                file {exchange_file} screen no title "# step T V lx ly lz exchanged_momentum"

fix             thermodat all print 200 &
                "$(step) $(temp) $(pe) $(ke) $(etotal) $(vol) $(lx) $(ly) $(lz)" &
                file {thermo_file} screen no title "# step temp pe ke etotal vol lx ly lz"

run             {run_steps}

unfix           f_nve
unfix           mp
unfix           vprof
unfix           mpprint
unfix           thermodat
write_restart   prod_final.restart
write_data      prod_final.data
{render_block}"""


def build_viscosity_einstein_input_text(
    *,
    temperature: float,
    replica_seed: int,
    decorrelation_steps: int,
    run_steps: int,
    timestep_fs: float,
    sample_every_steps: int,
    equil_restart_relpath: str = "../equil_nvt.restart",
    thermostat: str = "nve",
    tensor_file: str = "pressure_tensor.dat",
    thermo_file: str = "einstein_thermo.dat",
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    if thermostat not in {"nve", "nvt"}:
        raise ValueError("Einstein 粘度生产段 thermostat 必须是 nve 或 nvt")

    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )
    if thermostat == "nve":
        production_fix = "fix             f_prod all nve"
    else:
        production_fix = "fix             f_prod all nvt temp ${T} ${T} 100.0"

    return f"""# ================================================================
# LAMMPS Input: replica Einstein viscosity production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> equilibrium production
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}
variable        sample_every equal {sample_every_steps}
variable        step_v equal step
variable        temp_v equal temp
variable        press_v equal press
variable        density_v equal density
variable        vol_v equal vol
variable        lx_v equal lx
variable        ly_v equal ly
variable        lz_v equal lz
variable        etotal_v equal etotal

thermo          1000
thermo_style    custom step temp press density vol lx ly lz pe ke etotal pxy pxz pyz
thermo_modify   flush yes

fix             mom all momentum 100 linear 1 1 1
velocity        all create ${{T}} ${{seed}} rot yes dist gaussian

# ================================================================
# Stage E1: short NVT decorrelation
# ================================================================
fix             f_nvt_decor all nvt temp ${{T}} ${{T}} 100.0
run             {decorrelation_steps}
unfix           f_nvt_decor
write_restart   replica_init.restart
write_data      replica_init.data
reset_timestep  0

# ================================================================
# Stage E2: equilibrium production for Einstein viscosity
# ================================================================
unfix           mom
{production_fix}

variable        pxy_v equal pxy
variable        pxz_v equal pxz
variable        pyz_v equal pyz
fix             ptdat all print ${{sample_every}} &
                "${{step_v}} ${{temp_v}} ${{vol_v}} ${{pxy_v}} ${{pxz_v}} ${{pyz_v}}" &
                file {tensor_file} screen no title "# step T V pxy pxz pyz"

fix             thermodat all print 200 &
                "${{step_v}} ${{temp_v}} ${{press_v}} ${{density_v}} ${{vol_v}} ${{lx_v}} ${{ly_v}} ${{lz_v}} ${{etotal_v}}" &
                file {thermo_file} screen no title "# step temp press density vol lx ly lz etotal"

run             {run_steps}

unfix           f_prod
unfix           ptdat
unfix           thermodat
write_restart   prod_final.restart
write_data      prod_final.data
{render_block}"""


def build_viscosity_ppm_input_text(
    *,
    temperature: float,
    replica_seed: int,
    decorrelation_steps: int,
    run_steps: int,
    timestep_fs: float,
    accel_amplitude_A_fs2: float,
    pdim: str = "z",
    vdim: str = "x",
    nbin: int = 20,
    profile_every_steps: int = 10,
    profile_repeat: int = 100,
    profile_freq: int = 1000,
    profile_file: str = "velocity_profile.dat",
    thermo_file: str = "ppm_thermo.dat",
    equil_restart_relpath: str = "../equil_nvt.restart",
    render_cfg: dict[str, Any] | None = None,
    lammps_settings: dict[str, Any] | None = None,
) -> str:
    if pdim not in {"x", "y", "z"} or vdim not in {"x", "y", "z"}:
        raise ValueError("PPM 要求 pdim / vdim 只能是 x / y / z")
    if pdim == vdim:
        raise ValueError("PPM 要求驱动方向和响应速度方向不同")
    if nbin <= 0:
        raise ValueError("PPM 要求 nbin 为正整数")

    forcefield_block = _build_restart_forcefield_block(lammps_settings)
    render_block = build_render_last_frame_text(
        render_cfg=render_cfg or {},
        filename=str((render_cfg or {}).get("replica_image_file", "final_frame.png")),
    )

    coord_var = {"x": "x", "y": "y", "z": "z"}[pdim]
    lo_var = {"x": "xlo", "y": "ylo", "z": "zlo"}[pdim]
    box_var = {"x": "lx", "y": "ly", "z": "lz"}[pdim]
    force_var = {"x": "fx", "y": "fy", "z": "fz"}[vdim]
    velocity_component = {"x": "vx", "y": "vy", "z": "vz"}[vdim]

    return f"""# ================================================================
# LAMMPS Input: replica PPM viscosity production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> NVT + periodic body force
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

{forcefield_block}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        {timestep_fs:.6f}

variable        T equal {temperature:.6f}
variable        seed equal {replica_seed}
variable        Nbin equal {nbin}
variable        accel_amp equal {accel_amplitude_A_fs2:.12g}
variable        force_pref equal 2390.057361376673
variable        step_v equal step
variable        temp_v equal temp
variable        press_v equal press
variable        density_v equal density
variable        vol_v equal vol
variable        lx_v equal lx
variable        ly_v equal ly
variable        lz_v equal lz
variable        etotal_v equal etotal

thermo          1000
thermo_style    custom step temp press density vol lx ly lz pe ke etotal
thermo_modify   flush yes

fix             mom all momentum 100 linear 1 1 1
velocity        all create ${{T}} ${{seed}} rot yes dist gaussian

# ================================================================
# Stage P1: short NVT decorrelation
# ================================================================
fix             f_nvt_decor all nvt temp ${{T}} ${{T}} 100.0
run             {decorrelation_steps}
unfix           f_nvt_decor
write_restart   replica_init.restart
write_data      replica_init.data
reset_timestep  0

# ================================================================
# Stage P2: periodic perturbation method production
# ================================================================
unfix           mom
fix             f_nvt all nvt temp ${{T}} ${{T}} 100.0
variable        accel atom ${{accel_amp}}*cos(2.0*PI*(({coord_var}-{lo_var})/{box_var}))
variable        {force_var}_atom atom ${{force_pref}}*mass(type)*v_accel
fix             drive all addforce {"v_fx_atom" if vdim == "x" else "0.0"} {"v_fy_atom" if vdim == "y" else "0.0"} {"v_fz_atom" if vdim == "z" else "0.0"}

variable        binw equal l{pdim}/${{Nbin}}
compute         layers all chunk/atom bin/1d {pdim} lower ${{binw}} units box
fix             vprof all ave/chunk {profile_every_steps} {profile_repeat} {profile_freq} layers {velocity_component} file {profile_file}

fix             thermodat all print 200 &
                "${{step_v}} ${{temp_v}} ${{press_v}} ${{density_v}} ${{vol_v}} ${{lx_v}} ${{ly_v}} ${{lz_v}} ${{etotal_v}}" &
                file {thermo_file} screen no title "# step temp press density vol lx ly lz etotal"

run             {run_steps}

unfix           f_nvt
unfix           drive
unfix           vprof
unfix           thermodat
write_restart   prod_final.restart
write_data      prod_final.data
{render_block}"""


def build_input_text(**kwargs) -> str:
    """Backward-compatible alias for the old single-file workflow."""
    return build_equilibration_input_text(**kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write LAMMPS input files")
    parser.add_argument("--workdir", required=True, help="Case directory")
    parser.add_argument("--mode", choices=["equilibration", "replica"], default="equilibration")
    parser.add_argument("--temperature", type=float, default=300.0, help="Simulation temperature in K")
    parser.add_argument("--seed", type=int, default=20260316, help="Random seed")
    parser.add_argument("--npt_pre_50_steps", type=int, default=200000, help="50 atm NPT steps")
    parser.add_argument("--npt_pre_10_steps", type=int, default=300000, help="10 atm NPT steps")
    parser.add_argument("--npt_1atm_steps", type=int, default=2000000, help="1 atm NPT steps")
    parser.add_argument("--nvt_steps", type=int, default=1000000, help="NVT equilibration steps")
    parser.add_argument("--decorrelation_steps", type=int, default=500000, help="Replica decorrelation steps")
    parser.add_argument("--gk_steps", type=int, default=10000000, help="NVE+GK production steps")
    parser.add_argument("--timestep_fs", type=float, default=1.0, help="Time step in fs")
    parser.add_argument("--sample_every_steps", type=int, default=5, help="ACF sampling interval in steps")
    parser.add_argument("--corr_points", type=int, default=8000, help="Number of correlation points")
    parser.add_argument("--equil_restart_relpath", default="../equil_nvt.restart", help="Relative path to equilibration restart for replica mode")
    parser.add_argument("--out", default=None, help="Output input filename")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "equilibration":
        text = build_equilibration_input_text(
            temperature=args.temperature,
            seed=args.seed,
            npt_pre_50_steps=args.npt_pre_50_steps,
            npt_pre_10_steps=args.npt_pre_10_steps,
            npt_1atm_steps=args.npt_1atm_steps,
            nvt_steps=args.nvt_steps,
            timestep_fs=args.timestep_fs,
        )
        out_name = args.out or "in.equil.lammps"
    else:
        text = build_replica_input_text(
            temperature=args.temperature,
            replica_seed=args.seed,
            decorrelation_steps=args.decorrelation_steps,
            gk_steps=args.gk_steps,
            timestep_fs=args.timestep_fs,
            sample_every_steps=args.sample_every_steps,
            corr_points=args.corr_points,
            equil_restart_relpath=args.equil_restart_relpath,
        )
        out_name = args.out or "in.replica.lammps"

    out_file = workdir / out_name
    out_file.write_text(text)
    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
