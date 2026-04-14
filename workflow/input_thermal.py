"""Thermal-workflow LAMMPS input builders."""
from __future__ import annotations

from typing import Any

from .input_common import (
    _build_forcefield_block,
    _build_restart_forcefield_block,
    build_render_last_frame_text,
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
# Generated automatically by workflow/input_thermal.py
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
# Generated automatically by workflow/input_thermal.py
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


__all__ = [
    "build_equilibration_input_text",
    "build_rnemd_replica_input_text",
]
