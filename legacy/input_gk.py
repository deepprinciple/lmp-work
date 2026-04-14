"""Legacy GK replica input builder retained for manual HFACF generation."""

from __future__ import annotations

from typing import Any

from workflow.input_common import (
    _build_restart_forcefield_block,
    build_render_last_frame_text,
)


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
# Generated automatically by legacy/input_gk.py
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


__all__ = ["build_replica_input_text"]
