"""Generate LAMMPS input files for the Green-Kubo workflow."""
from __future__ import annotations

import argparse
from pathlib import Path


def build_equilibration_input_text(
    *,
    temperature: float,
    seed: int,
    npt_pre_50_steps: int,
    npt_pre_10_steps: int,
    npt_1atm_steps: int,
    nvt_steps: int,
    timestep_fs: float,
) -> str:
    return f"""# ================================================================
# LAMMPS Input: common equilibration
# Generated automatically by write_lammps_input.py
# Flow: system.data -> minimization -> NPT -> NVT
# ================================================================

units           real
atom_style      full
boundary        p p p

pair_style      lj/cut/coul/long 12.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  fourier
improper_style  cvff

read_data       system.data

pair_modify     mix arithmetic tail yes
kspace_style    pppm 1.0e-5
special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333

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
# Stage 2: NPT pre-compression at 50 atm
# ================================================================
fix             f_npt_50 all npt temp ${{T}} ${{T}} 100.0 iso 50.0 50.0 1000.0 drag 2.0
run             {npt_pre_50_steps}
unfix           f_npt_50
write_restart   stage2_npt_50atm.restart
write_data      stage2_npt_50atm.data

# ================================================================
# Stage 3: NPT relaxation at 10 atm
# ================================================================
fix             f_npt_10 all npt temp ${{T}} ${{T}} 100.0 iso 10.0 10.0 1000.0 drag 2.0
run             {npt_pre_10_steps}
unfix           f_npt_10
write_restart   stage3_npt_10atm.restart
write_data      stage3_npt_10atm.data

# ================================================================
# Stage 4: NPT equilibration at 1 atm
# ================================================================
fix             f_npt_1 all npt temp ${{T}} ${{T}} 100.0 iso 1.0 1.0 1000.0 drag 2.0
run             {npt_1atm_steps}
unfix           f_npt_1
write_restart   equil_npt.restart
write_data      equil_npt.data

# ================================================================
# Stage 5: NVT equilibration at fixed volume
# ================================================================
fix             nvtlog all ave/time 100 10 1000 v_t_inst v_p_inst v_rho_inst v_vol_inst v_pe_inst v_ke_inst v_et_inst file thermo_nvt.dat
fix             f_nvt_eq all nvt temp ${{T}} ${{T}} 100.0
run             {nvt_steps}
unfix           f_nvt_eq
unfix           nvtlog
write_restart   equil_nvt.restart
write_data      equil_nvt.data
unfix           nptlog
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
) -> str:
    return f"""# ================================================================
# LAMMPS Input: replica Green-Kubo production
# Generated automatically by write_lammps_input.py
# Flow: equil_nvt.restart -> short NVT decorrelation -> NVE + GK
# ================================================================

units           real
atom_style      full
boundary        p p p

read_restart    {equil_restart_relpath}

pair_modify     mix arithmetic tail yes
kspace_style    pppm 1.0e-5
special_bonds   lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333

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
"""


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
