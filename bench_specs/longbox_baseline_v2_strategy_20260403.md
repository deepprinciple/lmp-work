# Long-Box Baseline v2

This baseline keeps the original square-box EMD/GK workflow intact as the rollback path for viscosity fitting, while defining a shared elongated-box branch for split-property production.

## Current baseline

- Shared model: orthorhombic liquid box, default aspect ratio `1:1:3`.
- Shared equilibration: one elongated-box NPT/NVT preparation, saved as `equil_nvt.restart`.
- Thermal conductivity production: `fix thermal/conductivity` on the long box.
- Viscosity production: current accepted path is Green-Kubo viscosity.
- Fallback path: existing square-box GK viscosity pipeline remains untouched.

## Why keep the split production

- `fix thermal/conductivity` and `fix viscosity` both impose non-equilibrium fluxes.
- Running them together in one production trajectory is not a safe baseline because flux coupling can distort both gradients.
- Shared build and equilibration are still useful, so the baseline shares those stages and separates only the production/analysis stage.

## Current methanol baseline assets

- Thermal rNEMD branch: available under `methanol_500mol_thermal_rnemd_longbox_20260403/thermal_conductivity/method_rnemd/`.
- Long-box GK viscosity control: available under `methanol_500mol_thermal_rnemd_longbox_20260403/viscosity/method_gk/`.

## Archived exploratory branch

- Long-box `fix viscosity` control was tested on methanol and archived because it did not improve convergence relative to GK and remained noticeably below experiment within the 500 ps control run.

## Current decision

Use this as the active baseline for the 5-molecule validation set:

- viscosity via GK
- thermal conductivity via `fix thermal/conductivity`

Keep `fix viscosity` outside the baseline unless a later dedicated test clearly outperforms GK.
