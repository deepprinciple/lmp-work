# Green-Kubo Workflow

Main entrypoint:

```bash
python runlmp.py --input config.yaml
```

The workflow is now:

1. Build a liquid model from SMILES
2. Generate `system.data`
3. Write one common equilibration input
4. Run minimization + NPT + NVT once to get `equil_nvt.restart`
5. Generate multiple replica inputs
6. Run short decorrelation + NVE+GK for each replica
7. Average replica results and write `result.csv`

## Main files

- `runlmp.py`: one-shot workflow entrypoint
- `config.yaml`: active config file
- `config.example.yaml`: config template
- `build_system.py`: build `system.data` from SMILES
- `write_lammps_input.py`: generate equilibration and replica input files
- `analyze_hfacf_vacf.py`: standalone post-processing script

## Quick start

```bash
python runlmp.py --input config.yaml
```

## Replica behavior

The workflow does not rerun the whole system preparation for each replica.

Instead it does:

- one shared equilibration in `case.workdir`
- then `simulation.n_replicas` independent replica runs in:
  - `replica_01/`
  - `replica_02/`
  - `replica_03/`
  - ...

Each replica starts from `equil_nvt.restart`, gets a different velocity seed, runs a short NVT decorrelation, then runs NVE + Green-Kubo production.

## Output structure

Inside `case.workdir`, the important outputs are:

- `system.data`
- `in.equil.lammps`
- `equil_nvt.restart`
- `result.csv`
- `replica_01/`
- `replica_02/`
- ...

Inside each replica directory:

- `in.replica.lammps`
- `hfacf.dat`
- `vacf.dat`
- `gk_data.dat`
- `prod_final.restart`

## Result format

`result.csv` contains:

- one `summary` row with replica-averaged transport properties
- one `replica` row for each replica

Summary row fields include:

- `eta_cP`
- `eta_sem_cP` (replica-to-replica SEM)
- `kappa_W_mK`
- `kappa_sem_W_mK` (replica-to-replica SEM)

## Important config fields

### `simulation`

- `n_replicas`: number of replicas, default `3`
- `replica_seed_base`: first replica seed
- `replica_seed_step`: seed increment between replicas
- `replica_decorrelation_steps`: short NVT decorrelation length
- `gk_steps`: production length per replica

### `analysis`

- `output`: default `result.csv`
- `visc_t_max_ps`
- `kappa_t_max_ps`
- `smooth_window`
- `n_blocks`

## Requirement

`runlmp.py` requires `PyYAML`:

```bash
pip install pyyaml
```
