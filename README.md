# lmp-work — LAMMPS 输运性质计算工作流

基于 LAMMPS 的液体输运性质计算流水线，目前支持两条工作流：

| 工作流 | 方法 | 力场 | 入口 |
|--------|------|------|------|
| **热导率** | reverse-NEMD | OpenFF / LigParGen | `md_run.py` |
| **粘度** | Green-Kubo 压力张量 block ACF | ByteDance BAMBOO | `md_viscosity.py` |

---

## 环境搭建

推荐先创建独立 Python 环境：

```bash
conda create -n lmp-work python=3.11 -y
conda activate lmp-work
pip install -r requirements.txt
```

如果默认热导率路径使用 `OpenFF + AM1BCC` 电荷分配时报缺依赖，建议再补：

```bash
conda install -c conda-forge ambertools -y
```

```bash
conda install -c conda-forge rdkit openff-toolkit openff-interchange openff-units -y
```

除了 Python 包，还需要准备外部程序：

- `packmol`：两条工作流都需要，并且命令应在 `PATH` 中
- `lmp_mpi` 或等价 LAMMPS 可执行程序：两条工作流都需要
- BAMBOO-enabled LAMMPS + BAMBOO 模型文件：粘度工作流需要
- `BOSS`：仅当你选择 `LigParGen + BOSS` 路径时需要

## 热导率工作流

SMILES → OpenFF 参数化 → Packmol 建盒 → LAMMPS reverse-NEMD → κ

```bash
python md_run.py --config configs/config.yaml
```

短程 smoke test：

```bash
python md_run.py --config configs/config.smoke.yaml
```

### 核心文件（热导率）

- `md_run.py` — 主流程入口
- `workflow/input_thermal.py` — 平衡段和 rNEMD 输入生成
- `analysis/hfacf.py` — HFACF 热导率分析核心（也被粘度工作流复用）
- `analysis/rnemd.py` — rNEMD 热导率拟合与收敛分析
- `scripts/analyze_thermal_hfacf.py` — HFACF 后处理 CLI
- `configs/config.yaml` — 主配置模板

### 力场选择

- 默认：`forcefield.engine: openff`，`forcefield.charge_method: am1bcc`
- 可选：`forcefield.engine: ligpargen`（依赖 LigParGen/BOSS 运行环境）

### 授权与分发声明

- 本仓库不再包含 `BOSS` 程序及其兼容运行时副本。
- `BOSS` 为受限授权软件，使用者需自行向权利方申请授权并在本地安装。
- `LigParGen + BOSS` 路径仅保留接口能力；实际运行需通过本地 `wrapper_script`
  指向你自己的合法安装环境。

---

## 粘度工作流（BAMBOO Green-Kubo）

components/composition → RDKit/ion preset → 多组分 Packmol 建盒 → BAMBOO `atom_style full` data → NVT 预热 → NPT 压回目标密度 → 短 NVT 收尾 → NVE + 压力张量/应力 ACF → η

> **依赖**：需要使用 ByteDance BAMBOO 对应的 LAMMPS 构建，并准备 BAMBOO 模型文件。
> 默认配置会把 `/root/bamboo/pair/lammps/output` prepend 到 `PATH`，并给 LAMMPS 加
> `-k on g 1 -sf kk`；可在 `lammps:` 配置里关闭或覆盖。

```bash
python md_viscosity.py --config configs/viscosity.yaml
```

短程 smoke test：

```bash
python md_viscosity.py --config configs/viscosity.smoke.yaml
```

最小 replicate demo：

```yaml
replica:
  n_replicates: 3
  seed_stride:  1000
```

这会在 `case.workdir` 下生成 `replica_01/`、`replica_02/`、`replica_03/`，每个副本使用同一套
体系和时长参数，但 `simulation.seed` 会自动错开；顶层额外写出
`viscosity_replicas_summary.json`，汇总 `eta_replica_mean/std/sem`。

如果你要研究“模拟时长是否足够”，建议保持 `n_replicates` 固定不变，只修改
`simulation.prod_steps` 分多次运行比较，而不是让不同 replica 混用不同长度。

最小 length scan demo：

```yaml
replica:
  n_replicates: 3

length_scan:
  enabled: true
  prod_steps:
    - 1000000   # 1.0 ns @ 1.0 fs
    - 2000000   # 2.0 ns @ 1.0 fs
    - 4000000   # 4.0 ns @ 1.0 fs
```

这会在 `case.workdir` 下生成：

- `len_1000000/`
- `len_2000000/`
- `len_4000000/`

每个长度点内部再按 `replica_01/02/03` 组织；顶层额外写出
`viscosity_length_scan_summary.json`，把每个长度点的 `eta_mean/std/sem` 汇总到一起。

编辑 `configs/viscosity.yaml`，至少修改：

```yaml
case:
  name:    "LiPF6_EC_DMC"
  workdir: "./cases/lipf6_ec_dmc_viscosity"

components:
  - { name: EC,  kind: solvent, smiles: "C1COC(=O)O1", charge_source: openff_am1bcc }
  - { name: DMC, kind: solvent, smiles: "COC(=O)OC",   charge_source: openff_am1bcc }
  - { name: Li,  kind: cation,  preset: "Li+" }
  - { name: PF6, kind: anion,   preset: "PF6-" }

composition:
  mode: counts
  counts:
    EC:  182
    DMC: 182
    Li:  14
    PF6: 14

bamboo:
  model_file: "/path/to/bamboo_model.pt"   # 必填
```

### 四个可独立跳过的阶段

```
build_system → write_input → run_lammps → analyze
```

在配置文件的 `run:` 节中将对应项设为 `false` 即可跳过：

```yaml
run:
  build_system: false   # 已有 in.data + build_report.json 时跳过
  write_input:  true
  run_lammps:   true
  analyze:      true
```

如果你只想重跑平衡段或只想重跑 GK 生产段，可以继续细分：

```yaml
run:
  run_lammps: true
  run_equil:  false   # 复用已有 equil_nvt.restart
  run_gk:     true
```

同一 `workdir` 下重复执行时，粘度入口会自动判断哪些产物可以直接复用：

- `build_system`：已有 `in.data` + `build_report.json`
- `write_input`：已有 `in.equil.lammps` 和 `in.gk.lammps`
- `run_lammps`：已有 `equil_nvt.restart` + `npt_thermo.dat`，或已有 `stress_acf.dat` + `pressure_tensor.dat` + `gk_thermo.dat`
- `analyze`：已有 `viscosity_summary.json` + `viscosity_analysis.png`；`pressure_blocks` 方法还会复用 `viscosity_running.csv` 和 `viscosity_blocks.csv`

如果上游输入更新了，下游阶段会自动失效并重新执行。每次运行还会在 `workdir`
写出 `state.json` 和 `stage.done`，记录各阶段的 `done / reused / skipped / failed`
状态，方便续跑和排错。

### 核心文件（粘度）

- `md_viscosity.py` — 主流程入口
- `core/composition.py` — 多组分 counts / molality / molarity 解析
- `core/bamboo_structure.py` — SMILES → RDKit 3D 结构，或 ion preset 几何
- `core/bamboo_packing.py` — 多组分 Packmol 建盒
- `core/data_builder.py` — BAMBOO `atom_style full` data writer
- `workflow/bamboo_system.py` — YAML components/composition → `in.data` + `build_report.json`
- `workflow/input_viscosity.py` — BAMBOO NPT 平衡 + NVE Green-Kubo 输入生成
- `analysis/viscosity.py` — Green-Kubo η 后处理（ACF 兼容路径 + pressure-block SEM）
- `configs/viscosity.yaml` — 官方默认长程配置模板
- `configs/viscosity.smoke.yaml` — 短程 smoke test 模板

### 技术说明

- BAMBOO 使用 `atom_style full`，`Atoms` 行为 `id mol_id type charge x y z`
- `core/data_builder.py` 会按 atomic number 给元素分配 type，并把同一顺序写进 `pair_coeff`
- `pair_style bamboo` 的参数默认来自 BAMBOO 示例：`[5.0, 5.0, 10.0, 1]`
- 默认启用 PPPM，体系总电荷会在写 `in.data` 前检查为近似中性
- 默认通过 Kokkos 后缀运行：`lmp_mpi -k on g 1 -sf kk`
- 默认运行环境等价于先执行：`export PATH="/root/bamboo/pair/lammps/output:$PATH"`
- 平衡阶段会额外写出 `npt_thermo.dat`，并默认检查 NPT 尾段平均密度和最终 NVT 尾段平均温度
- GK 生产阶段默认额外写出 `pressure_tensor.dat`，包含 `step pxy pxz pyz`
- 长程配置默认使用 `analysis.method: pressure_blocks`：从压力张量时间序列分 block 计算 ACF、running η 和 block SEM，并用固定长度窗口选择 plateau
- smoke 配置默认保留 `analysis.method: acf`：直接分析 `stress_acf.dat`，适合快速验证流程，不应作为 production 粘度
- 单位换算：`ETA_CONV = atm² × Å³ × fs / k_B × 10³ ≈ 7.44×10⁻¹⁰ mPa·s·K/Å³ per atm²·fs`

---

## 公共基础设施

两条工作流共享以下模块：

| 模块 | 作用 |
|------|------|
| `core/structure.py` | SMILES → RDKit 3D 结构 |
| `core/packing.py` | Packmol 多分子建盒 |
| `core/bamboo_structure.py` / `core/bamboo_packing.py` | BAMBOO 粘度工作流的多组分建系 |
| `analysis/hfacf.py` | `parse_ave_correlate_detail`、`_integrate_acf`、收敛窗口检测 |
| `utils/constants.py` | 物理常数与单位换算 |
| `utils/io.py` | LAMMPS data 文件读写 |
| `workflow/config.py` | YAML 配置解析 |

## CLI 参数约定

- 统一使用 `kebab-case` 参数名，如 `--temperature-k`、`--timestep-fs`
- 需要显式单位时，直接写进参数名：`-k`、`-fs`、`-a`、`-a3`、`-g-cm3`

## 其他 CLI 示例

构建 `system.data`（热导率工作流）：

```bash
python scripts/build_system.py \
  --smiles CCO \
  --name ethanol \
  --workdir ./demo_ethanol \
  --n-molecules 200 \
  --density-g-cm3 0.789
```

HFACF 热导率后处理：

```bash
python scripts/analyze_thermal_hfacf.py \
  --hfacf-file ./hfacf.dat \
  --gk-data-file ./gk_data.dat \
  --temperature-k 298.15
```

启用 LAMMPS GPU 包（热导率工作流）时，将配置中的 `lammps.use_gpu` 改为 `true`
并设置 `lammps.gpu_count`。粘度工作流使用 BAMBOO/Kokkos，相关参数在 `lammps.kokkos`、
`lammps.kokkos_gpus` 和 `lammps.suffix` 中配置。
