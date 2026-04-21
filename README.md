# lmp-work — LAMMPS 输运性质计算工作流

基于 LAMMPS 的液体输运性质计算流水线，目前支持两条工作流：

| 工作流 | 方法 | 力场 | 入口 |
|--------|------|------|------|
| **热导率** | reverse-NEMD | OpenFF / LigParGen | `md_run.py` |
| **粘度** | Green-Kubo 应力 ACF | ANI 神经网络势 | `md_viscosity.py` |

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
- `lammps-ani` + ANI 模型文件：粘度工作流需要
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

## 粘度工作流（ANI Green-Kubo）

SMILES → Packmol 建盒 → ANI pair style → NVT 预热 → NPT 压回目标密度 → 短 NVT 收尾 → NVE + 应力 ACF → η

> **依赖**：需要编译安装 [lammps-ani](https://github.com/roitberg-group/lammps-ani)
> 并准备 ANI TorchScript 模型文件（`ani2x.pt`）。

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
    - 1000000   # 0.5 ns @ 0.5 fs
    - 2000000   # 1.0 ns @ 0.5 fs
    - 4000000   # 2.0 ns @ 0.5 fs
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
  smiles:  "CCO"          # 目标分子 SMILES
  name:    "ethanol"
  workdir: "./cases/ethanol_viscosity_300K"

ani:
  model_file: "/path/to/ani2x.pt"   # ← 必填
  device:     "cuda"
```

### 四个可独立跳过的阶段

```
build_system → write_input → run_lammps → analyze
```

在配置文件的 `run:` 节中将对应项设为 `false` 即可跳过：

```yaml
run:
  build_system: false   # 已有 system.data 时跳过
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

- `build_system`：已有 `system.data`
- `write_input`：已有 `in.equil.lammps` 和 `in.gk.lammps`
- `run_lammps`：已有 `equil_nvt.restart` + `npt_thermo.dat`，或已有 `stress_acf.dat` + `gk_thermo.dat`
- `analyze`：已有 `viscosity_summary.json` + `viscosity_analysis.png`

如果上游输入更新了，下游阶段会自动失效并重新执行。每次运行还会在 `workdir`
写出 `state.json` 和 `stage.done`，记录各阶段的 `done / reused / skipped / failed`
状态，方便续跑和排错。

### 核心文件（粘度）

- `md_viscosity.py` — 主流程入口
- `core/data_builder.py` — `AniDataBuilder`：Packmol XYZ → LAMMPS atomic data
- `workflow/input_ani_viscosity.py` — NVT 预热 + NPT 密度平衡 + 短 NVT 收尾 + NVE GK 输入生成
- `analysis/viscosity.py` — Green-Kubo η 积分（复用 `hfacf` 引擎）
- `configs/viscosity.yaml` — 官方默认长程配置模板
- `configs/viscosity.smoke.yaml` — 短程 smoke test 模板

### 技术说明

- ANI 通过原子质量识别元素，`pair_coeff` 只需 `* *`，无需写元素符号
- 必须使用 `pyaev full`（CUAEV 不支持 virial/stress 计算）
- 非 Kokkos 模式要求 `newton off`
- 平衡阶段会额外写出 `npt_thermo.dat`，并默认检查 NPT 尾段平均密度和最终 NVT 尾段平均温度
- 粘度入口会先尝试补齐保守的 ANI 运行环境；如需完全手动控制，可设 `lammps.auto_environment: false`
- 单位换算：`ETA_CONV = atm² × Å³ × fs / k_B × 10³ ≈ 7.44×10⁻¹⁰ mPa·s·K/Å³ per atm²·fs`

---

## 公共基础设施

两条工作流共享以下模块：

| 模块 | 作用 |
|------|------|
| `core/structure.py` | SMILES → RDKit 3D 结构 |
| `core/packing.py` | Packmol 多分子建盒 |
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
并设置 `lammps.gpu_count`。粘度工作流的 GPU 由 ANI/PyTorch 自行管理，无需此项。
