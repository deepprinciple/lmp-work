# lmp-work — LAMMPS 输运性质计算工作流

基于 LAMMPS 的液体输运性质计算流水线，目前支持两条工作流：

| 工作流 | 方法 | 力场 | 入口 |
|--------|------|------|------|
| **热导率** | reverse-NEMD | OpenFF / LigParGen | `md_run.py` |
| **粘度** | Green-Kubo 应力 ACF | ANI 神经网络势 | `md_viscosity.py` |

---

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

SMILES → Packmol 建盒 → ANI pair style → NVT 平衡 → NVE + 应力 ACF → η

> **依赖**：需要编译安装 [lammps-ani](https://github.com/roitberg-group/lammps-ani)
> 并准备 ANI TorchScript 模型文件（`ani2x.pt`）。

```bash
python md_viscosity.py --config configs/viscosity.yaml
```

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

### 核心文件（粘度）

- `md_viscosity.py` — 主流程入口
- `core/data_builder.py` — `AniDataBuilder`：Packmol XYZ → LAMMPS atomic data
- `workflow/input_ani_viscosity.py` — NVT 平衡段和 NVE GK 产出段输入生成
- `analysis/viscosity.py` — Green-Kubo η 积分（复用 `hfacf` 引擎）
- `configs/viscosity.yaml` — 粘度配置模板（水，300 K 示例）

### 技术说明

- ANI 通过原子质量识别元素，`pair_coeff` 只需 `* *`，无需写元素符号
- 必须使用 `pyaev full`（CUAEV 不支持 virial/stress 计算）
- 非 Kokkos 模式要求 `newton off`
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
