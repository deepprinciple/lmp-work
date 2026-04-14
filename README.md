# Thermal-Only Workflow

该目录已精简为热导率单线流程，核心能力如下：

1. 读取 `SMILES` 并生成 3D 分子结构
2. 力场参数化（默认 `OpenFF + AM1BCC`，可选 `LigParGen + BOSS`）
3. 使用 `Packmol` 建立液体体系
4. 生成并运行 LAMMPS reverse-NEMD 热导率模拟
5. 对热导率进行后处理拟合与可视化

## 主入口

```bash
python scripts/run_thermal_rnemd.py --config configs/config.yaml
```

短程 smoke test：

```bash
python scripts/run_thermal_rnemd.py --config configs/config.smoke.yaml
```

## CLI 参数约定

- 统一使用 `kebab-case` 参数名，如 `--temperature-k`、`--timestep-fs`
- 需要显式单位时，直接写进参数名：`-k`、`-fs`、`-a`、`-a3`、`-g-cm3`
- CLI 只保留这一套参数命名，不再兼容旧别名

## 核心文件

- `scripts/build_system.py`: `SMILES -> system.data`
- `scripts/run_thermal_rnemd.py`: 热导率主流程入口
- `scripts/write_lammps_input.py`: 手动生成平衡段/GK 输入的 CLI
- `analysis/hfacf.py`: HFACF 热导率分析与后处理核心
- `workflow/input_thermal.py`: 平衡段和 rNEMD 输入生成
- `analysis/rnemd.py`: rNEMD 热导率拟合与收敛分析
- `scripts/analyze_thermal_hfacf.py`: HFACF 热导率后处理 CLI

## 目录说明

- `scripts/`: 当前建议使用的 CLI 入口
- `scripts/env/`: 环境兼容 shell 脚本
- `configs/`: 当前主线配置文件
- `workflow/`: 热导率主线所需的运行与输入生成模块
- `legacy/`: 仅保留手动 GK 输入模板

## 力场选择

- 默认：`forcefield.engine: openff`，`forcefield.charge_method: am1bcc`
- 可选：`forcefield.engine: ligpargen`（依赖 LigParGen/BOSS 运行环境）

## 授权与分发声明

- 本仓库不再包含 `BOSS` 程序及其兼容运行时副本。
- `BOSS` 为受限授权软件，使用者需自行向权利方申请授权并在本地安装。
- `LigParGen + BOSS` 路径仅保留接口能力；实际运行需通过本地 `wrapper_script`
  指向你自己的合法安装环境。

## HFACF 后处理示例

```bash
python scripts/analyze_thermal_hfacf.py \
  --hfacf-file ./hfacf.dat \
  --gk-data-file ./gk_data.dat \
  --temperature-k 298.15
```

输出：

- `thermal_hfacf_summary.json`
- `thermal_hfacf_analysis.png`

## 其他 CLI 示例

构建 `system.data`：

```bash
python scripts/build_system.py \
  --smiles CCO \
  --name ethanol \
  --workdir ./demo_ethanol \
  --n-molecules 200 \
  --density-g-cm3 0.789
```

生成手动平衡/GK 输入：

```bash
python scripts/write_lammps_input.py \
  --workdir ./demo_case \
  --mode replica \
  --temperature-k 298.15 \
  --decorrelation-steps 500000 \
  --gk-steps 10000000
```
