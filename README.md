# Thermal-Only Workflow

该目录已精简为热导率单线流程，核心能力如下：

1. 读取 `SMILES` 并生成 3D 分子结构
2. 力场参数化（默认 `OpenFF + AM1BCC`，可选 `LigParGen + BOSS`）
3. 使用 `Packmol` 建立液体体系
4. 生成并运行 LAMMPS reverse-NEMD 热导率模拟
5. 对热导率进行后处理拟合与可视化

## 环境准备

完整环境安装请优先参考 `INSTALL.md`。

- `INSTALL.md` 覆盖了 `LAMMPS`、`Packmol`、`OpenFF`、`AmberTools`、MPI/GPU 相关依赖的安装方式
- 仓库根目录的 `requirements.txt` 只保留了补充安装的 Python 包，不负责完整环境引导

环境准备完成后，如果还需要补装仓库中的 Python 依赖，再执行：

```bash
pip install -r requirements.txt
```


## 主入口

```bash
python md_run.py --config configs/config.yaml
```

smoke test：

```bash
python md_run.py --config configs/config.smoke.yaml
```

## 结果在哪里看

主流程的所有中间文件和结果都会写到配置里的 `case.workdir` 目录。

- 默认主配置：`configs/config.yaml` -> `./cases/thf_thermal_only`
- smoke test 配置：`configs/config.smoke.yaml` -> `./cases/methanol_thermal_smoke`

通常跑完后，先看 `case.workdir` 根目录下的总汇总文件，再看各个 replica 子目录里的单次结果：

```text
case.workdir/
├─ build_system_report.json
├─ system.data
├─ in.equil.lammps
├─ run_equil.log
├─ thermal_branch_manifest.json
├─ thermal_replica1.log
├─ thermal_replica2.log
├─ ...
├─ thermal_results.json
└─ thermal_conductivity/method_rnemd/
   ├─ replica_01/
   │  ├─ in.thermal.lammps
   │  ├─ temp_profile.dat
   │  ├─ thermal_exchange.dat
   │  ├─ thermal_rnemd_thermo.dat
   │  ├─ thermal_summary.json
   │  ├─ thermal_profile.png
   │  └─ thermal_kappa.png
   └─ replica_02/
      └─ ...
```

重点看这几个文件：

- `thermal_results.json`: 整个案例的总汇总，包含所有 replica 的热导率结果和均值
- `thermal_conductivity/method_rnemd/replica_XX/thermal_summary.json`: 单个 replica 的详细分析结果
- `thermal_conductivity/method_rnemd/replica_XX/thermal_profile.png`: 最新温度剖面拟合图
- `thermal_conductivity/method_rnemd/replica_XX/thermal_kappa.png`: 热导率随时间的收敛图
- `run_equil.log` 和 `thermal_replica*.log`: LAMMPS 运行日志，排查报错时先看这里

如果你改了 `thermal_rnemd.summary_file`、`batch_summary_file`、`profile_plot_file`、`kappa_plot_file`
或 `replica_dir_prefix`，实际文件名和目录会随配置一起变化。

## CLI 参数约定

- 统一使用 `kebab-case` 参数名，如 `--temperature-k`、`--timestep-fs`
- 需要显式单位时，直接写进参数名：`-k`、`-fs`、`-a`、`-a3`、`-g-cm3`
- CLI 只保留这一套参数命名，不再兼容旧别名

## 核心文件

- `md_run.py`: 热导率主流程入口
- `scripts/build_system.py`: `SMILES -> system.data`
- `scripts/write_lammps_input.py`: 手动生成平衡段/GK 输入的 CLI
- `analysis/hfacf.py`: HFACF 热导率分析与后处理核心
- `workflow/input_thermal.py`: 平衡段和 rNEMD 输入生成
- `analysis/rnemd.py`: rNEMD 热导率拟合与收敛分析
- `scripts/analyze_thermal_hfacf.py`: HFACF 热导率后处理 CLI

## 目录说明

- `scripts/`: 辅助 CLI 与环境脚本
- `scripts/env/`: 环境兼容 shell 脚本
- `configs/config.yaml`: 当前主线默认配置，按 CPU 安全默认值提供
- `configs/config.smoke.yaml`: 短程 smoke test 配置
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

启用 GPU 时，将 `configs/config.yaml` 或你自己的配置中的 `lammps.use_gpu`
改为 `true`，并把 `lammps.gpu_count` 设为实际可用卡数。
