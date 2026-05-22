# lmp-work: BAMBOO viscosity GPU demo

这个分支是给 GPU 供应商验证 BAMBOO/LAMMPS 运行环境用的最小 demo。当前只保留
`LiPF6 + EC/DMC` 电解液的 Green-Kubo 剪切粘度流程：

```text
build_system -> write_input -> run_lammps -> analyze
```

工作流会生成 BAMBOO `atom_style full` 的 `in.data`、写出 LAMMPS 输入文件，
运行平衡和 NVT 生产段，并从压力张量后处理得到粘度。

## 1. 环境

建议用 conda-forge 安装 Python 依赖和建模工具：

```bash
conda create -n lmp-work -c conda-forge python=3.11 \
  numpy pyyaml "matplotlib>=3.9" rdkit openff-toolkit openff-units \
  ambertools packmol -y
conda activate lmp-work
```

还需要安装 ByteDance BAMBOO，并使用 BAMBOO 自己编译的 LAMMPS。安装方法请参考
[bytedance/bamboo](https://github.com/bytedance/bamboo)，模型请使用
`paper_new_disp.pt`。

开发机默认路径是：

```text
/root/bamboo/benchmark/paper_new_disp.pt
/root/bamboo/pair/lammps/output
```

供应商机器路径通常不同，运行前必须在 `configs/viscosity.yaml` 里改成实际位置：

```yaml
bamboo:
  model_file: "/path/to/paper_new_disp.pt"

lammps:
  executable: "lmp_mpi"
  env:
    PATH: "/path/to/bamboo/pair/lammps/output:${PATH}"
```

等价手动设置是：

```bash
export PATH="/path/to/bamboo/pair/lammps/output:$PATH"
```

## 2. 运行

主 demo：

```bash
python md_viscosity.py --config configs/viscosity.yaml
```

短程 smoke test：

```bash
python md_viscosity.py --config configs/viscosity.smoke.yaml
```

默认主配置只跑一个生产长度，不做多副本或 length scan。要改变 GPU 验证任务时长，
直接修改：

```yaml
simulation:
  timestep_fs: 1.0
  prod_steps: 2000000   # 2 ns at 1 fs
```

生产段默认参考 BAMBOO 示例，使用 NVT：

```text
fix nvt temp 300 300 10
fix ave/time 1 1 1 v_pxy v_pxz v_pyz file dump_pressure.out
```

## 3. 只重跑后处理

如果 LAMMPS 已经跑完，只想重新分析结果，把配置中的 `run` 改成：

```yaml
run:
  build_system: false
  write_input:  false
  run_lammps:   false
  run_equil:    false
  run_gk:       false
  analyze:      true
```

然后重新执行同一条命令即可：

```bash
python md_viscosity.py --config configs/viscosity.yaml
```

## 4. 关键输出

默认工作目录是 `cases/lipf6_ec_dmc_viscosity/`。主要文件：

- `in.data`: BAMBOO data file。
- `in.equil.lammps`: minimize + NVT + NPT + short NVT 平衡输入。
- `equil_nvt.restart`: 生产段初始 restart。
- `in.gk.lammps`: NVT pressure-trace 生产输入。
- `run_gk.log`: LAMMPS 生产段日志。
- `dump_pressure.out`: `step pxy pxz pyz` 原始压力张量。
- `gk_thermo.dat`: 生产段温度、压力和体积。
- `viscosity_summary.json`: 粘度结果和窗口诊断。
- `viscosity_running.csv`: running viscosity 曲线数据。
- `viscosity_blocks.csv`: block plateau 数据。
- `viscosity_analysis.png`: ACF 和 running viscosity 图。

主配置使用 `analysis.method: pressure_blocks`，从 `dump_pressure.out` 切 block
计算 ACF、running viscosity 和 block SEM。`viscosity_summary.json` 中
`window_quality: ok` 表示自动选择的积分窗口通过了当前稳定性检查。

## 5. 常用修改项

通常只需要改这些字段：

```yaml
case:
  workdir: "./cases/lipf6_ec_dmc_viscosity"

composition:
  counts:
    EC:  182
    DMC: 182
    Li:  14
    PF6: 14

bamboo:
  model_file: "/path/to/paper_new_disp.pt"

lammps:
  mpi_ranks: 1
  kokkos_gpus: 1
  env:
    PATH: "/path/to/bamboo/pair/lammps/output:${PATH}"
```

溶剂结构和 AM1-BCC 电荷由 RDKit/OpenFF 生成；`Li+` 和 `PF6-` 使用仓库内置
preset。BAMBOO 负责相互作用，生成的 `in.data` 不包含 bonds、angles 或
dihedrals。

## 6. 代码入口

- `md_viscosity.py`: 主入口。
- `configs/viscosity.yaml`: 默认 GPU 验证配置。
- `configs/viscosity.smoke.yaml`: 快速 smoke test 配置。
- `workflow/input_viscosity.py`: LAMMPS 输入生成。
- `workflow/bamboo_system.py`: 体系和 `in.data` 生成。
- `analysis/viscosity.py`: 粘度后处理。
