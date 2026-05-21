# lmp-work — BAMBOO 粘度 GPU 测试 demo

这个分支只保留一个用途：用 ByteDance BAMBOO 跑 Green-Kubo 剪切粘度，用来给 GPU
供应商做 MD 任务跑通性、性能和结果一致性验证。旧的其他输运性质 workflow
和生产向参数化路径已经移除。

默认 case 是 `LiPF6 + EC/DMC`。工作流会生成 BAMBOO `atom_style full`
的 `in.data`、LAMMPS 输入文件，运行 NPT/NVT，并对压力张量做粘度后处理。

## 环境搭建

建议先创建独立 Python 环境：

```bash
conda create -n lmp-work python=3.11 -y
conda activate lmp-work
pip install -r requirements.txt
```

还需要外部程序：

- `packmol`：用于生成初始盒子，命令需要在 `PATH` 中。
- `ambertools`：默认溶剂电荷为 OpenFF AM1-BCC，OpenFF toolkit 通常需要
  AmberTools 的 `sqm`。缺失时可用 `conda install -c conda-forge ambertools -y`。
- BAMBOO-enabled LAMMPS：必须使用 BAMBOO 自己编译过的 LAMMPS。
- BAMBOO 模型文件：请按 [bytedance/bamboo](https://github.com/bytedance/bamboo)
  官方说明安装 BAMBOO，并选用 `paper_new_disp.pt` 模型。

开发机上当前模型路径是：

```text
/root/bamboo/benchmark/paper_new_disp.pt
```

这只是本机路径。GPU 供应商机器上的安装目录很可能不同，运行前需要在
`configs/viscosity.yaml` 里改成目标机器上的实际路径：

```yaml
bamboo:
  model_file: "/path/to/bamboo/benchmark/paper_new_disp.pt"
```

同理，BAMBOO LAMMPS 的可执行程序目录也要指向目标机器的安装位置。开发机等价于先执行：

```bash
export PATH="/root/bamboo/pair/lammps/output:$PATH"
```

对应 YAML 配置是：

```yaml
lammps:
  executable: "lmp_mpi"
  env:
    PATH: "/root/bamboo/pair/lammps/output:${PATH}"
```

生成的 LAMMPS 输入中，`pair_coeff` 会由 `bamboo.model_file` 和体系元素顺序一起渲染，例如：

```text
pair_coeff      /root/bamboo/benchmark/paper_new_disp.pt   H LI C O F P
```

如果模型文件路径不对，LAMMPS 会在启动 BAMBOO pair style 时失败。

## 运行

主配置：

```bash
python md_viscosity.py --config configs/viscosity.yaml
```

短程 smoke test：

```bash
python md_viscosity.py --config configs/viscosity.smoke.yaml
```

默认只运行一个 case，不会自动展开多副本或多段模拟时长。需要调整 GPU 验证任务长度时，
直接改 `simulation.prod_steps`：

```yaml
simulation:
  timestep_fs: 1.0
  prod_ensemble: "nvt"
  prod_temp_damp_fs: 10.0
  prod_steps: 2000000   # 2.0 ns @ 1.0 fs
```

默认生产段参考 BAMBOO 示例输入，使用 `fix nvt temp 300 300 10`，并用
`fix ave/time 1 1 1` 输出压力张量。`prod_ensemble` 可以改成 `"nve"`，但 NVE
下必须额外检查能量漂移；1 fs 在 NVT 中能跑稳，不代表 NVE 中也一定守恒良好。

## 配置重点

`configs/viscosity.yaml` 里通常只需要改这些部分：

```yaml
case:
  name: "LiPF6_EC_DMC"
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
  model_file: "/path/to/paper_new_disp.pt"

lammps:
  env:
    PATH: "/path/to/bamboo/pair/lammps/output:${PATH}"
```

溶剂电荷默认仍用 OpenFF toolkit 的 AM1-BCC 生成；离子使用仓库内置 preset。
BAMBOO 本身负责短程相互作用，`in.data` 中不会写 Bonds/Angles/Dihedrals。

## 阶段控制

工作流阶段：

```text
build_system -> write_input -> run_lammps -> analyze
```

在配置文件的 `run:` 节中将对应项设为 `false` 可以跳过：

```yaml
run:
  build_system: false   # 已有 in.data + build_report.json 时跳过
  write_input:  true
  run_lammps:   true
  run_equil:    false   # 复用已有 equil_nvt.restart
  run_gk:       true
  analyze:      true
```

同一 `workdir` 下重复执行时，入口会复用仍然有效的产物；上游输入更新后，下游阶段会重新执行。
`state.json` 和 `stage.done` 会记录每个阶段的 `done / reused / skipped / failed` 状态。

## 输出文件

主要输入和运行产物：

- `in.data`：BAMBOO `atom_style full` data file。
- `in.equil.lammps`：minimize + NVT + NPT + short NVT 平衡输入。
- `equil_nvt.restart`：生产段初始 restart。
- `in.gk.lammps`：NVT Green-Kubo pressure-trace 生产输入。
- `npt_thermo.dat`：平衡段密度/温度 trace。
- `gk_thermo.dat`：生产段温度/压力/体积 trace。
- `dump_pressure.out`：`step pxy pxz pyz` 原始压力张量。
- `stress_acf.dat`：LAMMPS `fix ave/correlate` 的应力 ACF。
- `nvt.data`：生产段结束后的 data 文件。

后处理产物：

- `viscosity_summary.json`：最终粘度、窗口选择和诊断信息。
- `viscosity_running.csv`：running viscosity 曲线。
- `viscosity_blocks.csv`：pressure-block plateau 数据。
- `viscosity_analysis.png`：ACF 与 running η 图。

长程配置默认使用 `analysis.method: pressure_blocks`，会从 `dump_pressure.out`
切 block 做 ACF、running η 和 block SEM。smoke 配置默认使用 `analysis.method: acf`，
只用于快速验证流程，不建议把 smoke 结果当成正式粘度。

## 核心文件

- `md_viscosity.py`：主流程入口。
- `configs/viscosity.yaml`：默认 GPU 验证配置。
- `configs/viscosity.smoke.yaml`：短程 smoke test 配置。
- `workflow/bamboo_system.py`：components/composition 到 `in.data`。
- `workflow/input_viscosity.py`：BAMBOO LAMMPS 输入生成。
- `core/bamboo_structure.py`：SMILES/RDKit 结构和 ion preset。
- `core/bamboo_packing.py`：Packmol 建盒。
- `core/data_builder.py`：BAMBOO data writer。
- `analysis/viscosity.py`：粘度后处理。
- `analysis/correlation.py`：LAMMPS ACF 解析和 running integral 公共逻辑。
