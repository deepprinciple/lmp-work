# 热导率 workflow 整理说明

本文档只整理当前项目中与热导率 `rNEMD` 管线直接相关的入口、配置和输出约定，
不改变项目整体架构，也不引入新的独立子系统。

## 当前定位

当前项目的基线策略仍然是：

- 粘度：`EMD + GK`
- 热导率：长盒子 `reverse NEMD`，即 LAMMPS `fix thermal/conductivity`

热导率 workflow 保持在现有项目结构内，重点是让入口、配置和输出布局更清楚。

## 主要入口

- `run_thermal_rnemd.py`

它负责：

1. 读取 YAML 配置
2. 构建长盒子液体体系
3. 写共享平衡输入
4. 写热导率 replica 输入
5. 运行 `fix thermal/conductivity`
6. 汇总每个 replica 的热导率分析结果

## 关键模块

- `write_lammps_input.py`
  - `build_rnemd_replica_input_text()`
  - 负责生成热导率生产段的 LAMMPS 输入

- `analysis/rnemd.py`
  - `analyze_rnemd_replica()`
  - 负责读取 `thermal_exchange.dat` 和 `temp_profile.dat`
  - 输出 `thermal_rnemd_summary.json` 和收敛图

- `diagnostics/live_rnemd_monitor.py`
  - 负责对热导率生产段做实时监控

- `bench_runs/generate_5mol_baseline_configs.py`
  - 用于批量生成包含 `thermal_rnemd` 配置段的验证 YAML

## 当前目录约定

热导率分支在单个 case 下的默认布局为：

```text
case.workdir/
├── in.equil.lammps
├── equil_nvt.restart
├── thermal_branch_manifest.json
└── thermal_conductivity/
    └── method_rnemd/
        ├── replica_01/
        ├── replica_02/
        └── ...
```

每个 replica 目录下的关键文件包括：

- `in.thermal_rnemd.lammps`
- `thermal_exchange.dat`
- `temp_profile.dat`
- `thermal_rnemd_thermo.dat`
- `thermal_rnemd_summary.json`
- `thermal_profile_latest.png`
- `thermal_kappa_running.png`

## 配置重点

热导率 workflow 的核心配置位于 YAML 的 `thermal_rnemd:` 段，例如：

- `box_aspect_ratio`
- `replica_dir_prefix`
- `replica_input_filename`
- `replica_log_prefix`
- `run_steps`
- `edim`
- `nbin`
- `swap_every_steps`
- `profile_file`
- `exchange_file`
- `summary_file`
- `batch_summary_file`

这些字段的默认布局现在统一由 `run_thermal_rnemd.py` 内部解析，避免路径和文件名默认值散落在多个分支判断里。

## 使用建议

- 如果目标是正式热导率生产，优先从 `run_thermal_rnemd.py` 进入。
- 如果目标是粘度和热导率共同验证，仍然按当前 split strategy 使用各自入口，不要把两种生产段硬塞进同一条运行脚本。
- 如果只是做实时监控或问题定位，优先查看每个 replica 下的：
  - `thermal_exchange.dat`
  - `temp_profile.dat`
  - `thermal_rnemd_summary.json`

## 当前整理范围

这次整理的目标只是：

- 让热导率主入口更易读
- 让默认目录/文件布局更集中
- 让文档入口更清晰

不包括：

- 把 thermal workflow 抽成独立架构
- 改写整体项目目录
- 替换现有基线策略
