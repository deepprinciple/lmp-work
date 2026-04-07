# Thermal rNEMD Workflow

这个副本工程只保留一条主线：

- 从 `SMILES` 建模与参数化
- 生成长盒子液体体系
- 使用 LAMMPS `fix thermal/conductivity` 进行 `reverse NEMD`
- 对热导率结果做数值分析与可视化后处理

它不再以“一个项目同时追求粘度、介电常数、热导率都精准”为目标。

副本范围说明见 `docs/README_thermal_only.md`。

## 主入口

```bash
python run_thermal_rnemd.py --input config.yaml
```

短程 smoke test:

```bash
python run_thermal_rnemd.py --input config.smoke.yaml
```

## 工作流

1. `build_system.py` 从 `SMILES` 构建液体体系
2. `run_thermal_rnemd.py` 写共享平衡输入与 rNEMD replica 输入
3. 平衡段输出 `equil_nvt.restart`
4. 每个 replica 运行 `fix thermal/conductivity`
5. `analysis/rnemd.py` 基于温度剖面和交换热量计算热导率
6. `analyze_thermal_hfacf.py` 可选地对 `HFACF` 做 Green-Kubo 后处理和可视化

## 主要文件

- `run_thermal_rnemd.py`：热导率主入口
- `config.yaml`：正式热导率配置模板
- `config.smoke.yaml`：小体系 CPU smoke 配置
- `build_system.py`：从 `SMILES` 到 `system.data`
- `write_lammps_input.py`：平衡段与 rNEMD 输入生成
- `analysis/rnemd.py`：rNEMD 热导率分析
- `diagnostics/live_rnemd_monitor.py`：实时监控
- `analyze_thermal_hfacf.py`：HFACF 与热导率后处理、作图

## 结果目录

默认热导率 replica 位于：

```text
case.workdir/thermal_conductivity/method_rnemd/replica_01/
```

每个 replica 的关键输出包括：

- `thermal_exchange.dat`
- `temp_profile.dat`
- `thermal_rnemd_thermo.dat`
- `thermal_rnemd_summary.json`
- `thermal_profile_latest.png`
- `thermal_kappa_running.png`

## HFACF 后处理

如果你有 `hfacf.dat` 和 `gk_data.dat`，可以额外做 Green-Kubo 热导率后处理：

```bash
python analyze_thermal_hfacf.py \
  --hfacf ./hfacf.dat \
  --gk-data ./gk_data.dat \
  --temperature-k 298.15
```

它会输出：

- `thermal_hfacf_summary.json`
- `thermal_hfacf_analysis.png`

图中包含：

- 归一化 `HFACF`
- 热导率 running integral / 数值收敛曲线

## 力场参数化

当前仍支持：

- `OpenFF`
- `LigParGen/BOSS`

对应配置位于 `forcefield:` 段，例如：

- `engine: openff`
- `engine: ligpargen`

## 备注

- 这个目录是从更大的项目复制出来的 thermal-only 副本。
- 后续所有热导率专用重构都应在这个目录内继续，而不是回到原始项目里做。
