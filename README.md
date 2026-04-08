# Thermal-Only Workflow

该目录已精简为热导率单线流程，核心能力如下：

1. 读取 `SMILES` 并生成 3D 分子结构
2. 力场参数化（默认 `OpenFF + AM1BCC`，可选 `LigParGen + BOSS`）
3. 使用 `Packmol` 建立液体体系
4. 生成并运行 LAMMPS reverse-NEMD 热导率模拟
5. 对热导率进行后处理拟合与可视化

## 主入口

```bash
python run_thermal_rnemd.py --input config.yaml
```

短程 smoke test：

```bash
python run_thermal_rnemd.py --input config.smoke.yaml
```

## 核心文件

- `build_system.py`: `SMILES -> system.data`
- `run_thermal_rnemd.py`: 热导率主流程入口
- `write_lammps_input.py`: 平衡段和 rNEMD 输入生成
- `analysis/rnemd.py`: rNEMD 热导率拟合与收敛分析
- `analyze_thermal_hfacf.py`: HFACF 后处理和图形输出

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
python analyze_thermal_hfacf.py \
  --hfacf ./hfacf.dat \
  --gk-data ./gk_data.dat \
  --temperature-k 298.15
```

输出：

- `thermal_hfacf_summary.json`
- `thermal_hfacf_analysis.png`
