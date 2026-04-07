# Thermal-only scope

这个目录是从更大的多性质项目复制出来的热导率专用副本。

当前目标只保留：

- 建模与参数化
- 长盒子液体体系生成
- `fix thermal/conductivity` reverse NEMD 生产
- 热导率数值分析
- `HFACF` 后处理与可视化

不再继续维护：

- 粘度 workflow
- 介电常数 workflow
- 与方盒子粘度验证相关的批量脚本

## 主入口

- `run_thermal_rnemd.py`
- `analyze_thermal_hfacf.py`
- `diagnostics/live_rnemd_monitor.py`

## 基本原则

- 后续功能开发优先围绕热导率一条主线
- 如果需要新方法探索，优先作为热导率专用后处理脚本加入
- 不再把“多种性质共用一条生产流程”作为设计目标
