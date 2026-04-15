# LAMMPS 环境搭建手册（GPU + OpenMPI + FFTW3 + KOKKOS）

目标能力：

- 使用 `mpirun` 并行（OpenMPI）
- 使用 GPU package 加速（`-sf gpu -pk gpu ...`）
- 使用 FFTW3（`kspace_style pppm` 依赖）
- 保留 `KOKKOS` 支持，便于后续切换到 KOKKOS 路径

---

## 1. 输入文件所需功能清单

当前 `in.lammps` 需要的关键样式：

- `pair_style lj/cut/coul/long`
- `kspace_style pppm`
- `bond_style harmonic`
- `angle_style harmonic`
- `dihedral_style fourier`
- `improper_style cvff`
- `fix npt / nvt / nve / ave/correlate`
- `compute heat/flux / stress/atom`

对应建议开启的 LAMMPS 包：

- `KSPACE`
- `MOLECULE`
- `EXTRA-MOLECULE`
- `EXTRA-FIX`
- `EXTRA-COMPUTE`
- `GPU`
- `OPENMP`（可选但推荐）

---

## 2. 系统依赖安装（Ubuntu/Debian）

如果你是普通用户，请加 `sudo`；如果你已是 root，可去掉 `sudo`。

```bash
apt-get update
apt-get install -y \
  build-essential cmake git pkg-config \
  openmpi-bin libopenmpi-dev \
  libfftw3-dev libfftw3-mpi-dev \
  libjpeg-dev libpng-dev zlib1g-dev
```

确认 CUDA 工具链可用（GPU package 需要）：

```bash
nvidia-smi
nvcc --version
```

如果 `nvcc` 不存在，请先安装 CUDA Toolkit（版本需与驱动兼容）。

---

## 2.1 安装 Packmol

优先推荐系统包安装：

```bash
apt-get install -y packmol
```

如果系统源没有，可用 conda-forge：

```bash
conda install -c conda-forge packmol
```

验收：

```bash
which packmol
packmol < /dev/null || true
```

---

## 2.2 安装 OpenFF（含 Python 环境）

当前工作流（`build_system.py` / `forcefields/openff.py`）需要：

- `rdkit`
- `openff-toolkit`
- `openff-interchange`

推荐单独 Conda 环境（避免污染系统 Python）：

```bash
conda create -n lmp python=3.11 -y
conda activate lmp
conda install -c conda-forge rdkit openff-toolkit openff-interchange -y
```

如果你要稳定使用 `am1bcc` 电荷，建议再装：

```bash
conda install -c conda-forge ambertools -y
```

> 说明：本仓库 `OpenFF` 已支持电荷回退（`gasteiger`），即使不装 AmberTools 通常也能继续。

验收：

```bash
python - <<'PY'
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
print("RDKit OK:", Chem.MolFromSmiles("CCO") is not None)
print("OpenFF OK:", Molecule, ForceField)
PY
```

---

## 3. 获取 LAMMPS 源码

```bash
cd ~
git clone https://github.com/lammps/lammps.git
cd lammps
```

可选：固定到某个稳定版本（避免开发版波动）：

```bash
git checkout stable
```

---

## 4. CMake 配置（启用 KOKKOS）

> 下面以安装到 `/opt/lammps` 为例。
> `GPU_ARCH` 要按你的显卡改：A30/A100 常见是 `sm_80`。

```bash
cd ~/lammps
cmake -S cmake -B build_gpu_mpi \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX=/opt/lammps \
  -D BUILD_MPI=ON \
  -D BUILD_OMP=ON \
  -D BUILD_SHARED_LIBS=ON \
  -D PKG_OPENMP=ON \
  -D PKG_GPU=ON \
  -D GPU_API=cuda \
  -D GPU_ARCH=sm_80 \
  -D GPU_PREC=mixed \
  -D PKG_KSPACE=ON \
  -D PKG_MOLECULE=ON \
  -D PKG_EXTRA-MOLECULE=ON \
  -D PKG_EXTRA-COMPUTE=ON \
  -D PKG_EXTRA-FIX=ON \
  -D FFT=FFTW3 \
  -D PKG_KOKKOS=ON
```

> 如果 CMake 找不到 FFTW3，可追加：
>
> `-D CMAKE_PREFIX_PATH=/usr`

---

## 5. 编译与安装

```bash
cmake --build ~/lammps/build_gpu_mpi -j"$(nproc)"
cmake --install ~/lammps/build_gpu_mpi
```

配置运行时库路径（避免 `liblammps.so` 找不到）：

```bash
echo "/opt/lammps/lib" > /etc/ld.so.conf.d/lammps.conf
ldconfig
```

---

## 6. 安装验收

### 6.1 查看编译能力

```bash
/opt/lammps/bin/lmp -h | egrep "MPI|GPU package|FFT library|KSPACE|MOLECULE|EXTRA-MOLECULE|EXTRA-COMPUTE|EXTRA-FIX|KOKKOS"
```

你应看到：

- `MPI ... Open MPI`
- `GPU package API: CUDA`
- `FFT library = FFTW3 ...`
- 包列表包含 `KSPACE/MOLECULE/EXTRA-*/KOKKOS`

并且包列表里应能看到 `KOKKOS`，这样后续需要切到 KOKKOS 方案时无需重新编译。

### 6.1.1 Packmol / OpenFF 联合验收

```bash
which packmol
python - <<'PY'
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
mol = Molecule.from_smiles("Cc1ccccc1")
ff = ForceField("openff-2.0.0.offxml")
print("OpenFF test OK:", mol.n_atoms, ff.__class__.__name__)
PY
```

### 6.2 运行你的输入文件

在 `in.lammps` 所在目录执行：

```bash
cd /Users/heisenberg/work/gk_workflow
export OMP_NUM_THREADS=1
mpirun --allow-run-as-root -np 8 \
  /opt/lammps/bin/lmp -sf gpu -pk gpu 1 \
  -in in.lammps -log run_in_lammps.log
```

说明：

- `-sf gpu`：自动给可支持样式加 `/gpu` 后缀（主要是 pair/kspace）
- `-pk gpu 1`：使用 1 张 GPU
- `OMP_NUM_THREADS=1`：避免 MPI x OMP 过度竞争（通常更稳）

---

## 7. 与当前流程匹配的建议

- 先用 `-np 6~10` 做基准，选最快的核数（不是越大越快）。
- 保持 `in.lammps` 中力场样式不变，不要手写 `bond_style .../gpu`（多数并不存在）。
- 如果你后续还要介电常数并在 LAMMPS 端用 `compute dipole`，可额外开启 `PKG_DIPOLE=ON`。

---

## 8. 常见问题排查

### 8.1 `Could not find CMAKE_ROOT`

通常是 Conda 里的 CMake 损坏。改用系统 CMake：

```bash
which cmake
/usr/bin/cmake --version
```

### 8.2 `liblammps.so.0: cannot open shared object file`

没做 `ldconfig` 或 `LD_LIBRARY_PATH` 未包含 `/opt/lammps/lib`。

### 8.3 `mpirun ... not enough slots`

- 降低 `-np`
- 或加 `--oversubscribe`（仅测试用）

### 8.4 root 运行 OpenMPI 被拦截

可用：

```bash
mpirun --allow-run-as-root ...
```

---

## 9. 推荐最小可用命令汇总

```bash
# configure
cmake -S cmake -B build_gpu_mpi \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX=/opt/lammps \
  -D BUILD_MPI=ON -D BUILD_OMP=ON -D BUILD_SHARED_LIBS=ON \
  -D PKG_OPENMP=ON -D PKG_GPU=ON -D GPU_API=cuda -D GPU_ARCH=sm_80 -D GPU_PREC=mixed \
  -D PKG_KSPACE=ON -D PKG_MOLECULE=ON -D PKG_EXTRA-MOLECULE=ON \
  -D PKG_EXTRA-COMPUTE=ON -D PKG_EXTRA-FIX=ON \
  -D FFT=FFTW3 -D PKG_KOKKOS=ON

# build + install
cmake --build build_gpu_mpi -j"$(nproc)"
cmake --install build_gpu_mpi

# run
export OMP_NUM_THREADS=1
mpirun --allow-run-as-root -np 8 \
  /opt/lammps/bin/lmp -sf gpu -pk gpu 1 \
  -in in.lammps -log run_in_lammps.log
```
