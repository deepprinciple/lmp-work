#!/usr/bin/env bash
set -euo pipefail

CUDA11_HOME="/usr/lib/cuda"
LMP_BIN="/opt/lammps/bin/lmp"

# Keep the driver-facing libraries untouched and only force the CUDA 11.5
# user-space toolchain ahead of the incompatible /usr/local/cuda 12.8 install.
export CUDA_HOME="${CUDA11_HOME}"
export PATH="/usr/bin:/usr/sbin:/bin:/sbin:/usr/local/sbin:/usr/local/bin"
export LD_LIBRARY_PATH="${CUDA11_HOME}/lib64:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

exec "${LMP_BIN}" "$@"
