#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPAT_DIR="${SCRIPT_DIR}/.gpu11_compat_libs"

mkdir -p "${COMPAT_DIR}"

ln -sfn "/lib/x86_64-linux-gnu/libnvidia-ml.so.470.129.06" "${COMPAT_DIR}/libnvidia-ml.so.1"

export LD_LIBRARY_PATH="${COMPAT_DIR}:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu"

exec /usr/bin/nvidia-smi "$@"
