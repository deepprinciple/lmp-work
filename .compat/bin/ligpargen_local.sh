#!/usr/bin/env bash
set -euo pipefail

ROOT="/root/lmp-work"
CONDA_BIN="/root/miniconda3/bin/conda"
COMPAT_BIN="$ROOT/.compat/bin"
BOSS_SRC="$ROOT/boss"
BOSS_RUNTIME="$ROOT/.compat/boss_runtime"
LIGPARGEN_SRC="$ROOT/.vendor/ligpargen"

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "Missing conda executable: $CONDA_BIN" >&2
  exit 1
fi

if [[ ! -d "$LIGPARGEN_SRC" ]]; then
  echo "Missing LigParGen source tree: $LIGPARGEN_SRC" >&2
  exit 1
fi

if [[ ! -d "$BOSS_SRC" ]]; then
  echo "Missing BOSS source tree: $BOSS_SRC" >&2
  exit 1
fi

mkdir -p "$BOSS_RUNTIME"

for name in \
  00README.txt Copyright.txt Licenses autozmat bossman miscexec molecules \
  notes oplsaa.par oplsaa.sb oplsua.par oplsua.sb org1box org2box par \
  par.form sb scripts solbox testjobs watbox xx yy
do
  if [[ -e "$BOSS_SRC/$name" && ! -e "$BOSS_RUNTIME/$name" ]]; then
    ln -s "$BOSS_SRC/$name" "$BOSS_RUNTIME/$name"
  fi
done

ln -sf "$COMPAT_BIN/boss32_wrapper.sh" "$BOSS_RUNTIME/BOSS"

export PATH="$COMPAT_BIN:$PATH"
export BOSSdir="$BOSS_RUNTIME"
export PYTHONPATH="$LIGPARGEN_SRC${PYTHONPATH:+:$PYTHONPATH}"

exec "$CONDA_BIN" run -n HTMD python -m ligpargen.ligpargen "$@"
