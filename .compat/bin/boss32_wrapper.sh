#!/usr/bin/env bash
set -euo pipefail

ROOT="/root/lmp-work"
REAL_BOSS_DIR="$ROOT/boss"
QEMU_BIN="$ROOT/.compat/qemu/usr/bin/qemu-i386-static"
LOADER_BIN="$ROOT/.compat/boss32/lib32/ld-linux.so.2"
LIB_PATH="$ROOT/.compat/boss32/lib32:$ROOT/.compat/boss32/usr/lib32"

if [[ ! -x "$QEMU_BIN" ]]; then
  echo "Missing qemu-i386-static: $QEMU_BIN" >&2
  exit 1
fi

if [[ ! -x "$REAL_BOSS_DIR/BOSS" ]]; then
  echo "Missing BOSS binary: $REAL_BOSS_DIR/BOSS" >&2
  exit 1
fi

exec "$QEMU_BIN" "$LOADER_BIN" \
  --library-path "$LIB_PATH" \
  "$REAL_BOSS_DIR/BOSS" "$@"
