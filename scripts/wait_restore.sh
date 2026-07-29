#!/usr/bin/env bash
# After the design-A restore finishes, refresh the R11 integrity manifest so it
# indexes the complete 768-trajectory design-A tree.
set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
until grep -q "^\[r11\] finished" rev11_designA_restore.log 2>/dev/null; do sleep 60; done
echo "[wait] $(grep '^\[r11\] finished' rev11_designA_restore.log | tail -1)"
"$PY" rev11_finalize_manifest.py
