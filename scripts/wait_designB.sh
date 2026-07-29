#!/usr/bin/env bash
# Block until design-B finishes, then refresh tables and restore design-A.
set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
until grep -q "^\[r11\] finished" rev11_designB_run4.log 2>/dev/null; do sleep 60; done
echo "[wait] $(grep '^\[r11\] finished' rev11_designB_run4.log | tail -1)"
"$PY" rev11_manuscript_tables.py
echo "[wait] tables refreshed"
"$PY" rev11_full_convergence.py run --workers 5 > rev11_designA_restore.log 2>&1
echo "[wait] design-A restore done: $(tail -1 rev11_designA_restore.log)"
