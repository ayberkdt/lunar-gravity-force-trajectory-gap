#!/usr/bin/env bash
# Extra run for the extended (~10:00) window: after BOTH Atallah benchmarks
# (design A and design B) finish, run the accuracy-tolerance sweep with a 09:45
# guard. Separate from the overnight orchestrator, so it never disturbs running
# work; it only starts once the benchmarks are done and the cores are free.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
log() { echo "[r12-sweep-after $(date '+%H:%M:%S')] $*"; }

log "waiting for both Atallah benchmarks to finish"
until grep -q "all R12 overnight stages complete" run_r12_overnight.log 2>/dev/null; do
  sleep 60
done
log "benchmarks done; starting the accuracy-tolerance sweep"

"$PY" rev12_atallah_sweep.py run --workers 5 \
  --deadline 2026-07-25T09:45:00+03:00 > rev12_sweep.log 2>&1
log "sweep exited: $(tail -1 rev12_sweep.log)"
log "extended-window work complete"
