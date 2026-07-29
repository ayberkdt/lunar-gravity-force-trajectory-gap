#!/usr/bin/env bash
# R11 overnight orchestration.
#
# Stage A (64-orbit vector-tolerance convergence) and the geometric
# verification are already running when this starts; this script waits for
# Stage A to release its five workers and then keeps all six cores busy with
# the remaining campaigns, ordered so the long pole starts first.
#
#   1. corrected_blend/tighter   (long pole, ~6 h)      1 core
#   2. corrected_blend/baseline  (~4 h)                 1 core
#   3. phase sweep, 24 phases    (~2 h)                 4 cores
#   4. once the sweep is done: the four light blend cases
#   5. collect the blend index
#
# Progress is detected from the artifacts each stage writes, not from process
# tables: pgrep cannot see native Windows processes from Git Bash.  Every stage
# caches on a config hash, so re-running this script after an interruption
# re-uses finished work.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
CASES="../metrics/r11_cases/blend_lro_vector"
log() { echo "[orchestrator $(date '+%H:%M:%S')] $*"; }

log "waiting for stage A (rev11_full_convergence) to finish"
until grep -q "^\[r11\] finished" rev11_full_convergence.log 2>/dev/null; do
  sleep 60
done
log "stage A finished: $(grep '^\[r11\] finished' rev11_full_convergence.log)"

log "launching corrected_blend tighter + baseline (long poles)"
nohup "$PY" rev11_blend_lro_vector.py one --policy corrected_blend --level tighter \
  > rev11_blend_ct.log 2>&1 &
nohup "$PY" rev11_blend_lro_vector.py one --policy corrected_blend --level baseline \
  > rev11_blend_cb.log 2>&1 &

log "launching 24-phase sweep with 4 workers"
"$PY" rev11_phase_sweep.py run --workers 4 > rev11_phase_sweep.log 2>&1
log "phase sweep exited with status $?"

log "launching the four light blend cases"
for spec in "truth_N600 tighter" "truth_N600 baseline" \
            "fixed_N120 tighter" "fixed_N120 baseline"; do
  set -- $spec
  nohup "$PY" rev11_blend_lro_vector.py one --policy "$1" --level "$2" \
    > "rev11_blend_$1_$2.log" 2>&1 &
done

log "waiting for all six blend sidecars"
while [ "$(ls "$CASES"/*.json 2>/dev/null | grep -vc smoke)" -lt 6 ]; do
  sleep 60
done

log "collecting blend index"
"$PY" rev11_blend_lro_vector.py collect
log "all R11 stages finished"
