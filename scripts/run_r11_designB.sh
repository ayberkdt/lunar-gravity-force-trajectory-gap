#!/usr/bin/env bash
# design-B independent confirmation: prepass -> full convergence.
#
# The prepass (n_work / n_critical / sub-50 km degree-adequacy) is already
# running when this starts.  This waits for it to write the rows file, then
# launches the 64-orbit x 6-policy x 2-level vector-tolerance convergence with
# a deadline guard well inside 17:00.  Both phases cache on config hashes, so
# re-running after an interruption resumes.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
log() { echo "[designB $(date '+%H:%M:%S')] $*"; }

log "waiting for prepass to finish"
until grep -q "^\[prepass\] wrote" rev11_designB_prepass.log 2>/dev/null; do
  sleep 30
done
log "prepass done: $(grep '^\[prepass\] wrote' rev11_designB_prepass.log)"

# Abort if the prepass did not produce a full 64-row population.
if ! grep -q "wrote 64/64 rows" rev11_designB_prepass.log; then
  log "prepass incomplete; not starting the convergence run"
  exit 1
fi

log "launching design-B convergence (5 workers, deadline 16:30 local)"
"$PY" rev11_designB_convergence.py run --workers 5 \
  --deadline 2026-07-24T16:30:00+03:00 > rev11_designB_convergence.log 2>&1
status=$?
log "design-B convergence exited with status $status"
"$PY" rev11_designB_convergence.py status
log "finished"
