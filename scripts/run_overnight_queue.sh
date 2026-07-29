#!/usr/bin/env bash
# Overnight continuation. Each stage caches by config hash, so a stage that is
# already complete costs seconds and a stage cut by its deadline resumes.
set -u
cd "$(dirname "$0")"
log() { echo "[night $(date +%H:%M:%S)] $*"; }

log "R20 sixty-day span check (8 orbits x 5 k)"
python -u rev20_span_longarc.py run --workers 11 --deadline-min 260 \
    > rev20_span_longarc.log 2>&1
log "R20 done (exit $?)"

log "R18 budget robustness: A beta=1.50"
python -u rev18_span_sweep.py run --design A --beta 1.50 --workers 11 \
    --deadline-min 150 > rev18_span_A_beta_1.50.log 2>&1
log "A beta=1.50 done (exit $?)"

log "overnight queue complete"
