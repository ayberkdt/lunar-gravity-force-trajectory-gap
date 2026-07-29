#!/usr/bin/env bash
# R19: interior member vs a constant degree matched on REALIZED total work.
# Then, with whatever time remains, the R18 budget-robustness stages.
set -u
cd "$(dirname "$0")"
log() { echo "[queue $(date +%H:%M:%S)] $*"; }

log "R19 design A"
python -u rev19_equal_total_work.py run --design A --workers 11 \
    --deadline-min 100 > rev19_A.log 2>&1
log "R19 design A done (exit $?)"

log "R19 design B"
python -u rev19_equal_total_work.py run --design B --workers 11 \
    --deadline-min 100 > rev19_B.log 2>&1
log "R19 design B done (exit $?)"

log "R18 budget robustness: A beta=0.50"
python -u rev18_span_sweep.py run --design A --beta 0.50 --workers 11 \
    --deadline-min 110 > rev18_span_A_beta_0.50.log 2>&1
log "A beta=0.50 done (exit $?)"

log "R18 budget robustness: B beta=0.50"
python -u rev18_span_sweep.py run --design B --beta 0.50 --workers 11 \
    --deadline-min 110 > rev18_span_B_beta_0.50.log 2>&1
log "B beta=0.50 done (exit $?)"

log "queue complete"
