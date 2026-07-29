#!/usr/bin/env bash
# R18 budget robustness: does the interior optimum survive a change of budget?
#
# beta = 1.00 is already archived for both designs. This queue adds the
# compute-starved and compute-rich ends, using the R14 records that exist:
#   beta = 0.50  designs A and B   (radial rule loses catastrophically there)
#   beta = 1.50  design A only     (no R14 design-B record at that budget)
#
# Runs are sequential so the pool is never oversubscribed. Each stage has its
# own deadline; a stage that overruns is cut and leaves the earlier stages
# intact, because every stage writes its own summary file.
set -u
cd "$(dirname "$0")"

log() { echo "[queue $(date +%H:%M:%S)] $*"; }

stage() {   # design beta deadline_min
  local d=$1 b=$2 dl=$3
  log "start design=$d beta=$b deadline=${dl}min"
  python -u rev18_span_sweep.py run --design "$d" --beta "$b" \
      --workers 11 --deadline-min "$dl" \
      > "rev18_span_${d}_beta_${b}.log" 2>&1
  log "done  design=$d beta=$b (exit $?)"
}

stage A 0.50 110
stage B 0.50 110
stage A 1.50 220

log "queue complete"
