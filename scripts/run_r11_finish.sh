#!/usr/bin/env bash
# Complete the R11 campaign now that the compute window runs to 23:00.
#
# Sequential, so no two stages compete for cores and there is no race with the
# earlier (now-dormant) postprocess orchestrator, which watches a different log.
#   1. finish design-B: resumes from the 40 cached orbits and runs the
#      remaining 24 (the cheapest, since heavy-first already did the N900 ones);
#      fc.run's CancelledError bug is fixed, so the index is written on exit.
#   2. regenerate every manuscript table + descriptives JSON (now full design-B).
#   3. restore the four design-A raw orbit trees (003/028/035/060) deleted during
#      the earlier bug cleanup: rev11 caches on a config hash, so 60 orbits are
#      skipped and only these re-run.
# All stages are idempotent and resumable.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
log() { echo "[finish $(date '+%H:%M:%S')] $*"; }

log "stage 1: completing design-B (resumes 40 cached + 24 remaining)"
"$PY" rev11_designB_convergence.py run --workers 5 \
  --deadline 2026-07-24T22:30:00+03:00 > rev11_designB_run2.log 2>&1
log "design-B exited $? -- $(tail -1 rev11_designB_run2.log)"

log "stage 2: regenerating manuscript tables (full design-B)"
"$PY" rev11_manuscript_tables.py
log "tables regenerated"

log "stage 3: restoring design-A raw sidecars (orbits 003/028/035/060)"
"$PY" rev11_full_convergence.py run --workers 5 > rev11_designA_restore.log 2>&1
log "design-A restore exited $? -- $(tail -1 rev11_designA_restore.log)"

# integrity: no design-B-population leftovers in the design-A tree
"$PY" - <<'PYEOF'
import json, glob
m = "../metrics"
dA = {r['sobol_index']: r['design_point']['initial_state_si']
      for r in json.load(open(m + "/r10_sobolA_baseline_truth_corrected.json"))['rows']}
dB = {r['sobol_index']: r['design_point']['initial_state_si']
      for r in json.load(open(m + "/r11_designB_rows.json"))['rows']}
leak = sum(1 for f in glob.glob(m + "/r11_cases/convergence/*/*.json")
           if 'smoke' not in f
           for d in [json.load(open(f))]
           if d['config']['initial_state_si'] == dB.get(d['config']['sobol_index'])
           and d['config']['initial_state_si'] != dA.get(d['config']['sobol_index']))
print(f"[finish] design-A tree design-B leftovers: {leak}")
PYEOF

log "all R11 stages complete; artifacts ready"
