#!/usr/bin/env bash
# Post-design-B mechanical preparation for the coordinated manuscript pass.
#
# Waits for the design-B convergence run to finish, then does the two
# mechanical steps so that the (model-driven) prose update starts from ready,
# consistent artifacts:
#   1. regenerate every manuscript LaTeX table + the descriptives JSON,
#      now including the design-B replication table;
#   2. restore the four design-A raw orbit trees (003/028/035/060) that were
#      deleted while cleaning up the earlier path-routing bug -- rev11 caches on
#      a config hash, so the other 60 orbits are skipped and only these re-run.
#
# Detection is by log grep (pgrep cannot see native Windows python from Git
# Bash).  Everything is idempotent and resumable.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
log() { echo "[postproc $(date '+%H:%M:%S')] $*"; }

log "waiting for design-B convergence to finish"
until grep -q "^\[r11\] finished" rev11_designB_convergence.log 2>/dev/null; do
  sleep 60
done
log "design-B finished: $(grep '^\[r11\] finished' rev11_designB_convergence.log | tail -1)"

log "regenerating manuscript tables (now including design-B)"
"$PY" rev11_manuscript_tables.py
log "tables regenerated"

log "restoring design-A raw sidecars for orbits 003/028/035/060"
"$PY" rev11_full_convergence.py run --workers 5 \
  > rev11_designA_restore.log 2>&1
log "design-A restore exited with status $? -- $(tail -1 rev11_designA_restore.log)"

# final integrity: design-A tree should have no design-B-population leftovers
"$PY" - <<'PYEOF'
import json, glob, os
m = os.path.join(os.path.dirname(os.getcwd()), "metrics") if os.path.basename(os.getcwd())=="python_codes" else "../metrics"
m = "../metrics"
dA = {r['sobol_index']: r['design_point']['initial_state_si']
      for r in json.load(open(m + "/r10_sobolA_baseline_truth_corrected.json"))['rows']}
dB = {r['sobol_index']: r['design_point']['initial_state_si']
      for r in json.load(open(m + "/r11_designB_rows.json"))['rows']}
leak = 0
for f in glob.glob(m + "/r11_cases/convergence/*/*.json"):
    if 'smoke' in f: continue
    try: d = json.load(open(f))
    except: continue
    i = d['config']['sobol_index']; st = d['config']['initial_state_si']
    if st == dB.get(i) and st != dA.get(i): leak += 1
print(f"[postproc] design-A tree design-B leftovers: {leak}")
PYEOF

log "post-processing complete; artifacts ready for the prose pass"
