#!/usr/bin/env bash
# R14 overnight finish, ordered to fit a 09:00 deadline and to respect the
# pre-registered priority order rather than the order the first launcher used.
#
#   1. A beta = 1.50            (Priority 4) -- extends the trajectory trend
#   2. wait for the variational (Priority 6) -- already running in parallel
#   3. serial measured-time panel (Priority 3, MANDATORY) -- needs an idle
#      machine, so nothing else may run during it
#   4. A beta = 3.00            (Priority 4) -- most expensive, deadline-guarded
#   5. finalize: tables, figure, manifest, deliverables, both PDFs
#
# Stages 1 and 4 carry --deadline, so an overrun stops cleanly and writes a
# partial payload with stopped_for_deadline = true rather than being killed.
set -u
LOCK=.r14_overnight.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another R14 overnight launcher holds $LOCK -- refusing to start"; exit 9
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

W=5
D150="2026-07-26T03:46:00+03:00"
D300="2026-07-26T08:35:00+03:00"
VAR_JSON=../metrics/r14_variational_budget.json

stamp () { date "+%H:%M:%S"; }

echo "=== [$(stamp)] stage 1/5: design A beta 1.50 (deadline $D150) ==="
python rev14_budget_trajectory.py run --design A --beta 1.50 --workers $W --deadline "$D150"
echo "--- [$(stamp)] stage 1 exit $? ---"

echo "=== [$(stamp)] stage 2/5: waiting for the variational panel to finish ==="
while [ ! -f "$VAR_JSON" ]; do
  if ! pgrep -f rev14_variational_budget >/dev/null 2>&1; then
    echo "  variational process is gone and produced no output; continuing"
    break
  fi
  sleep 60
done
echo "--- [$(stamp)] variational done or absent ---"

echo "=== [$(stamp)] stage 3/5: serial measured-time panel (machine must be idle) ==="
# guard: the timing claim is only valid with no competing propagation
for i in $(seq 1 30); do
  n=$(pgrep -f 'rev14_budget_trajectory|rev14_variational_budget|rev14_oracle' | wc -l)
  [ "$n" -eq 0 ] && break
  echo "  waiting for $n competing python process(es) to exit"; sleep 30
done
python rev14_timing_budget.py select
python rev14_timing_budget.py run --repeats 3
python rev14_timing_budget.py aggregate
echo "--- [$(stamp)] stage 3 exit $? ---"

echo "=== [$(stamp)] stage 4/5: design A beta 3.00 (deadline $D300) ==="
python rev14_budget_trajectory.py run --design A --beta 3.00 --workers $W --deadline "$D300"
echo "--- [$(stamp)] stage 4 exit $? ---"

echo "=== [$(stamp)] stage 5/5: finalize ==="
python rev14_tables.py
python make_figures_r14.py
python rev14_finalize_manifest.py
python rev14_deliverables.py
cd ..
latexmk -pdf -interaction=nonstopmode supplement.tex >/dev/null 2>&1
latexmk -pdf -interaction=nonstopmode main.tex      >/dev/null 2>&1
latexmk -pdf -interaction=nonstopmode main.tex      >/dev/null 2>&1
echo "main.pdf   : $(grep -oE 'Output written on main.pdf \([0-9]+ pages' main.log)"
echo "supplement : $(grep -oE 'Output written on supplement.pdf \([0-9]+ pages' supplement.log)"
echo "latex errors: $(grep -cE '^! ' main.log supplement.log | paste -sd' ')"
echo "undefined refs: $(grep -coE \"Reference \\\`[^']+' on page [0-9]+ undefined\" main.log)"
echo "=== [$(stamp)] R14_OVERNIGHT_DONE ==="
