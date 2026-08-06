#!/usr/bin/env bash
# R14 overnight, second pass. Two faults in the first pass are corrected here:
#
#   * the idle-machine guard used `pgrep`, which does not exist in this Git-Bash
#     environment, so it passed instantly and the serial timing panel would have
#     been measured against a busy machine. It now uses the repository's own
#     psutil-based scan (wait_for_idle.py), which fails loudly on timeout.
#   * the timing panel itself crashed on a KeyError: the R13 high-degree cost
#     curve stores its table under `combined_rows`, not `rows`. Fixed in
#     rev14_timing_budget.py.
#
# Order is corrected too. The first pass had the MANDATORY measured-time panel
# (Priority 3) queued behind beta = 3.00 (Priority 4). It now runs first, and
# beta = 3.00 takes whatever time is left under its deadline guard.
set -u
LOCK=.r14_overnight2.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another launcher holds $LOCK -- refusing to start"; exit 9
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

W=5
D300="2026-07-26T08:20:00+03:00"
stamp () { date "+%H:%M:%S"; }

echo "=== [$(stamp)] stage 1/4: wait for the variational panel (Priority 6, in flight) ==="
python wait_for_idle.py --timeout-min 200 --poll-s 60
IDLE=$?
echo "--- [$(stamp)] idle guard exit $IDLE ---"

if [ "$IDLE" -ne 0 ]; then
  echo "!! machine never went idle; SKIPPING the timing panel rather than"
  echo "!! recording contended kernel times. Priority 3 stays 'not run'."
else
  echo "=== [$(stamp)] stage 2/4: serial measured-time panel (Priority 3) ==="
  python rev14_timing_budget.py select     || echo "!! select failed"
  python rev14_timing_budget.py run --repeats 3 || echo "!! run failed"
  python rev14_timing_budget.py aggregate   || echo "!! aggregate failed"
  echo "--- [$(stamp)] stage 2 exit $? ---"
fi

echo "=== [$(stamp)] stage 3/4: design A beta 3.00 (deadline $D300) ==="
python rev14_budget_trajectory.py run --design A --beta 3.00 --workers $W --deadline "$D300"
echo "--- [$(stamp)] stage 3 exit $? ---"

echo "=== [$(stamp)] stage 4/4: finalize ==="
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
echo "latex errors  : $(grep -hcE '^! ' main.log supplement.log | paste -sd'+' | bc 2>/dev/null || echo '?')"
echo "=== [$(stamp)] R14_OVERNIGHT2_DONE ==="
