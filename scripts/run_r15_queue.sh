#!/usr/bin/env bash
# R15 queue: the two P0 items from the audit-response plan, plus the cheap
# grid-convergence fix, run after the cadence check releases the machine.
#
#   1. R15-B  deployable (truth-free) budget calibration, design A then B
#   2. R15-A  budget-saturating and best-under-budget fixed comparators
#   3. finalize: tables, figure, deliverables, both PDFs
#
# Each stage is guarded: a failure is reported and the queue continues, so one
# broken stage cannot silently take the others down with it.
set -u
LOCK=.r15_queue.lock
if ! mkdir "$LOCK" 2>/dev/null; then echo "another R15 queue holds $LOCK"; exit 9; fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT
stamp () { date "+%H:%M:%S"; }

echo "=== [$(stamp)] waiting for the cadence check to release the machine ==="
python wait_for_idle.py --timeout-min 90 --poll-s 30
echo "--- [$(stamp)] idle guard exit $? ---"

echo "=== [$(stamp)] stage 1/3: deployable budget calibration, design A ==="
python rev15_deployable_calibration.py run --design A --workers 5 || echo "!! A failed"
echo "=== [$(stamp)] stage 1b: deployable budget calibration, design B ==="
python rev15_deployable_calibration.py run --design B --workers 5 || echo "!! B failed"

echo "=== [$(stamp)] stage 2/3: fixed-degree saturating and oracle comparators ==="
python rev15_fixed_oracle.py run --orbits 16 --workers 5 || echo "!! oracle failed"

echo "=== [$(stamp)] stage 3/3: finalize ==="
python rev14_tables.py
python make_figures_r14.py
python rev14_finalize_manifest.py
python rev14_deliverables.py
python check_assets.py
cd ..
latexmk -pdf -interaction=nonstopmode supplement.tex >/dev/null 2>&1
latexmk -pdf -interaction=nonstopmode main.tex      >/dev/null 2>&1
echo "main : $(grep -oE 'Output written on main.pdf \([0-9]+ pages' main.log)"
echo "supp : $(grep -oE 'Output written on supplement.pdf \([0-9]+ pages' supplement.log)"
echo "=== [$(stamp)] R15_QUEUE_DONE ==="
