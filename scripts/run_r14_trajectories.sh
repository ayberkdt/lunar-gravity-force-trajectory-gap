#!/usr/bin/env bash
# R14 trajectory campaign, staged per the pre-registered priority order.
#   Priority 2 (mandatory): beta = 1 on both populations.
#   Priority 4: pre-declared design-A budgets 0.75, 1.50, 3.00.
#   Priority 5: design B at beta = 1 always, plus beta = 0.50 -- the grid value
#               nearest the force-level crossover the Phase-A sweep placed
#               between 0.50 and 0.75, selected by the pre-registered adaptive
#               extension rule, which also requires the bracketing pair
#               {0.50, 0.75} on design A.
set -u
LOCK=.r14_traj.lock
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another R14 trajectory launcher holds $LOCK -- refusing to start"; exit 9
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT
W=5
run () {
  echo "=== design $1 beta $2 ==="
  python rev14_budget_trajectory.py run --design "$1" --beta "$2" --workers $W
  echo "--- stage exit $? ---"
}
run A 1.00
run B 1.00
run A 0.75
run A 0.50
run B 0.50
run A 1.50
run A 3.00
echo "ALL_R14_TRAJECTORIES_DONE"
