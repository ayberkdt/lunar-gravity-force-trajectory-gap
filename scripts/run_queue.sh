#!/bin/sh
# The population queue after design C, in one process so nothing idles between
# campaigns and nothing overlaps them: two campaigns on eleven workers each
# would halve both rather than finish either.
#
# Order: the out-of-box operational elliptical population first (it is what the
# widest geometry question needs, and it is cheap -- every perilune above 80 km,
# so every truth runs at degree 300), then the five geometry strata cheapest
# first. Each supervisor refuses to start a stage it cannot finish by the stop.
#
# The wait condition is a file on disk, not a process lookup. pgrep does not see
# Windows process command lines from this shell, so a process test reports the
# design C supervisor as gone the moment it is asked, which is exactly the
# mistake that started this queue on top of a running campaign once already.
cd "$(dirname "$0")" || exit 1

DONE_MARK="../metrics/r19_equal_total_work_C_beta_1.00.json"

while true; do
    if grep -q "=== campaign done" r29_campaign.log 2>/dev/null; then
        echo "[queue] design C supervisor logged completion"
        break
    fi
    if [ -f "$DONE_MARK" ]; then
        echo "[queue] design C reached its last budget"
        break
    fi
    sleep 60
done

echo "[queue] starting R31 operational elliptical at $(date)"
python rev30_campaign.py --registry r31 --stop-at "2026-08-02 21:30" \
    --workers 11 >> r31_campaign.out 2>&1
echo "[queue] R31 finished rc=$? at $(date)"

echo "[queue] starting R30 strata at $(date)"
python rev30_campaign.py --stop-at "2026-08-02 21:30" \
    --workers 11 >> r30_campaign.out 2>&1
echo "[queue] R30 finished rc=$? at $(date)"
