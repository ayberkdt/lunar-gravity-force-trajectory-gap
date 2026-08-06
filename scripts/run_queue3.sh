#!/bin/sh
# Fill the window to 09:20 and stop there, in decreasing value per hour.
#
# The order is not the frozen order; it is what a 7-hour window can actually
# finish. high_apolune holds no orbit below 80 km perilune, so it is a whole
# population for the price of half of one. equatorial and frozen_like cannot fit
# a ladder in what is left after it, so they get their operating point and
# calibration only -- an hour of work that a later window would otherwise have
# to spend before it could start propagating.
#
# The tail of the window goes to the low_perilune base, which is the one job
# here that is resumable at trajectory granularity: every finished trajectory is
# a sidecar and a raw array on disk, and the next run skips what is already
# there. Stopping it mid-population costs nothing, which is exactly what makes
# it the right thing to put last.
#
# Nothing starts unless the clock says it can finish. HARD is the wall.
cd "$(dirname "$0")" || exit 1

HARD="2026-08-03 09:20"
HARD_S=$(date -d "$HARD" +%s)

left_min() { echo $(( (HARD_S - $(date +%s)) / 60 )); }

can_start() {   # can_start <minutes needed> <label>
    L=$(left_min)
    if [ "$L" -lt "$1" ]; then
        echo "[queue3] skip $2: $L min left, needs $1"
        return 1
    fi
    echo "[queue3] start $2 with $L min left"
    return 0
}

while ! grep -q "\[queue2\] done" run_queue2.log 2>/dev/null; do
    [ "$(left_min)" -lt 5 ] && { echo "[queue3] wall reached while waiting"; exit 0; }
    sleep 60
done

if can_start 95 "high_apolune (full ladder)"; then
    python rev30_campaign.py --strata high_apolune --stop-at "$HARD" \
        --workers 11 >> r30_campaign.out 2>&1
fi

for s in equatorial frozen_like; do
    if can_start 40 "$s operating point"; then
        python rev30_stratum_ops.py --stratum "$s" op --workers 11 \
            >> r30_campaign.out 2>&1
    fi
    if can_start 15 "$s calibration"; then
        python rev30_stratum_ops.py --stratum "$s" cal --workers 11 \
            >> r30_campaign.out 2>&1
    fi
done

if can_start 20 "low_perilune base (resumable)"; then
    python rev30_stratum_base.py --stratum low_perilune run --workers 11 \
        --deadline "$(date -d "$HARD" +%Y-%m-%dT%H:%M:%S%:z)" \
        >> r30_campaign.out 2>&1
fi

echo "[queue3] propagation stopped at $(date)"
sh finalize_all.sh
echo "[queue3] done at $(date)"
