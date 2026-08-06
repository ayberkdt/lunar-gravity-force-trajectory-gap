#!/bin/sh
# Recovery queue: the two populations chosen for this window, in order.
#
# Both bases are already on disk and indexed; what runs here is the operating
# point, the calibration and the declared-budget ladder. The operational
# elliptical population goes first because every one of its orbits sits at the
# lower truth degree, so it is the one certain to finish inside the window;
# polar follows and its supervisor will decline the ladder rather than start one
# it cannot complete.
cd "$(dirname "$0")" || exit 1

STOP="2026-08-03 07:00"

echo "[queue2] operational elliptical at $(date)"
python rev30_campaign.py --registry r31 --stop-at "$STOP" --workers 11 \
    >> r31_campaign.out 2>&1
echo "[queue2] operational elliptical rc=$? at $(date)"

echo "[queue2] polar at $(date)"
python rev30_campaign.py --strata polar --stop-at "$STOP" --workers 11 \
    >> r30_campaign.out 2>&1
echo "[queue2] polar rc=$? at $(date)"

sh finalize_all.sh
echo "[queue2] done at $(date)"
