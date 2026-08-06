#!/usr/bin/env bash
# Run the corrected R20 once the overnight queue reports done. Waits on the log
# rather than the process table, which this shell cannot see.
set -u
cd "$(dirname "$0")"
until grep -q "overnight queue complete" chain_after_r19.log 2>/dev/null; do sleep 60; done
echo "[r20chain $(date +%H:%M:%S)] starting corrected R20"
python -u rev20_span_longarc.py run --workers 11 --deadline-min 230 \
    > rev20_span_longarc.log 2>&1
echo "[r20chain $(date +%H:%M:%S)] R20 done (exit $?)"
