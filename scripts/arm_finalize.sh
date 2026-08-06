#!/bin/sh
# Run finalize_all.sh once the population queue is done, or once the campaign's
# stop time has passed, whichever comes first. The second condition is the one
# that matters: if the queue dies for a reason nobody is awake to see, the
# results still get assembled from whatever reached the disk.
cd "$(dirname "$0")" || exit 1

DEADLINE=$(date -d "2026-08-02 21:45" +%s 2>/dev/null || echo 0)

while true; do
    if grep -q "R30 finished" run_queue.log 2>/dev/null; then
        echo "[arm] queue reported finished"
        break
    fi
    if [ "$DEADLINE" -gt 0 ] && [ "$(date +%s)" -ge "$DEADLINE" ]; then
        echo "[arm] campaign stop time passed; finalizing whatever is on disk"
        break
    fi
    sleep 120
done

sh finalize_all.sh
