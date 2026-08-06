#!/usr/bin/env bash
# Wait for the R19 queue to write its completion marker, then run the overnight
# continuation. The wait is on the log rather than on a process table, because
# pgrep does not see Windows processes from this shell.
set -u
cd "$(dirname "$0")"
until grep -q "queue complete" r19_queue.log 2>/dev/null; do sleep 60; done
echo "[chain $(date +%H:%M:%S)] R19 queue finished; starting overnight stages"
bash run_overnight_queue.sh
