#!/usr/bin/env bash
# Overnight Atallah benchmark: design-A primary (already running) -> design-B
# independent replication -> aggregate tables for both.
#
# The design-A campaign was launched separately and writes r12_atallah_campaign.json
# on exit. This waits for it, then runs the identical benchmark on design B (env
# R12_POPULATION=B points the campaign at design B's rows and reused truth/critical
# trees), then builds the tables for both. Everything caches on a config hash and
# is deadline-guarded, so the machine is never left idle and never overruns.

set -u
cd "$(dirname "$0")"
PY="D:/Masaustu/LUNAR_SIMULATION/.venv/Scripts/python.exe"
log() { echo "[r12-overnight $(date '+%H:%M:%S')] $*"; }

log "waiting for the design-A Atallah campaign to finish"
until grep -q "^\[atallah\] done" rev12_campaign.log 2>/dev/null; do
  sleep 60
done
log "design-A done: $(grep '^\[atallah\] done' rev12_campaign.log | tail -1)"

log "building design-A benchmark tables"
"$PY" rev12_atallah_tables.py || log "design-A table build reported an issue"

log "launching design-B Atallah benchmark (independent replication)"
R12_POPULATION=B "$PY" rev12_atallah_campaign.py run --workers 5 \
  --deadline 2026-07-25T06:45:00+03:00 > rev12_campaign_designB.log 2>&1
log "design-B exited: $(tail -1 rev12_campaign_designB.log)"

log "building design-B benchmark tables"
R12_POPULATION=B "$PY" rev12_atallah_tables.py || log "design-B table build reported an issue"

log "all R12 overnight stages complete"
