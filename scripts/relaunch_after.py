#!/usr/bin/env python3
"""Wait for a running campaign supervisor to exit, then start another one.

A supervisor's wall clock cannot be extended from outside. It reads --stop-at
once, at start, and the value lives only in its own memory; the ladder stage
is handed the remaining minutes as a deadline when it starts, so a stop time
set too early truncates the ladder rather than merely ending the session.

Killing the supervisor to pass a new one is the obvious move and the wrong one
here. The operating-point stage keeps every finished orbit in memory and
writes its record only after the last one lands, so a kill at 24 of 64 orbits
costs all 24 -- about eighty minutes at eight workers.

So this waits instead. The supervisor is idempotent: it skips a base, an
operating point, a calibration or a ladder that is already complete on disk,
and it refuses to start a stage it cannot finish. A second run therefore picks
up exactly where the first stopped, with the later clock.

Usage:
    python relaunch_after.py <pid> <stop-at> [--workers N] [--strata S]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import psutil

HERE = Path(__file__).resolve().parent
LOG = HERE / "relaunch_after.log"
POLL_SECONDS = 60
# If the watched process is still alive well past its own stop time it is
# finishing a stage it started in time, which is allowed; this is only a guard
# against waiting on a pid that will never exit.
GIVE_UP_HOURS = 24.0


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pid", type=int)
    ap.add_argument("stop_at")
    ap.add_argument("--workers", default="8")
    ap.add_argument("--strata", default="")
    ap.add_argument("--registry", default="r30")
    a = ap.parse_args()

    log(f"watching pid {a.pid}; will relaunch with stop {a.stop_at}")
    t0 = time.time()
    while True:
        if not psutil.pid_exists(a.pid):
            break
        try:
            proc = psutil.Process(a.pid)
            if proc.status() == psutil.STATUS_ZOMBIE:
                break
        except psutil.NoSuchProcess:
            break
        if (time.time() - t0) / 3600.0 > GIVE_UP_HOURS:
            log(f"pid {a.pid} still alive after {GIVE_UP_HOURS} h; giving up")
            return 1
        time.sleep(POLL_SECONDS)

    log(f"pid {a.pid} gone after {(time.time()-t0)/60:.0f} min; relaunching")
    cmd = [sys.executable, "rev30_campaign.py", "--registry", a.registry,
           "--stop-at", a.stop_at, "--workers", a.workers]
    if a.strata:
        cmd += ["--strata", a.strata]
    out = HERE / "r30_relaunch.log"
    with out.open("ab", buffering=0) as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL)
    log(f"relaunched supervisor exited rc={rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
