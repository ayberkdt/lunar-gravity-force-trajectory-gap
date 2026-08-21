"""Wait for the running R53 supervisor to exit, then relaunch it for SF.

The supervisor now on the machine was started with a 16:00 cutoff and will not
start a cell it cannot fit before then. SE is expected to land close enough to
that cutoff that SF would be skipped as declared-and-not-run, which is a
bookkeeping outcome rather than a real one: nothing is wrong with SF, there is
simply no clock left under that supervisor.

So this watcher does one thing. It waits for that supervisor to exit, checks
whether the frozen-like ladder actually landed, and if it did not, starts a new
supervisor with a later cutoff. The new one skips every cell already on disk, so
it costs nothing to run and picks up exactly the work that was left.

It refuses to start a second supervisor while one is alive, and it exits without
launching if SF is already complete or if the disk is below the campaign's own
floor, so a relaunch cannot be what fills the volume.

Usage:
    python rev53_relaunch_watch.py --stop-at "2026-08-14 20:00" --workers 11
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"

TARGET_KEY = "SF"
BETA = "0.62"
POLL_S = 30.0
DISK_FLOOR_MB = 500.0
LOG = HERE / "r53_relaunch_watch.log"


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def ladder_done(key: str) -> bool:
    p = METRICS / f"r19_equal_total_work_{key}_beta_{BETA}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def pid_alive(pid: int) -> bool:
    """Is that process id still running?

    tasklist rather than wmic, which Windows 11 no longer ships, and rather than
    os.kill(pid, 0), which does not mean on Windows what it means elsewhere.
    """
    out = subprocess.run(
        ["tasklist", "/FI", f"PID eq {pid}", "/NH"],
        capture_output=True, text=True).stdout
    return str(pid) in out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--pid", type=int, required=True,
                    help="the supervisor to wait on")
    a = ap.parse_args()

    log(f"watching pid {a.pid}; when it exits, {TARGET_KEY} runs with "
        f"cutoff {a.stop_at}")

    while pid_alive(a.pid):
        time.sleep(POLL_S)

    log("no supervisor is running any more")

    if ladder_done(TARGET_KEY):
        log(f"{TARGET_KEY} is already on disk; nothing to relaunch")
        return 0

    free = shutil.disk_usage(METRICS).free / (1024.0 * 1024.0)
    if free < DISK_FLOOR_MB:
        log(f"{free:.0f} MB free, below the campaign floor of "
            f"{DISK_FLOOR_MB:.0f}; not relaunching")
        return 1

    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    out = ROOT / "output" / f"r53_campaign_{stamp}_sf.stdout.log"
    cmd = [sys.executable, "launch_detached.py", str(out),
           "rev53_campaign.py", "--stop-at", a.stop_at,
           "--workers", str(a.workers)]
    log(f"relaunching: {' '.join(cmd[1:])}")
    rc = subprocess.call(cmd, cwd=str(HERE))
    log(f"launcher returned {rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
