"""Block until no competing propagation is running, then exit 0.

The measured-time panel is only valid on an idle machine: R13 measured the same
configurations running 1.49 to 2.56 times slower under a five-worker pool than
serially. The overnight orchestrator's first attempt at this guard used
``pgrep``, which does not exist in this Git-Bash environment, so the check
silently passed. This uses the repository's own psutil-based process scan
instead.

Usage:  python wait_for_idle.py [--timeout-min 240] [--poll-s 60]
"""

from __future__ import annotations

import argparse
import sys
import time

import rev10_sobol_confirmatory as base

# any script that propagates or sweeps the field; the list is a prefix match on
# the command line, so new campaign scripts must be added here or the guard will
# pass while they are still running
WATCH = ("rev14_budget_trajectory", "rev14_variational_budget", "rev14_oracle",
         "rev14_budget_pareto", "rev14_timing_budget",
         "rev15_cadence_check", "rev15_deployable_calibration",
         "rev15_fixed_oracle", "rev15_atallah_grid")


def competing() -> list[str]:
    out = []
    for p in base.other_python_processes():
        cmd = " ".join(p.get("cmdline") or [])
        if any(w in cmd for w in WATCH):
            out.append(f"{p['pid']}:{next(w for w in WATCH if w in cmd)}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout-min", type=float, default=240.0)
    ap.add_argument("--poll-s", type=float, default=60.0)
    a = ap.parse_args()
    deadline = time.time() + a.timeout_min * 60.0
    while True:
        busy = competing()
        if not busy:
            print("[idle] no competing propagation; safe to time", flush=True)
            return 0
        if time.time() > deadline:
            print(f"[idle] TIMEOUT with {len(busy)} still running: {busy}; "
                  "refusing to declare the machine idle", flush=True)
            return 3
        print(f"[idle] waiting on {len(busy)}: {busy}", flush=True)
        time.sleep(a.poll_s)


if __name__ == "__main__":
    sys.exit(main())
