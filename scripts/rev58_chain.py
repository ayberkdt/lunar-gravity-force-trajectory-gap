"""Runs the R58 endpoint work-matched control over the wide-elliptic budgets.

Sequential on purpose: two campaigns already hold eight and three workers, so
this takes a small fixed slice rather than competing for the machine. Worker
count is the throttle; process priority is not touched.

BLAS threads are pinned in the child environment. Without that, nested
threading inside the pool workers has broken the pool outright on this machine.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (design, beta) in the order the argument needs them: the declared budget
# first, because that is the cell the wide-elliptic claim rests on, then the
# budgets below it that locate the work-matched crossing.
CELLS = [
    ("OEU", 1.00), ("OE", 1.00),
    ("OEU", 0.75), ("OE", 0.75),
    ("OEU", 0.62), ("OE", 0.62),
    ("OEU", 0.50), ("OE", 0.50),
]

# The low-perilune base is paused for the duration of this chain, so the cores
# it held are available here. R57 keeps its three.
WORKERS = 8

# Resumed when the chain finishes, so the machine does not idle behind it. The
# stop time is the one already set for that campaign, not a new one: it simply
# uses whatever is left of its own window.
RESUME_LOW_PERILUNE = [
    "rev30_campaign.py", "--stop-at", "2026-08-16 17:00",
    "--workers", "8", "--strata", "low_perilune",
]

CHILD_ENV = dict(os.environ)
CHILD_ENV.update({
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONUNBUFFERED": "1",
})

NO_WINDOW = 0x08000000 | 0x00000200      # CREATE_NO_WINDOW | NEW_PROCESS_GROUP


def main() -> int:
    t0 = time.time()
    for design, beta in CELLS:
        print(f"\n=== R58 {design} beta={beta:.2f} "
              f"(t+{(time.time()-t0)/60:.0f} min) ===", flush=True)
        cmd = [sys.executable, "rev58_endpoint_equal_work.py", "run",
               "--design", design, "--beta", f"{beta:.2f}",
               "--workers", str(WORKERS), "--deadline-min", "90"]
        rc = subprocess.call(cmd, cwd=str(HERE), env=CHILD_ENV,
                             creationflags=NO_WINDOW)
        print(f"=== R58 {design} beta={beta:.2f} exit={rc} ===", flush=True)
    print(f"\n[r58-chain] done, wall={(time.time()-t0)/60:.1f} min", flush=True)

    print("[r58-chain] resuming the low-perilune base", flush=True)
    subprocess.Popen([sys.executable] + RESUME_LOW_PERILUNE,
                     cwd=str(HERE), env=CHILD_ENV, creationflags=NO_WINDOW,
                     stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
