"""Run both R68 arms in one exclusive session, endpoint first.

The arms are timed, so they cannot share the machine with each other or with
anything else: they run in this one process, one cell at a time, in the
registered order. The endpoint arm goes first because it is the arm the
manuscript's principal claim rests on; if the window closes early, the
campaign has the arm that matters and the interior arm resumes from its state
file next time.

Each arm's pipeline writes its state after every cell, so a kill at any point
loses at most the cell in flight.

Usage:
    python rev68_chain.py [--cutoff 2026-08-22T09:00]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rev68_timing_full as r68            # noqa: E402


class _Args:
    def __init__(self, arm, cutoff):
        self.arm = arm
        self.cutoff = cutoff


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cutoff", default=None,
                    help="ISO 8601; each arm stops cleanly before the next "
                         "cell once passed")
    args = ap.parse_args()
    t0 = time.time()
    for arm in ("endpoint", "interior"):
        print(f"=== r68 arm {arm} starting "
              f"({(time.time()-t0)/3600:.2f} h into the session)", flush=True)
        rc = r68.pipeline(_Args(arm, args.cutoff))
        print(f"=== r68 arm {arm} returned {rc} "
              f"({(time.time()-t0)/3600:.2f} h)", flush=True)
        if rc not in (0,):
            print("[r68-chain] stopping: the arm did not complete cleanly",
                  flush=True)
            return rc
    print(f"[r68-chain] both arms done in {(time.time()-t0)/3600:.2f} h",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
