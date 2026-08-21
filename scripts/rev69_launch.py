"""Wait for a quiet machine, then run R69 (O61) to completion.

The timing protocol needs an idle machine, and `idle_or_die` only sees other
python processes: the load that lost R65's first attempt does not have to be
python at all. probe_kernel_quiet.py measures the thing that actually matters,
the kernel's own throughput against the archived idle-machine curve, so this
launcher gates on the probe rather than on the process table.

It polls until the probe passes, then starts the pipeline in this process. If
the machine never quiets down inside the deadline it exits without
propagating anything, because a contaminated timing campaign is worse than no
campaign: R65's discarded attempt cost a night and every one of its cells
reported an in-band ratio while the degree it had selected was wrong.

Usage:
    python rev69_launch.py [--deadline-hours 6] [--poll-minutes 10]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent


def probe() -> bool:
    r = subprocess.run([sys.executable, str(HERE / "probe_kernel_quiet.py")],
                       cwd=HERE, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    tail = [ln for ln in r.stdout.splitlines() if ln.strip()][-1:]
    print(f"[probe] rc={r.returncode} {tail[0] if tail else ''}", flush=True)
    return r.returncode == 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--deadline-hours", type=float, default=6.0)
    ap.add_argument("--poll-minutes", type=float, default=10.0)
    ap.add_argument("--cutoff", default=None,
                    help="ISO 8601 passed through to the pipeline")
    args = ap.parse_args()

    deadline = datetime.now() + timedelta(hours=args.deadline_hours)
    while datetime.now() < deadline:
        if probe():
            print(f"[r69-launch] machine quiet at "
                  f"{datetime.now().isoformat(timespec='seconds')}; starting",
                  flush=True)
            cmd = [sys.executable, str(HERE / "rev69_reference_free.py"),
                   "pipeline"]
            if args.cutoff:
                cmd += ["--cutoff", args.cutoff]
            return subprocess.run(cmd, cwd=HERE).returncode
        print(f"[r69-launch] not quiet; waiting {args.poll_minutes:.0f} min "
              f"(deadline {deadline.isoformat(timespec='minutes')})",
              flush=True)
        time.sleep(args.poll_minutes * 60.0)

    print("[r69-launch] deadline passed without a quiet machine; nothing "
          "was propagated", flush=True)
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
