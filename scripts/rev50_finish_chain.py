"""Run everything that is still owed after R50's registered grid, in order.

The machine sat idle for seven hours last night because a run died and the next
step was waiting on a person. This chain removes the gaps between the four
things that are already decided:

  1. the amended budgets, 1.25 and 1.50, on both R50 blocks;
  2. the ceiling probe, re-run on an idle machine so its provenance check reads
     the parent inputs by hash rather than racing a campaign that is writing its
     neighbours;
  3. the R51 registration, which refuses to write itself unless the probe
     supports a control at degree 600;
  4. the R51 campaign.

Each step is skipped when its output is already on disk, and a step that fails
stops the chain rather than letting the next one run on a missing input. Step 3
failing is not necessarily an error: the freeze script is written to refuse when
the probe says the ceiling never binds, and a refusal there is a finding, not a
fault, so it is logged as one.

The chain can be cut short at any of its four steps. Stopping after the probe
leaves the registration and the control to be started by hand, which is the
right default when the person who decides whether eight hours of machine time
buys anything wants to read the probe first.

Usage:
    python rev50_finish_chain.py --stop-at "2026-08-13 12:00" --workers 11
    python rev50_finish_chain.py --stop-at "..." --stop-after probe
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"
LOG = HERE / "r50_finish_chain.log"

AMENDED = [1.25, 1.50]
BLOCK_KEYS = {"span_ladder_a": "RS1", "span_ladder_b": "RS2"}


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def run(name: str, cmd: list[str]) -> int:
    log(f"START {name}: {' '.join(cmd[1:])}")
    t0 = time.time()
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT)
    log(f"END   {name}: rc={rc} in {(time.time()-t0)/60:.1f} min")
    return rc


def ladder_done(key: str, beta: float) -> bool:
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--stop-after",
                    choices=("amended", "probe", "freeze", "campaign"),
                    default="campaign",
                    help="last step to run; the chain stops cleanly after it")
    a = ap.parse_args()
    py = sys.executable
    steps = ["amended", "probe", "freeze", "campaign"]
    last = steps.index(a.stop_after)
    log(f"=== finish chain: stop at {a.stop_at}, workers {a.workers}, "
        f"running {', '.join(steps[:last + 1])} ===")

    # 1. the amended budgets, on whichever blocks carry the registered grid.
    #    rev50_campaign refuses the amendment for a block whose registered
    #    budgets are not all present, so an incomplete block is skipped by the
    #    gate rather than by a decision taken here.
    todo = [b for b, k in BLOCK_KEYS.items()
            if any(not ladder_done(k, beta) for beta in AMENDED)]
    if todo:
        rc = run("amended_budgets",
                 [py, "rev50_campaign.py", "--stop-at", a.stop_at,
                  "--workers", str(a.workers), "--betas", "1.25,1.50",
                  "--blocks", ",".join(todo)])
        if rc != 0:
            log("the amended budgets did not complete; the chain stops here "
                "rather than registering a control on a moving parent")
            return 1
    else:
        log("amended budgets already on disk")

    if last < steps.index("probe"):
        log("stopping after the amended budgets, as asked")
        return 0

    # 2. the probe, on an idle machine.
    rc = run("ceiling_probe",
             [py, "rev51_ceiling_probe.py", "--workers", str(a.workers)])
    if rc != 0:
        log("the ceiling probe failed; R51 cannot be registered without it")
        return 1

    if last < steps.index("freeze"):
        log("stopping after the ceiling probe, as asked. The registration and "
            "the control are left to be started by hand: read "
            "metrics/r51_ceiling_probe.json first, then "
            "rev51_uncapped_freeze.py and rev51_campaign.py.")
        return 0

    # 3. the registration, which is allowed to refuse.
    rc = run("r51_freeze", [py, "rev51_uncapped_freeze.py"])
    if rc != 0:
        log("the R51 registration refused. Read r51_ceiling_probe.json: the "
            "freeze script refuses when 600 does not clear the demand and when "
            "the ceiling never binds, and the second of those is a result "
            "about the ladder rather than a failure of the chain.")
        return 0

    if last < steps.index("campaign"):
        log("stopping after the registration, as asked")
        return 0

    # 4. the control.
    rc = run("r51_campaign",
             [py, "rev51_campaign.py", "--stop-at", a.stop_at,
              "--workers", str(a.workers)])
    log(f"=== finish chain done, rc={rc} ===")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
