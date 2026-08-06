"""R25 bisection stage: the full chain at beta = 0.62, design A then design B.

The beta = 0.75 result put the sign change between 0.50 and 0.75. This halves
that bracket. The midpoint is taken as 0.62 rather than 0.625 because every
driver tags records with f"beta_{beta:.2f}", so 0.625 would write files named
beta_0.62 while carrying 0.625 inside them. No design has an archived budget
record at this value, so both run the whole chain rather than the two stages
design A needed at 0.75:

    rev14_budget_trajectory  ->  rev18_span_sweep  ->  rev19_equal_total_work

The chain is written once here and called per design, because the three stages
differ only in which design they name.

Ordering is design A first. A is where the bracket was established, so a
complete A chain localizes the crossing on the same design that produced the
bracket; B, if the window allows, tests whether the localization holds on an
independently scrambled design.

A design is only started if the window can plausibly hold its whole chain. A
half-finished budget record would leave the two downstream stages running on a
partial panel, and a bracket localized on a truncated panel is worth less than
no bracket at all.

Usage:
    python rev25_bisect_queue.py --stop-at "2026-07-31 08:30"
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "python_codes"
METRICS = ROOT / "metrics"
LOG = CODE / "r25_queue.log"

BETA = 0.62
TAG = f"beta_{BETA:.2f}"
RESERVE_MIN = 20.0
# a chain needs roughly two hours at the rates the earlier stages measured;
# below this there is no point starting one
CHAIN_MIN = 110.0

# every way the design-B queue can end
B_DONE = re.compile(
    r"stage B complete|stage B skipped|B2 and B3 cannot run|"
    r"no useful window for B1|window closed after B")


def say(msg: str) -> None:
    line = f"[r25C {datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def minutes_left(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def wait_for_stage_B(stop_at: datetime) -> bool:
    say("waiting for the design-B queue to finish before starting")
    while True:
        try:
            if B_DONE.search(LOG.read_text(encoding="utf-8", errors="replace")):
                say("design-B queue has ended")
                return True
        except OSError:
            pass
        if minutes_left(stop_at) - RESERVE_MIN <= CHAIN_MIN:
            say("window can no longer hold a chain; bisection step skipped")
            return False
        time.sleep(60)


def run(cmd: list[str], label: str) -> int:
    say(f"{label}: $ {' '.join(cmd[1:])}")
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    for stream in (proc.stdout, proc.stderr):
        if stream:
            with LOG.open("a", encoding="utf-8") as fh:
                fh.write(stream)
            tail = [ln for ln in stream.strip().splitlines() if ln.strip()]
            for ln in tail[-6:]:
                print("   " + ln, flush=True)
    say(f"{label}: exit {proc.returncode} after {(time.time()-t0)/60:.1f} min")
    return proc.returncode


def chain(design: str, stop_at: datetime, workers: int) -> int:
    """Budget campaign, span sweep, realized-work comparison, in that order."""
    py = sys.executable
    worst = 0

    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= CHAIN_MIN:
        say(f"design {design}: {left:.0f} min left, not enough for a chain; "
            f"skipped rather than truncated")
        return 0
    say(f"design {design}: chain opens with {left:.0f} min")

    # R14 parses a naive deadline as UTC, so the local stop is made explicit.
    b1_stop = datetime.now() + timedelta(minutes=min(left * 0.45, left))
    worst = max(worst, run(
        [py, str(CODE / "rev14_budget_trajectory.py"), "run",
         "--design", design, "--beta", str(BETA), "--workers", str(workers),
         "--deadline", b1_stop.astimezone().isoformat()],
        f"{design}1 budget"))

    record = METRICS / f"r14_trajectory_{design}_{TAG}.json"
    if not record.exists():
        say(f"design {design}: {record.name} not produced; chain stops here")
        return max(worst, 1)

    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 15:
        say(f"design {design}: window closed after the budget stage")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev18_span_sweep.py"), "run",
         "--design", design, "--beta", str(BETA), "--workers", str(workers),
         "--deadline-min", f"{min(left * 0.65, left):.0f}"],
        f"{design}2 span"))
    run([py, str(CODE / "rev18_span_sweep.py"), "summarize",
         "--design", design, "--beta", str(BETA), "--from-disk"],
        f"{design}2 summary")

    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 10:
        say(f"design {design}: window closed after the span stage")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev19_equal_total_work.py"), "run",
         "--design", design, "--beta", str(BETA), "--workers", str(workers),
         "--deadline-min", f"{left:.0f}"], f"{design}3 equal-work"))
    run([py, str(CODE / "rev19_equal_total_work.py"), "summarize",
         "--design", design, "--beta", str(BETA)], f"{design}3 summary")
    return worst


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=11)
    args = ap.parse_args()

    stop_at = datetime.strptime(args.stop_at, "%Y-%m-%d %H:%M")
    if not wait_for_stage_B(stop_at):
        return 0

    worst = 0
    for design in ("A", "B"):
        worst = max(worst, chain(design, stop_at, args.workers))
    say(f"bisection step done (worst exit {worst}); "
        f"{minutes_left(stop_at):.0f} min to spare")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
