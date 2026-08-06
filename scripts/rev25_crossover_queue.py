"""R25 (O33): locate the budget at which the interior-member result changes sign.

The constructive claim holds at beta = 1 and fails at beta = 0.5, and nothing in
the paper says where between them it turns over. That gap is the first thing a
reader will press on, because a claim that survives at one budget and dies at
half of it is only as useful as the width of the interval it is known on.

Design A carries an archived R14 trajectory record at beta = 0.75, so the
midpoint can be measured with the two protocols already in the paper and no new
machinery:

  stage 1  the R18 span sweep at beta = 0.75, run exactly as it was run at every
           other budget, which produces the interior member and its realized
           work
  stage 2  the R19 realized-work comparison at beta = 0.75, run verbatim, which
           scores that member against a constant degree matched on realized
           total quadratic work

Design B has no archived R14 record at beta = 0.75, so this is a single-design
result and is reported as one. Building the missing B record would be a new
budget campaign, not a control, and is out of scope here.

Both stages are resumable and both take a deadline. This runner converts a wall
clock stop time into per-stage deadlines rather than trusting a duration
estimate, so the queue stops when it is told to and not when a stage happens to
finish. Whatever is complete at that point is summarized and reported with its
completion fraction; the panel order is the archived Sobol index order, which
does not correlate with cost, so an early stop shortens the panel without
tilting it.

Usage:
    python rev25_crossover_queue.py --stop-at "2026-07-31 08:30"
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "python_codes"
METRICS = ROOT / "metrics"
LOG = CODE / "r25_queue.log"

BETA = 0.75
DESIGN = "A"
# leave the last slice of the window for the summaries and tables, which are
# cheap but must not be the thing that gets cut off
RESERVE_MIN = 20.0


def say(msg: str) -> None:
    line = f"[r25 {datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def minutes_left(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def run(cmd: list[str], budget_min: float) -> int:
    say(f"$ {' '.join(cmd[1:])}   (budget {budget_min:.0f} min)")
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    for stream in (proc.stdout, proc.stderr):
        if stream:
            with LOG.open("a", encoding="utf-8") as fh:
                fh.write(stream)
            tail = [ln for ln in stream.strip().splitlines() if ln.strip()]
            for ln in tail[-6:]:
                print("   " + ln, flush=True)
    say(f"exit {proc.returncode} after {(time.time()-t0)/60:.1f} min")
    return proc.returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True,
                    help='wall clock stop, e.g. "2026-07-31 08:30"')
    ap.add_argument("--workers", type=int, default=11)
    args = ap.parse_args()

    stop_at = datetime.strptime(args.stop_at, "%Y-%m-%d %H:%M")
    total = minutes_left(stop_at)
    if total <= RESERVE_MIN:
        say(f"stop time {stop_at:%Y-%m-%d %H:%M} is already here; nothing run")
        return 1

    say(f"window opens; stop at {stop_at:%Y-%m-%d %H:%M} "
        f"({total/60:.1f} h), design {DESIGN}, beta {BETA}")
    py = sys.executable
    worst = 0

    # ---- stage 1: the span sweep at the midpoint budget
    budget = max(0.0, minutes_left(stop_at) - RESERVE_MIN)
    # keep a share for stage 2, which cannot start until this one has a record
    stage1 = min(budget * 0.6, budget)
    rc = run([py, str(CODE / "rev18_span_sweep.py"), "run",
              "--design", DESIGN, "--beta", str(BETA),
              "--workers", str(args.workers),
              "--deadline-min", f"{stage1:.0f}"], stage1)
    worst = max(worst, rc)

    # rebuild the record from whatever landed, so a truncated stage still
    # leaves a readable record for stage 2
    run([py, str(CODE / "rev18_span_sweep.py"), "summarize",
         "--design", DESIGN, "--beta", str(BETA), "--from-disk"], 5)

    span = METRICS / f"r18_span_sweep_{DESIGN}_beta_{BETA:.2f}.json"
    if not span.exists():
        say(f"{span.name} was not produced; stage 2 cannot run")
        return 1

    # ---- stage 2: the realized-work comparison on that member
    budget = max(0.0, minutes_left(stop_at) - RESERVE_MIN)
    if budget <= 5:
        say("no time left for stage 2; stage 1 record is on disk and "
            "the queue is resumable")
        return worst
    rc = run([py, str(CODE / "rev19_equal_total_work.py"), "run",
              "--design", DESIGN, "--beta", str(BETA),
              "--workers", str(args.workers),
              "--deadline-min", f"{budget:.0f}"], budget)
    worst = max(worst, rc)

    run([py, str(CODE / "rev19_equal_total_work.py"), "summarize",
         "--design", DESIGN, "--beta", str(BETA)], 5)

    say(f"queue complete (worst exit {worst}); "
        f"{minutes_left(stop_at):.0f} min to spare")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
