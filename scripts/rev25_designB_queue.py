"""R25 stage B: extend the beta = 0.75 crossover to design B.

Design B has no archived R14 trajectory record at beta = 0.75, which is why the
first queue is design A only. Building that record is a third stage rather than
a control, so it runs after design A has finished and only with whatever window
is left. All 128 design-B truth arrays are already on disk, so R14 here
propagates the two budget policies and reuses every truth.

    stage B1  rev14_budget_trajectory  the budget campaign at beta = 0.75
    stage B2  rev18_span_sweep         the interpolation family on that record
    stage B3  rev19_equal_total_work   the realized-work comparison

A truncated stage is not a wasted one: each driver writes its record with the
orbits it finished, the panel is walked in archived Sobol index order, and that
order does not correlate with cost, so an early stop shortens the panel without
tilting it. The completion fraction travels with the result.

Two things this runner is careful about. It waits for the design-A queue to
finish rather than starting alongside it, because eleven workers are already
busy and two campaigns would only slow each other. And it hands R14 an
offset-aware deadline: that driver parses a naive timestamp as UTC, so a bare
"08:30" would be read as 08:30 UTC and run three hours past the intended stop.

Usage:
    python rev25_designB_queue.py --stop-at "2026-07-31 08:30"
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "python_codes"
METRICS = ROOT / "metrics"
LOG = CODE / "r25_queue.log"
A_MARKER = "queue complete"

BETA = 0.75
DESIGN = "B"
RESERVE_MIN = 20.0
# if design A has not finished by this point there is no useful window left for
# a three-stage chain, so stage B is skipped rather than started half-cocked
MIN_USEFUL_MIN = 45.0


def say(msg: str) -> None:
    line = f"[r25B {datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def minutes_left(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def wait_for_design_A(stop_at: datetime) -> bool:
    """Block until the design-A queue reports completion."""
    say("waiting for the design-A queue to finish before starting")
    while True:
        try:
            if A_MARKER in LOG.read_text(encoding="utf-8", errors="replace"):
                say("design A reported complete")
                return True
        except OSError:
            pass
        left = minutes_left(stop_at) - RESERVE_MIN
        if left <= MIN_USEFUL_MIN:
            say(f"only {left:.0f} min left and design A is still running; "
                f"stage B skipped rather than started on a stub")
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=11)
    args = ap.parse_args()

    stop_at = datetime.strptime(args.stop_at, "%Y-%m-%d %H:%M")
    if not wait_for_design_A(stop_at):
        return 0

    py = sys.executable
    worst = 0
    say(f"stage B opens with {minutes_left(stop_at):.0f} min to the stop")

    # ---- B1: the budget campaign. R14 wants an ISO deadline and reads a naive
    # one as UTC, so the local stop is converted explicitly.
    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= MIN_USEFUL_MIN:
        say("no useful window for B1; stopping")
        return 0
    b1_stop = datetime.now() + timedelta(minutes=min(left * 0.55, left))
    worst = max(worst, run(
        [py, str(CODE / "rev14_budget_trajectory.py"), "run",
         "--design", DESIGN, "--beta", str(BETA),
         "--workers", str(args.workers),
         "--deadline", b1_stop.astimezone().isoformat()], "B1 budget"))

    record = METRICS / f"r14_trajectory_{DESIGN}_beta_{BETA:.2f}.json"
    if not record.exists():
        say(f"{record.name} not produced; B2 and B3 cannot run")
        return 1

    # ---- B2: the span sweep on that record
    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 10:
        say("window closed after B1; the budget record is on disk and the "
            "chain is resumable")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev18_span_sweep.py"), "run",
         "--design", DESIGN, "--beta", str(BETA),
         "--workers", str(args.workers),
         "--deadline-min", f"{min(left * 0.6, left):.0f}"], "B2 span"))
    run([py, str(CODE / "rev18_span_sweep.py"), "summarize",
         "--design", DESIGN, "--beta", str(BETA), "--from-disk"], "B2 summary")

    # ---- B3: the realized-work comparison
    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 10:
        say("window closed after B2; span record is on disk, chain resumable")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev19_equal_total_work.py"), "run",
         "--design", DESIGN, "--beta", str(BETA),
         "--workers", str(args.workers),
         "--deadline-min", f"{left:.0f}"], "B3 equal-work"))
    run([py, str(CODE / "rev19_equal_total_work.py"), "summarize",
         "--design", DESIGN, "--beta", str(BETA)], "B3 summary")

    say(f"stage B complete (worst exit {worst}); "
        f"{minutes_left(stop_at):.0f} min to spare")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
