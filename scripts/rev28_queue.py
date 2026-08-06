"""Drive the beta = 0.62 bisection chain on both designs.

Six stages, in order, each a subprocess so that a stage that dies takes only
itself down:

    A: R14 trajectory -> R18 span sweep -> R19 realized-work comparison
    B: R14 trajectory -> R18 span sweep -> R19 realized-work comparison

The order is A first because the R25 amendment says so, and because if the
window closes early a complete design A is worth more than two half designs.
R18 depends on R14's record for its tolerances and endpoints; R19 depends on
R18's k = 0.50 and k = 0.00 entries for realized work and the constant degree.
Neither dependency is checked by those scripts, so the queue checks it here:
a stage whose input is missing or truncated does not start.

Wall clock, not duration estimates. --stop-at is a local timestamp and every
stage gets a deadline derived from what is actually left, so a slow stage
shortens the next one instead of overrunning the wall.

The R14 clock trap, which is real and cost this campaign a run once already:
rev14_budget_trajectory.parse_deadline reads a naive ISO string as UTC. The
machine is UTC+3, so a bare "23:00" would be taken as 02:00 local the next day
and overshoot by three hours. This queue always passes an offset-aware string.
R18 and R19 take --deadline-min instead and are not affected.

Costs are measured, not guessed, from the R25 campaign on the same machine:
R14 ~47 min, R18 ~60-70 min, R19 ~18-21 min, so a full chain is ~140 min. A
stage is skipped rather than started if what remains cannot finish it, because
starting a chain that cannot complete is the one failure this queue exists to
prevent.

Usage:
    python rev28_queue.py --stop-at "2026-08-01 06:00" --workers 11
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
BETA = 0.62
TAG = f"beta_{BETA:.2f}"

# measured on this machine during R25; the guard is deliberately above them
COST_MIN = {"r14": 60.0, "r18": 85.0, "r19": 30.0}
RESERVE_MIN = 5.0          # never hand a stage the last of the window


def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def minutes_left(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def r14_out(design: str) -> Path:
    return METRICS / f"r14_trajectory_{design}_{TAG}.json"


def r18_out(design: str) -> Path:
    return METRICS / f"r18_span_sweep_{design}_{TAG}.json"


def r19_out(design: str) -> Path:
    return METRICS / f"r19_equal_total_work_{design}_{TAG}.json"


def rows_in(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return len(json.loads(path.read_text(encoding="utf-8")).get("rows", []))
    except Exception:
        return 0


def call(argv: list[str], cwd: Path = HERE) -> int:
    log(f"$ {' '.join(argv[1:])}")
    t0 = time.time()
    rc = subprocess.call(argv, cwd=str(cwd))
    log(f"  -> rc={rc} in {(time.time() - t0) / 60:.1f} min")
    return rc


def stage_r14(design: str, workers: int, stop_at: datetime) -> bool:
    out = r14_out(design)
    if rows_in(out) >= 64:
        log(f"R14 {design}: {out.name} already complete, skipping")
        return True
    budget = minutes_left(stop_at) - RESERVE_MIN
    if budget < COST_MIN["r14"]:
        log(f"R14 {design}: {budget:.0f} min left, needs {COST_MIN['r14']:.0f}; "
            f"not starting a chain that cannot finish")
        return False
    # offset-aware on purpose: R14 reads a naive string as UTC (see module docstring)
    deadline = datetime.fromtimestamp(
        time.time() + budget * 60.0).astimezone().replace(microsecond=0)
    return call([sys.executable, "rev28_r14_beta062.py", "--design", design,
                 "--workers", str(workers),
                 "--deadline", deadline.isoformat()]) == 0


def stage_r18(design: str, workers: int, stop_at: datetime) -> bool:
    if rows_in(r14_out(design)) < 64:
        log(f"R18 {design}: R14 record incomplete, refusing to start")
        return False
    budget = minutes_left(stop_at) - RESERVE_MIN
    if budget < COST_MIN["r18"]:
        log(f"R18 {design}: {budget:.0f} min left, needs {COST_MIN['r18']:.0f}")
        return False
    rc = call([sys.executable, "rev18_span_sweep.py", "run", "--design", design,
               "--beta", str(BETA), "--workers", str(workers),
               "--deadline-min", f"{budget:.0f}"])
    if rc != 0:
        return False
    return call([sys.executable, "rev18_span_sweep.py", "summarize",
                 "--design", design, "--beta", str(BETA)]) == 0


def stage_r19(design: str, workers: int, stop_at: datetime) -> bool:
    span = r18_out(design)
    if not span.exists():
        log(f"R19 {design}: {span.name} missing, refusing to start")
        return False
    try:
        entries = json.loads(span.read_text(encoding="utf-8"))["rows"]
    except Exception as exc:
        log(f"R19 {design}: cannot read {span.name} ({exc})")
        return False
    usable = sum(1 for r in entries
                 if "0.50" in r.get("entries", {})
                 and "0.00" in r.get("entries", {}))
    if usable == 0:
        log(f"R19 {design}: no orbit carries both k=0.50 and k=0.00 entries")
        return False
    log(f"R19 {design}: {usable}/{len(entries)} orbits carry both endpoints")
    budget = minutes_left(stop_at) - RESERVE_MIN
    if budget < COST_MIN["r19"]:
        log(f"R19 {design}: {budget:.0f} min left, needs {COST_MIN['r19']:.0f}")
        return False
    rc = call([sys.executable, "rev19_equal_total_work.py", "run",
               "--design", design, "--beta", str(BETA),
               "--workers", str(workers), "--deadline-min", f"{budget:.0f}"])
    if rc != 0:
        return False
    return call([sys.executable, "rev19_equal_total_work.py", "summarize",
                 "--design", design, "--beta", str(BETA)]) == 0


def report(design: str) -> None:
    path = r19_out(design)
    if not path.exists():
        log(f"design {design}: no R19 record")
        return
    d = json.loads(path.read_text(encoding="utf-8"))
    s = d.get("summary", d)
    log(f"design {design} beta {BETA}: {json.dumps(s)[:400]}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True,
                    help='local wall clock, e.g. "2026-08-01 06:00"')
    ap.add_argument("--workers", type=int, default=11)
    a = ap.parse_args()
    stop_at = datetime.fromisoformat(a.stop_at)

    cal = METRICS / "r28_budget_pareto_beta_0.62.json"
    if not cal.exists():
        log(f"abort: {cal.name} missing")
        return 2
    log(f"beta = {BETA} bisection chain, stop at {stop_at:%Y-%m-%d %H:%M} "
        f"({minutes_left(stop_at):.0f} min)")
    log(f"declared post hoc by r28_calibration_amendment.json; "
        f"outcomes E/F/G fixed in r25_preregistration_amendment.json")

    for design in ("A", "B"):
        log(f"===== design {design} =====")
        for name, fn in (("R14", stage_r14), ("R18", stage_r18),
                         ("R19", stage_r19)):
            if not fn(design, a.workers, stop_at):
                log(f"design {design}: stopped at {name}")
                break
        else:
            log(f"design {design}: chain complete")
            report(design)

    log("queue complete")
    for design in ("A", "B"):
        for path in (r14_out(design), r18_out(design), r19_out(design)):
            log(f"  {path.name}: {'present' if path.exists() else 'absent'}"
                f"{f' ({rows_in(path)} rows)' if path.exists() else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
