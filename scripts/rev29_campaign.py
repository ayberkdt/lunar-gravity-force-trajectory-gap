"""R29 supervisor: run the design-C campaign to a wall clock, in order.

The stages depend on each other in one direction only -- base, operating point,
calibration, then one ladder per budget -- so the supervisor is a queue rather
than a scheduler. What it has to get right is the same thing the R25 supervisor
had to get right: never start a stage it cannot finish, and never idle when
there is a stage it can.

Each stage runs as a subprocess. A ladder that dies takes its budget down and
nothing else, model caches do not leak from one budget into the next, and the
log keeps the stage boundaries visible.

Budgets are propagated in the pre-registered order: the two that bracket the
crossing on designs A and B first, so that a campaign that runs out of time has
still answered the question it was opened for.

Chain cost is measured rather than assumed. The first completed ladder sets the
estimate for the rest; before that the estimate comes from the design A and B
chains of R25 (140 min), rounded up, because starting a chain that cannot
finish is the one mistake this file exists to prevent.

Usage:
    python rev29_campaign.py --stop-at "2026-08-02 23:00" --workers 11
    python rev29_campaign.py --status
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

LOG = HERE / "r29_campaign.log"
PROGRESS = METRICS / "r29_campaign_progress.json"

BUDGET_ORDER = [0.75, 0.50, 1.00, 0.62, 1.25, 1.50, 2.00]

# minutes; the floor under "can this still finish", not a prediction
CHAIN_MIN_DEFAULT = 165.0
BASE_MIN = 45.0
OPPOINT_MIN = 25.0
PARETO_MIN = 60.0
RESERVE_MIN = 20.0          # left for summaries and tables at the end


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def remaining_min(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def read_progress() -> dict:
    if PROGRESS.exists():
        return json.loads(PROGRESS.read_text(encoding="utf-8"))
    return {"schema": "r29_campaign_progress_v1", "stages": []}


def write_progress(p: dict) -> None:
    PROGRESS.parent.mkdir(parents=True, exist_ok=True)
    PROGRESS.write_text(json.dumps(p, indent=2), encoding="utf-8")


def record(name: str, rc: int, minutes: float, note: str = "") -> None:
    p = read_progress()
    p["stages"].append({"stage": name, "rc": rc,
                        "minutes": round(minutes, 1),
                        "finished_local": datetime.now().isoformat(
                            timespec="seconds"),
                        "note": note})
    p["updated_local"] = datetime.now().isoformat(timespec="seconds")
    write_progress(p)


def run_stage(name: str, cmd: list[str]) -> tuple[int, float]:
    log(f"START {name}: {' '.join(cmd[1:])}")
    t0 = time.time()
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT)
    dt = (time.time() - t0) / 60.0
    log(f"END   {name}: rc={rc} in {dt:.1f} min")
    record(name, rc, dt)
    return rc, dt


def base_complete() -> bool:
    p = METRICS / "r26_designC_convergence.json"
    if not p.exists():
        return False
    return bool(json.loads(p.read_text(encoding="utf-8")).get("complete"))


def ladder_done(beta: float) -> bool:
    return (METRICS / f"r19_equal_total_work_C_beta_{beta:.2f}.json").exists()


def status() -> int:
    print(f"base complete: {base_complete()}")
    for f, label in ((METRICS / "r29_designC_operating_point.json", "operating point"),
                     (METRICS / "r29_budget_pareto_designC.json", "calibration")):
        print(f"{label}: {'yes' if f.exists() else 'no'}")
    for b in BUDGET_ORDER:
        print(f"  beta {b:.2f}: {'done' if ladder_done(b) else '-'}")
    if PROGRESS.exists():
        p = read_progress()
        for s in p["stages"]:
            print(f"  {s['finished_local']}  {s['stage']:<18} rc={s['rc']} "
                  f"{s['minutes']:>6.1f} min")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", help="local wall clock, e.g. '2026-08-02 23:00'")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--budgets", default="",
                    help="comma-separated subset of the pre-registered order; "
                         "empty means all of it")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.status:
        return status()
    order = ([float(x) for x in a.budgets.split(",")] if a.budgets
             else list(BUDGET_ORDER))
    if any(b not in BUDGET_ORDER for b in order):
        ap.error(f"budgets outside the pre-registered order {BUDGET_ORDER}")
    if not a.stop_at:
        ap.error("--stop-at is required")

    stop_at = datetime.fromisoformat(a.stop_at)
    py = sys.executable
    log(f"=== R29 campaign: stop at {stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h), workers={a.workers} ===")

    # ---- stage 1: the base
    if not base_complete():
        if remaining_min(stop_at) < BASE_MIN + RESERVE_MIN:
            log("no time for the design-C base; nothing else can run")
            return 4
        deadline = stop_at.astimezone().isoformat()
        rc, _ = run_stage("base", [py, "rev26_designC_base.py", "run",
                                   "--workers", str(a.workers),
                                   "--deadline", deadline])
        if not base_complete():
            log("base incomplete after its run; the campaign stops here "
                "(no budget is quoted from a partial base)")
            return 1

    # ---- stage 2: the accuracy-target operating point
    op = METRICS / "r29_designC_operating_point.json"
    if not op.exists():
        if remaining_min(stop_at) < OPPOINT_MIN + RESERVE_MIN:
            log("no time for the operating point")
            return 4
        rc, _ = run_stage("operating_point",
                          [py, "rev29_designC_operating_point.py",
                           "--workers", str(a.workers)])
        if rc != 0 or not op.exists():
            log("operating point failed; the calibration cannot run")
            return 1

    # ---- stage 3: the calibration
    cal = METRICS / "r29_budget_pareto_designC.json"
    if not cal.exists():
        if remaining_min(stop_at) < PARETO_MIN + RESERVE_MIN:
            log("no time for the calibration")
            return 4
        rc, _ = run_stage("calibration", [py, "rev29_designC_pareto.py",
                                          "--workers", str(a.workers)])
        if rc != 0 or not cal.exists():
            log("calibration failed or was refused; no ladder can run")
            return 1

    # ---- stage 4: the ladders
    chain_est = CHAIN_MIN_DEFAULT
    measured: list[float] = []
    for beta in order:
        if ladder_done(beta):
            log(f"beta {beta:.2f}: already on disk, skipping")
            continue
        left = remaining_min(stop_at) - RESERVE_MIN
        if left < chain_est:
            log(f"beta {beta:.2f}: {left:.0f} min left, chain needs "
                f"{chain_est:.0f}; stopping rather than starting a chain that "
                f"cannot finish")
            break
        rc, dt = run_stage(f"ladder_beta_{beta:.2f}",
                           [py, "rev29_designC_ladder.py",
                            "--beta", f"{beta:.2f}",
                            "--workers", str(a.workers),
                            "--deadline-min", f"{left - 5.0:.0f}"])
        if rc == 0 and ladder_done(beta):
            measured.append(dt)
            chain_est = max(max(measured) * 1.15, 60.0)
            log(f"chain estimate now {chain_est:.0f} min "
                f"(measured {', '.join(f'{m:.0f}' for m in measured)})")
        else:
            log(f"beta {beta:.2f}: rc={rc}, ladder incomplete; moving on")

    log(f"=== campaign done, {remaining_min(stop_at):.0f} min before the stop ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
