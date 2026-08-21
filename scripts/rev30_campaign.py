"""R30 supervisor: run the four strata to a wall clock, cheapest first.

Order is by cost, not by interest. high_apolune holds no orbit below 80 km
perilune, so every one of its truths runs at degree 300 rather than 900; it
finishes in a fraction of the time the other three take and it exercises the
whole stratum pipeline end to end before the expensive populations commit to it.
The remaining three are equally expensive and run in the order they were frozen.

A stratum is only useful complete: a partial base cannot be calibrated and a
partial ladder cannot be quoted. So the supervisor refuses to start a stage it
cannot finish, and when the clock runs short it stops rather than leaving a
half-propagated population behind that nothing may be said about.

The declared budget beta = 1 is propagated on every stratum before any stratum
gets a second budget, which is what the R30 registration commits to.

Usage:
    python rev30_campaign.py --stop-at "2026-08-02 21:00" --workers 11
    python rev30_campaign.py --status
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"

import population_registry as registry

REG = "r30"          # rebound by main() when another registration is driven
LOG = HERE / "r30_campaign.log"
PROGRESS = METRICS / "r30_campaign_progress.json"

# cheapest first, then in frozen order
ORDER = ["high_apolune", "polar", "equatorial", "frozen_like",
         "low_perilune"]
PRIMARY_BETA = 1.00
SECOND_PASS_BETAS = [0.75, 0.50]

BASE_MIN_CHEAP = 75.0        # high_apolune: no 900-degree truths
BASE_MIN_FULL = 200.0        # prepass ~100 min + two-policy convergence ~70 min
# low_perilune: every truth at the highest degree. The archived attempt's own
# session_wall_ns gives 6 orbits in 108.9 min at eleven workers, so the 58 that
# remain are about 17.5 h there and 21-29 h at four. The old 300 min would wave
# through a base that could not possibly finish; a base that stops short is not
# wasted, because it resumes from its sidecars, but the stratum is then still
# uncalibrated and nothing downstream may start.
BASE_MIN_DEEP = 300.0
OP_MIN = 60.0
CAL_MIN = 40.0
LADDER_MIN = 150.0
RESERVE_MIN = 25.0


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def remaining_min(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def record(name: str, rc: int, minutes: float) -> None:
    p = (json.loads(PROGRESS.read_text(encoding="utf-8"))
         if PROGRESS.exists()
         else {"schema": "r30_campaign_progress_v1", "stages": []})
    p["stages"].append({"stage": name, "rc": rc, "minutes": round(minutes, 1),
                        "finished_local": datetime.now().isoformat(
                            timespec="seconds")})
    p["updated_local"] = datetime.now().isoformat(timespec="seconds")
    PROGRESS.write_text(json.dumps(p, indent=2), encoding="utf-8")


# One BLAS thread per worker, and no console window. Both were learned on
# 15 August: a pool whose workers each spin up a full BLAS thread pool dies
# with BrokenProcessPool within a minute on this host, at one worker as
# readily as at four, while the same trajectory runs clean in-process; and a
# child started DETACHED has no console for its own workers to inherit, so
# Windows draws a new console window for each of them.
CREATE_NO_WINDOW = 0x08000000
CREATE_NEW_PROCESS_GROUP = 0x00000200
CHILD_ENV = {"PYTHONUNBUFFERED": "1", "OMP_NUM_THREADS": "1",
             "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
             "NUMEXPR_NUM_THREADS": "1"}


# A stratum's op, calibration and three ladders write about half a gigabyte
# between them. The floor is set well above that because the volume this runs
# on has been filled twice mid-campaign, and a stage that dies on ENOSPC costs
# its whole elapsed time: frozen_like's beta = 0.50 ladder wrote a zero-byte
# summary on 8 August and the supervisor then skipped the retry.
DISK_FLOOR_GB = 2.0


def free_gb() -> float:
    return shutil.disk_usage(str(ROOT)).free / 2**30


def run_stage(name: str, cmd: list[str]) -> tuple[int, float]:
    free = free_gb()
    if free < DISK_FLOOR_GB:
        log(f"REFUSE {name}: {free:.1f} GB free, floor is {DISK_FLOOR_GB} GB")
        record(name, -1, 0.0)
        return -1, 0.0
    log(f"START {name}: {' '.join(cmd[1:])} ({free:.1f} GB free)")
    t0 = time.time()
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL,
                             env=dict(os.environ, **CHILD_ENV),
                             creationflags=CREATE_NO_WINDOW
                             | CREATE_NEW_PROCESS_GROUP)
    dt = (time.time() - t0) / 60.0
    log(f"END   {name}: rc={rc} in {dt:.1f} min")
    record(name, rc, dt)
    return rc, dt


def key_of(stratum: str) -> str:
    return registry.spec(REG, stratum)["design_key"]


def base_complete(stratum: str) -> bool:
    p = METRICS / f"{REG}_{stratum}_convergence.json"
    return (p.exists()
            and bool(json.loads(p.read_text(encoding="utf-8")).get("complete")))


def op_done(stratum: str) -> bool:
    return (METRICS / f"{REG}_{stratum}_operating_point.json").exists()


def cal_done(stratum: str) -> bool:
    return (METRICS / f"{REG}_budget_pareto_{stratum}.json").exists()


def ladder_done(stratum: str, beta: float) -> bool:
    """A ladder counts as done only if its record is readable and scored.

    Existence is not enough. On 2026-08-08 the frozen_like beta = 0.50 stage
    ran out of disk while writing, leaving a zero-byte file; an existence test
    read that as success, so the supervisor skipped the retry and returned with
    348 minutes of its window unused, and --status reported the ladder present.
    Requiring the summary block to parse makes a truncated or empty write count
    as what it is.
    """
    p = (METRICS / f"r19_equal_total_work_{key_of(stratum)}"
         f"_beta_{beta:.2f}.json")
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def status() -> int:
    for s in ORDER:
        done = [f"{b:.2f}" for b in [PRIMARY_BETA] + SECOND_PASS_BETAS
                if ladder_done(s, b)]
        print(f"{s:<14} base={'y' if base_complete(s) else '-'} "
              f"op={'y' if op_done(s) else '-'} "
              f"cal={'y' if cal_done(s) else '-'} "
              f"ladders={','.join(done) or '-'}")
    if PROGRESS.exists():
        for st in json.loads(PROGRESS.read_text(encoding="utf-8"))["stages"]:
            print(f"  {st['finished_local']}  {st['stage']:<28} "
                  f"rc={st['rc']} {st['minutes']:>6.1f} min")
    return 0


def prepare(stratum: str, stop_at: datetime, workers: int, py: str) -> bool:
    """Base, operating point and calibration. True when the stratum is ready
    for a ladder."""
    if not base_complete(stratum):
        need = {"high_apolune": BASE_MIN_CHEAP,
                "low_perilune": BASE_MIN_DEEP}.get(stratum, BASE_MIN_FULL)
        if remaining_min(stop_at) < need + OP_MIN + CAL_MIN + RESERVE_MIN:
            log(f"{stratum}: not enough time for a base that could then be "
                f"calibrated and run; skipping the stratum entirely")
            return False
        rc, _ = run_stage(f"{stratum}/base",
                          [py, "rev30_stratum_base.py", "--registry", REG,
                           "--stratum", stratum,
                           "run", "--workers", str(workers),
                           "--deadline", stop_at.astimezone().isoformat()])
        if not base_complete(stratum):
            log(f"{stratum}: base incomplete; nothing is quoted from it")
            return False
    if not op_done(stratum):
        if remaining_min(stop_at) < OP_MIN + CAL_MIN + RESERVE_MIN:
            log(f"{stratum}: no time for the operating point")
            return False
        rc, _ = run_stage(f"{stratum}/op",
                          [py, "rev30_stratum_ops.py", "--registry", REG,
                           "--stratum", stratum,
                           "op", "--workers", str(workers)])
        if rc != 0:
            log(f"{stratum}: operating point failed")
            return False
    if not cal_done(stratum):
        if remaining_min(stop_at) < CAL_MIN + RESERVE_MIN:
            log(f"{stratum}: no time for the calibration")
            return False
        rc, _ = run_stage(f"{stratum}/cal",
                          [py, "rev30_stratum_ops.py", "--registry", REG,
                           "--stratum", stratum,
                           "cal", "--workers", str(workers)])
        if rc != 0:
            log(f"{stratum}: calibration failed or was refused")
            return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--strata", default="", help="comma-separated subset")
    ap.add_argument("--registry", default="r30",
                    help="which registration declares the populations to run")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    global REG, ORDER, LOG, PROGRESS
    REG = a.registry
    if REG != "r30":
        ORDER = sorted(registry.populations(REG))
        LOG = HERE / f"{REG}_campaign.log"
        PROGRESS = METRICS / f"{REG}_campaign_progress.json"
    if a.status:
        return status()
    if not a.stop_at:
        ap.error("--stop-at is required")

    order = [s.strip() for s in a.strata.split(",") if s.strip()] or ORDER
    unknown = [s for s in order if s not in ORDER]
    if unknown:
        ap.error(f"unknown strata {unknown}; known: {ORDER}")

    stop_at = datetime.fromisoformat(a.stop_at)
    py = sys.executable
    log(f"=== R30 strata campaign: stop at "
        f"{stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h), order={order} ===")

    ladder_est = LADDER_MIN
    measured: list[float] = []
    for stratum in order:
        if ladder_done(stratum, PRIMARY_BETA):
            log(f"{stratum}: beta {PRIMARY_BETA:.2f} already on disk")
            continue
        if not prepare(stratum, stop_at, a.workers, py):
            continue
        left = remaining_min(stop_at) - RESERVE_MIN
        if left < ladder_est:
            log(f"{stratum}: {left:.0f} min left, a ladder needs "
                f"{ladder_est:.0f}; stopping rather than half-running one")
            break
        rc, dt = run_stage(f"{stratum}/ladder_{PRIMARY_BETA:.2f}",
                           [py, "rev30_stratum_ops.py", "--registry", REG,
                            "--stratum", stratum,
                            "ladder", "--beta", f"{PRIMARY_BETA:.2f}",
                            "--workers", str(a.workers),
                            "--deadline-min", f"{left - 5.0:.0f}"])
        if rc == 0 and ladder_done(stratum, PRIMARY_BETA):
            measured.append(dt)
            ladder_est = max(max(measured) * 1.15, 60.0)
            log(f"ladder estimate now {ladder_est:.0f} min")
        else:
            log(f"{stratum}: ladder rc={rc}, incomplete")

    # second pass: only once every stratum has the declared budget
    for beta in SECOND_PASS_BETAS:
        if not all(ladder_done(s, PRIMARY_BETA) for s in order):
            log("second pass skipped: not every stratum carries the declared "
                "budget yet")
            break
        for stratum in order:
            if ladder_done(stratum, beta):
                continue
            left = remaining_min(stop_at) - RESERVE_MIN
            if left < ladder_est:
                log(f"second pass: {left:.0f} min left, stopping")
                return 0
            run_stage(f"{stratum}/ladder_{beta:.2f}",
                      [py, "rev30_stratum_ops.py", "--registry", REG,
                       "--stratum", stratum,
                       "ladder", "--beta", f"{beta:.2f}",
                       "--workers", str(a.workers),
                       "--deadline-min", f"{left - 5.0:.0f}"])

    log(f"=== strata campaign done, {remaining_min(stop_at):.0f} min before "
        f"the stop ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
