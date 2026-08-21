"""R51 supervisor: the cap-lifted control on block A of the span ladder.

The shape of rev38_campaign, which ran the same control on the operational
elliptical population, with the priors moved to what a degree-600 reference
actually cost there: a base of 288 minutes against 72 at degree 300, and ladders
of 43, 40 and 36 minutes, because lifting the ceiling redistributes degrees at a
fixed mean square rather than buying more of them.

The prepass is not re-run and cannot be: the rows are the R50 rows with one
field changed, written by rev51_uncapped_freeze.py, so that n_critical and
n_work are the parent's and the comparator is the same comparator. The operating
point runs before the base, because it reads those rows and nothing else, and
five hours of propagation should not be what discovers a misconfiguration. The
reindex step is inside the chain rather than left to be remembered.

Usage:
    python rev51_campaign.py --stop-at "2026-08-13 12:00" --workers 11
    python rev51_campaign.py --status
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"

REG = "r51"
STRATUM = "span_ladder_a_uncapped"
KEY = "RS1U"
BETAS = [1.00, 0.75, 0.62, 0.50]

LOG = HERE / "r51_campaign.log"
PROGRESS = METRICS / "r51_campaign_progress.json"

# measured on the degree-600 control of the operational elliptical population
OP_MIN = 15.0
BASE_MIN = 320.0
CAL_MIN = 15.0
LADDER_PRIOR_MIN = 60.0
RESERVE_MIN = 20.0

CONV = METRICS / f"{REG}_{STRATUM}_convergence.json"
OP = METRICS / f"{REG}_{STRATUM}_operating_point.json"
CAL = METRICS / f"{REG}_budget_pareto_{STRATUM}.json"
ROWS = METRICS / f"{REG}_{STRATUM}_rows.json"


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
         else {"schema": "r51_campaign_progress_v1", "stages": []})
    p["stages"].append({"stage": name, "rc": rc, "minutes": round(minutes, 1),
                        "finished_local": datetime.now().isoformat(
                            timespec="seconds")})
    p["updated_local"] = datetime.now().isoformat(timespec="seconds")
    PROGRESS.write_text(json.dumps(p, indent=2), encoding="utf-8")


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
    return (CONV.exists()
            and bool(json.loads(CONV.read_text(encoding="utf-8"))
                     .get("complete")))


def base_propagated() -> bool:
    raw = METRICS / "r11_raw" / f"stratum_{STRATUM}_convergence"
    if not raw.exists():
        return False
    want = {f"{p}_{l}.npz" for p in ("truth", "fixed_critical")
            for l in ("tight", "tighter")}
    dirs = [d for d in raw.iterdir() if d.is_dir()]
    if len(dirs) != 64:
        return False
    return all(want <= {f.name for f in d.iterdir()} for d in dirs)


def ladder_done(beta: float) -> bool:
    p = METRICS / f"r19_equal_total_work_{KEY}_beta_{beta:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def status() -> int:
    print(f"{STRATUM} ({KEY}) rows={'y' if ROWS.exists() else '-'} "
          f"op={'y' if OP.exists() else '-'} "
          f"propagated={'y' if base_propagated() else '-'} "
          f"indexed={'y' if base_complete() else '-'} "
          f"cal={'y' if CAL.exists() else '-'} "
          f"ladders={','.join(f'{b:.2f}' for b in BETAS if ladder_done(b)) or '-'}")
    if PROGRESS.exists():
        for st in json.loads(PROGRESS.read_text(encoding="utf-8"))["stages"]:
            print(f"  {st['finished_local']}  {st['stage']:<24} "
                  f"rc={st['rc']} {st['minutes']:>6.1f} min")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.status:
        return status()
    if not a.stop_at:
        ap.error("--stop-at is required")
    if not ROWS.exists():
        log("rows are missing; run rev51_uncapped_freeze.py first")
        return 1

    stop_at = datetime.fromisoformat(a.stop_at)
    py = sys.executable
    log(f"=== R51 cap-lifted span-ladder control: stop at "
        f"{stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h) ===")

    if not OP.exists():
        rc, _ = run_stage("op", [py, "rev30_stratum_ops.py", "--registry", REG,
                                 "--stratum", STRATUM, "op",
                                 "--workers", str(a.workers)])
        if rc != 0:
            log("operating point failed")
            return 1

    if not base_complete():
        if not base_propagated():
            base_stop = stop_at - timedelta(
                minutes=CAL_MIN + LADDER_PRIOR_MIN + RESERVE_MIN)
            if remaining_min(base_stop) < 30.0:
                log("no room for a base that could then be calibrated and "
                    "laddered")
                return 1
            log(f"base deadline {base_stop.isoformat(timespec='minutes')}, "
                f"leaving {CAL_MIN + LADDER_PRIOR_MIN:.0f} min for a "
                f"calibration and the declared budget")
            run_stage("base",
                      [py, "rev30_stratum_base.py", "--registry", REG,
                       "--stratum", STRATUM, "run", "--workers",
                       str(a.workers), "--deadline",
                       base_stop.astimezone().isoformat()])
        if not base_propagated():
            log("base did not finish propagating; nothing is quoted from it. "
                "The sidecars on disk are valid and a later run resumes.")
            return 1
        rc, _ = run_stage("reindex",
                          [py, "rev32_reindex_two_policy_base.py",
                           "--registry", REG, "--stratum", STRATUM])
        if rc != 0 or not base_complete():
            log("reindex failed; the base is propagated but not indexed")
            return 1

    if not CAL.exists():
        if remaining_min(stop_at) < CAL_MIN + RESERVE_MIN:
            log("no time for the calibration")
            return 1
        rc, _ = run_stage("cal", [py, "rev30_stratum_ops.py", "--registry", REG,
                                  "--stratum", STRATUM, "cal",
                                  "--workers", str(a.workers)])
        if rc != 0:
            log("calibration failed or was refused")
            return 1

    est = LADDER_PRIOR_MIN
    measured: list[float] = []
    for beta in BETAS:
        if ladder_done(beta):
            log(f"beta {beta:.2f} already on disk")
            continue
        left = remaining_min(stop_at) - RESERVE_MIN
        if left < est:
            log(f"{left:.0f} min left, a ladder needs {est:.0f}; stopping "
                f"rather than half-running one")
            break
        rc, dt = run_stage(f"ladder_{beta:.2f}",
                           [py, "rev30_stratum_ops.py", "--registry", REG,
                            "--stratum", STRATUM, "ladder", "--beta",
                            f"{beta:.2f}", "--workers", str(a.workers),
                            "--deadline-min", f"{left - 5.0:.0f}"])
        if rc == 0 and ladder_done(beta):
            measured.append(dt)
            est = max(max(measured) * 1.15, 40.0)
            log(f"ladder estimate now {est:.0f} min")
        else:
            log(f"beta {beta:.2f}: rc={rc}, incomplete")
            break

    log("=== R51 done ===")
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
