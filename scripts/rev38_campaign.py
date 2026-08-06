"""R38 supervisor: base, reindex, operating point, calibration, ladders.

Same shape as rev30_campaign, with two differences that are not cosmetic.

First, the reindex step is in the chain rather than left to be remembered. A
two-policy base is propagated correctly and then indexed incorrectly, because
the pinned convergence driver assembles its summary from the four policies this
scope switches off; rev32_reindex_two_policy_base rebuilds the index from the
sidecars. R30 lost thirty-four hours of ladder time to that gap being outside
the supervisor. It is inside this one.

Second, the ladder time estimate is measured rather than assumed. R30 carried a
150-minute prior because it did not know what a stratum ladder cost. This
population's parent measured 46, 44 and 35 minutes for the three budgets on the
same sixty-four orbits, and the ladder here propagates the same policies at the
same realized work -- the calibration holds mean N squared fixed, so lifting the
ceiling redistributes degrees without buying more of them. A 150-minute prior
would refuse a ladder there is time for.

Every stage is skipped when its output already exists, and the base is
checkpointed per trajectory, so an interrupted run resumes where it stopped.
Stages run in the registered budget order -- beta = 1.00 first, because it is
the budget that carries the published number this control exists to check --
and the supervisor stops rather than starting a stage it cannot finish.

Usage:
    python rev38_campaign.py --stop-at "2026-08-05 08:30" --workers 11
    python rev38_campaign.py --status
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

REG = "r38"
STRATUM = "operational_elliptical_uncapped"
KEY = "OEU"

LOG = HERE / "r38_campaign.log"
PROGRESS = METRICS / "r38_campaign_progress.json"

BETAS = [1.00, 0.75, 0.50]

# measured on the capped parent, same orbits, same policies, same realized work
LADDER_PRIOR_MIN = 75.0
OP_MIN = 35.0
CAL_MIN = 25.0
RESERVE_MIN = 15.0


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


CONV = METRICS / f"{REG}_{STRATUM}_convergence.json"
OP = METRICS / f"{REG}_{STRATUM}_operating_point.json"
CAL = METRICS / f"{REG}_budget_pareto_{STRATUM}.json"


def base_complete() -> bool:
    return (CONV.exists()
            and bool(json.loads(CONV.read_text(encoding="utf-8")).get("complete")))


def base_propagated() -> bool:
    """Every trajectory on disk, whether or not the index says so."""
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
    return (METRICS / f"r19_equal_total_work_{KEY}_beta_{beta:.2f}.json").exists()


def status() -> int:
    print(f"{STRATUM}: base={'y' if base_complete() else '-'} "
          f"propagated={'y' if base_propagated() else '-'} "
          f"op={'y' if OP.exists() else '-'} cal={'y' if CAL.exists() else '-'} "
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

    stop_at = datetime.fromisoformat(a.stop_at)
    py = sys.executable
    log(f"=== R38 cap-lifted control: stop at "
        f"{stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h) ===")

    # The operating point reads the prepass rows and nothing else -- no truth
    # trajectory, no index -- so it runs before the base rather than after it.
    # Four hours of propagation should not be what discovers that a later stage
    # is misconfigured.
    if not OP.exists():
        rc, _ = run_stage("op", [py, "rev30_stratum_ops.py", "--registry", REG,
                                 "--stratum", STRATUM, "op",
                                 "--workers", str(a.workers)])
        if rc != 0:
            log("operating point failed")
            return 1

    if not base_complete():
        if not base_propagated():
            # The base gets the window minus what one calibration and one
            # ladder need. A base that ran to the wall would leave sixty-four
            # truths on disk and not one comparison to quote from them.
            base_stop = stop_at - timedelta(minutes=CAL_MIN + LADDER_PRIOR_MIN
                                            + RESERVE_MIN)
            log(f"base deadline {base_stop.isoformat(timespec='minutes')} "
                f"({remaining_min(base_stop)/60:.1f} h), leaving "
                f"{CAL_MIN + LADDER_PRIOR_MIN:.0f} min for a calibration and "
                f"the declared budget")
            run_stage("base",
                      [py, "rev30_stratum_base.py", "--registry", REG,
                       "--stratum", STRATUM, "run", "--workers",
                       str(a.workers), "--deadline",
                       base_stop.astimezone().isoformat()])
        if not base_propagated():
            log("base did not finish propagating; nothing is quoted from it. "
                "The sidecars on disk are valid and the run resumes from them.")
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

    log("=== done ===")
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
