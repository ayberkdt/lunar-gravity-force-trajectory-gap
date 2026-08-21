"""R50 supervisor: the paired radial-span ladder, block A then block B.

The same shape as rev38_campaign, with the stages in the order that hurt least
when the clock runs out, and three differences.

The prepass runs first and the operating point second, before any trajectory is
propagated, because both read the design and the prepass rows and nothing else.
Ninety minutes of propagation should not be what discovers that a later stage is
misconfigured. The reindex step is inside the chain rather than left to be
remembered: a two-policy base is propagated correctly and then indexed
incorrectly by the pinned driver, and R30 lost thirty-four hours to that gap
sitting outside a supervisor.

Block B is conditional, and the condition is measured rather than assumed. Its
whole chain has to fit in what is left, using block A's own measured stage times
on the same sixty-four orbits, plus a margin. A block that starts and cannot
finish leaves a panel that carries some budgets and not others, and the
registration says what is quoted from a partial ladder: nothing.

Every stage is skipped when its output already exists and the base is
checkpointed per trajectory, so an interrupted run resumes where it stopped.

Usage:
    python rev50_campaign.py --stop-at "2026-08-12 11:30" --workers 11
    python rev50_campaign.py --status
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

REG = "r50"
BLOCKS = [("span_ladder_a", "RS1"), ("span_ladder_b", "RS2")]
BETAS = [1.00, 0.75, 0.62, 0.50]
# r50_budget_extension_amendment.json adds these, declared before any ladder
# ran and conditional on the clock. They are not run unless asked for by name,
# so a default invocation carries the registered grid and nothing else.
AMENDMENT_BETAS = [1.25, 1.50]

LOG = HERE / "r50_campaign.log"
PROGRESS = METRICS / "r50_campaign_progress.json"

# priors, from the R31 population: same perilune band, same adopted degree 300,
# same sixty-four orbits, base 72 min and ladders 46, 44 and 35 min.
PREPASS_MIN = 12.0
OP_MIN = 10.0
BASE_MIN = 95.0
CAL_MIN = 10.0
LADDER_PRIOR_MIN = 55.0
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
         else {"schema": "r50_campaign_progress_v1", "stages": []})
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


def rows_path(name: str) -> Path:
    return METRICS / f"{REG}_{name}_rows.json"


def conv_path(name: str) -> Path:
    return METRICS / f"{REG}_{name}_convergence.json"


def op_path(name: str) -> Path:
    return METRICS / f"{REG}_{name}_operating_point.json"


def cal_path(name: str) -> Path:
    return METRICS / f"{REG}_budget_pareto_{name}.json"


def base_complete(name: str) -> bool:
    p = conv_path(name)
    return (p.exists()
            and bool(json.loads(p.read_text(encoding="utf-8")).get("complete")))


def base_propagated(name: str) -> bool:
    """Every trajectory on disk, whether or not the index says so."""
    raw = METRICS / "r11_raw" / f"stratum_{name}_convergence"
    if not raw.exists():
        return False
    want = {f"{p}_{l}.npz" for p in ("truth", "fixed_critical")
            for l in ("tight", "tighter")}
    dirs = [d for d in raw.iterdir() if d.is_dir()]
    if len(dirs) != 64:
        return False
    return all(want <= {f.name for f in d.iterdir()} for d in dirs)


def ladder_done(key: str, beta: float) -> bool:
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def status() -> int:
    for name, key in BLOCKS:
        print(f"{name:<15} ({key}) rows={'y' if rows_path(name).exists() else '-'} "
              f"op={'y' if op_path(name).exists() else '-'} "
              f"propagated={'y' if base_propagated(name) else '-'} "
              f"indexed={'y' if base_complete(name) else '-'} "
              f"cal={'y' if cal_path(name).exists() else '-'} "
              f"ladders="
              f"{','.join(f'{b:.2f}' for b in BETAS if ladder_done(key, b)) or '-'}")
    if PROGRESS.exists():
        for st in json.loads(PROGRESS.read_text(encoding="utf-8"))["stages"]:
            print(f"  {st['finished_local']}  {st['stage']:<28} "
                  f"rc={st['rc']} {st['minutes']:>6.1f} min")
    return 0


def run_block(name: str, key: str, stop_at: datetime, workers: int, py: str,
              ladder_est: float, betas: list[float]) -> tuple[bool, float, float]:
    """Returns (block finished every budget, spent minutes, ladder estimate)."""
    spent = 0.0

    if not rows_path(name).exists():
        if remaining_min(stop_at) < PREPASS_MIN + OP_MIN + BASE_MIN + CAL_MIN \
                + ladder_est + RESERVE_MIN:
            log(f"{name}: not enough time for a chain that could be quoted; "
                f"the block is not started")
            return False, spent, ladder_est
        rc, dt = run_stage(f"{name}/prepass",
                           [py, "rev30_stratum_base.py", "--registry", REG,
                            "--stratum", name, "prepass",
                            "--workers", str(workers)])
        spent += dt
        if rc != 0 or not rows_path(name).exists():
            log(f"{name}: prepass failed")
            return False, spent, ladder_est

    if not op_path(name).exists():
        rc, dt = run_stage(f"{name}/op",
                           [py, "rev30_stratum_ops.py", "--registry", REG,
                            "--stratum", name, "op", "--workers", str(workers)])
        spent += dt
        if rc != 0:
            log(f"{name}: operating point failed")
            return False, spent, ladder_est

    if not base_complete(name):
        if not base_propagated(name):
            # the base gets the window minus what one calibration and one
            # ladder need: sixty-four truths on disk and not one comparison to
            # quote from them is the worst way to spend a night.
            base_stop = stop_at - timedelta(
                minutes=CAL_MIN + ladder_est + RESERVE_MIN)
            if remaining_min(base_stop) < 20.0:
                log(f"{name}: no room for a base that could then be calibrated "
                    f"and laddered")
                return False, spent, ladder_est
            log(f"{name}: base deadline "
                f"{base_stop.isoformat(timespec='minutes')}, leaving "
                f"{CAL_MIN + ladder_est:.0f} min for a calibration and one "
                f"ladder")
            _, dt = run_stage(f"{name}/base",
                              [py, "rev30_stratum_base.py", "--registry", REG,
                               "--stratum", name, "run", "--workers",
                               str(workers), "--deadline",
                               base_stop.astimezone().isoformat()])
            spent += dt
        if not base_propagated(name):
            log(f"{name}: base did not finish propagating; nothing is quoted "
                f"from it. The sidecars on disk are valid and a later run "
                f"resumes from them.")
            return False, spent, ladder_est
        rc, dt = run_stage(f"{name}/reindex",
                           [py, "rev32_reindex_two_policy_base.py",
                            "--registry", REG, "--stratum", name])
        spent += dt
        if rc != 0 or not base_complete(name):
            log(f"{name}: reindex failed; the base is propagated but not "
                f"indexed")
            return False, spent, ladder_est

    if not cal_path(name).exists():
        if remaining_min(stop_at) < CAL_MIN + RESERVE_MIN:
            log(f"{name}: no time for the calibration")
            return False, spent, ladder_est
        rc, dt = run_stage(f"{name}/cal",
                           [py, "rev30_stratum_ops.py", "--registry", REG,
                            "--stratum", name, "cal", "--workers",
                            str(workers)])
        spent += dt
        if rc != 0:
            log(f"{name}: calibration failed or was refused")
            return False, spent, ladder_est

    measured: list[float] = []
    for beta in betas:
        if ladder_done(key, beta):
            log(f"{name}: beta {beta:.2f} already on disk")
            continue
        left = remaining_min(stop_at) - RESERVE_MIN
        if left < ladder_est:
            log(f"{name}: {left:.0f} min left, a ladder needs "
                f"{ladder_est:.0f}; stopping rather than half-running one")
            return False, spent, ladder_est
        rc, dt = run_stage(f"{name}/ladder_{beta:.2f}",
                           [py, "rev30_stratum_ops.py", "--registry", REG,
                            "--stratum", name, "ladder", "--beta",
                            f"{beta:.2f}", "--workers", str(workers),
                            "--deadline-min", f"{left - 5.0:.0f}"])
        spent += dt
        if rc == 0 and ladder_done(key, beta):
            measured.append(dt)
            ladder_est = max(max(measured) * 1.15, 30.0)
            log(f"{name}: ladder estimate now {ladder_est:.0f} min")
        else:
            log(f"{name}: beta {beta:.2f} rc={rc}, incomplete")
            return False, spent, ladder_est

    return True, spent, ladder_est


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--blocks", default="",
                    help="comma-separated subset of the declared blocks")
    ap.add_argument("--betas", default="",
                    help=("comma-separated budgets; the registered four by "
                          "default. The amendment budgets are run only when "
                          "named here, and only once the registered four are "
                          "on disk for the block."))
    a = ap.parse_args()
    if a.status:
        return status()
    if not a.stop_at:
        ap.error("--stop-at is required")

    stop_at = datetime.fromisoformat(a.stop_at)
    py = sys.executable
    wanted = [s.strip() for s in a.blocks.split(",") if s.strip()]
    blocks = [b for b in BLOCKS if not wanted or b[0] in wanted]
    betas = [float(s) for s in a.betas.split(",") if s.strip()] or list(BETAS)
    unknown = [b for b in betas if b not in BETAS + AMENDMENT_BETAS]
    if unknown:
        ap.error(f"{unknown} is not a registered or amended budget; the "
                 f"registered grid is {BETAS} and the amendment adds "
                 f"{AMENDMENT_BETAS}")
    if any(b in AMENDMENT_BETAS for b in betas):
        if not (METRICS / "r50_budget_extension_amendment.json").exists():
            ap.error("an amendment budget was asked for and "
                     "r50_budget_extension_amendment.json does not exist")
        # the amendment is conditional on the registered grid being carried
        # first; an amended cell standing where a registered one is missing
        # would be a grid chosen after the fact.
        for name, key in blocks:
            missing = [f"{b:.2f}" for b in BETAS if not ladder_done(key, b)]
            if missing:
                ap.error(f"{name} is missing registered budgets {missing}; the "
                         f"amendment runs after the registered grid, not "
                         f"instead of part of it")
    log(f"=== R50 paired radial-span ladder: stop at "
        f"{stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h), blocks="
        f"{[b[0] for b in blocks]}, betas={betas} ===")

    ladder_est = LADDER_PRIOR_MIN
    first_spent = 0.0
    for n, (name, key) in enumerate(blocks):
        if n > 0:
            # block B is conditional on the clock, measured on block A rather
            # than on a prior. A block that starts and cannot finish leaves a
            # panel carrying some budgets and not others.
            need = first_spent * 1.10 if first_spent else (
                PREPASS_MIN + OP_MIN + BASE_MIN + CAL_MIN
                + len(betas) * ladder_est)
            left = remaining_min(stop_at) - RESERVE_MIN
            if left < need:
                log(f"{name}: {left:.0f} min left and the measured chain needs "
                    f"{need:.0f}; the block is reported as declared and not "
                    f"run")
                break
            log(f"{name}: {left:.0f} min left against a measured chain of "
                f"{need:.0f}; starting")
        done, spent, ladder_est = run_block(name, key, stop_at, a.workers, py,
                                            ladder_est, betas)
        if n == 0:
            first_spent = spent
        log(f"{name}: {'complete' if done else 'incomplete'} after "
            f"{spent:.0f} min")
        if not done:
            break

    log(f"=== R50 done, {remaining_min(stop_at):.0f} min before the stop ===")
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
