"""Does the geometry crossing survive realistic lunar dynamics?

The force--trajectory reversal was replicated under DE440/MOON_PA with Earth
and Sun third-body gravity and eclipsed SRP on the two confirmatory designs,
and it survived. The geometry result was not: the wide-elliptic population is
where the endpoint ordering reverses, where the paper's headline
"radial wins at the declared budget while the constant degree still wins at
half of it" comes from, and it has only ever been propagated gravity-only. The
title now says "budget--geometry crossings", so the crossing itself is what a
reviewer will ask about, and this campaign asks whether it is a property of the
dynamics or of the isolated system.

Population: OEU, the ceiling-free wide-elliptic population, because that is the
one the printed numbers are read from. The plan for each budget was frozen by
revJ2_fullforce.py's own select command before any trajectory ran, under the
same nested "first N by Sobol index" rule the archived J2 campaign used, and
the policies are read from the frozen budget records rather than recalibrated.

Budgets run in the order 1.00, 0.75, 0.62. That is deliberate: 1.00 carries the
headline result, 0.75 and 0.62 bracket the crossing, so a run that stops early
leaves the most valuable cell complete rather than three partial ones.

Resumption is per trajectory. revJ2_fullforce.command_run skips any trajectory
whose case files already exist, so a stop costs at most the one in flight and a
later chain continues where this one left off.

Four workers, not eleven: on 15 August the pool broke abruptly at eleven, eight
and six on this host and held at four, and four measured only about a fifth
slower on a whole ladder. Children run unbuffered so the progress log is
written as it happens rather than eight kilobytes at a time, which is what the
stall guard reads.

Usage:
    python launch_detached.py ../output/r56_night.log rev56_night_chain.py
"""

from __future__ import annotations

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
OUTPUT = ROOT / "output"

DESIGN = "OEU"
BUDGETS = ("1.00", "0.75", "0.62")

STOP_AT = "2026-08-16 09:00"

# Tried in order, per budget, dropping on a broken pool rather than abandoning
# the budget. Measured on this host: four workers with the threads pinned ran
# 8.6 min and four trajectories and then broke; one worker with the same
# pinning ran its whole window and exited cleanly on the deadline. So
# concurrency is a real factor here even after the thread fix, and the ladder
# ends at the setting that has actually been observed to survive.
WORKER_LADDER = (2, 1)

# Measured: design A ran 144 trajectories in 64.6 min at eleven workers with
# three of its twenty-four orbits above reference degree 300. Every orbit here
# is above it, which roughly doubles the per-budget cost, and four workers add
# about a fifth again. A budget is not started unless this much time remains.
PRIOR_MIN = 260.0

# Raw arrays are offloaded to the D: volume by revJ2_fullforce; the case JSON
# stays on the metrics volume. Both are checked, because either filling stops
# the run and the disk has killed two campaigns already.
FLOOR_METRICS_MB = 700.0
FLOOR_RAW_MB = 2000.0

STALL_MIN = 45.0

# CREATE_NO_WINDOW, not DETACHED_PROCESS. Both keep the child out of this
# shell's console so a stray Ctrl event cannot reach it, but they differ in
# what the child's own pool workers inherit: a DETACHED parent has no console
# at all, so Windows hands each spawned worker a brand new one and four console
# windows appear on the desktop. CREATE_NO_WINDOW gives the parent a console
# that is never displayed, the workers inherit that, and nothing is drawn.
# CREATE_NEW_PROCESS_GROUP still keeps the whole tree out of this shell's group.
CREATE_NO_WINDOW = 0x08000000
CREATE_NEW_PROCESS_GROUP = 0x00000200

LOG = OUTPUT / "r56_night_chain.log"
REPORT = OUTPUT / "r56_morning_report.md"


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    try:
        with LOG.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except OSError:
        pass


def raw_root() -> Path:
    return Path(os.environ.get("JCAMP_RAW_ROOT",
                               r"D:\makale_raw_offload\jgcd"))


def free_mb(p: Path) -> float:
    try:
        return shutil.disk_usage(p).free / (1024.0 * 1024.0)
    except OSError:
        return -1.0


def remaining_min(stop: datetime) -> float:
    return (stop - datetime.now()).total_seconds() / 60.0


def score_path(beta: str) -> Path:
    suffix = "" if beta == "1.00" else f"_beta_{beta}"
    return METRICS / f"rJ2_score{suffix}_design{DESIGN}.json"


def trajectories_done(beta: str) -> int:
    """How many of this budget's trajectories already have both case files.

    The pool breaking is not the same as the pass achieving nothing, and the
    difference decides whether re-entering at the same worker count is worth a
    retry or is a loop. Counting what landed answers that from the disk rather
    than from the exit code.
    """
    plan_suffix = "" if beta == "1.00" else f"_beta_{beta}"
    plan = METRICS / f"rJ2_plan{plan_suffix}_design{DESIGN}.json"
    if not plan.exists():
        return 0
    try:
        rows = json.loads(plan.read_text(encoding="utf-8"))["rows"]
    except (ValueError, OSError, KeyError):
        return 0
    root = METRICS / "rJ2_cases"
    n = 0
    for row in rows:
        case = root / f"J2_design{DESIGN}_{int(row['sobol_index']):03d}"
        if not case.is_dir():
            continue
        n += sum(1 for p in case.glob("*.json") if p.stat().st_size > 0)
    return n


def done(beta: str) -> bool:
    p = score_path(beta)
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        # the score record names its result "verdict"; checking for "outcome"
        # here made every scored budget read as unscored and put a wrong state
        # string in the report
        return "verdict" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def spawn(args: list[str], beta: str, log_path: Path,
          wall_limit_min: float) -> int:
    # One BLAS thread per worker. Without this the pool dies on this host with
    # BrokenProcessPool within half a minute, at four workers and at one alike,
    # while the identical trajectory runs clean in-process in 19 s: the worker
    # is being killed by thread oversubscription inside numpy, not by anything
    # about the trajectory. Pinning the threads was the difference between
    # "0 ok, pool broken" and "1 ok, 0 failed" on the same command.
    env = dict(os.environ, PYTHONUNBUFFERED="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               JCAMP_J2_DESIGN=DESIGN, JCAMP_BETA=beta)
    log(f"start beta={beta}: {' '.join(args)} "
        f"(wall limit {wall_limit_min:.0f} min)")
    t0 = time.time()
    with log_path.open("a", encoding="utf-8") as fh:
        proc = subprocess.Popen(
            [sys.executable] + args, cwd=str(HERE), env=env,
            stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
            creationflags=CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP)
        last_size = log_path.stat().st_size
        last_change = time.time()
        while True:
            try:
                rc = proc.wait(timeout=60)
                break
            except subprocess.TimeoutExpired:
                size = log_path.stat().st_size
                if size != last_size:
                    last_size, last_change = size, time.time()
                elif (time.time() - last_change) / 60.0 > STALL_MIN:
                    log(f"beta={beta}: no output for {STALL_MIN:.0f} min; "
                        f"killing the stage and its workers")
                    subprocess.call(["taskkill", "/PID", str(proc.pid),
                                     "/T", "/F"],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
                    rc = proc.wait()
                    break
                if (time.time() - t0) / 60.0 > wall_limit_min:
                    log(f"beta={beta}: wall limit passed; killing the stage")
                    subprocess.call(["taskkill", "/PID", str(proc.pid),
                                     "/T", "/F"],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
                    rc = proc.wait()
                    break
    log(f"end   beta={beta}: rc={rc} in {(time.time()-t0)/60.0:.1f} min")
    return rc


def run_budget(beta: str, stop: datetime) -> str:
    if done(beta):
        log(f"beta={beta}: already scored, skipped")
        return "already on disk"
    m, r = free_mb(METRICS), free_mb(raw_root())
    if m < FLOOR_METRICS_MB or (r >= 0 and r < FLOOR_RAW_MB):
        log(f"beta={beta}: {m:.0f} MB free on metrics and {r:.0f} MB on the "
            f"raw volume; below the floor, not starting")
        return "not started: disk floor"
    left = remaining_min(stop)
    if left < PRIOR_MIN:
        log(f"beta={beta}: {left:.0f} min left, a budget needs about "
            f"{PRIOR_MIN:.0f}; declared and not run")
        return "declared and not run: clock"

    # A broken pool loses only the trajectories in flight: command_run skips
    # every case already on disk, so dropping a worker and re-entering resumes
    # rather than restarts. The loop therefore keeps going at the same budget
    # until a pass returns cleanly or the ladder and the clock are exhausted.
    rc = 1
    for workers in WORKER_LADDER:
        while rc != 0:
            left = remaining_min(stop)
            if left < 20.0:
                return "out of clock partway through"
            deadline = datetime.now().timestamp() + (left - 15.0) * 60.0
            dl = datetime.fromtimestamp(deadline).isoformat(timespec="seconds")
            before = trajectories_done(beta)
            rc = spawn(["revJ2_fullforce.py", "run",
                        "--workers", str(workers), "--deadline", dl],
                       beta, OUTPUT / f"r56_run_{beta}.log", left + 30.0)
            after = trajectories_done(beta)
            log(f"beta={beta}: pass at {workers} workers rc={rc}, "
                f"{after - before} trajectories added, {after} on disk")
            if rc == 0:
                break
            if after == before:
                log(f"beta={beta}: that pass added nothing; dropping a worker")
                break
    if rc != 0:
        return f"run returned {rc} at every worker count tried"
    rc = spawn(["revJ2_fullforce.py", "score"], beta,
               OUTPUT / f"r56_run_{beta}.log", 30.0)
    if rc != 0:
        return f"score returned {rc}"
    return "complete" if done(beta) else "ran but produced no scored outcome"


def summarise() -> list[str]:
    out = []
    for beta in BUDGETS:
        p = score_path(beta)
        if not p.exists():
            out.append(f"  beta={beta}  not scored")
            continue
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            out.append(f"  beta={beta}  score file unreadable")
            continue
        out.append(f"  beta={beta}  outcome={d.get('outcome')}  "
                   f"{json.dumps({k: v for k, v in d.items() if k in ('resolved', 'sign_agreement', 'reversal', 'n_orbits')})}")
    return out


def main() -> int:
    started = datetime.now()
    OUTPUT.mkdir(exist_ok=True)
    log(f"=== R56: does the geometry crossing survive expanded dynamics? "
        f"population {DESIGN}, budgets {', '.join(BUDGETS)} ===")
    log(f"stop at {STOP_AT}; {free_mb(METRICS):.0f} MB free on metrics, "
        f"{free_mb(raw_root()):.0f} MB on the raw volume")
    states = {}
    try:
        stop = datetime.fromisoformat(STOP_AT)
        for beta in BUDGETS:
            states[beta] = run_budget(beta, stop)
            log(f"beta={beta}: {states[beta]}")
    except Exception as exc:                       # noqa: BLE001
        log(f"the campaign raised {type(exc).__name__}: {exc}")
    lines = [
        "# R56 night chain, wide-elliptic population under expanded dynamics",
        "",
        f"Started {started.isoformat(timespec='seconds')}, "
        f"finished {datetime.now().isoformat(timespec='seconds')}.",
        f"Population {DESIGN}, 24 orbits, policies read from the frozen budget "
        f"records and not recalibrated.",
        "",
        "## Budgets",
        "",
        *[f"- beta={b}: {states.get(b, 'not reached')}" for b in BUDGETS],
        "",
        "## Scores",
        "",
        *summarise(),
        "",
        "## If a budget is missing",
        "",
        "Resumption is per trajectory: rerun this chain and it continues where",
        "it stopped. Nothing is half-written, because a trajectory is only",
        "counted once both of its case files exist.",
        "",
        "## Before quoting any of this",
        "",
        "1. The plans were frozen before propagation; their digests are in",
        "   `metrics/rJ2_plan*_designOEU.json`.",
        "2. `revJ2_fullforce.py` was generalised to name populations beyond the",
        "   two confirmatory designs, so the rJ manifest needs re-sealing and",
        "   the archived A and B results should be confirmed unchanged.",
        "3. The wide-elliptic reversal is currently reported gravity-only in",
        "   Section VIII and the Discussion; whatever this returns belongs",
        "   beside it, not in place of it.",
        "",
    ]
    tmp = REPORT.with_suffix(".md.tmp")
    tmp.write_text("\n".join(lines), encoding="utf-8")
    tmp.replace(REPORT)
    log(f"report written to {REPORT}")
    log("=== R56 done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
