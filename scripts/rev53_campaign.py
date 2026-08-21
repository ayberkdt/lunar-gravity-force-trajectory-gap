"""R53 supervisor: the post-hoc budget beta = 0.62 on the seven populations of
the regime map that do not have it.

The shape of rev52_campaign.py with one difference that matters. R51 and R52
each drove one population through a whole chain -- operating point, base,
calibration, ladders -- and the base was the expensive stage. Here every
population already has all of that on disk, verified before this file was
written: seven calibration records whose budgets_computed already contains 0.62,
seven complete bases, seven operating points. Nothing is calibrated, nothing is
propagated twice, and no stage of any earlier campaign is re-entered. What runs
is one ladder per population, the same ladder those populations already ran at
0.50, 0.75 and 1.00, at a fourth budget.

So the loop is over populations rather than over budgets, and the two drivers
differ by population: design C has its own ladder driver because it reaches the
archived chain through import-time redirections, and the other six go through
the generic stratum driver with their own registry and stratum names. Both are
sha256-pinned in sealed manifests and neither is edited.

The cell order is in r53_preregistration.json and is read from there rather than
repeated here, so that the order which runs is the order that was registered. A
cell whose ladder record is already on disk is skipped and said to be skipped; a
cell the clock cannot fit is left unrun and reported as declared and not run,
which the registration names as outcome Z. Nothing is half-run: the estimate
carries a margin over the measured prior of that population's own earlier
ladders, and a cell is not started unless the whole of it fits.

Usage:
    python rev53_campaign.py --stop-at "2026-08-14 10:45" --workers 11
    python rev53_campaign.py --status
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

REG = "r53"
BETA = 0.62

PREREG = METRICS / f"{REG}_preregistration.json"
LOG = HERE / f"{REG}_campaign.log"
PROGRESS = METRICS / f"{REG}_campaign_progress.json"

# a cell is not started unless this much time remains, and the reserve is kept
# clear of the stop so that a finishing ladder is never racing the clock.
RESERVE_MIN = 10.0

# A cell is not started unless this much disk remains either. The first attempt
# at the design-C cell died on OSError 28 with the volume at zero, thirty-four
# minutes in, having written most of a tree it could not finish; the atomic JSON
# write is what kept a truncated record off the disk, but the trees were left to
# be cleaned by hand. The two cells that had completed measured 181 MB each, so
# the floor below is that footprint with better than a factor of two over it.
# The same failure took a ladder of another campaign on 8 August. A clock guard
# without a disk guard only covers the way a run is expected to end.
CELL_MB = 181.0
DISK_FLOOR_MB = 500.0


def free_mb(path: Path) -> float:
    return shutil.disk_usage(path).free / (1024.0 * 1024.0)


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def cells() -> list[dict]:
    if not PREREG.exists():
        raise SystemExit(f"{PREREG.name} is missing; the registration is written "
                         f"before the campaign, never after it")
    reg = json.loads(PREREG.read_text(encoding="utf-8"))
    if reg["budget"] != BETA:
        raise SystemExit(f"the registration is for beta={reg['budget']}, not {BETA}")
    return sorted(reg["cells"], key=lambda c: c["order"])


def remaining_min(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def read_progress() -> dict:
    """The progress file, or a fresh one if it is missing or unreadable.

    It was written straight to its own path, so when the volume filled during
    the design-C cell the write left a zero-byte file and the next --status
    raised out of json.loads. A supervisor's own bookkeeping should not be the
    thing that cannot survive the failure it is there to record.
    """
    blank = {"schema": "r53_campaign_progress_v1", "beta": BETA, "stages": []}
    if not PROGRESS.exists() or PROGRESS.stat().st_size == 0:
        return blank
    try:
        return json.loads(PROGRESS.read_text(encoding="utf-8"))
    except ValueError:
        return blank


def record(name: str, rc: int, minutes: float, note: str = "") -> None:
    p = read_progress()
    p["stages"].append({"stage": name, "rc": rc, "minutes": round(minutes, 1),
                        "note": note,
                        "finished_local": datetime.now().isoformat(
                            timespec="seconds")})
    p["updated_local"] = datetime.now().isoformat(timespec="seconds")
    # written beside the target and moved into place, so a full volume leaves
    # the previous record intact rather than truncating it to nothing
    tmp = PROGRESS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(p, indent=2), encoding="utf-8")
    tmp.replace(PROGRESS)


def ladder_done(key: str) -> bool:
    p = METRICS / f"r19_equal_total_work_{key}_beta_{BETA:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def traj_complete(key: str) -> bool:
    """True when this cell's R14 trajectory record is already finished.

    rev30_stratum_ops.stage_ladder always begins at R14 and rev14.run has no
    resume: it rebuilds every task and overwrites the record. The equatorial
    cell finished its 64-orbit R14 on 14 August and then lost its span stage, so
    restarting the cell would re-propagate two hours of trajectories that are
    already on disk and overwrite a good record to do it. Where R14 is finished
    the cell is handed to the resume runner instead, which refuses anything less
    than a complete, failure-free 64-row record.
    """
    p = METRICS / f"r14_trajectory_{key}_beta_{BETA:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        rec = json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False
    return bool(rec.get("complete")) and len(rec.get("rows", [])) == 64 \
        and not rec.get("failures")


def command(cell: dict, workers: int, deadline_min: float) -> list[str]:
    py = sys.executable
    key = cell["design_key"]
    if traj_complete(key) and cell["driver"] != "rev29_designC_ladder.py":
        return [py, "rev53_resume_stages.py", "--registry", cell["registry"],
                "--stratum", cell["population"], "--beta", f"{BETA:.2f}",
                "--workers", str(workers),
                "--deadline-min", f"{deadline_min:.0f}"]
    if cell["driver"] == "rev29_designC_ladder.py":
        return [py, "rev29_designC_ladder.py", "--beta", f"{BETA:.2f}",
                "--workers", str(workers),
                "--deadline-min", f"{deadline_min:.0f}"]
    return [py, "rev30_stratum_ops.py", "--registry", cell["registry"],
            "--stratum", cell["population"], "ladder", "--beta", f"{BETA:.2f}",
            "--workers", str(workers),
            "--deadline-min", f"{deadline_min:.0f}"]


# The first launch of this campaign lost its ladder three seconds in with exit
# code 0xC000013A, STATUS_CONTROL_C_EXIT: a console control event reached the
# child. launch_detached.py protects the supervisor and the supervisor survived,
# but a child started with a plain subprocess call does not inherit that
# protection by itself. Every cell is therefore started with no console of its
# own and in its own process group, which is the same protection the supervisor
# has and the reason its own worker pool was never the thing that died.
DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200


# The SE cell hung for two hours and forty-five minutes without writing a line.
# Its span-sweep stage never printed even its own header, so it stopped inside
# pool construction rather than part way through the work, and both the cell's
# clock guard and this supervisor's were useless: a process blocked on a pool
# that will never answer cannot check a deadline, and this supervisor was
# blocked in turn waiting for it. A healthy stage is never quiet for long -- the
# span sweep prints its header within seconds and one line per trajectory after
# that, and the longest healthy stage in this campaign ran 21.7 minutes end to
# end. Silence past this threshold therefore means stopped, not slow.
STALL_MIN = 30.0


def run_cell(cell: dict, workers: int, deadline_min: float) -> tuple[int, float]:
    name = f"ladder_{cell['design_key']}_{BETA:.2f}"
    cmd = command(cell, workers, deadline_min)
    log(f"START {name}: {' '.join(cmd[1:])}")
    t0 = time.time()
    # The stall guard below measures progress by the growth of this log, and a
    # child whose stdout is a file block-buffers at 8 KB. The span stage prints
    # about ninety bytes per trajectory and none of its prints ask to be
    # flushed, so a healthy stage can stay silent on disk for far longer than
    # the guard allows and be killed for finishing nothing it had not in fact
    # finished. Unbuffering the child is what makes the guard measure work
    # rather than buffering.
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    with LOG.open("a", encoding="utf-8") as fh:
        proc = subprocess.Popen(
            cmd, cwd=str(HERE), stdout=fh, stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL, env=env,
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)
        rc = None
        last_size, last_change = LOG.stat().st_size, time.time()
        while rc is None:
            try:
                rc = proc.wait(timeout=60)
            except subprocess.TimeoutExpired:
                size = LOG.stat().st_size
                if size != last_size:
                    last_size, last_change = size, time.time()
                elif (time.time() - last_change) / 60.0 > STALL_MIN:
                    quiet = (time.time() - last_change) / 60.0
                    log(f"{name}: no output for {quiet:.0f} min; killing the "
                        f"stage and its workers rather than waiting on a pool "
                        f"that has stopped answering")
                    # /T because the workers are children of the stage, and
                    # killing the stage alone would orphan them onto the CPU
                    subprocess.call(["taskkill", "/PID", str(proc.pid), "/T",
                                     "/F"], stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
                    rc = proc.wait()
                    break
    dt = (time.time() - t0) / 60.0
    log(f"END   {name}: rc={rc} in {dt:.1f} min")
    return rc, dt


def status() -> int:
    print(f"R53: beta={BETA:.2f} on the seven populations that lack it")
    for c in cells():
        state = "done" if ladder_done(c["design_key"]) else "-"
        print(f"  {c['order']}. {c['design_key']:<4} "
              f"{c['population']:<36} {state}")
    if PROGRESS.exists() and read_progress()["stages"]:
        print()
        for st in read_progress()["stages"]:
            print(f"  {st['finished_local']}  {st['stage']:<24} "
                  f"rc={st['rc']} {st['minutes']:>6.1f} min  {st.get('note','')}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at")
    ap.add_argument("--workers", type=int, default=11)
    # Every prior_min in the registration was measured at eleven workers. The
    # five cells that ran did so at eleven; the two that did not each broke a
    # pool worker abruptly at eleven, eight and six, and only four workers held.
    # A cell run at four workers therefore takes several times its registered
    # prior, and the guard that refuses to start a cell it cannot finish is only
    # honest if the prior it compares against is the one for the worker count
    # actually being used. This scales the prior; it does not change the cells,
    # their order, or what any of them runs.
    ap.add_argument("--prior-scale", type=float, default=1.0)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.prior_scale <= 0:
        ap.error("--prior-scale must be positive")
    if a.status:
        return status()
    if not a.stop_at:
        ap.error("--stop-at is required")

    stop_at = datetime.fromisoformat(a.stop_at)
    plan = cells()
    log(f"=== R53 post-hoc budget column, beta={BETA:.2f}: stop at "
        f"{stop_at.isoformat(timespec='minutes')} "
        f"({remaining_min(stop_at)/60:.1f} h), {len(plan)} cells, "
        f"{a.workers} workers, prior x{a.prior_scale:g} ===")

    unrun: list[str] = []
    for cell in plan:
        key = cell["design_key"]
        if ladder_done(key):
            log(f"{key}: already on disk, skipped")
            record(f"ladder_{key}_{BETA:.2f}", 0, 0.0, "already on disk")
            continue
        disk = free_mb(METRICS)
        if disk < DISK_FLOOR_MB:
            log(f"{key}: {disk:.0f} MB free on the metrics volume, a cell "
                f"needs about {CELL_MB:.0f} and this campaign will not start "
                f"one below {DISK_FLOOR_MB:.0f}; stopping rather than filling "
                f"the disk mid-tree")
            record(f"ladder_{key}_{BETA:.2f}", -2, 0.0,
                   f"not started: {disk:.0f} MB free")
            unrun.extend(c["design_key"] for c in plan
                         if c["order"] >= cell["order"]
                         and not ladder_done(c["design_key"]))
            break
        left = remaining_min(stop_at) - RESERVE_MIN
        need = cell["prior_min"] * a.prior_scale
        if left < need:
            log(f"{key}: {left:.0f} min left, this ladder needs {need:.0f} "
                f"({cell['measured_prior']}); declared and not run")
            record(f"ladder_{key}_{BETA:.2f}", -1, 0.0,
                   "declared and not run: outcome Z")
            unrun.append(key)
            continue
        rc, dt = run_cell(cell, a.workers, left - 5.0)
        ok = rc == 0 and ladder_done(key)
        record(f"ladder_{key}_{BETA:.2f}", rc, dt,
               "complete" if ok else "incomplete, nothing is quoted from it")
        if not ok:
            log(f"{key}: rc={rc} and the record is incomplete; nothing is "
                f"quoted from it and the campaign stops rather than starting "
                f"another cell on an unexplained failure")
            unrun.extend(c["design_key"] for c in plan
                         if c["order"] > cell["order"])
            break

    if unrun:
        log(f"not run: {', '.join(unrun)} -- declared and not run (outcome Z)")
    log("=== R53 done ===")
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
