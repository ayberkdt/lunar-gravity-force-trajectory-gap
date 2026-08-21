"""Supervisor for R61 (O42-ext): walk the registered cells inside one window.

It carries the four things earlier campaigns lost a night to:

  * a disk floor checked before every cell, because a full volume has killed
    two runs and left a zero-byte record behind that looked finished;
  * a completion test that requires size and a parseable summary, not
    ``exists()`` -- the ``ladder_done()`` defect of 2026-08-08;
  * atomic progress writes, so a kill mid-write cannot corrupt the record of
    what has run;
  * a stop-at that is checked against each cell's measured cost before the
    cell starts, so the window ends on a finished cell rather than a partial
    one.

The cell order comes from the registration and is never reordered here.

Usage:
    python rev61_campaign.py --stop-at 2026-08-19T11:00 --workers 10
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
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
PREREG = METRICS / "r61_preregistration.json"
PROGRESS = METRICS / "r61_campaign_progress.json"
DRIVER = HERE / "rev61_equal_work_tighter_ext.py"

DISK_FLOOR_GB = 6.0
MIN_RECORD_BYTES = 2000
EXPECTED_ORBITS = 64


def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def record_path(key: str, beta: float) -> Path:
    return METRICS / f"r61_equal_work_tighter_{key}_{beta_tag(beta)}.json"


def cell_done(key: str, beta: float) -> bool:
    """Size, a parseable summary, and the full orbit count -- never bare
    existence.

    A cell cut by the deadline still writes a record: run() calls summarize()
    whatever happened, and that record summarises the orbits that did finish.
    Testing for the record, or for a truthy orbit count, would read a partial
    cell as a finished one, which is the 2026-08-08 ladder defect wearing a
    different hat. Every registered cell planned 64 orbits with none censored,
    so 64 is the completion count; a cell short of it is resumed, and the
    worker skips the sidecars it already wrote.
    """
    p = record_path(key, beta)
    try:
        if not p.exists() or p.stat().st_size < MIN_RECORD_BYTES:
            return False
        d = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    s = d.get("summary")
    return isinstance(s, dict) and s.get("orbits") == EXPECTED_ORBITS


def free_gb(path: Path) -> float:
    return shutil.disk_usage(str(path)).free / 1024 ** 3


def write_progress(payload: dict) -> None:
    tmp = PROGRESS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, PROGRESS)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True,
                    help="local ISO time, e.g. 2026-08-19T11:00")
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()

    stop_at = datetime.fromisoformat(a.stop_at).timestamp()
    reg = json.loads(PREREG.read_text(encoding="utf-8"))
    cells = reg["cells"]

    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

    state = {"schema": "r61_campaign_progress_v1",
             "started_local": datetime.now().isoformat(timespec="seconds"),
             "stop_at_local": a.stop_at, "workers": a.workers, "cells": []}
    write_progress(state)

    durations = []
    for i, cell in enumerate(cells, 1):
        pop, key, beta = cell["population"], cell["design_key"], cell["beta"]
        label = f"{pop}@{beta_tag(beta)}"

        if cell_done(key, beta):
            print(f"[r61] {label}: already complete, skipping", flush=True)
            state["cells"].append({"cell": label, "status": "skipped_done"})
            write_progress(state)
            continue

        now = time.time()
        if now >= stop_at:
            print(f"[r61] window closed before {label}; stopping", flush=True)
            state["cells"].append({"cell": label, "status": "not_reached"})
            write_progress(state)
            break

        # Do not start a cell the window cannot finish: use the slowest cell
        # measured so far, and R44's slowest cell (45 min) until one exists.
        projected = max(durations) if durations else 45 * 60.0
        if now + projected > stop_at:
            left = (stop_at - now) / 60
            print(f"[r61] {left:.0f} min left, slowest cell so far "
                  f"{projected / 60:.0f} min; not starting {label}",
                  flush=True)
            state["cells"].append({"cell": label,
                                   "status": "not_started_window"})
            write_progress(state)
            break

        gb = free_gb(METRICS)
        if gb < DISK_FLOOR_GB:
            print(f"[r61] ABORT before {label}: {gb:.1f} GB free on the "
                  f"metrics volume, floor is {DISK_FLOOR_GB} GB", flush=True)
            state["cells"].append({"cell": label, "status": "abort_disk",
                                   "free_gb": round(gb, 1)})
            write_progress(state)
            break

        deadline_min = max(1.0, (stop_at - now) / 60.0)
        print(f"=== [{i}/{len(cells)}] R61 {label} start "
              f"{datetime.now():%H:%M:%S} ({gb:.0f} GB free) ===", flush=True)
        t0 = time.time()
        env["R61_POP"] = pop
        rc = subprocess.call(
            [sys.executable, str(DRIVER), "run", "--beta", f"{beta:.2f}",
             "--workers", str(a.workers),
             "--deadline-min", f"{deadline_min:.1f}"],
            cwd=str(HERE), env=env)
        wall = time.time() - t0

        done = cell_done(key, beta)
        if done:
            durations.append(wall)
        status = "complete" if done else ("partial" if rc == 0 else "failed")
        print(f"=== R61 {label} {status} rc={rc} "
              f"{wall / 60:.1f} min ===", flush=True)
        state["cells"].append({"cell": label, "status": status, "rc": rc,
                               "minutes": round(wall / 60, 1),
                               "finished_local": datetime.now().isoformat(
                                   timespec="seconds")})
        state["updated_local"] = datetime.now().isoformat(timespec="seconds")
        write_progress(state)

        if status == "failed":
            print(f"[r61] {label} failed; continuing to the next cell",
                  flush=True)

    state["ended_local"] = datetime.now().isoformat(timespec="seconds")
    write_progress(state)
    n_done = sum(1 for c in state["cells"] if c["status"] == "complete")
    print(f"[r61] window over: {n_done} cells completed this run", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
