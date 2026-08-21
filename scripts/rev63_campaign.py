"""Supervisor for R63 (O55): walk the registered ladder cells in one window.

Same discipline as rev61_campaign.py, and for the same reasons: a disk floor
checked before each cell, a completion test that requires the full identity
count rather than the presence of a record, atomic progress writes, and a
stop-at that refuses to start a cell the window cannot finish.

Usage:
    python rev63_campaign.py --stop-at 2026-08-19T18:00 --workers 10
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
PREREG = METRICS / "r63_preregistration.json"
PROGRESS = METRICS / "r63_campaign_progress.json"
DRIVER = HERE / "rev63_ladder_uncapped_rematch.py"

DISK_FLOOR_GB = 6.0
MIN_RECORD_BYTES = 1000
EXPECTED_ORBITS = 64          # 16 identities at each of the four levels


def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def record_path(key: str, beta: float) -> Path:
    return METRICS / f"r63_ladder_uncapped_{key}_{beta_tag(beta)}.json"


def cell_done(key: str, beta: float) -> bool:
    """Size, a parseable summary, and the full identity count. A cell cut by
    the deadline still writes a record summarising whatever finished, so a
    truthy orbit count would read a partial cell as a complete one."""
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
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()

    stop_at = datetime.fromisoformat(a.stop_at).timestamp()
    cells = json.loads(PREREG.read_text(encoding="utf-8"))["cells"]

    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

    state = {"schema": "r63_campaign_progress_v1",
             "started_local": datetime.now().isoformat(timespec="seconds"),
             "stop_at_local": a.stop_at, "workers": a.workers, "cells": []}
    write_progress(state)

    durations = []
    for i, cell in enumerate(cells, 1):
        pop, key, beta = cell["population"], cell["design_key"], cell["beta"]
        label = f"{pop}@{beta_tag(beta)}"

        if cell_done(key, beta):
            print(f"[r63] {label}: already complete, skipping", flush=True)
            state["cells"].append({"cell": label, "status": "skipped_done"})
            write_progress(state)
            continue

        now = time.time()
        projected = max(durations) if durations else 25 * 60.0
        if now >= stop_at or now + projected > stop_at:
            left = max(0.0, (stop_at - now) / 60)
            print(f"[r63] {left:.0f} min left, slowest cell so far "
                  f"{projected / 60:.0f} min; not starting {label}",
                  flush=True)
            state["cells"].append({"cell": label,
                                   "status": "not_started_window"})
            write_progress(state)
            break

        gb = free_gb(METRICS)
        if gb < DISK_FLOOR_GB:
            print(f"[r63] ABORT before {label}: {gb:.1f} GB free, floor is "
                  f"{DISK_FLOOR_GB} GB", flush=True)
            state["cells"].append({"cell": label, "status": "abort_disk",
                                   "free_gb": round(gb, 1)})
            write_progress(state)
            break

        deadline_min = max(1.0, (stop_at - now) / 60.0)
        print(f"=== [{i}/{len(cells)}] R63 {label} start "
              f"{datetime.now():%H:%M:%S} ({gb:.0f} GB free) ===", flush=True)
        t0 = time.time()
        env["R63_POP"] = pop
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
        print(f"=== R63 {label} {status} rc={rc} {wall / 60:.1f} min ===",
              flush=True)
        state["cells"].append({"cell": label, "status": status, "rc": rc,
                               "minutes": round(wall / 60, 1),
                               "finished_local": datetime.now().isoformat(
                                   timespec="seconds")})
        state["updated_local"] = datetime.now().isoformat(timespec="seconds")
        write_progress(state)

    state["ended_local"] = datetime.now().isoformat(timespec="seconds")
    write_progress(state)
    n = sum(1 for c in state["cells"] if c["status"] == "complete")
    print(f"[r63] window over: {n} cells completed this run", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
