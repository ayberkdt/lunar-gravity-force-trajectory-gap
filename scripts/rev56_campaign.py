"""Supervisor for R56 (O56): member, then comparator, then summarize.

The two propagation stages are ordered by a real dependency, not by taste: the
comparator's degree is sized from the member's own tighter-level telemetry, so
the member must finish first. Each stage carries the disk floor and the
stop-at the other campaigns use, and progress is written atomically after
every stage.

Usage:
    python rev56_campaign.py --stop-at 2026-08-20T11:00 --workers 8
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
PROGRESS = METRICS / "r56_campaign_progress.json"
DRIVER = HERE / "rev56_longarc_interior.py"
RESULT = METRICS / "r56_longarc_interior.json"

DISK_FLOOR_GB = 6.0
PANEL_SIZE = 8


def free_gb(path: Path) -> float:
    return shutil.disk_usage(str(path)).free / 1024 ** 3


def write_progress(payload: dict) -> None:
    tmp = PROGRESS.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, PROGRESS)


def stage_done(which: str) -> int:
    """How many orbits already carry both levels of one stage."""
    root = METRICS / "r56_cases" / which
    if not root.is_dir():
        return 0
    n = 0
    for d in sorted(root.glob("sobolA_*")):
        if (d / "arc_tight.json").exists() and (d / "arc_tighter.json").exists():
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    stop_at = datetime.fromisoformat(a.stop_at).timestamp()

    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

    state = {"schema": "r56_campaign_progress_v1",
             "started_local": datetime.now().isoformat(timespec="seconds"),
             "stop_at_local": a.stop_at, "workers": a.workers, "stages": []}
    write_progress(state)

    for which in ("member", "comparator"):
        have = stage_done(which)
        if have >= PANEL_SIZE:
            print(f"[r56] {which}: all {have} orbits complete, skipping",
                  flush=True)
            state["stages"].append({"stage": which, "status": "skipped_done"})
            write_progress(state)
            continue
        now = time.time()
        if now >= stop_at:
            print(f"[r56] window closed before {which}", flush=True)
            state["stages"].append({"stage": which, "status": "not_reached"})
            write_progress(state)
            break
        gb = free_gb(METRICS)
        if gb < DISK_FLOOR_GB:
            print(f"[r56] ABORT before {which}: {gb:.1f} GB free", flush=True)
            state["stages"].append({"stage": which, "status": "abort_disk"})
            write_progress(state)
            break
        deadline_min = max(1.0, (stop_at - now) / 60.0)
        print(f"=== R56 {which} start {datetime.now():%H:%M:%S} "
              f"({have}/{PANEL_SIZE} done, {gb:.0f} GB free) ===", flush=True)
        t0 = time.time()
        rc = subprocess.call(
            [sys.executable, str(DRIVER), which, "--workers", str(a.workers),
             "--deadline-min", f"{deadline_min:.1f}"], cwd=str(HERE), env=env)
        wall = (time.time() - t0) / 60
        got = stage_done(which)
        status = "complete" if got >= PANEL_SIZE else (
            "partial" if rc == 0 else "failed")
        print(f"=== R56 {which} {status} rc={rc} {got}/{PANEL_SIZE} "
              f"{wall:.1f} min ===", flush=True)
        state["stages"].append({"stage": which, "status": status, "rc": rc,
                                "orbits": got, "minutes": round(wall, 1),
                                "finished_local": datetime.now().isoformat(
                                    timespec="seconds")})
        write_progress(state)
        if status == "failed":
            break

    rc = subprocess.call([sys.executable, str(DRIVER), "summarize"],
                         cwd=str(HERE), env=env)
    state["summarize_rc"] = rc
    state["ended_local"] = datetime.now().isoformat(timespec="seconds")
    write_progress(state)
    print(f"[r56] done; summarize rc={rc}, result "
          f"{'written' if RESULT.exists() else 'absent'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
