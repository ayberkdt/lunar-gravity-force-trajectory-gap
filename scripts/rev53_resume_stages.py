"""Finish a ladder cell whose trajectory stage is already on disk.

The equatorial cell wrote a complete R14 record at 13:09 -- 64 rows, no
failures -- and then its span stage stopped answering. Every attempt to restart
the cell re-enters rev30_stratum_ops.stage_ladder, which begins by re-running
R14, and that re-run now kills a pool worker abruptly within a minute. The
trajectory stage is therefore both finished and unrepeatable today, which is a
bad reason to abandon a cell that only needs its two remaining stages.

This runner calls those two stages exactly as stage_ladder calls them: the same
sealed drivers, the same design key, the same budget, the same namespace
fields. It skips nothing that has not already been produced and it writes
nothing R14 owns. If the crash is a property of the machine rather than of the
R14 worker, these stages will fail the same way and say so.

It refuses to start unless the R14 record for the cell is present, complete and
free of failures, because a span sweep over a partial trajectory set would be a
partial cell wearing a whole cell's file name.

Usage:
    python rev53_resume_stages.py --registry r30 --stratum equatorial \
        --beta 0.62 --workers 6 --deadline-min 200
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", required=True)
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--deadline-min", type=float, default=200.0)
    a = ap.parse_args()

    # rev30_stratum_ops binds its registry, stratum and design key from argv at
    # import time, so argv has to carry them before it is imported.
    sys.argv = [sys.argv[0], "--registry", a.registry,
                "--stratum", a.stratum, "status"]
    import rev30_stratum_ops as ops

    key, tag = ops.KEY, f"beta_{a.beta:.2f}"
    traj = METRICS / f"r14_trajectory_{key}_{tag}.json"
    if not traj.exists():
        print(f"[abort] {traj.name} is missing; this runner resumes a cell, it "
              f"does not start one")
        return 2
    rec = json.loads(traj.read_text(encoding="utf-8"))
    if not rec.get("complete") or len(rec.get("rows", [])) != 64 \
            or rec.get("failures"):
        print(f"[abort] {traj.name} is not a complete 64-row record; refusing "
              f"to sweep over it")
        return 2
    print(f"[resume] {key} {tag}: R14 record complete, {len(rec['rows'])} rows, "
          f"written {rec.get('ended_utc')}")

    t0 = time.time()

    def left(reserve: float) -> float:
        return a.deadline_min - (time.time() - t0) / 60.0 - reserve

    span = METRICS / f"r18_span_sweep_{key}_{tag}.json"
    if span.exists():
        print(f"[resume] {span.name} already on disk, skipping R18")
    else:
        print(f"[resume] R18 {a.stratum} beta={a.beta:.2f}", flush=True)
        rc = ops.r18.run(argparse.Namespace(
            design=key, beta=a.beta, workers=a.workers,
            deadline_min=left(10.0), limit=0))
        if rc != 0:
            print(f"[resume] R18 returned {rc}")
            return rc

    print(f"[resume] R19 {a.stratum} beta={a.beta:.2f}", flush=True)
    rc = ops.r19.run(argparse.Namespace(
        design=key, beta=a.beta, workers=a.workers, deadline_min=left(5.0)))
    if rc != 0:
        print(f"[resume] R19 returned {rc}")
        return rc

    print(f"[resume] {key} {tag} done in {(time.time()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
