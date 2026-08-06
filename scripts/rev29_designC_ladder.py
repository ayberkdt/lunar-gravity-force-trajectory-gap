"""R29 stage C: the archived budget ladder, run on design C at one budget.

The three drivers of the ladder are sha256-pinned in sealed manifests and are
not edited. They are called here through redirections applied at module import:

  * design C is registered in rev14_budget_trajectory.DESIGNS, which is what
    R14 reads for the design source and the reused R11 truth trees, and what
    R18 and R19 reach through in turn;

  * rev14_budget_trajectory.PARETO is pointed at the design-C calibration
    record, because build_specs reads a frozen calibration and the archived one
    has no design C in it. R28 established that this redirection is safe:
    build_specs runs once in the parent and its output travels in the task
    payload, so no worker ever opens a calibration file;

  * rev19_equal_total_work.ANCHOR_BETA is moved off 1.0. R19 writes the anchor
    budget with no budget suffix, and its manifest claims exactly the subtrees
    that carry no suffix. A design-C run at beta = 1 would therefore land inside
    the R19 inventory and break the partition the integrity check enforces.
    Moving the anchor gives design C an explicit beta_1.00 suffix like every
    other budget, which no existing manifest claims.

The last redirection has to survive process spawn, since R19 builds its pool
after the patch, so R19's worker is replaced by a wrapper defined in this module:
unpickling it in a child imports this module, which applies every redirection
above before the wrapper runs the original worker.

Usage:
    python rev29_designC_ladder.py --beta 0.75 --workers 11 --deadline-min 200
    python rev29_designC_ladder.py --beta 0.75 --stage r18 --workers 11
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import rev14_budget_trajectory as r14
import rev18_span_sweep as r18
import rev19_equal_total_work as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PREREG = METRICS / "r29_preregistration.json"
CALIBRATION = METRICS / "r29_budget_pareto_designC.json"
CONVERGENCE = METRICS / "r26_designC_convergence.json"

# --- redirections, at import so that spawned children see them too ----------
r14.DESIGNS["C"] = {
    "rows": METRICS / "r26_designC_rows.json",
    "reuse_case": METRICS / "r11_cases" / "designC_convergence",
    "reuse_raw": METRICS / "r11_raw" / "designC_convergence",
}
r14.PARETO = CALIBRATION

# Any value no budget equals; 1.0 must stop being the unsuffixed special case.
r19.ANCHOR_BETA = -1.0
_R19_WORKER = r19.worker


def r19_worker(task: dict) -> dict:
    """R19's worker, reached through this module so the child process applies
    the redirections above -- above all the anchor move, which decides where
    the child writes."""
    return _R19_WORKER(task)


r19.worker = r19_worker


# ---------------------------------------------------------------------------
def guard(beta: float) -> None:
    if not PREREG.exists():
        raise SystemExit(f"{PREREG.name} missing; run rev29_preregister.py")
    if not CALIBRATION.exists():
        raise SystemExit(f"{CALIBRATION.name} missing; run "
                         f"rev29_designC_pareto.py")
    cal = json.loads(CALIBRATION.read_text(encoding="utf-8"))
    key = r14.tag_of(beta)
    rows = cal["designs"]["C"]["rows"]
    if not rows or key not in rows[0]["budgets"]:
        raise SystemExit(f"design C has no calibration at {key}")
    if not CONVERGENCE.exists():
        raise SystemExit(f"{CONVERGENCE.name} missing; the design-C base has "
                         f"not been generated")
    conv = json.loads(CONVERGENCE.read_text(encoding="utf-8"))
    if not conv.get("complete"):
        raise SystemExit("the design-C base is incomplete; R29 does not quote a "
                         "budget from a partial base")


def stage_r14(beta: float, workers: int, deadline_min: float, limit: int) -> int:
    out = METRICS / f"r14_trajectory_C_{r14.tag_of(beta)}.json"
    deadline = datetime.now(timezone.utc) + timedelta(minutes=deadline_min)
    print(f"[r29-ladder] R14 design C beta={beta:.2f} "
          f"deadline={deadline.astimezone().isoformat(timespec='minutes')}",
          flush=True)
    rc = r14.run("C", beta, workers, deadline, limit)
    print(f"[r29-ladder] R14 rc={rc} -> {out.name}", flush=True)
    return rc


def stage_r18(beta: float, workers: int, deadline_min: float, limit: int) -> int:
    args = argparse.Namespace(design="C", beta=beta, workers=workers,
                              deadline_min=deadline_min, limit=limit)
    print(f"[r29-ladder] R18 design C beta={beta:.2f} "
          f"budget={deadline_min:.0f} min", flush=True)
    rc = r18.run(args)
    print(f"[r29-ladder] R18 rc={rc}", flush=True)
    return rc


def stage_r19(beta: float, workers: int, deadline_min: float) -> int:
    args = argparse.Namespace(design="C", beta=beta, workers=workers,
                              deadline_min=deadline_min)
    print(f"[r29-ladder] R19 design C beta={beta:.2f} "
          f"budget={deadline_min:.0f} min", flush=True)
    rc = r19.run(args)
    print(f"[r29-ladder] R19 rc={rc} -> {r19.out_path('C', beta).name}",
          flush=True)
    return rc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline-min", type=float, default=240.0,
                    help="wall-clock budget for the whole chain")
    ap.add_argument("--stage", choices=("all", "r14", "r18", "r19"),
                    default="all")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    guard(a.beta)
    t0 = time.time()

    def left(reserve: float) -> float:
        return a.deadline_min - (time.time() - t0) / 60.0 - reserve

    if a.stage in ("all", "r14"):
        rc = stage_r14(a.beta, a.workers, left(0.0), a.limit)
        if rc not in (0, 3):
            return rc
    if a.stage in ("all", "r18"):
        if left(10.0) <= 5.0:
            print("[r29-ladder] no time left for R18; stopping here", flush=True)
            return 4
        rc = stage_r18(a.beta, a.workers, left(10.0), a.limit)
        if rc != 0:
            return rc
    if a.stage in ("all", "r19"):
        if left(5.0) <= 5.0:
            print("[r29-ladder] no time left for R19; stopping here", flush=True)
            return 4
        rc = stage_r19(a.beta, a.workers, left(5.0))
        if rc != 0:
            return rc
    print(f"[r29-ladder] beta={a.beta:.2f} chain done in "
          f"{(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
