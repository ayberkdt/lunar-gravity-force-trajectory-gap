"""R30: operating point, calibration and budget ladder for one stratum.

The design-C pipeline generalized to a stratum, and generalized by reuse rather
than by copying: the per-orbit operating-point worker is R29's, the calibration
worker is rev14_budget_pareto's, and the ladder is the three pinned drivers.
What this file adds is a stratum's identity -- where its rows live, where its
truth trees are, where its configs and records go -- and the redirections that
put that identity in front of code which knows only about designs A and B.

Every redirection is applied at module import, not inside a function, because
ProcessPoolExecutor on Windows re-imports this module in each child to unpickle
the callables below. A redirection applied in main() would exist in the parent
and nowhere else, and the children would either fail on a missing design key or,
worse, resolve one that exists and write into another population's tree.

Two differences from designs A and B are deliberate and are recorded in the
stratum's own files rather than assumed:

  * the accuracy-target operating point is regenerated under the R12 rule,
    because a stratum has no R12 campaign. R29 checked that regeneration
    against the archived R12 configurations of designs A and B; this file
    points at that check rather than repeating it per stratum;

  * beta = 1 is written with an explicit budget suffix. R19's manifest claims
    exactly the subtrees carrying no suffix, so an unsuffixed stratum record
    would land inside another campaign's inventory.

Usage:
    python rev30_stratum_ops.py --stratum polar op --workers 11
    python rev30_stratum_ops.py --stratum polar cal --workers 11
    python rev30_stratum_ops.py --stratum polar ladder --beta 1.00 --workers 11
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
import population_registry as registry


def _from_argv(flag: str, default: str | None = None) -> str:
    for i, arg in enumerate(sys.argv):
        if arg == flag and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if arg.startswith(flag + "="):
            return arg.split("=", 1)[1]
    if default is None:
        raise SystemExit(f"{flag} is required")
    return default


REGISTRY = _from_argv("--registry", "r30")
STRATUM = _from_argv("--stratum")
PREREG = METRICS / f"{REGISTRY}_preregistration.json"
_PREREG = registry.registration(REGISTRY)
SPEC = registry.spec(REGISTRY, STRATUM)
KEY = SPEC["design_key"]

DESIGN = METRICS / SPEC["file"]
ROWS = METRICS / f"{REGISTRY}_{STRATUM}_rows.json"
CONVERGENCE = METRICS / f"{REGISTRY}_{STRATUM}_convergence.json"
CASE_ROOT_OP = METRICS / f"{REGISTRY}_cases" / f"atallah_{STRATUM}"
OP_RECORD = METRICS / f"{REGISTRY}_{STRATUM}_operating_point.json"
CALIBRATION = METRICS / f"{REGISTRY}_budget_pareto_{STRATUM}.json"
ARCHIVE = METRICS / "r14_budget_pareto.json"
R29_OP_CHECK = METRICS / "r29_designC_operating_point.json"

import rev10_sobol_confirmatory as base                     # noqa: E402
import rev14_budget_pareto as pareto                        # noqa: E402
import rev14_budget_trajectory as r14                       # noqa: E402
import rev18_span_sweep as r18                              # noqa: E402
import rev19_equal_total_work as r19                        # noqa: E402
import rev29_designC_operating_point as op29                # noqa: E402

# --- identity, registered at import so spawned children share it ------------
pareto.DESIGNS[KEY] = {
    "rows": ROWS,
    "r12_case": CASE_ROOT_OP,
    "r11_raw": METRICS / "r11_raw" / f"stratum_{STRATUM}_convergence",
}
r14.DESIGNS[KEY] = {
    "rows": ROWS,
    "reuse_case": METRICS / "r11_cases" / f"stratum_{STRATUM}_convergence",
    "reuse_raw": METRICS / "r11_raw" / f"stratum_{STRATUM}_convergence",
}
r14.PARETO = CALIBRATION
r19.ANCHOR_BETA = -1.0
_R19_WORKER = r19.worker


def r19_worker(task: dict) -> dict:
    return _R19_WORKER(task)


r19.worker = r19_worker


def pareto_worker(task: dict) -> dict:
    return pareto.worker(task)


# ---------------------------------------------------------------- operating point
def write_configs(rec: dict, prereg_sha: str) -> None:
    """The two files rev14_budget_pareto.worker opens, and nothing more."""
    d = CASE_ROOT_OP / f"sobolA_{rec['sobol_index']:03d}"
    common = {
        "sobol_index": rec["sobol_index"],
        "adopted_truth_degree": rec["adopted_truth_degree"],
        "original_truth_degree": rec["original_truth_degree"],
        "n_critical": rec["n_critical"],
        "initial_state_si": rec["initial_state_si"],
        "atallah_tol_accel_m_s2": rec["atallah_tol_accel_m_s2"],
        "atallah_reference": ("Atallah et al. 2022, J.Astronaut.Sci. "
                              "69(3):745-766, Eq.28/20"),
        "perilune_match": ("tol = actual worst-case accel truncation error of "
                           "N_crit vs adopted truth at perilune; Atallah N_req "
                           "binned on 10-km altitude bins, capped at adopted"),
        "atallah_degree_table": rec["atallah_degree_table"],
        "stratum": STRATUM, "design_key": KEY,
        "source": base.provenance(),
    }
    note = (f"regenerated by rev30_stratum_ops.py under the R12 rule; the "
            f"{STRATUM} stratum has no R12 campaign. Configuration only: no "
            f"trajectory is stored for this policy, because the Phase-A "
            f"calibration reads the tolerance, the degree table and the "
            f"work-matched degree and nothing else. The regeneration rule was "
            f"checked against the archived R12 configurations of designs A and "
            f"B in {R29_OP_CHECK.name}.")
    base.atomic_json(d / "atallah_tight.json", {
        "schema": f"{REGISTRY}_operating_point_config_v1",
        "created_utc": base.utc_now(),
        "config": {**common, "policy": "atallah", "level": "tight",
                   "policy_spec": {"kind": "atallah_radial_adaptive",
                                   "tol": rec["atallah_tol_accel_m_s2"]}},
        "config_only": True, "note": note,
        "preregistration_sha256": prereg_sha,
        "tight_atallah_telemetry": rec["telemetry"],
        "propagation_status": rec["propagation_status"]})
    base.atomic_json(d / "fixed_work_atallah_tight.json", {
        "schema": f"{REGISTRY}_operating_point_config_v1",
        "created_utc": base.utc_now(),
        "config": {**common, "policy": "fixed_work_atallah", "level": "tight",
                   "policy_spec": {"kind": "fixed_work", "degree": rec["n_work"],
                                   "source": ("sqrt(<N^2>) of tight Atallah RHS "
                                              "history")}},
        "config_only": True, "note": note,
        "preregistration_sha256": prereg_sha})


def stage_op(workers: int, limit: int) -> int:
    if not ROWS.exists():
        print(f"[abort] {ROWS.name} missing; the prepass has not run")
        return 2
    if not R29_OP_CHECK.exists():
        print(f"[abort] {R29_OP_CHECK.name} missing: the regeneration rule is "
              f"checked against designs A and B once, in R29, and a stratum "
              f"does not get to skip that check by running first")
        return 2
    check = json.loads(R29_OP_CHECK.read_text(encoding="utf-8"))
    bad = [m for c in check.get("reproduction_check", {}).values()
           for m in c.get("pure_mismatches", [])]
    if bad:
        print(f"[abort] the R29 reproduction check recorded {len(bad)} "
              f"mismatches; the rule is not the archived one")
        return 2

    rows = json.loads(ROWS.read_text(encoding="utf-8"))["rows"]
    if limit:
        rows = rows[:limit]
    done, fails = op29.run_pool(rows, workers, f"{STRATUM} operating point")
    for rec in done:
        write_configs(rec, _PREREG["preregistration_sha256"])
    tols = [r["atallah_tol_accel_m_s2"] for r in done]
    works = [r["n_work"] for r in done]
    base.atomic_json(OP_RECORD, {
        "schema": f"{REGISTRY}_operating_point_v1",
        "created_utc": base.utc_now(),
        "stratum": STRATUM, "design_key": KEY,
        "preregistration": PREREG.name,
        "preregistration_sha256": _PREREG["preregistration_sha256"],
        "rule": ("R12 accuracy-target operating point, regenerated; the rule "
                 "was validated against designs A and B in "
                 f"{R29_OP_CHECK.name}"),
        "orbits": len(done), "failures": fails, "rows": done,
        "case_root": str(CASE_ROOT_OP.relative_to(ROOT)),
        "summary": {
            "orbits": len(done),
            "tol_accel_m_s2": {"min": min(tols), "median": float(np.median(tols)),
                               "max": max(tols)} if tols else None,
            "n_work": {"min": min(works), "median": float(np.median(works)),
                       "max": max(works)} if works else None,
        },
        "source": base.provenance()})
    print(f"[written] {OP_RECORD.name}: {len(done)} orbits, "
          f"{len(fails)} failures")
    return 1 if fails else 0


# ------------------------------------------------------------------ calibration
def stage_cal(workers: int, limit: int) -> int:
    if not CONVERGENCE.exists():
        print(f"[abort] {CONVERGENCE.name} missing; the base has not run")
        return 2
    conv = json.loads(CONVERGENCE.read_text(encoding="utf-8"))
    if not conv.get("complete"):
        print(f"[abort] the {STRATUM} base is incomplete; no budget is "
              f"calibrated from a partial base")
        return 2
    if not CASE_ROOT_OP.exists():
        print(f"[abort] operating point missing; run the op stage first")
        return 2

    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    archive_sha = base.file_hash(ARCHIVE)
    grid = list(archived["budget_grid"])
    betas = sorted(set(grid + [0.62]))
    keys = [f"beta_{b:.2f}" for b in betas] + ["original"]

    rows = json.loads(ROWS.read_text(encoding="utf-8"))["rows"]
    if limit:
        rows = rows[:limit]
    tasks = [{"design": KEY, "row": r, "betas": betas} for r in rows]
    print(f"[r30-cal] {STRATUM}: {len(tasks)} orbits x {len(betas)} budgets "
          f"+ operating point", flush=True)
    done, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(pareto_worker, t) for t in tasks]
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['sobol_index']:03d} {rec.get('message')}",
                      flush=True)
            else:
                done.append(rec)
            if n % 8 == 0 or n == len(tasks):
                print(f"  [{n:3d}/{len(tasks)}] "
                      f"elapsed={(time.time()-t0)/60:.1f} min", flush=True)
    done.sort(key=lambda r: r["sobol_index"])
    if fails:
        print(f"[abort] {len(fails)} worker failures; no partial calibration "
              f"is written")
        return 1

    payload = {
        "schema": f"{REGISTRY}_budget_pareto_population_v1",
        "created_utc": base.utc_now(),
        "population": STRATUM, "design_key": KEY,
        "preregistration": PREREG.name,
        "preregistration_sha256": _PREREG["preregistration_sha256"],
        "design_sha256": SPEC["design_sha256"],
        "reference_level": pareto.LEVEL,
        "frozen_budget_grid": grid,
        "budgets_computed": betas,
        "post_hoc_budget": {"beta": 0.62,
                            "status": "declared post hoc, not pre-registered",
                            "declared_in": "r28_calibration_amendment.json"},
        "computed_by": ("rev14_budget_pareto.worker, called verbatim through "
                        "rev30_stratum_ops.pareto_worker"),
        "parent_record": ARCHIVE.name,
        "parent_record_sha256": archive_sha,
        "archive_untouched": ("this stratum is never merged into the archived "
                              "calibration, which is sha256-pinned in three "
                              "sealed manifests"),
        "admissibility_check": {
            "delegated_to": "r29_budget_pareto_designC.json",
            "what": ("the exact reproduction of the archived calibration on "
                     "every orbit and budget of designs A and B was measured "
                     "once, in R29, with this same worker and these same "
                     "inputs; it is not repeated per stratum"),
        },
        "designs": {KEY: {
            "rows": done, "failures": fails,
            "summary": {k: pareto.summarize_budget(done, k) for k in keys},
            "beta_original": pareto.stat(
                [r["budgets"]["original"]["atallah"]["beta_achieved"]
                 for r in done]),
        }},
        "source": base.provenance(),
    }
    base.atomic_json(CALIBRATION, payload)
    if base.file_hash(ARCHIVE) != archive_sha:
        print(f"[FAIL] {ARCHIVE.name} changed during the run")
        return 1
    for k in keys:
        e = payload["designs"][KEY]["summary"][k]
        if not e.get("orbits"):
            continue
        print(f"  {k}: n={e['orbits']} censored={e['censored']} "
              f"R_a median={e['R_a_defect_rms']['median']:.3g}", flush=True)
    print(f"[written] {CALIBRATION.name}: {len(done)} orbits x {len(betas)} "
          f"budgets; archive untouched {archive_sha[:16]}")
    return 0


# ---------------------------------------------------------------------- ladder
def stage_ladder(beta: float, workers: int, deadline_min: float,
                 limit: int) -> int:
    if not CALIBRATION.exists():
        print(f"[abort] {CALIBRATION.name} missing")
        return 2
    cal = json.loads(CALIBRATION.read_text(encoding="utf-8"))
    key = r14.tag_of(beta)
    rows = cal["designs"][KEY]["rows"]
    if not rows or key not in rows[0]["budgets"]:
        print(f"[abort] {STRATUM} has no calibration at {key}")
        return 2

    t0 = time.time()

    def left(reserve: float) -> float:
        return deadline_min - (time.time() - t0) / 60.0 - reserve

    deadline = datetime.now(timezone.utc) + timedelta(minutes=left(0.0))
    print(f"[r30-ladder] R14 {STRATUM} ({KEY}) beta={beta:.2f}", flush=True)
    rc = r14.run(KEY, beta, workers, deadline, limit)
    if rc not in (0, 3):
        return rc
    if left(10.0) <= 5.0:
        print("[r30-ladder] no time left for R18", flush=True)
        return 4
    print(f"[r30-ladder] R18 {STRATUM} beta={beta:.2f}", flush=True)
    rc = r18.run(argparse.Namespace(design=KEY, beta=beta, workers=workers,
                                    deadline_min=left(10.0), limit=limit))
    if rc != 0:
        return rc
    if left(5.0) <= 5.0:
        print("[r30-ladder] no time left for R19", flush=True)
        return 4
    print(f"[r30-ladder] R19 {STRATUM} beta={beta:.2f}", flush=True)
    rc = r19.run(argparse.Namespace(design=KEY, beta=beta, workers=workers,
                                    deadline_min=left(5.0)))
    if rc != 0:
        return rc
    print(f"[r30-ladder] {STRATUM} beta={beta:.2f} chain done in "
          f"{(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--registry", default="r30")
    ap.add_argument("stage", choices=("op", "cal", "ladder", "status"))
    ap.add_argument("--beta", type=float, default=1.00)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline-min", type=float, default=240.0)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    if a.stage == "status":
        print(json.dumps({
            "registry": REGISTRY, "population": STRATUM, "design_key": KEY,
            "rows": ROWS.exists(),
            "base_complete": (CONVERGENCE.exists() and json.loads(
                CONVERGENCE.read_text(encoding="utf-8")).get("complete")),
            "operating_point": OP_RECORD.exists(),
            "calibration": CALIBRATION.exists(),
            "ladders": sorted(
                p.name for p in METRICS.glob(
                    f"r19_equal_total_work_{KEY}_beta_*.json")),
        }, indent=2))
        return 0
    if a.stage == "op":
        return stage_op(a.workers, a.limit)
    if a.stage == "cal":
        return stage_cal(a.workers, a.limit)
    return stage_ladder(a.beta, a.workers, a.deadline_min, a.limit)


if __name__ == "__main__":
    raise SystemExit(main())
