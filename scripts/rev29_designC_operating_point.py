"""R29 stage A: the accuracy-target operating point of design C.

rev14_budget_pareto.worker reads two archived R12 configurations for every
orbit -- the accuracy-target Atallah tolerance with its binned degree table, and
the work-matched fixed degree -- and carries that point through the calibration
as a non-grid budget. Design C has no R12 campaign, so that input does not
exist and the calibration cannot run without it.

Rather than change what the pinned calibration reads, this regenerates the same
quantity by the same rule, into a tree of its own:

    tol      = actual worst-case acceleration truncation error of N_crit
               against the adopted truth degree at perilune, on the same
               25 x 48 lat/lon grid;
    table    = the 10-km binned Atallah schedule from that tol, floor 2,
               capped at the adopted truth degree;
    N_work   = round(sqrt(<N^2>)) of the tight Atallah RHS degree history.

The first two are pure functions of the design point, so they must reproduce
the archived R12 configurations of designs A and B exactly; that is the gate.
The third comes from propagation telemetry, and the archived R12 campaign ran
under a different interpreter, so it is compared and reported per orbit rather
than gated -- a disagreement there is information, not a reason to stop.

No R12 artifact is read for writing, created, or moved. The design-C configs go
to metrics/r29_cases/atallah_designC, which no manifest partition claims.

The fixed-work sidecar holds a configuration and no trajectory: the calibration
reads policy_spec.degree from it and nothing else, and the file says so rather
than looking like a trajectory record that lost its arrays.

Usage:
    python rev29_designC_operating_point.py --workers 11
    python rev29_designC_operating_point.py --workers 11 --check-orbits 8
    python rev29_designC_operating_point.py --limit 2 --check-orbits 2   # smoke
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at
import rev12_atallah_campaign as camp
import rev14_budget_pareto as pareto

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PREREG = METRICS / "r29_preregistration.json"
ROWS_C = METRICS / "r26_designC_rows.json"
CASE_ROOT_C = METRICS / "r29_cases" / "atallah_designC"
OUTPUT = METRICS / "r29_designC_operating_point.json"

N_LAT, N_LON = 25, 48
BIN_KM = 10.0
FLOOR = 2


def worker(task: dict) -> dict:
    """One orbit: the accuracy-target tolerance, its table, and N_work."""
    row = task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        original = int(row["original_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)

        model, args = camp._model(adopted)
        g = camp._g(adopted)
        r_p = model.r_ref + hp_km * 1e3
        tol = at.actual_truncation_error_max(model, args, r_p, n_crit, adopted,
                                             n_lat=N_LAT, n_lon=N_LON)
        degree_fn, table = at.atallah_binned_schedule(
            model, g, tol, hp_km, ha_km, floor=FLOOR, cap=adopted, bin_km=BIN_KM)

        t, y, status, event, failure, telemetry = camp._propagate(
            model, args, y0, degree_fn, "tight")
        if status == "numerical_failure":
            return {"sobol_index": index, "status": "numerical_failure",
                    "message": failure}
        n_work = int(round(math.sqrt(telemetry["mean_degree_sq"])))

        return {
            "sobol_index": index,
            "status": "complete",
            "adopted_truth_degree": adopted,
            "original_truth_degree": original,
            "n_critical": n_crit,
            "initial_state_si": [float(v) for v in y0],
            "atallah_tol_accel_m_s2": float(tol),
            "atallah_degree_table": {str(k): int(v) for k, v in table.items()},
            "n_work": n_work,
            "telemetry": {"mean_degree": telemetry.get("mean_degree"),
                          "mean_degree_sq": telemetry.get("mean_degree_sq"),
                          "degree_range": telemetry.get("degree_range"),
                          "n_rhs": telemetry.get("n_rhs")},
            "propagation_status": status,
            "impacted": bool(event),
        }
    except Exception as exc:
        return {"sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def run_pool(rows: list, workers: int, label: str) -> tuple[list, list]:
    tasks = [{"row": r} for r in rows]
    print(f"[r29-op] {label}: {len(tasks)} orbits", flush=True)
    done, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(worker, t) for t in tasks]
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['sobol_index']:03d} {rec.get('message')}",
                      flush=True)
            else:
                done.append(rec)
            if n % 8 == 0 or n == len(tasks):
                el = time.time() - t0
                print(f"  [{n:3d}/{len(tasks)}] elapsed={el/60:5.1f} min "
                      f"eta={(len(tasks)-n)*el/n/60:5.1f} min", flush=True)
    done.sort(key=lambda r: r["sobol_index"])
    return done, fails


def archived_config(design: str, index: int) -> tuple[dict, dict]:
    d = pareto.DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
    at_cfg = json.loads((d / "atallah_tight.json").read_text())["config"]
    fw_cfg = json.loads((d / "fixed_work_atallah_tight.json").read_text())["config"]
    return at_cfg, fw_cfg


def check_against_archive(design: str, rows: list, workers: int) -> dict:
    """Reproduce the R12 operating point on a subsample of an archived design."""
    done, fails = run_pool(rows, workers, f"reproduction check, design {design}")
    pure_bad, work_bad = [], []
    for rec in done:
        i = rec["sobol_index"]
        at_cfg, fw_cfg = archived_config(design, i)
        if rec["atallah_tol_accel_m_s2"] != at_cfg["atallah_tol_accel_m_s2"]:
            pure_bad.append(f"{design}/{i:03d} tol "
                            f"{rec['atallah_tol_accel_m_s2']!r} != "
                            f"{at_cfg['atallah_tol_accel_m_s2']!r}")
        if rec["atallah_degree_table"] != at_cfg["atallah_degree_table"]:
            pure_bad.append(f"{design}/{i:03d} binned degree table differs")
        archived_work = int(fw_cfg["policy_spec"]["degree"])
        if rec["n_work"] != archived_work:
            work_bad.append({"design": design, "sobol_index": i,
                             "regenerated": rec["n_work"],
                             "archived": archived_work,
                             "rel_diff": rec["n_work"] / archived_work - 1.0})
    return {"orbits_checked": len(done), "failures": fails,
            "pure_mismatches": pure_bad, "n_work_differences": work_bad}


def write_configs(rec: dict, prereg_sha: str) -> None:
    """The two files rev14_budget_pareto.worker opens, and nothing more."""
    d = CASE_ROOT_C / f"sobolA_{rec['sobol_index']:03d}"
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
        "design": "C",
        "source": base.provenance(),
    }
    note = ("regenerated by rev29_designC_operating_point.py under the R12 rule; "
            "design C has no R12 campaign. This file carries a configuration "
            "only: no trajectory is stored for this policy, because the Phase-A "
            "calibration reads the tolerance, the degree table and the "
            "work-matched degree and nothing else.")
    base.atomic_json(d / "atallah_tight.json", {
        "schema": "r29_operating_point_config_v1",
        "created_utc": base.utc_now(),
        "config": {**common, "policy": "atallah", "level": "tight",
                   "policy_spec": {"kind": "atallah_radial_adaptive",
                                   "tol": rec["atallah_tol_accel_m_s2"]}},
        "config_only": True, "note": note,
        "preregistration_sha256": prereg_sha,
        "tight_atallah_telemetry": rec["telemetry"],
        "propagation_status": rec["propagation_status"]})
    base.atomic_json(d / "fixed_work_atallah_tight.json", {
        "schema": "r29_operating_point_config_v1",
        "created_utc": base.utc_now(),
        "config": {**common, "policy": "fixed_work_atallah", "level": "tight",
                   "policy_spec": {"kind": "fixed_work", "degree": rec["n_work"],
                                   "source": ("sqrt(<N^2>) of tight Atallah RHS "
                                              "history")}},
        "config_only": True, "note": note,
        "preregistration_sha256": prereg_sha})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--check-orbits", type=int, default=8,
                    help="orbits per archived design in the reproduction check")
    ap.add_argument("--limit", type=int, default=0, help="smoke: design-C orbits")
    a = ap.parse_args()

    if not PREREG.exists():
        print(f"[abort] {PREREG.name} missing; run rev29_preregister.py first")
        return 2
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    if not ROWS_C.exists():
        print(f"[abort] {ROWS_C.name} missing; the design-C prepass has not run")
        return 2

    payload = {
        "schema": "r29_designC_operating_point_v1",
        "created_utc": base.utc_now(),
        "design": "C",
        "preregistration": PREREG.name,
        "preregistration_sha256": prereg["preregistration_sha256"],
        "rule": ("R12 accuracy-target operating point: tol = actual worst-case "
                 "accel truncation error of N_crit vs adopted truth at perilune "
                 f"({N_LAT}x{N_LON} lat/lon grid), Atallah N_req binned on "
                 f"{BIN_KM:.0f}-km bins with floor {FLOOR} and cap = adopted, "
                 "N_work = round(sqrt(<N^2>)) of the tight Atallah RHS history"),
        "source": base.provenance(),
        "reproduction_check": {},
    }

    # ---- gate first: the pure quantities must reproduce the archive
    checks = {}
    if a.check_orbits:
        for d in ("A", "B"):
            rows = json.loads(
                pareto.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]
            checks[d] = check_against_archive(d, rows[:a.check_orbits], a.workers)
            c = checks[d]
            print(f"  design {d}: {c['orbits_checked']} orbits, "
                  f"tol+table {'exact' if not c['pure_mismatches'] else 'MISMATCH'}, "
                  f"N_work differs on {len(c['n_work_differences'])}", flush=True)
            for line in c["pure_mismatches"][:5]:
                print(f"    {line}", flush=True)
            for w in c["n_work_differences"][:5]:
                print(f"    N_work {d}/{w['sobol_index']:03d} "
                      f"{w['regenerated']} vs archived {w['archived']} "
                      f"({w['rel_diff']:+.3%})", flush=True)
    payload["reproduction_check"] = checks
    bad = [m for c in checks.values() for m in c["pure_mismatches"]]
    if bad:
        print(f"\n[abort] {len(bad)} mismatches in quantities that are pure "
              f"functions of the design point: the regenerated operating point "
              f"is not the archived rule, so design C gets none")
        return 1
    if any(c["failures"] for c in checks.values()):
        print("\n[abort] reproduction check had worker failures")
        return 1

    # ---- design C
    rows = json.loads(ROWS_C.read_text(encoding="utf-8"))["rows"]
    if a.limit:
        rows = rows[:a.limit]
    done, fails = run_pool(rows, a.workers, "design C operating point")
    for rec in done:
        write_configs(rec, prereg["preregistration_sha256"])

    payload["orbits"] = len(done)
    payload["failures"] = fails
    payload["rows"] = done
    payload["case_root"] = str(CASE_ROOT_C.relative_to(ROOT))
    tols = [r["atallah_tol_accel_m_s2"] for r in done]
    works = [r["n_work"] for r in done]
    payload["summary"] = {
        "orbits": len(done),
        "tol_accel_m_s2": {"min": min(tols), "median": float(np.median(tols)),
                           "max": max(tols)} if tols else None,
        "n_work": {"min": min(works), "median": float(np.median(works)),
                   "max": max(works)} if works else None,
    }
    base.atomic_json(OUTPUT, payload)
    print(f"\n[written] {OUTPUT.name}: {len(done)} orbits, {len(fails)} failures")
    print(f"  configs under {payload['case_root']}")
    if tols:
        print(f"  tol median {payload['summary']['tol_accel_m_s2']['median']:.3g} "
              f"m/s^2, N_work median "
              f"{payload['summary']['n_work']['median']:.0f}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
