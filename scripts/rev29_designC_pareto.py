"""R29 stage B: the Phase-A budget calibration of design C, in its own record.

rev14_budget_trajectory.build_specs does not compute a budget; it reads one from
a frozen calibration. Design C therefore cannot be propagated at any budget
until that calibration exists for it, and r14_budget_pareto.json cannot be the
place it goes: that file is sha256-pinned in three sealed manifests, and merging
a third design into it would break all three at once. Design C gets a separate
record, the same pattern R28 used for beta = 0.62.

Nothing is reimplemented. rev14_budget_pareto.worker is called verbatim, with
design C registered in its design table at module import so that every spawned
worker resolves the same paths as the parent -- a parent-only patch would leave
the children with a design table that has no C in it.

A separate record has to be the same object the archive is, and that is asserted
by measurement rather than by assumption: the same pass recomputes the entire
frozen grid for designs A and B and requires exact equality with the archived
calibration on every orbit, every budget, and the archived operating point. If
one scalar differs, no design-C calibration is written.

The budget grid is the frozen one. beta = 0.62 is carried as well, and is
flagged in the record as what it is: a post-hoc extension declared for designs A
and B in r28_calibration_amendment.json, inherited here without being promoted.

Usage:
    python rev29_designC_pareto.py --workers 11
    python rev29_designC_pareto.py --workers 11 --limit 2 --no-check   # smoke
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import rev10_sobol_confirmatory as base
import rev14_budget_pareto as pareto

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PREREG = METRICS / "r29_preregistration.json"
R28_AMENDMENT = METRICS / "r28_calibration_amendment.json"
ARCHIVE = METRICS / "r14_budget_pareto.json"
R28_RECORD = METRICS / "r28_budget_pareto_beta_0.62.json"
OUTPUT = METRICS / "r29_budget_pareto_designC.json"

POST_HOC_BETA = 0.62

# Registered at import, not in main(): ProcessPoolExecutor on Windows spawns a
# fresh interpreter per worker, which re-imports this module to unpickle the
# callable below. Patching in main() would leave every child without design C.
pareto.DESIGNS["C"] = {
    "rows": METRICS / "r26_designC_rows.json",
    "r12_case": METRICS / "r29_cases" / "atallah_designC",
    "r11_raw": METRICS / "r11_raw" / "designC_convergence",
}

CHECKED = [
    ("atallah", "tol_accel_m_s2"),
    ("atallah", "achieved_work"),
    ("atallah", "beta_achieved"),
    ("atallah", "work_mismatch"),
    ("atallah", "attainable"),
    ("atallah", "limit"),
    ("fixed", "degree"),
    ("fixed", "achieved_work"),
    ("fixed", "beta_achieved"),
    ("fixed", "work_mismatch"),
    ("fixed", "attainable"),
    ("fixed", "limit"),
]


def worker(task: dict) -> dict:
    """rev14_budget_pareto.worker, reached through this module so that the
    child process imports the design-C registration above before running it."""
    return pareto.worker(task)


def tag_of(beta: float) -> str:
    return f"beta_{beta:.2f}"


def compare(fresh: dict, archived: dict, design: str, index: int,
            key: str) -> list[str]:
    bad = []
    if bool(fresh["censored"]) != bool(archived["censored"]):
        bad.append(f"{design}/{index:03d}/{key} censored "
                   f"{fresh['censored']} != {archived['censored']}")
    for kind, field in CHECKED:
        a, b = fresh[kind][field], archived[kind][field]
        if a != b:
            bad.append(f"{design}/{index:03d}/{key} {kind}.{field} "
                       f"{a!r} != {b!r}")
    return bad


def run_design(design: str, rows: list, betas: list, workers: int) -> tuple:
    tasks = [{"design": design, "row": r, "betas": betas} for r in rows]
    print(f"[r29-cal] design {design}: {len(tasks)} orbits x "
          f"{len(betas)} budgets + operating point", flush=True)
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-check", action="store_true",
                    help="smoke only: skip the A/B reproduction gate")
    a = ap.parse_args()

    if not PREREG.exists():
        print(f"[abort] {PREREG.name} missing; run rev29_preregister.py first")
        return 2
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    rows_c_path = pareto.DESIGNS["C"]["rows"]
    if not rows_c_path.exists():
        print(f"[abort] {rows_c_path.name} missing; the design-C prepass has "
              f"not run")
        return 2
    op_root = pareto.DESIGNS["C"]["r12_case"]
    if not op_root.exists():
        print(f"[abort] {op_root.name} missing; run "
              f"rev29_designC_operating_point.py first")
        return 2

    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    archive_sha_before = base.file_hash(ARCHIVE)
    r28 = json.loads(R28_AMENDMENT.read_text(encoding="utf-8"))
    if archive_sha_before != r28["archive_integrity"]["archived_sha256"]:
        print(f"[abort] {ARCHIVE.name} does not match the sha256 on record")
        return 2

    grid = list(archived["budget_grid"])
    betas = sorted(set(grid + [POST_HOC_BETA]))
    keys = [tag_of(b) for b in betas] + ["original"]

    payload = {
        "schema": "r29_budget_pareto_designC_v1",
        "created_utc": base.utc_now(),
        "design": "C",
        "preregistration": PREREG.name,
        "preregistration_sha256": prereg["preregistration_sha256"],
        "parent_design_preregistration_sha256":
            prereg["parent"]["preregistration_sha256"],
        "design_sha256": prereg["parent"]["design_sha256"],
        "reference_level": pareto.LEVEL,
        "frozen_budget_grid": grid,
        "budgets_computed": betas,
        "post_hoc_budget": {
            "beta": POST_HOC_BETA,
            "status": "declared post hoc, not pre-registered",
            "declared_in": R28_AMENDMENT.name,
            "declared_sha256": r28["amendment_sha256"],
            "carries": ("the same flag on design C as on designs A and B; it is "
                        "reported beside the grid, never inside it"),
        },
        "operating_point_source": {
            "record": "r29_designC_operating_point.json",
            "note": ("regenerated under the R12 rule because design C has no "
                     "R12 campaign; see that record for the reproduction check "
                     "against the archived design A and B configurations"),
        },
        "computed_by": ("rev14_budget_pareto.worker, called verbatim through "
                        "rev29_designC_pareto.worker"),
        "parent_record": ARCHIVE.name,
        "parent_record_sha256": archive_sha_before,
        "source": base.provenance(),
        "admissibility_check": {
            "recomputed": ("the entire frozen grid and the archived operating "
                           "point, on every orbit of designs A and B"),
            "requires": "exact equality with the archived calibration",
        },
        "designs": {},
    }

    # ---- gate: the archived designs, recomputed
    all_bad: list[str] = []
    if not a.no_check:
        for d in ("A", "B"):
            rows = json.loads(
                pareto.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]
            if a.limit:
                rows = rows[:a.limit]
            done, fails = run_design(d, rows, grid, a.workers)
            if fails:
                print(f"[abort] design {d}: {len(fails)} worker failures")
                return 1
            arch_by_index = {int(r["sobol_index"]): r
                             for r in archived["designs"][d]["rows"]}
            bad = []
            for rec in done:
                i = int(rec["sobol_index"])
                for k in [tag_of(b) for b in grid] + ["original"]:
                    bad += compare(rec["budgets"][k],
                                   arch_by_index[i]["budgets"][k], d, i, k)
            all_bad += bad
            print(f"  reproduction of {len(grid)} archived budgets + operating "
                  f"point on {len(done)} orbits: "
                  f"{'exact' if not bad else str(len(bad)) + ' MISMATCHES'}",
                  flush=True)
            for line in bad[:10]:
                print(f"    {line}", flush=True)
            payload["admissibility_check"][d] = {
                "orbits": len(done), "budgets": len(grid),
                "mismatches": bad,
            }
        if all_bad:
            print(f"\n[abort] {len(all_bad)} reproduction mismatches: a separate "
                  f"record computed by the same code on the same inputs is not "
                  f"the same object, so design C is not calibrated")
            return 1
        payload["admissibility_check"]["result"] = \
            "exact_on_every_orbit_and_budget_of_designs_A_and_B"
    else:
        payload["admissibility_check"]["result"] = "SKIPPED (smoke run)"

    # ---- design C
    rows = json.loads(rows_c_path.read_text(encoding="utf-8"))["rows"]
    if a.limit:
        rows = rows[:a.limit]
    done, fails = run_design("C", rows, betas, a.workers)
    payload["designs"]["C"] = {
        "rows": done, "failures": fails,
        "summary": {k: pareto.summarize_budget(done, k) for k in keys},
        "beta_original": pareto.stat(
            [r["budgets"]["original"]["atallah"]["beta_achieved"]
             for r in done]) if done else None,
    }
    for k in keys:
        e = payload["designs"]["C"]["summary"][k]
        if not e.get("orbits"):
            print(f"  {k}: all {e.get('censored')} orbits censored", flush=True)
            continue
        flag = "  [post hoc]" if k == tag_of(POST_HOC_BETA) else ""
        print(f"  {k}: n={e['orbits']} censored={e['censored']} "
              f"R_a median={e['R_a_defect_rms']['median']:.3g} "
              f"At wins {e['atallah_wins_defect']}/{e['orbits']}{flag}",
              flush=True)

    if fails:
        print(f"\n[abort] design C: {len(fails)} worker failures; no partial "
              f"calibration is written")
        return 1

    base.atomic_json(OUTPUT, payload)
    archive_sha_after = base.file_hash(ARCHIVE)
    if archive_sha_after != archive_sha_before:
        print(f"[FAIL] {ARCHIVE.name} changed during the run")
        return 1
    print(f"\n[written] {OUTPUT.name}: {len(done)} orbits x {len(betas)} budgets")
    print(f"  archive untouched: {archive_sha_after[:16]}")
    print(f"  admissibility: {payload['admissibility_check']['result']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
