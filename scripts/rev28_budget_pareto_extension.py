"""R28: the Phase-A calibration point at beta = 0.62, in its own record.

Declared post hoc by r28_calibration_amendment.json, which explains why a budget
outside the frozen grid needs a registration act and what it does and does not
license. This script is the mechanical half of that: it produces the calibration
that rev14_budget_trajectory.build_specs needs and cannot compute for itself.

Nothing here reimplements the calibration. rev14_budget_pareto.worker is called
verbatim, on the same archived truth epochs, with the same binning, floor, cap,
bisection target and censoring rule. The only new input is the value of beta.

Two properties this file has to earn:

  * it must not touch the archive. r14_budget_pareto.json is sha256-pinned in
    three sealed manifests, so the new budget goes in a separate record and the
    archived file is verified byte-identical before and after.

  * a separate record must be the same object the archive is. The worker is
    asked for beta = 0.75 in the same pass, and every calibrated tolerance,
    comparator degree, achieved work and censoring flag it returns must equal
    the archived value exactly, on all 128 orbits. Reproduction is asserted by
    construction -- the defect evaluation is deterministic and integration-free
    -- so the check is cheap, and cheap checks of things assumed true are the
    ones worth running. If it fails anywhere, no 0.62 record is written.

The archived accuracy-target point ("original") comes along because the worker
computes it unconditionally; it is carried into the record unused, and it gives
the reproduction check a third independent quantity per orbit for free.

Usage:
    python rev28_budget_pareto_extension.py --workers 8
    python rev28_budget_pareto_extension.py --workers 8 --limit 4    # smoke
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
ARCHIVE = METRICS / "r14_budget_pareto.json"
AMENDMENT = METRICS / "r28_calibration_amendment.json"
OUTPUT = METRICS / "r28_budget_pareto_beta_0.62.json"

NEW_BETA = 0.62
CHECK_BETA = 0.75          # archived budget recomputed as an admissibility check

# every scalar build_specs reads, plus the two the check can compare for free
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


def tag_of(beta: float) -> str:
    """The key rev14_budget_pareto.worker writes and build_specs reads."""
    return f"beta_{beta:.2f}"


def compare(fresh: dict, archived: dict, design: str, index: int) -> list[str]:
    """Exact equality, not tolerance: same code, same inputs, same floats."""
    bad = []
    if bool(fresh["censored"]) != bool(archived["censored"]):
        bad.append(f"{design}/{index:03d} censored "
                   f"{fresh['censored']} != {archived['censored']}")
    for kind, field in CHECKED:
        a, b = fresh[kind][field], archived[kind][field]
        if a != b:
            bad.append(f"{design}/{index:03d} {kind}.{field} {a!r} != {b!r}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    if not AMENDMENT.exists():
        print(f"[abort] {AMENDMENT.name} missing: the grid extension is a "
              f"registration act and runs only once it is on record")
        return 2
    amend = json.loads(AMENDMENT.read_text(encoding="utf-8"))
    archive_sha_before = base.file_hash(ARCHIVE)
    if archive_sha_before != amend["archive_integrity"]["archived_sha256"]:
        print(f"[abort] {ARCHIVE.name} does not match the sha256 the amendment "
              f"recorded; the archive moved under this campaign")
        return 2
    print(f"[r28] amendment {amend['amendment_sha256'][:16]}, "
          f"archive {archive_sha_before[:16]} verified")

    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    check_key, new_key = tag_of(CHECK_BETA), tag_of(NEW_BETA)
    if new_key in archived["designs"]["A"]["rows"][0]["budgets"]:
        print(f"[abort] {new_key} already present in the archive")
        return 2

    betas = [NEW_BETA, CHECK_BETA]
    payload = {
        "schema": "r28_budget_pareto_extension_v1",
        "created_utc": base.utc_now(),
        "beta": NEW_BETA,
        "status": "post_hoc_declared_extension_of_the_frozen_phase_a_grid",
        "amendment": AMENDMENT.name,
        "amendment_sha256": amend["amendment_sha256"],
        "parent_record": ARCHIVE.name,
        "parent_record_sha256": archive_sha_before,
        "parent_preregistration_sha256": archived["preregistration_sha256"],
        "frozen_budget_grid": archived["budget_grid"],
        "reference_level": pareto.LEVEL,
        "computed_by": ("rev14_budget_pareto.worker, called verbatim; the only "
                        "new input is the value of beta"),
        "admissibility_check": {
            "recomputed_budget": CHECK_BETA,
            "requires": ("exact equality with the archived calibration on every "
                         "orbit of both designs"),
        },
        "source": base.provenance(),
        "designs": {},
    }

    all_bad: list[str] = []
    for d in ("A", "B"):
        rows = json.loads(
            pareto.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]
        if a.limit:
            rows = rows[:a.limit]
        tasks = [{"design": d, "row": r, "betas": betas} for r in rows]
        print(f"[r28] design {d}: {len(tasks)} orbits x beta "
              f"{{{NEW_BETA}, {CHECK_BETA}}} + original", flush=True)

        t0 = time.time()
        done, fails = [], []
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(pareto.worker, t): t for t in tasks}
            for n, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                if rec["status"] != "complete":
                    fails.append(rec)
                    print(f"  !! {rec['sobol_index']:03d} {rec['message']}",
                          flush=True)
                    continue
                done.append(rec)
                if n % 16 == 0 or n == len(tasks):
                    print(f"  [{n:3d}/{len(tasks)}] "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
        if fails:
            print(f"[abort] design {d}: {len(fails)} worker failures")
            return 1
        done.sort(key=lambda r: r["sobol_index"])

        arch_by_index = {int(r["sobol_index"]): r
                         for r in archived["designs"][d]["rows"]}
        bad = []
        for rec in done:
            i = int(rec["sobol_index"])
            bad += compare(rec["budgets"][check_key],
                           arch_by_index[i]["budgets"][check_key], d, i)
            bad += compare(rec["budgets"]["original"],
                           arch_by_index[i]["budgets"]["original"], d, i)
        all_bad += bad
        print(f"  reproduction of beta={CHECK_BETA} and the archived operating "
              f"point on {len(done)} orbits: "
              f"{'exact' if not bad else str(len(bad)) + ' MISMATCHES'}",
              flush=True)
        for line in bad[:10]:
            print(f"    {line}", flush=True)

        payload["designs"][d] = {
            "rows": done,
            "failures": fails,
            "summary": {k: pareto.summarize_budget(done, k)
                        for k in (new_key, check_key, "original")},
            "reproduction_mismatches": bad,
        }
        s = payload["designs"][d]["summary"][new_key]
        if s.get("orbits"):
            print(f"  {new_key}: n={s['orbits']} censored={s['censored']} "
                  f"R_a median={s['R_a_defect_rms']['median']:.3g} "
                  f"At wins {s['atallah_wins_defect']}/{s['orbits']} "
                  f"displ wins {s['atallah_wins_displacement_proxy']}/{s['orbits']}",
                  flush=True)

    if all_bad:
        print(f"\n[abort] {len(all_bad)} reproduction mismatches: a separate "
              f"record computed by the same code on the same inputs is not the "
              f"same object, so no beta = {NEW_BETA} record is written")
        return 1

    payload["admissibility_check"]["result"] = "exact_on_all_orbits_both_designs"
    base.atomic_json(OUTPUT, payload)

    archive_sha_after = base.file_hash(ARCHIVE)
    if archive_sha_after != archive_sha_before:
        print(f"[FAIL] {ARCHIVE.name} changed during the run "
              f"({archive_sha_before[:16]} -> {archive_sha_after[:16]})")
        return 1
    print(f"\n[written] {OUTPUT.name}")
    print(f"  archive untouched: {archive_sha_after[:16]}")
    print(f"  admissibility: exact reproduction of beta={CHECK_BETA} "
          f"and the archived operating point, both designs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
