"""R33: the variational panel again, with the gradient at the reference degree.

The forced variational reconstruction evaluates the reference gravity gradient
at degree 120 while the reference trajectory and the differenced accelerations
are at the orbit's adopted truth degree. The manuscript argues, and R21
measures, that the neglected part of the gradient enters the prediction at
second order and is small --- except on the two 31-km-perilune orbits of the
panel, where the pessimistic bound is not small at all.

Arguing that a term is negligible is weaker than not neglecting it. This runs
the same panel with the gradient evaluated at the reference degree, so the
approximation is removed rather than bounded, and the sign agreement can be
read without a caveat attached.

Nothing is reimplemented: rev13_variational_check runs verbatim, with its
gradient degree raised past every adopted truth degree so that its own
``min(GRADIENT_DEGREE, adopted)`` resolves to the reference degree on every
orbit. The patch is applied at module import, because ProcessPoolExecutor on
Windows re-imports this module in each worker to unpickle the callable below,
and a patch applied in main() would leave the children computing the degree-120
gradient while the parent reported otherwise --- the failure mode that would be
hardest to notice and worst to publish.

The archived degree-120 record is not touched; this writes its own, so the two
can be compared orbit by orbit.

Cost, measured from the archived run: the gradient contributes six field
evaluations per right-hand-side call, so raising it from 120 to 300 roughly
doubles the per-call work and to 900 roughly multiplies it by two and a half.
The two 900-degree orbits dominate the wall clock and are the two the exercise
exists for.

Usage:
    python rev33_variational_reference_gradient.py --workers 8
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev13_variational_check as vc

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

ARCHIVE = METRICS / "r13_variational_check.json"
OUTPUT = METRICS / "r33_variational_reference_gradient.json"

# Above every adopted truth degree in the archive (300, 600, 900), so that
# rev13's own cap resolves the gradient to the reference degree on every orbit.
vc.GRADIENT_DEGREE = 10_000


def worker(task: dict) -> dict:
    """rev13's worker, reached through this module so the child applies the
    gradient-degree patch above before running it."""
    return vc.worker(task)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--design", default="A")
    a = ap.parse_args()

    if not ARCHIVE.exists():
        print(f"[abort] {ARCHIVE.name} missing; this compares against it")
        return 2
    arch = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    if int(arch["gradient_degree"]) != 120:
        print(f"[abort] the archived record was not the degree-120 run")
        return 2
    tasks = [{"design": r["design"], "index": int(r["sobol_index"])}
             for r in arch["rows"]]
    print(f"[r33] {len(tasks)} orbits, gradient at the reference degree "
          f"(archive used 120)", flush=True)

    results, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec.get('index')} {rec.get('message')}", flush=True)
                continue
            rec["measured_tight_m"] = vc.measured(rec["design"],
                                                  rec["sobol_index"])
            m, p = rec["measured_tight_m"], rec["policies"]
            rec["calibration"] = {
                "critical_predicted_over_measured": (
                    p["fixed_critical"]["predicted_pos_rms_m"]
                    / m["fixed_critical"]),
                "predicted_ratio_fixed_over_atallah": (
                    p["fixed_work"]["predicted_pos_rms_m"]
                    / p["atallah"]["predicted_pos_rms_m"]),
                "predicted_gap_m": abs(p["fixed_work"]["predicted_pos_rms_m"]
                                       - p["atallah"]["predicted_pos_rms_m"]),
                "measured_threshold_m": (m["envelope_atallah"]
                                         + m["envelope_fixed_work"]),
            }
            rec["calibration"]["predicted_gap_over_threshold"] = (
                rec["calibration"]["predicted_gap_m"]
                / rec["calibration"]["measured_threshold_m"])
            results.append(rec)
            print(f"  [{n}/{len(tasks)}] idx={rec['sobol_index']:03d} "
                  f"N_ref={rec['adopted_truth_degree']} "
                  f"ratio={rec['calibration']['predicted_ratio_fixed_over_atallah']:.3f} "
                  f"elapsed={(time.time()-t0)/60:.1f} min", flush=True)
    results.sort(key=lambda r: (r["design"], r["sobol_index"]))

    # orbit-by-orbit against the archived degree-120 prediction
    by_index = {(r["design"], int(r["sobol_index"])): r for r in arch["rows"]}
    comparison, flips = [], 0
    for r in results:
        old = by_index.get((r["design"], int(r["sobol_index"])))
        if not old:
            continue
        new_ratio = r["calibration"]["predicted_ratio_fixed_over_atallah"]
        old_ratio = old["calibration"]["predicted_ratio_fixed_over_atallah"]
        flip = (new_ratio > 1.0) != (old_ratio > 1.0)
        flips += bool(flip)
        comparison.append({
            "design": r["design"], "sobol_index": r["sobol_index"],
            "adopted_truth_degree": r["adopted_truth_degree"],
            "perilune_km": None,
            "ratio_gradient_120": old_ratio,
            "ratio_gradient_reference": new_ratio,
            "relative_change": new_ratio / old_ratio - 1.0,
            "sign_flipped": bool(flip),
        })

    payload = {
        "schema": "r33_variational_reference_gradient_v1",
        "created_utc": base.utc_now(),
        "gradient_degree": "adopted truth degree of each orbit",
        "archive_compared": ARCHIVE.name,
        "archive_gradient_degree": arch["gradient_degree"],
        "purpose": ("removes the degree-120 gradient approximation from the "
                    "forced variational reconstruction rather than bounding it"),
        "rows": results, "failures": fails,
        "comparison_with_degree_120": comparison,
        "sign_flips": flips,
        "source": base.provenance(),
    }
    if comparison:
        ch = [abs(c["relative_change"]) for c in comparison]
        payload["summary"] = {
            "orbits": len(comparison),
            "sign_flips": flips,
            "abs_relative_change": {
                "median": float(np.median(ch)),
                "max": float(max(ch)),
            },
            "predicted_ratio_above_one": sum(
                c["ratio_gradient_reference"] > 1.0 for c in comparison),
        }
    base.atomic_json(OUTPUT, payload)
    print(f"\n[written] {OUTPUT.name}: {len(results)} orbits, "
          f"{flips} sign flips against the degree-120 run")
    for c in comparison:
        print(f"  {c['design']}{c['sobol_index']:03d} N_ref="
              f"{c['adopted_truth_degree']:3d}  "
              f"{c['ratio_gradient_120']:.3f} -> "
              f"{c['ratio_gradient_reference']:.3f}  "
              f"({c['relative_change']:+.1%})"
              f"{'  SIGN FLIP' if c['sign_flipped'] else ''}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
