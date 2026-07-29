"""Extension of the measured kernel cost curve above degree 300 (R13).

The time-matched comparator of rev13_timing_match.py needs c(N) up to the
degrees the Atallah rule actually uses at low perilune, which exceed the degree
300 ceiling of the archived curve (``r12_kernel_cost_curve.json``). This script
measures the same quantity, single-threaded and after per-degree warm-up, on the
degree ladder 320--900, and writes a combined curve. It must run on an idle
machine; the archived low-degree points are re-measured as a consistency check.

Usage:
    python rev13_cost_curve_high.py --repeats 200
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from rev3_common import load_model, kernel_args, warmup
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

METRICS = Path(__file__).resolve().parents[1] / "metrics"
BASE_CURVE = METRICS / "r12_kernel_cost_curve.json"
OUTPUT = METRICS / "r13_kernel_cost_curve_high.json"

HIGH_DEGREES = [320, 360, 400, 450, 500, 550, 600, 700, 800, 900]
RECHECK = [100, 200, 300]


def measure(args, degree: int, repeats: int, r_m: float) -> dict:
    x, y, z = r_m, 0.0, 0.0
    for _ in range(50):
        sh_accel_fixed_numba(x, y, z, degree, *args)
    samples = np.empty(repeats)
    for i in range(repeats):
        t0 = time.perf_counter_ns()
        sh_accel_fixed_numba(x, y, z, degree, *args)
        samples[i] = time.perf_counter_ns() - t0
    return {"degree": degree,
            "per_call_ns_median": float(np.median(samples)),
            "per_call_ns_p10": float(np.percentile(samples, 10)),
            "per_call_ns_p90": float(np.percentile(samples, 90)),
            "per_call_ns_min": float(samples.min())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=200)
    a = ap.parse_args()
    model = load_model(900)
    args = kernel_args(model)
    warmup(model, args)
    r_m = model.r_ref + 100e3

    rows = []
    for n in RECHECK + HIGH_DEGREES:
        row = measure(args, n, a.repeats, r_m)
        rows.append(row)
        print(f"  N={n:4d}  {row['per_call_ns_median'] / 1000:8.1f} us", flush=True)

    base = json.loads(BASE_CURVE.read_text())
    old = {r["degree"]: r["per_call_ns_median"] for r in base["rows"]}
    check = [{"degree": n,
              "archived_ns": old[n],
              "remeasured_ns": next(r["per_call_ns_median"] for r in rows
                                    if r["degree"] == n),
              "ratio": next(r["per_call_ns_median"] for r in rows
                            if r["degree"] == n) / old[n]}
             for n in RECHECK if n in old]
    combined = sorted(
        {r["degree"]: r for r in base["rows"] + [r for r in rows
                                                 if r["degree"] not in RECHECK]}.values(),
        key=lambda r: r["degree"])
    deg = np.array([r["degree"] for r in combined], float)
    ns = np.array([r["per_call_ns_median"] for r in combined], float)
    k = float(np.sum(ns * deg ** 2) / np.sum(deg ** 4))
    resid = float(np.sqrt(np.mean((ns - k * deg ** 2) ** 2)) / np.mean(ns))
    payload = {"schema": "r13_kernel_cost_curve_high_v1",
               "machine_single_threaded": True, "repeats": a.repeats,
               "eval_radius_m": r_m, "model_degree": 900,
               "high_rows": [r for r in rows if r["degree"] not in RECHECK],
               "repeatability_check": check,
               "combined_rows": combined,
               "pure_N2_fit_combined": {"k": k, "rel_rms_residual": resid},
               "ratio_900_over_300": float(
                   next(r["per_call_ns_median"] for r in rows if r["degree"] == 900)
                   / next(r["per_call_ns_median"] for r in rows if r["degree"] == 300)),
               "quadratic_expectation_900_over_300": 9.0}
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUTPUT.name}")
    for c in check:
        print(f"  recheck N={c['degree']}: archived {c['archived_ns']/1000:.1f} us, "
              f"remeasured {c['remeasured_ns']/1000:.1f} us, ratio {c['ratio']:.3f}")
    print(f"  c(900)/c(300) = {payload['ratio_900_over_300']:.2f} "
          f"(pure quadratic would be 9.00)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
