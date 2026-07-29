"""Measured gravity-kernel cost curve c(N) (R12).

No-regret calibration that supports both candidate overnight tasks:
  * the direct-benchmark timing panel (measured, not proxy, seconds), and
  * the measured-gravity-time-matched fixed comparator, which replaces the
    quadratic ``<N^2>`` work proxy the manuscript currently admits as a
    limitation.

It measures the wall time of a single fixed-degree gravity evaluation at each
degree N on this machine, single-threaded, with warm Numba kernels and many
repeats, and reports the median per-call cost plus a fit of c(N) against the
N^2 proxy so the proxy can be validated or corrected.

Run single-threaded and with no other Python load for meaningful numbers:
    python rev12_kernel_cost_curve.py --repeats 400
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from rev3_common import load_model, kernel_args, warmup, OUT
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

DEGREES = list(range(20, 261, 10)) + [280, 300]


def measure(model, args, degree: int, repeats: int, r_m: float) -> dict:
    # a fixed body-fixed evaluation point at the requested radius
    x, y, z = r_m, 0.0, 0.0
    # warm the specific degree
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
    ap.add_argument("--repeats", type=int, default=400)
    a = ap.parse_args()

    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    r_m = model.r_ref + 100e3

    rows = []
    for n in DEGREES:
        row = measure(model, args, n, a.repeats, r_m)
        rows.append(row)
        print(f"  N={n:4d}  {row['per_call_ns_median']/1e3:8.2f} us/call "
              f"(min {row['per_call_ns_min']/1e3:.2f})", flush=True)

    # Fit c(N) = a*N^2 + b*N + c and compare with the pure-N^2 proxy.
    N = np.array([r["degree"] for r in rows], dtype=float)
    c = np.array([r["per_call_ns_median"] for r in rows]) / 1e9  # seconds
    A = np.vstack([N * N, N, np.ones_like(N)]).T
    coef, *_ = np.linalg.lstsq(A, c, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((c - pred) ** 2))
    ss_tot = float(np.sum((c - c.mean()) ** 2))
    # pure quadratic-through-scale proxy c ~ k*N^2
    k = float(np.sum(c * N * N) / np.sum(N ** 4))
    quad_pred = k * N * N
    quad_rel_resid = float(np.sqrt(np.mean(((c - quad_pred) / c) ** 2)))

    payload = {
        "schema": "r12_kernel_cost_curve_v1",
        "machine_single_threaded": True,
        "repeats": a.repeats,
        "eval_radius_m": r_m,
        "degrees": DEGREES,
        "rows": rows,
        "fit_quadratic_full": {"a_N2": coef[0], "b_N": coef[1], "c_0": coef[2],
                               "r2": 1.0 - ss_res / ss_tot},
        "fit_pure_N2_proxy": {"k": k,
                              "rel_rms_residual": quad_rel_resid,
                              "interpretation": (
                                  "how well the <N^2> work proxy tracks the "
                                  "measured per-call cost; small => proxy fair")},
    }
    (OUT / "r12_kernel_cost_curve.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] r12_kernel_cost_curve.json  "
          f"pure-N^2 proxy rel-RMS residual = {quad_rel_resid*100:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
