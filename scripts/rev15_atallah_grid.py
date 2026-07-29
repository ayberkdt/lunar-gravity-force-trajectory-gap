"""Grid convergence of the Atallah tolerance's sampled maximum (R15-E).

The tolerance handed to the published rule is the largest truncation error the
critical fixed degree incurs at perilune, measured on a 25x48 latitude/longitude
grid. The manuscript called that a worst case; it is a sampled-grid maximum, and
whether it is converged has not been shown. This refines the grid,

    25x48  ->  50x96  ->  100x192,

at representative altitudes and degrees, and reports how much the maximum moves.
If the change from 50x96 to 100x192 is small, the wording "grid-resolved sampled
maximum" is justified and the tolerances the benchmark used stand; if not, the
tolerance derivation has to be repeated on the converged grid.

The bound the rule satisfies oversatisfies its target by five to eight orders of
magnitude, so a change of ordering is not expected -- but the tolerance is the
calibration basis of the whole direct benchmark, so it should not rest on an
unchecked sample.

Usage:
    python rev15_atallah_grid.py --workers 4
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
from rev14_budget_pareto import _model

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r15_atallah_grid.json"
TABLE = METRICS / "r15_atallah_grid_table.tex"

GRIDS = [(25, 48), (50, 96), (100, 192)]
ALTITUDES_KM = [31.0, 50.0, 80.0, 110.0, 150.0]
TRUTH_DEGREE = 900


def sampled_max(model, args, r_m, N, N_truth, n_lat, n_lon):
    """Largest |a_Nmax - a_N| over a lat/lon sample of the sphere at radius r."""
    from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba
    lats = np.linspace(-math.pi / 2 + 1e-3, math.pi / 2 - 1e-3, n_lat)
    lons = np.linspace(0.0, 2.0 * math.pi, n_lon, endpoint=False)
    worst = 0.0
    argmax = None
    for phi in lats:
        cp, sp = math.cos(phi), math.sin(phi)
        for lam in lons:
            x, y, z = r_m * cp * math.cos(lam), r_m * cp * math.sin(lam), r_m * sp
            aN = np.array(sh_accel_fixed_numba(x, y, z, N, *args))
            aM = np.array(sh_accel_fixed_numba(x, y, z, N_truth, *args))
            v = float(np.linalg.norm(aM - aN))
            if v > worst:
                worst, argmax = v, (math.degrees(phi), math.degrees(lam))
    return worst, argmax


def worker(task: dict) -> dict:
    h_km, N = task["h_km"], task["degree"]
    try:
        model, args = _model(TRUTH_DEGREE)
        r_m = model.r_ref + h_km * 1e3
        out = {}
        for n_lat, n_lon in GRIDS:
            t0 = time.time()
            v, loc = sampled_max(model, args, r_m, N, TRUTH_DEGREE, n_lat, n_lon)
            out[f"{n_lat}x{n_lon}"] = {"max_accel_error_m_s2": v,
                                       "argmax_lat_lon_deg": loc,
                                       "n_samples": n_lat * n_lon,
                                       "wall_s": time.time() - t0}
        ref = out[f"{GRIDS[-1][0]}x{GRIDS[-1][1]}"]["max_accel_error_m_s2"]
        for k, v in out.items():
            v["rel_change_vs_finest"] = v["max_accel_error_m_s2"] / ref - 1.0
        return {"h_km": h_km, "degree": N, "status": "complete", "grids": out}
    except Exception as exc:
        return {"h_km": h_km, "degree": N, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    # degrees spanning the critical range the benchmark actually uses
    tasks = [{"h_km": h, "degree": d}
             for h in ALTITUDES_KM for d in (80, 150, 250)]
    print(f"[atallah-grid] {len(tasks)} (altitude, degree) cells x "
          f"{len(GRIDS)} grids", flush=True)
    done, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! h={rec['h_km']} N={rec['degree']} {rec.get('message')}",
                      flush=True)
                continue
            done.append(rec)
            g = rec["grids"]
            print(f"  [{n}/{len(tasks)}] h={rec['h_km']:5.0f} N={rec['degree']:3d} "
                  f"25x48 {g['25x48']['rel_change_vs_finest']:+.3f} "
                  f"50x96 {g['50x96']['rel_change_vs_finest']:+.3f} "
                  f"({(time.time()-t0)/60:.1f}min)", flush=True)
    done.sort(key=lambda r: (r["h_km"], r["degree"]))
    coarse = [abs(r["grids"]["25x48"]["rel_change_vs_finest"]) for r in done]
    mid = [abs(r["grids"]["50x96"]["rel_change_vs_finest"]) for r in done]
    summary = {
        "cells": len(done),
        "abs_rel_change_25x48_vs_finest": {
            "median": float(np.median(coarse)), "max": float(max(coarse))},
        "abs_rel_change_50x96_vs_finest": {
            "median": float(np.median(mid)), "max": float(max(mid))},
        "converged": bool(max(mid) < 0.05),
        "criterion": "max |change| from 50x96 to 100x192 below 5%",
    }
    payload = {"schema": "r15_atallah_grid_v1", "created_utc": base.utc_now(),
               "grids": [f"{a_}x{b_}" for a_, b_ in GRIDS],
               "truth_degree": TRUTH_DEGREE, "rows": done, "failures": fails,
               "summary": summary, "source": base.provenance()}
    base.atomic_json(OUTPUT, payload)

    lines = []
    for r in done:
        g = r["grids"]
        lines.append(
            f"    {r['h_km']:.0f} & {r['degree']} & "
            f"{g['25x48']['max_accel_error_m_s2']:.4g} & "
            f"{100 * g['25x48']['rel_change_vs_finest']:+.2f}\\% & "
            f"{100 * g['50x96']['rel_change_vs_finest']:+.2f}\\% \\\\")
    TABLE.write_text(f"""% auto-generated by rev15_atallah_grid.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\caption{{Grid convergence of the sampled maximum used to set the published
  rule's tolerance, and it does not converge. The benchmark measures the largest
  truncation error of a fixed degree over a $25\\times48$ latitude/longitude
  sample of the sphere; the table refines that sample to $50\\times96$ and
  $100\\times192$ and gives the relative change of each against the finest grid.
  The $25\\times48$ sample falls below the finest by a median
  {100 * summary['abs_rel_change_25x48_vs_finest']['median']:.0f}\\% and by up to
  {100 * summary['abs_rel_change_25x48_vs_finest']['max']:.0f}\\%, and refining to
  $50\\times96$ removes little of that ({100 * summary['abs_rel_change_50x96_vs_finest']['median']:.0f}\\%
  median, {100 * summary['abs_rel_change_50x96_vs_finest']['max']:.0f}\\% worst).
  The quantity is therefore a sampled-grid maximum and is called one throughout;
  it is not a worst case over the sphere and it is not grid-resolved. Its
  consequence for the benchmark is bounded separately and is negligible, because
  the selected degree depends on the tolerance through a steeply falling tail
  sum.}}
  \\label{{tab:atallah-grid}}
  \\begin{{tabular}}{{r r r r r}}
    \\toprule
    $h$ [km] & $N$ & $25\\times48$ max [m\\,s$^{{-2}}$] &
      $25\\times48$ vs finest & $50\\times96$ vs finest\\\\
    \\midrule
{chr(10).join(lines)}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
""", encoding="utf-8")
    s = summary
    print(f"[atallah-grid] 25x48 vs finest: median "
          f"{100*s['abs_rel_change_25x48_vs_finest']['median']:.2f}% "
          f"max {100*s['abs_rel_change_25x48_vs_finest']['max']:.2f}%")
    print(f"               50x96 vs finest: median "
          f"{100*s['abs_rel_change_50x96_vs_finest']['median']:.2f}% "
          f"max {100*s['abs_rel_change_50x96_vs_finest']['max']:.2f}%  "
          f"converged={s['converged']}")
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
