"""Reproducible verification of the Atallah (2022) implementation on JGGRX (R12).

The source paper reports its own acceptance evidence on GL0660B and EGM2008,
neither of which is used here, so the implementation is verified against the
paper's own criteria on the field this study actually propagates:

  V1  Eq.(29) bounding property: the analytical per-degree term a_hat_n(r)
      dominates the measured degree-n acceleration contribution of the real
      field at every tested radius and degree.
  V2  Eq.(31)/Fig.4 selection property: at N_req(r,tol) the measured worst-case
      truncation error of the real field is below the requested tolerance, over
      a radius x tolerance matrix.
  V3  monotonicity: N_req decreases with radius and increases with tightening
      tolerance.
  V4  the selection curve N_req(r,tol) itself, tabulated for reporting.

Worst cases are taken over a latitude/longitude grid on the sphere of the given
radius, which is the quantity the analytical bound is written for.

Usage:
    python rev12_atallah_verification.py --degree 300
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at

METRICS = Path(__file__).resolve().parents[1] / "metrics"
OUTPUT = METRICS / "r12_atallah_verification.json"
TABLE = METRICS / "r12_atallah_verification_table.tex"

ALT_KM = [30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0]
TOLS = [1e-6, 1e-8, 1e-10, 1e-12]
DEGREE_LADDER = [2, 5, 10, 20, 40, 60, 90, 120, 160, 200, 240, 280, 300,
                 400, 500, 600]
GRID = {"n_lat": 13, "n_lon": 24}


def _sci(v, digits=1):
    if v == 0.0:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    return f"{v / 10.0**e:.{digits}f}\\times 10^{{{e}}}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--degree", type=int, default=300,
                    help="model degree used for the verification field")
    a = ap.parse_args()
    model = base.load_model(a.degree)
    args = base.kernel_args(model)
    base.warmup(model, args)
    g = at.precompute_Sn(model, a.degree)
    payload = {"schema": "r12_atallah_verification_v1",
               "created_utc": base.utc_now(),
               "field": "JGGRX_1800F", "model_degree": a.degree,
               "grid": GRID, "source": base.provenance()}

    # V1 -- Eq.(29) bounding property
    v1 = []
    for h in (50.0, 100.0, 300.0):
        r = model.r_ref + h * 1e3
        a_hat = at.a_hat_series(r, model, g)
        for n in DEGREE_LADDER:
            if n < 2 or n > a.degree:
                continue
            actual = at.actual_degree_contribution_max(model, args, r, n, **GRID)
            v1.append({"altitude_km": h, "degree": n,
                       "measured_max": actual, "bound": float(a_hat[n]),
                       "ratio": actual / float(a_hat[n])})
            print(f"  V1 h={h:5.0f} n={n:3d} actual/bound={v1[-1]['ratio']:.4f}",
                  flush=True)
    payload["bounding"] = {
        "checks": v1, "violations": sum(x["ratio"] > 1.0 for x in v1),
        "tightest_ratio": max(x["ratio"] for x in v1)}

    # V2 -- selection property over a radius x tolerance matrix
    v2 = []
    for h in ALT_KM:
        r = model.r_ref + h * 1e3
        for tol in TOLS:
            n = at.n_req(r, tol, model, g, floor=2, cap=a.degree)
            err = at.actual_truncation_error_max(model, args, r, n, a.degree, **GRID)
            v2.append({"altitude_km": h, "tol": tol, "n_req": int(n),
                       "measured_error": err, "satisfied": bool(err < tol),
                       "saturated": bool(n >= a.degree)})
            print(f"  V2 h={h:5.0f} tol={tol:.0e} N={n:3d} err={err:.2e} "
                  f"{'OK' if err < tol else 'FAIL'}", flush=True)
    payload["selection"] = {
        "matrix": v2, "failures": sum(not x["satisfied"] for x in v2),
        "note": ("cells where N_req reaches the model degree are saturated: the "
                 "measured error is then zero because the comparison field is "
                 "the model itself")}

    # V3 -- monotonicity
    r_list = [model.r_ref + h * 1e3 for h in ALT_KM]
    n_vs_r = [at.n_req(r, 1e-9, model, g, floor=2, cap=a.degree) for r in r_list]
    n_vs_tol = [at.n_req(model.r_ref + 100e3, t, model, g, floor=2, cap=a.degree)
                for t in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)]
    payload["monotonicity"] = {
        "altitudes_km": ALT_KM, "n_req_at_tol_1e-9": n_vs_r,
        "monotone_decreasing_in_radius": bool(
            all(x >= y for x, y in zip(n_vs_r, n_vs_r[1:]))),
        "tolerances": [1e-4, 1e-6, 1e-8, 1e-10, 1e-12],
        "n_req_at_100km": n_vs_tol,
        "monotone_increasing_with_tightening": bool(
            all(x <= y for x, y in zip(n_vs_tol, n_vs_tol[1:])))}

    # V4 -- selection curve
    alt_curve = list(np.arange(20.0, 505.0, 20.0))
    payload["selection_curve"] = {
        "altitudes_km": alt_curve,
        "n_req": {f"{t:.0e}": [at.n_req(model.r_ref + h * 1e3, t, model, g,
                                        floor=2, cap=a.degree)
                               for h in alt_curve] for t in TOLS}}

    base.atomic_json(OUTPUT, payload)

    # LaTeX: the radius x tolerance matrix, which is the compact statement of V2
    head = " & ".join(f"$10^{{{int(math.log10(t))}}}$" for t in TOLS)
    body = []
    for h in ALT_KM:
        cells = []
        for tol in TOLS:
            c = next(x for x in v2 if x["altitude_km"] == h and x["tol"] == tol)
            mark = "$^{\\dagger}$" if c["saturated"] else ""
            err = "0" if c["measured_error"] == 0.0 else f"${_sci(c['measured_error'])}$"
            cells.append(f"{c['n_req']}{mark} / {err}")
        body.append(f"    {h:.0f} & " + " & ".join(cells) + "\\\\")
    tightest = payload["bounding"]["tightest_ratio"]
    tex = f"""% auto-generated by rev12_atallah_verification.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\caption{{Verification of the Atallah selection rule on JGGRX\\_1800F
  (degree {a.degree} verification field), against the acceptance criteria of the
  source paper. Each cell gives the selected degree $N_{{\\mathrm{{req}}}}(r,
  \\varepsilon_a)$ and the measured worst-case acceleration truncation error of
  the real field at that degree, over a $13\\times24$ latitude/longitude grid on
  the sphere of that altitude. The requirement is that the measured error is
  below the requested tolerance in every cell; it is met in
  {len(v2) - payload['selection']['failures']} of {len(v2)} cells and violated in
  {payload['selection']['failures']}. Cells marked $\\dagger$ are saturated at the
  verification-field degree, where the measured error is identically zero. The
  Eq.~(29) bounding property holds with
  {payload['bounding']['violations']} violations over
  {len(v1)} radius/degree checks, the tightest measured-to-bound ratio being
  {tightest:.3f}.}}
  \\label{{tab:atallah-verification}}
  \\begin{{tabular}}{{r cccc}}
    \\toprule
    & \\multicolumn{{4}}{{c}}{{$\\varepsilon_a$ [m\\,s$^{{-2}}$]: $N_{{\\mathrm{{req}}}}$ / measured error [m\\,s$^{{-2}}$]}}\\\\
    \\cmidrule(lr){{2-5}}
    $h$ [km] & {head}\\\\
    \\midrule
{chr(10).join(body)}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""
    TABLE.write_text(tex, encoding="utf-8")
    print(f"[written] {OUTPUT.name}, {TABLE.name}")
    print(json.dumps({k: payload[k] for k in ("monotonicity",)}, indent=2))
    print(f"  bounding violations: {payload['bounding']['violations']}, "
          f"selection failures: {payload['selection']['failures']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
