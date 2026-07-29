"""Sensitivity of the forced-variational prediction to the gradient degree (R21).

The forced variational system of the manuscript,

    d/dt dr = dv,
    d/dt dv = G(t) dr + Delta_a_P(t),

evaluates the reference gravity gradient G at degree 120 rather than at the
adopted truth degree. The manuscript argues that the resulting error enters at
second order, because the omitted part of the gradient, dG, multiplies dr,
which is itself first-order small. That argument is only complete if dG is
shown to be small, which is what this script measures.

For every orbit of the eight-orbit variational panel it samples the archived
truth trajectory, evaluates G at degree 120 and at the orbit's adopted truth
degree by the same central differences the variational solve uses, and reports

  * the relative gradient truncation ||dG|| / ||G||, and
  * the neglected forcing ||dG||_2 * dr_rms measured against the first-order
    forcing Delta_a_rms that drives the same channel.

The second ratio is the quantity that decides whether the prediction can be
read quantitatively: it is the fraction of the forcing that the degree-120
gradient misrepresents.

Usage:
    python rev21_gradient_sensitivity.py [--epochs 241]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev10_sobol_confirmatory as base          # noqa: E402
import rev14_budget_trajectory as r14            # noqa: E402
from rev13_variational_check import _model, gradient, GRADIENT_DEGREE  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PANEL = METRICS / "r14_variational_budget.json"
PARETO = METRICS / "r14_budget_pareto.json"
OUT = METRICS / "r21_gradient_sensitivity.json"
TABLE = METRICS / "r21_gradient_sensitivity_table.tex"
BETA_KEY = "beta_1.00"


def defect_rms(pareto: dict, design: str, index: int) -> dict:
    for row in pareto["designs"][design]["rows"]:
        if int(row["sobol_index"]) != index:
            continue
        b = row["budgets"][BETA_KEY]
        out = {}
        for pol, key in (("atallah_budget", "atallah"), ("fixed_budget", "fixed")):
            d = b[key]["defect"]
            out[pol] = float(d.get("rms") or d.get("defect_rms_m_s2") or
                             d.get("rms_m_s2"))
        return out
    return {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=241)
    a = ap.parse_args()

    panel = json.loads(PANEL.read_text(encoding="utf-8"))
    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    rows = []
    t_start = time.time()

    for rec in panel["rows"]:
        design = rec["design"]
        index = int(rec["sobol_index"])
        adopted = int(rec["adopted_truth_degree"])
        gdeg = min(GRADIENT_DEGREE, adopted)
        _, raw = r14.reuse_paths(design, index, "truth", "tight")
        if not raw.exists():
            print(f"  !! missing truth npz for {design}{index:03d}", flush=True)
            continue
        d = np.load(raw)
        t_s, state = d["t_s"], d["state_si"]
        if state.shape[0] == 6:          # stored as (6, N)
            state = state.T
        step = max(1, len(t_s) // a.epochs)
        idx = np.arange(0, len(t_s), step)

        _, g_args = _model(gdeg)
        _, a_args = _model(adopted)

        rel_f, dg_2 = [], []
        for i in idx:
            r_vec = np.asarray(state[i][:3], dtype=float)
            t = float(t_s[i])
            g_lo = gradient(r_vec, t, gdeg, g_args)
            g_hi = gradient(r_vec, t, adopted, a_args)
            dg = g_hi - g_lo
            rel_f.append(np.linalg.norm(dg, "fro") / np.linalg.norm(g_lo, "fro"))
            dg_2.append(np.linalg.norm(dg, 2))

        rel_f = np.asarray(rel_f)
        dg_2 = np.asarray(dg_2)
        defects = defect_rms(pareto, design, index)
        pol_out = {}
        for pol in ("atallah_budget", "fixed_budget"):
            dr = float(rec["policies"][pol]["predicted_pos_rms_m"])
            da = defects.get(pol)
            pol_out[pol] = {
                "predicted_pos_rms_m": dr,
                "defect_rms_m_s2": da,
                # typical contribution: median gradient truncation over the arc
                "neglected_over_forcing_median": (float(np.median(dg_2)) * dr / da)
                                                 if da else None,
                # pessimistic pairing: worst gradient truncation anywhere on the
                # arc against the RMS displacement, which never co-occur
                "neglected_over_forcing_max": (float(dg_2.max()) * dr / da)
                                              if da else None}

        rows.append({
            "design": design, "sobol_index": index, "hp_km": rec["hp_km"],
            "adopted_truth_degree": adopted, "gradient_degree": gdeg,
            "epochs": int(len(idx)),
            "rel_frobenius": {"median": float(np.median(rel_f)),
                              "max": float(rel_f.max())},
            "dG_spectral_median_s2": float(np.median(dg_2)),
            "dG_spectral_max_s2": float(dg_2.max()),
            "policies": pol_out})
        print(f"  [{len(rows)}] {design}{index:03d} hp={rec['hp_km']:.0f} km "
              f"rel|dG|/|G| median={np.median(rel_f):.2e} max={rel_f.max():.2e} "
              f"elapsed={(time.time()-t_start)/60:.1f} min", flush=True)

    allrel = [r["rel_frobenius"]["max"] for r in rows]
    allratio = [p["neglected_over_forcing_max"]
                for r in rows for p in r["policies"].values()
                if p["neglected_over_forcing_max"] is not None]
    allratio_med = [p["neglected_over_forcing_median"]
                    for r in rows for p in r["policies"].values()
                    if p["neglected_over_forcing_median"] is not None]
    high = [p["neglected_over_forcing_max"]
            for r in rows if r["hp_km"] >= 50.0
            for p in r["policies"].values()
            if p["neglected_over_forcing_max"] is not None]
    payload = {
        "schema": "r21_gradient_sensitivity_v1",
        "created_utc": base.utc_now(),
        "description": "relative gravity-gradient truncation at the variational "
                       "gradient degree, and the forcing it misrepresents",
        "gradient_degree": GRADIENT_DEGREE,
        "fd_step_m": 1.0,
        "rows": rows,
        "summary": {
            "orbits": len(rows),
            "max_rel_frobenius": max(allrel) if allrel else None,
            "max_neglected_over_forcing": max(allratio) if allratio else None,
            "max_neglected_over_forcing_perilune_ge_50km": (max(high)
                                                           if high else None),
            "median_neglected_over_forcing_worstcase_pairing": (
                float(np.median(allratio)) if allratio else None),
            "median_neglected_over_forcing_typical": (
                float(np.median(allratio_med)) if allratio_med else None)}}
    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    def sci(x: float) -> str:
        s = f"{x:.1e}"
        mant, exp = s.split("e")
        return f"${mant}\\times 10^{{{int(exp)}}}$"

    lines = []
    for r in rows:
        typ = [p["neglected_over_forcing_median"] for p in r["policies"].values()
               if p["neglected_over_forcing_median"] is not None]
        wor = [p["neglected_over_forcing_max"] for p in r["policies"].values()
               if p["neglected_over_forcing_max"] is not None]
        lines.append(
            f"    {r['design']}{r['sobol_index']:03d} & {r['hp_km']:.0f} & "
            f"{r['adopted_truth_degree']} & "
            f"{sci(r['rel_frobenius']['median'])} & {sci(r['rel_frobenius']['max'])} & "
            f"{sci(max(typ))} & {sci(max(wor))} \\\\")
    TABLE.write_text(
        "% auto-generated by rev21_gradient_sensitivity.py -- do not edit\n"
        "\\begin{tabular}{@{}l r r r r r r@{}}\n"
        "\\toprule\n"
        "Orbit & $h_p$ [km] & $N_{\\mathrm{truth}}$ & "
        "\\multicolumn{2}{c}{$\\|\\delta G\\|_F/\\|G\\|_F$} & "
        "\\multicolumn{2}{c}{$\\|\\delta G\\|_2\\,\\delta r_{\\mathrm{rms}}/"
        "\\|\\Delta a\\|_{\\mathrm{rms}}$} \\\\\n"
        "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n"
        " & & & median & max & typical & worst pairing \\\\\n"
        "\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n",
        encoding="utf-8")
    print(f"[gradient-sensitivity] {len(rows)} orbits -> {OUT.name}", flush=True)
    print(json.dumps(payload["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
