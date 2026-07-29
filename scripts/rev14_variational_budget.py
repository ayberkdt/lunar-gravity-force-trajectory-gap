"""O25 Phase E: forced-variational mechanism check at equal budget (R14).

The force-level sweep says the budget-calibrated radial history has the smaller
truncation defect at beta = 1. Whether that survives dynamical filtering is a
separate question, because long-arc displacement is the state-transition-weighted
integral of the defect and not its norm. This script answers it with the
manuscript's own forced variational system,

    d/dt dr = dv,
    d/dt dv = G(t) dr + Delta_a_P(t),      Delta_a_P = a_{N_P} - a_{N_truth},

integrated along one shared reference for both policies at once, so a single
augmented integration gives both linear predictions under an identical gradient.

At beta = 1 the equal-budget comparator degree equals the critical-altitude
degree on every orbit, so the comparator is also the calibration channel: its
measured seven-day error is large and tolerance-stable, and agreement between its
linear prediction and its measured error certifies the prediction for the radial
policy, whose measured error is closer to the numerical floor.

Usage:
    python rev14_variational_budget.py --orbits 8 --workers 4
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
from scipy.integrate import solve_ivp

import rev10_sobol_confirmatory as base
import rev12_atallah as at
from rev13_variational_check import (accel_inertial, gradient, binned_lookup,
                                     ric_axes, _model, GRADIENT_DEGREE,
                                     RTOL, ATOL, MAX_STEP_S, DURATION,
                                     OUTPUT_STEP)
from rev14_budget_pareto import BIN_KM, FLOOR

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PARETO = METRICS / "r14_budget_pareto.json"
OUTPUT = METRICS / "r14_variational_budget.json"
TABLE = METRICS / "r14_variational_budget_table.tex"
BETA = 1.00
ROWS = {"A": METRICS / "r10_sobolA_baseline_truth_corrected.json",
        "B": METRICS / "r11_designB_rows.json"}


def worker(task: dict) -> dict:
    design, row, spec = task["design"], task["row"], task["spec"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        n_fixed = int(spec["fixed"]["degree"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        gdeg = min(GRADIENT_DEGREE, adopted)
        _, gargs = _model(gdeg)
        g = at.precompute_Sn(model, adopted)
        _, table = at.atallah_binned_schedule(
            model, g, float(spec["atallah"]["tol_accel_m_s2"]), hp_km, ha_km,
            floor=FLOOR, cap=adopted, bin_km=BIN_KM)
        deg_at = binned_lookup({str(k): int(v) for k, v in table.items()},
                               bin_km=BIN_KM)
        r_ref = model.r_ref
        policies = ("atallah_budget", "fixed_budget")

        def rhs(t, Y):
            r, v = Y[0:3], Y[3:6]
            a_ref = accel_inertial(r, t, adopted, args)
            G = gradient(r, t, gdeg, gargs)
            alt = float(np.linalg.norm(r)) - r_ref
            degs = (deg_at(alt), n_fixed)
            dY = np.empty_like(Y)
            dY[0:3] = v
            dY[3:6] = a_ref
            for k, n in enumerate(degs):
                off = 6 + 6 * k
                da = accel_inertial(r, t, n, args) - a_ref
                dY[off:off + 3] = Y[off + 3:off + 6]
                dY[off + 3:off + 6] = G @ Y[off:off + 3] + da
            return dY

        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        Y0 = np.zeros(6 + 6 * len(policies))
        Y0[0:6] = y0
        t0 = time.time()
        sol = solve_ivp(rhs, (0.0, DURATION), Y0, method="DOP853", t_eval=grid,
                        rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S)
        if not sol.success:
            return {"design": design, "sobol_index": index,
                    "status": "integration_failed", "message": sol.message}
        ref = sol.y[0:6]
        out = {"design": design, "sobol_index": index, "status": "complete",
               "adopted_truth_degree": adopted, "n_critical": n_crit,
               "n_fixed": n_fixed, "hp_km": hp_km,
               "comparator_is_critical": bool(n_fixed == n_crit),
               "wall_s": time.time() - t0, "n_rhs": int(sol.nfev), "policies": {}}
        for k, name in enumerate(policies):
            off = 6 + 6 * k
            d = sol.y[off:off + 3]
            ric = np.empty((d.shape[1], 3))
            for j in range(d.shape[1]):
                axes = ric_axes(ref[0:3, j], ref[3:6, j])
                ric[j] = [ax @ d[:, j] for ax in axes]
            mag = np.linalg.norm(d, axis=0)
            out["policies"][name] = {
                "predicted_pos_rms_m": float(np.sqrt(np.mean(mag ** 2))),
                "predicted_pos_final_m": float(mag[-1]),
                "predicted_radial_final_m": float(ric[-1, 0]),
                "predicted_in_track_final_m": float(ric[-1, 1]),
                "predicted_cross_track_final_m": float(ric[-1, 2]),
                "predicted_in_track_rms_m": float(np.sqrt(np.mean(ric[:, 1] ** 2))),
                "predicted_radial_rms_m": float(np.sqrt(np.mean(ric[:, 0] ** 2))),
            }
        pa = out["policies"]["atallah_budget"]["predicted_pos_rms_m"]
        pf = out["policies"]["fixed_budget"]["predicted_pos_rms_m"]
        out["predicted_ratio_fixed_over_atallah"] = (pf / pa) if pa > 0 else None
        return out
    except Exception as exc:
        return {"design": design, "sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def measured(design: str, index: int) -> dict:
    p = METRICS / f"r14_trajectory_{design}_beta_{BETA:.2f}.json"
    if not p.exists():
        return {}
    for r in json.loads(p.read_text(encoding="utf-8"))["rows"]:
        if int(r["sobol_index"]) == index:
            return {
                "atallah_budget": r["policies"]["atallah_budget"]["error_tight"]["pos_rms_m"],
                "fixed_budget": r["policies"]["fixed_budget"]["error_tight"]["pos_rms_m"],
                "envelope_atallah": r["policies"]["atallah_budget"]["truth_inclusive_envelope_m"],
                "envelope_fixed": r["policies"]["fixed_budget"]["truth_inclusive_envelope_m"],
                "rho_measured": r["comparison"]["rho_budget"],
                "resolved": r["comparison"]["resolved"]}
    return {}


def stat(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)], float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max())}


def sig2(v) -> str:
    """Two significant figures in LaTeX math, never rounded to a bare 0.00.

    A ratio column that prints 0.00 against 0.00 cannot be checked by a reader,
    and the whole point of the measured column is that it can be.
    """
    if v is None or not np.isfinite(v):
        return "--"
    if v == 0:
        return "$0$"
    exp = int(np.floor(np.log10(abs(v))))
    if -2 <= exp <= 2:
        return f"${v:.{max(0, 1 - exp)}f}$"
    return f"${v / 10.0 ** exp:.1f}\\times 10^{{{exp}}}$"


def build_table(rows) -> str:
    lines = []
    for r in rows:
        m = r.get("measured", {})
        pa = r["policies"]["atallah_budget"]
        pf = r["policies"]["fixed_budget"]
        cal = r.get("calibration_ratio_fixed")
        lines.append(
            f"    {r['design']} & {r['sobol_index']:03d} & {r['hp_km']:.0f} & "
            f"{r['n_fixed']} & {pa['predicted_pos_rms_m']:.3f} & "
            f"{pf['predicted_pos_rms_m']:.3f} & "
            f"{sig2(r['predicted_ratio_fixed_over_atallah'])} & "
            f"{m.get('atallah_budget', float('nan')):.3f} & "
            f"{m.get('fixed_budget', float('nan')):.3f} & "
            f"{cal:.3f}".replace("nan", "--") + " & "
            f"{sig2(m.get('rho_measured'))}"
            + "\\\\")
    body = "\n".join(lines)
    return f"""% auto-generated by rev14_variational_budget.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Forced-variational mechanism check at equal budget. Both policies are
  carried along one shared reference under one gravity gradient, so the
  predictions differ only through their truncation defects. At $\\beta = 1$ the
  comparator degree equals the critical-altitude degree, whose measured error is
  large and tolerance-stable, so it doubles as the calibration channel: the
  calibration ratio is its predicted over its measured seven-day RMS error.
  $\\rho$ columns are fixed over radial, so values above unity favor the radial
  allocation; both are given to two significant figures so that the predicted
  and measured ratios can be compared where they are small. Errors are in
  metres.}}
  \\label{{tab:budget-variational}}
  \\begin{{tabular}}{{l r r r r r r r r r r}}
    \\toprule
    Des. & idx & $h_p$ [km] & $N_F$ & pred.\\ $E_{{\\mathrm{{At}}}}$ &
      pred.\\ $E_{{\\mathrm{{fix}}}}$ & pred.\\ $\\rho$ &
      meas.\\ $E_{{\\mathrm{{At}}}}$ &
      meas.\\ $E_{{\\mathrm{{fix}}}}$ & calib. & meas.\\ $\\rho$\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orbits", type=int, default=8)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    key = f"beta_{BETA:.2f}"
    tasks = []
    for design in ("A", "B"):
        rows = {int(r["sobol_index"]): r
                for r in json.loads(ROWS[design].read_text())["rows"]}
        pr = [r for r in pareto["designs"][design]["rows"]
              if not r["budgets"][key]["censored"]]
        pr.sort(key=lambda r: r["hp_km"])
        picks = [pr[int(i)] for i in
                 np.linspace(0, len(pr) - 1, a.orbits // 2).round()]
        for p in picks:
            idx = int(p["sobol_index"])
            tasks.append({"design": design, "row": rows[idx],
                          "spec": p["budgets"][key]})
    print(f"[variational] {len(tasks)} orbits at beta={BETA}", flush=True)
    done, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['sobol_index']:03d} {rec.get('message')}", flush=True)
                continue
            rec["measured"] = measured(rec["design"], rec["sobol_index"])
            if rec["measured"]:
                mf = rec["measured"]["fixed_budget"]
                pf = rec["policies"]["fixed_budget"]["predicted_pos_rms_m"]
                rec["calibration_ratio_fixed"] = (pf / mf) if mf > 0 else None
            done.append(rec)
            print(f"  [{n}/{len(tasks)}] {rec['design']}{rec['sobol_index']:03d} "
                  f"pred rho={rec['predicted_ratio_fixed_over_atallah']:.3g} "
                  f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
    done.sort(key=lambda r: (r["design"], r["sobol_index"]))
    payload = {"schema": "r14_variational_budget_v1", "created_utc": base.utc_now(),
               "beta": BETA, "gradient_degree": GRADIENT_DEGREE,
               "rows": done, "failures": fails,
               "summary": {
                   "orbits": len(done),
                   "calibration_ratio_fixed": stat(
                       [r.get("calibration_ratio_fixed") for r in done]),
                   "predicted_ratio": stat(
                       [r["predicted_ratio_fixed_over_atallah"] for r in done]),
                   "predicted_favors_atallah": int(sum(
                       (r["predicted_ratio_fixed_over_atallah"] or 0) > 1.0
                       for r in done)),
                   "comparator_is_critical": int(sum(
                       r["comparator_is_critical"] for r in done))},
               "source": base.provenance()}
    base.atomic_json(OUTPUT, payload)
    if done:
        TABLE.write_text(build_table(done), encoding="utf-8")
    s = payload["summary"]
    print(f"[variational] orbits={s['orbits']} "
          f"calibration median={s['calibration_ratio_fixed']['median']:.3f} "
          f"predicted rho median={s['predicted_ratio']['median']:.3g} "
          f"favoring radial on {s['predicted_favors_atallah']}/{s['orbits']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
