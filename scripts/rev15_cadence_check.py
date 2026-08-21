"""Output-cadence convergence of the force-defect measurement (R15-D).

Every force-level number in O25/O26 is a statistic over the archived truth
trajectory sampled on a uniform 120-s grid. On an eccentric arc the perilune
passage -- which dominates the defect RMS -- lasts only minutes, so 120 s may or
may not resolve it. This tests both directions:

  * coarser, by decimating the archived grid (240 s, 480 s), which is free;
  * finer, by re-propagating the reference with 30-s and 10-s output grids. The
    solver takes the same steps either way; only the number of dense-output
    evaluations changes, so the extra cost is small compared with the
    integration itself.

For each orbit and cadence we recompute the defect ratio between the
budget-calibrated radial history and its budget-saturating fixed comparator, and
report whether the median ratio and the per-orbit win ordering move.

Usage:
    python rev15_cadence_check.py --orbits 8 --workers 4
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
from rev14_budget_pareto import (DESIGNS, LEVEL, BIN_KM, FLOOR, accel_inertial,
                                 degrees_from_table, calibrate_tolerance,
                                 fixed_degree_for, _model, _g)

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r15_cadence_check.json"
TABLE = METRICS / "r15_cadence_check_table.tex"

BETAS = [0.50, 1.00, 2.00]
DURATION = 7.0 * base.DAY
MAX_STEP = 60.0
FINE_STEPS = [30.0, 10.0]        # re-propagated
COARSE_FACTORS = [2, 4]          # decimated from the archived 120-s grid
RTOL = 1.0e-13
ATOL = np.array([1.0e-6] * 3 + [1.0e-9] * 3)   # the archived "tighter" level


def defect_rms(t, y, model, args, adopted, degrees_fn, n_fixed):
    """RMS truncation defect of a degree history and of a fixed degree."""
    r_all = y[0:3, :].T
    h_km = (np.linalg.norm(r_all, axis=1) - model.r_ref) / 1e3
    deg_a = degrees_fn(h_km)
    sa = sf = 0.0
    n = len(t)
    for j in range(n):
        rj, tj = r_all[j], float(t[j])
        a_truth = accel_inertial(rj, tj, adopted, args)
        da = accel_inertial(rj, tj, int(deg_a[j]), args) - a_truth
        df = accel_inertial(rj, tj, n_fixed, args) - a_truth
        sa += float(da @ da)
        sf += float(df @ df)
    return math.sqrt(sa / n), math.sqrt(sf / n), float(np.mean(deg_a.astype(float) ** 2))


def worker(task: dict) -> dict:
    design, row = task["design"], task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        g = _g(adopted)

        # archived 120-s reference states
        t120, y120 = base.load_raw(DESIGNS[design]["r11_raw"]
                                   / f"sobolA_{index:03d}" / f"truth_{LEVEL}.npz")
        h120 = (np.linalg.norm(y120[0:3, :].T, axis=1) - model.r_ref) / 1e3

        # calibrate once on the archived grid, exactly as O25 does, and freeze
        # the resulting policies so every cadence measures the same histories
        policies = {}
        for beta in BETAS:
            cal = calibrate_tolerance(model, g, hp_km, ha_km, adopted, h120,
                                      beta * n_crit ** 2)
            n_f, cens = fixed_degree_for(beta, n_crit, adopted)
            _, table = at.atallah_binned_schedule(model, g, cal["tol"], hp_km,
                                                  ha_km, floor=FLOOR, cap=adopted,
                                                  bin_km=BIN_KM)
            policies[beta] = {"tol": cal["tol"], "n_fixed": n_f,
                              "censored": bool(cens or not cal["attainable"]),
                              "table": {float(k): int(v) for k, v in table.items()}}

        grids = {}
        grids[120.0] = (t120, y120)
        for f in COARSE_FACTORS:
            grids[120.0 * f] = (t120[::f], y120[:, ::f])
        for step in FINE_STEPS:
            grid = np.arange(0.0, DURATION + 0.5 * step, step)
            tt, yy, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, lambda _t, _h, n=adopted: n, args,
                RTOL, ATOL, max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"design": design, "sobol_index": index,
                        "status": "truth_failure", "message": fail}
            grids[step] = (tt, yy)

        out = {}
        for beta, p in policies.items():
            if p["censored"]:
                continue
            per_cadence = {}
            for step in sorted(grids):
                t, y = grids[step]
                ea, ef, work = defect_rms(
                    t, y, model, args, adopted,
                    lambda h, tab=p["table"]: degrees_from_table(tab, h),
                    p["n_fixed"])
                per_cadence[f"{step:.0f}"] = {
                    "n_epochs": int(len(t)),
                    "atallah_defect_rms": ea, "fixed_defect_rms": ef,
                    "R_a": (ef / ea) if ea > 0 else None,
                    "atallah_smaller": bool(ea < ef),
                    "sampled_work": work,
                    "beta_achieved": work / (n_crit ** 2)}
            ref = per_cadence["120"]
            for k, v in per_cadence.items():
                v["R_a_rel_change_vs_120s"] = (
                    None if (v["R_a"] is None or ref["R_a"] is None)
                    else v["R_a"] / ref["R_a"] - 1.0)
                v["beta_rel_change_vs_120s"] = (
                    v["beta_achieved"] / ref["beta_achieved"] - 1.0)
            out[f"beta_{beta:.2f}"] = {"n_fixed": p["n_fixed"],
                                       "cadences": per_cadence}
        return {"design": design, "sobol_index": index, "hp_km": hp_km,
                "adopted_truth_degree": adopted, "status": "complete",
                "budgets": out}
    except Exception as exc:
        return {"design": design, "sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def summarize(rows):
    out = {}
    for key in (f"beta_{b:.2f}" for b in BETAS):
        entries = [r["budgets"][key] for r in rows if key in r["budgets"]]
        if not entries:
            continue
        cads = sorted({c for e in entries for c in e["cadences"]}, key=float)
        per = {}
        for c in cads:
            R = [e["cadences"][c]["R_a"] for e in entries if c in e["cadences"]]
            rel = [abs(e["cadences"][c]["R_a_rel_change_vs_120s"])
                   for e in entries if c in e["cadences"]
                   and e["cadences"][c]["R_a_rel_change_vs_120s"] is not None]
            wins = sum(e["cadences"][c]["automatic"] if False else
                       e["cadences"][c]["atallah_smaller"]
                       for e in entries if c in e["cadences"])
            per[c] = {"orbits": len(R), "median_R_a": float(np.median(R)),
                      "atallah_smaller": int(wins),
                      "max_abs_rel_change_vs_120s": (max(rel) if rel else 0.0),
                      "median_abs_rel_change_vs_120s": (float(np.median(rel))
                                                        if rel else 0.0)}
        # ordering stability: does any orbit flip its winner against 120 s?
        flips = 0
        for e in entries:
            ref = e["cadences"]["120"]["atallah_smaller"]
            for c, v in e["cadences"].items():
                if v["atallah_smaller"] != ref:
                    flips += 1
                    break
        out[key] = {"orbits": len(entries), "per_cadence": per,
                    "orbits_with_any_winner_flip": flips}
    return out


def build_table(payload) -> str:
    lines = []
    for key, s in payload["summary"].items():
        beta = key.replace("beta_", "")
        for c in sorted(s["per_cadence"], key=float):
            v = s["per_cadence"][c]
            lines.append(
                f"    {beta} & {c} & {v['orbits']} & {v['median_R_a']:.3g} & "
                f"{v['atallah_smaller']}/{v['orbits']} & "
                f"{100 * v['median_abs_rel_change_vs_120s']:.1f}\\% & "
                f"{100 * v['max_abs_rel_change_vs_120s']:.1f}\\%\\\\")
    body = "\n".join(lines)
    flips = sum(s["orbits_with_any_winner_flip"] for s in payload["summary"].values())
    return f"""% auto-generated by rev15_cadence_check.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Output-cadence convergence of the force-defect measurement. Every
  force-level statistic in this paper is taken on the archived 120-s reference grid;
  on an eccentric arc the perilune passage that dominates the defect lasts only
  minutes, so the grid is tested in both directions: decimated to 240 and 480~s,
  and re-propagated at 30 and 10~s. The last two columns give the median and
  worst per-orbit change of the defect ratio against the 120-s value. Across all
  cadences and budgets, {flips} of the tested orbits change which policy has the
  smaller defect.}}
  \\label{{tab:cadence-check}}
  \\begin{{tabular}}{{l r r r c r r}}
    \\toprule
    $\\beta$ & cadence [s] & orbits & median $R_a$ & At.\\ smaller &
      median $|\\Delta R_a|$ & max $|\\Delta R_a|$\\\\
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
    tasks = []
    for design in ("A", "B"):
        rows = json.loads(DESIGNS[design]["rows"].read_text())["rows"]
        rows.sort(key=lambda r: r["design_point"]["hp_km"])
        picks = [rows[int(i)] for i in
                 np.linspace(0, len(rows) - 1, a.orbits // 2).round()]
        tasks += [{"design": design, "row": r} for r in picks]
    print(f"[cadence] {len(tasks)} orbits x {len(BETAS)} budgets x 5 cadences",
          flush=True)
    t0 = time.time()
    done, fails = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['design']}{rec['sobol_index']:03d} "
                      f"{rec.get('message')}", flush=True)
                continue
            done.append(rec)
            print(f"  [{n}/{len(tasks)}] {rec['design']}{rec['sobol_index']:03d} "
                  f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
    done.sort(key=lambda r: (r["design"], r["sobol_index"]))
    payload = {"schema": "r15_cadence_check_v1", "created_utc": base.utc_now(),
               "betas": BETAS, "coarse_factors": COARSE_FACTORS,
               "fine_output_steps_s": FINE_STEPS,
               "rows": done, "failures": fails, "summary": summarize(done),
               "source": base.provenance()}
    base.atomic_json(OUTPUT, payload)
    if done:
        TABLE.write_text(build_table(payload), encoding="utf-8")
    for k, s in payload["summary"].items():
        print(f"  {k}: " + ", ".join(
            f"{c}s R_a={v['median_R_a']:.3g} (dmed={100*v['median_abs_rel_change_vs_120s']:.1f}%)"
            for c, v in sorted(s["per_cadence"].items(), key=lambda kv: float(kv[0])))
            + f" | winner flips: {s['orbits_with_any_winner_flip']}/{s['orbits']}",
            flush=True)
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
