"""Numerical-floor check for the designed-matrix and month-long results.

Reviewer concern: the 24-orbit population and the 28-day extension include
policies whose errors are only tens of meters, so the integration floor must be
shown to sit well below them. Rather than re-run all 24 orbits at three
tolerances, this selects orbits spanning the policy-error distribution and
re-integrates each at a tighter tolerance; the run-to-run difference is the
integration envelope E_num, and the report gives E_policy / E_num.

Two blocks:
  * 24-orbit matrix: pick min / p25 / median / p75 / max orbits by
    best-schedule RMS (Stage-2), re-run truth, fixed_crit, fixed_work and the
    best schedule at baseline (rtol=1e-12) and tight (rtol=1e-13).
  * 28-day: polar and LRO-like, re-run at the Stage-3 tolerance (rtol=1e-11)
    and a tighter one (rtol=1e-12); best schedule and fixed_crit.

Run: ``.venv\\Scripts\\python.exe robustness_numerical_floor_check.py``
     add ``--smoke`` for a 1-day / 2-orbit shakedown.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from rev3_common import DAY, dump, err_stats, kernel_args, load_model, \
    propagate, warmup, degree_power
from rev7_doe_screening import (CAP, alt_sched, emp_nmin_exact, emp_table,
                                initial_state, kaula_table)

METRICS = Path(__file__).resolve().parents[1] / "metrics"
STAGE2 = METRICS / "r7_doe_matrix_stage2.json"


def _rms_diff(a, b):
    d = np.linalg.norm(a[:3] - b[:3], axis=0)
    return float(np.sqrt(np.mean(d * d)))


def select_orbits(rows, k=5):
    scheds = ("sched_emp", "sched_down", "sched_up")
    best = []
    for r in rows:
        e = min(r["policies"][s]["pos_rms_m"] for s in scheds)
        best.append((e, r))
    best.sort(key=lambda x: x[0])
    n = len(best)
    idx = sorted(set(int(round(q * (n - 1))) for q in (0.0, 0.25, 0.5, 0.75, 1.0)))
    return [best[i][1] for i in idx][:k]


def floor_for_orbit(orb, dur, models, tables, powers, tol_pairs):
    truth_deg = 600 if orb["hp_km"] < 50.0 else 300
    model, args = models[truth_deg]
    y0 = initial_state(model, orb)
    grid = np.arange(0.0, dur + 1.0, 300.0)
    n_crit = min(CAP, emp_nmin_exact(powers[truth_deg], model.r_ref,
                                     orb["hp_km"] * 1e3))
    scheds = {"sched_emp": alt_sched(tables[("emp", truth_deg)]),
              "sched_down": alt_sched(tables[("down", truth_deg)]),
              "sched_up": alt_sched(tables[("up", truth_deg)])}

    # work-matched degree from the empirical schedule mean N^2
    def run(degfun, rtol, atol):
        sol, rhs, _ = propagate(model, y0, dur, grid, degfun, args, rtol, atol,
                                max_step=120.0)
        return sol.y, rhs

    (r_base, a_base), (r_tight, a_tight) = tol_pairs
    truth_base, _ = run(lambda t, h: truth_deg, r_base, a_base)
    truth_tight, _ = run(lambda t, h: truth_deg, r_tight, a_tight)
    _, emp_rhs = run(scheds["sched_emp"], r_base, a_base)
    n_work = int(round(math.sqrt(emp_rhs.sum_deg_sq / max(emp_rhs.n_calls, 1))))

    policies = {"fixed_crit": (lambda t, h, n=n_crit: n),
                "fixed_work": (lambda t, h, n=n_work: n)}
    # pick the best schedule by baseline error
    sched_err = {}
    sched_base = {}
    for name, fn in scheds.items():
        yb, _ = run(fn, r_base, a_base)
        sched_base[name] = yb
        sched_err[name] = _rms_diff(yb, truth_base)
    best_sched = min(sched_err, key=sched_err.get)
    policies[best_sched] = scheds[best_sched]

    out = {"name": orb["name"], "truth_degree": truth_deg,
           "n_critical": int(n_crit), "n_work": n_work,
           "best_schedule": best_sched,
           "truth_num_envelope_rms_m": _rms_diff(truth_base, truth_tight),
           "policies": {}}
    for name, fn in policies.items():
        yb, _ = run(fn, r_base, a_base)
        yt, _ = run(fn, r_tight, a_tight)
        e_policy = _rms_diff(yb, truth_base)
        e_num = _rms_diff(yb, yt)
        out["policies"][name] = {
            "E_policy_rms_m": e_policy,
            "E_numerical_envelope_rms_m": e_num,
            "ratio_policy_over_envelope": e_policy / max(e_num, 1e-12)}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()

    rows = json.loads(STAGE2.read_text())["rows"]
    picks = select_orbits(rows, k=2 if a.smoke else 5)

    need = sorted({600 if o["hp_km"] < 50.0 else 300 for o in picks}
                  | {300, 600})
    models, powers, tables = {}, {}, {}
    for d in need:
        m = load_model(d); ar = kernel_args(m); warmup(m, ar)
        models[d] = (m, ar); powers[d] = degree_power(m)
        tables[("down", d)] = kaula_table(m, "down")
        tables[("up", d)] = kaula_table(m, "up")
        tables[("emp", d)] = emp_table(m, powers[d])
        print(f"[model] degree {d} loaded", flush=True)

    dur7 = 1.0 * DAY if a.smoke else 7.0 * DAY
    dur28 = 1.0 * DAY if a.smoke else 28.0 * DAY

    payload = {"schema": "robustness_numerical_floor_check_v1",
               "smoke": a.smoke, "matrix_7day": [], "longarc_28day": []}

    # 24-orbit matrix subset at 7 days: baseline 1e-12 vs tight 1e-13
    print("== 24-orbit matrix numerical floor (7 d) ==", flush=True)
    for orb in picks:
        r = floor_for_orbit(orb, dur7, models, tables, powers,
                            ((1e-12, 1e-5), (1e-13, 1e-6)))
        payload["matrix_7day"].append(r)
        worst = min(p["ratio_policy_over_envelope"] for p in r["policies"].values())
        print(f"  {orb['name']:16s} min ratio E_policy/E_num = {worst:.0f}",
              flush=True)

    # 28-day polar + LRO: baseline 1e-11 vs tight 1e-12
    print("== 28-day numerical floor ==", flush=True)
    canon = {o["name"]: o for o in rows}
    for nm in ("c2_50x300_polar", "c6_lro_30x216"):
        orb = canon.get(nm) or {"name": nm}
        # reconstruct canonical geometry if not in rows
        if "hp_km" not in orb:
            from rev7_doe_screening import CANONICAL
            for n2, hp, ha, i, w, om in CANONICAL:
                if n2 == nm:
                    orb = {"name": nm, "hp_km": hp, "ha_km": ha,
                           "incl_deg": i, "argp_deg": w, "raan_deg": om}
        r = floor_for_orbit(orb, dur28, models, tables, powers,
                            ((1e-11, 1e-6), (1e-12, 1e-7)))
        payload["longarc_28day"].append(r)
        worst = min(p["ratio_policy_over_envelope"] for p in r["policies"].values())
        print(f"  {nm:16s} min ratio E_policy/E_num = {worst:.0f}", flush=True)

    dump("robustness_numerical_floor_check.json", payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
