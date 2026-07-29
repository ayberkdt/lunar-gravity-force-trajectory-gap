"""Test 9 / Stage 3: 28-day long-arc validation on four representative orbits.

Stage 2 (7-day, 24-orbit matrix) is complete; Stage 3 checks whether the
policy ranking is preserved over a full lunar month and whether the scheduling
in-track error keeps growing secularly or saturates.  Per plan section 9.4 the
four representative geometries are chosen from the canonical set by regime:

    c2_50x300_polar : polar case where scheduling loses badly
    c6_lro_30x216   : LRO-like case where scheduling beats work-matched fixed
    c3_50x300_i60   : intermediate-inclination regime
    c4_50x300_i5    : near-equatorial / non-polar regime

Each orbit is propagated for 28 days against its truth field (N=600 for the
30 km perilune LRO case, N=300 otherwise; the N=600->900 adequacy of the LRO
truth is established separately in Test 3).  In-track and total-RMS errors are
recorded at weekly checkpoints (7, 14, 21, 28 d) so secular vs. saturating
growth can be read off directly, and the rho_crit / rho_work ratios are
compared with the Stage 2 seven-day values.

Full run:  ``.venv\\Scripts\\python.exe rev7_doe_stage3_longarc.py``
Smoke run: add ``--smoke`` (two representative orbits, 6 hours).
"""

from __future__ import annotations

import argparse
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from rev3_common import (DAY, Rhs, dump, err_stats, kernel_args, load_model,
                         propagate, warmup, degree_power)
from rev7_doe_screening import (CANONICAL, CAP, alt_sched,
                                emp_nmin_exact, emp_table, initial_state,
                                kaula_table, perilune_geography)

# Long-arc tolerance family (matches the Test 3 seven-day settings, not the
# tighter rev7 stage-1/2 rtol=1e-12): a 28-day arc needs a step-size cap and a
# slightly looser rtol to avoid adaptive step collapse.  The same tolerance is
# used for truth and every policy within this stage, so model-relative
# differences remain internally consistent.
L_RTOL, L_ATOL, L_MAX_STEP = 1.0e-11, 1.0e-6, 120.0

OUT_STEP = 300.0
STAGE3_NAMES = ("c2_50x300_polar", "c3_50x300_i60", "c4_50x300_i5",
                "c6_lro_30x216")
CHECKPOINT_DAYS = (7, 14, 21, 28)


def _orbit_dict(name: str) -> dict:
    for n, hp, ha, i, w, om in CANONICAL:
        if n == name:
            return {"name": n, "family": "canonical", "hp_km": hp, "ha_km": ha,
                    "incl_deg": i, "argp_deg": w, "raan_deg": om}
    raise KeyError(name)


def _checkpoint_stats(sol_y, truth_y, t_grid, duration):
    """In-track and RMS error at each weekly checkpoint (cumulative RMS up to
    the checkpoint; instantaneous signed in-track at the checkpoint)."""
    out = {}
    for d in CHECKPOINT_DAYS:
        t_cp = d * DAY
        if t_cp > duration + 1.0:
            continue
        idx = int(np.searchsorted(t_grid, t_cp - 1e-6))
        idx = min(idx, t_grid.size - 1)
        st = err_stats(sol_y[:, :idx + 1], truth_y[:, :idx + 1])
        out[f"d{d}"] = {
            "pos_rms_m": st["pos_rms_m"],
            "pos_at_checkpoint_m": st["pos_final_m"],
            "in_track_rms_m": st["ric_rms_m"]["in_track"],
            "in_track_at_checkpoint_m": st["ric_final_m"]["in_track"],
        }
    return out


def _impact_event(model):
    """Terminal event when the trajectory reaches the reference sphere."""
    def event(t, y):
        return float(np.linalg.norm(y[:3]) - model.r_ref)
    event.terminal = True
    event.direction = -1
    return event


def _truth_with_impact(model, y0, dur, t_grid, truth_deg, args):
    """Propagate the truth degree with a surface-impact terminal event.
    Returns (sol_y_on_grid_or_None, impacted, impact_t)."""
    rhs = Rhs(model, lambda t, h: truth_deg, args)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853", t_eval=t_grid,
                    rtol=L_RTOL, atol=L_ATOL, max_step=L_MAX_STEP,
                    events=_impact_event(model))
    wall = time.perf_counter() - t0
    if sol.t_events[0].size > 0:
        return None, True, float(sol.t_events[0][0]), rhs, wall
    if not sol.success:
        return None, True, float(sol.t[-1]), rhs, wall
    return sol.y, False, None, rhs, wall


def run_orbit_long(orb: dict, models: dict, tables: dict, powers: dict,
                   dur: float) -> dict:
    truth_deg = 600 if orb["hp_km"] < 50.0 else 300
    model, args = models[truth_deg]
    y0 = initial_state(model, orb)
    t_grid = np.arange(0.0, dur + 1.0, OUT_STEP)

    n_crit = min(CAP, emp_nmin_exact(powers[truth_deg], model.r_ref,
                                     orb["hp_km"] * 1e3))
    row = dict(orb)
    row.update(perilune_geography(y0))
    rp = model.r_ref + orb["hp_km"] * 1e3
    ra = model.r_ref + orb["ha_km"] * 1e3
    row["eccentricity"] = (ra - rp) / (ra + rp)
    row["dwell_ratio_ra_over_rp"] = ra / rp
    row["truth_degree"] = truth_deg
    row["n_critical_empirical"] = int(n_crit)

    t0 = time.perf_counter()
    truth_y, impacted, impact_t, truth_rhs, truth_wall = _truth_with_impact(
        model, y0, dur, t_grid, truth_deg, args)
    if impacted:
        row["impacted"] = True
        row["impact_day"] = impact_t / DAY
        row["note"] = ("reference orbit reaches the surface before the arc "
                       "ends; not a viable month-long geometry, policy "
                       "comparison omitted")
        row["orbit_wall_s"] = time.perf_counter() - t0
        print(f"  truth N={truth_deg}: IMPACT at day "
              f"{impact_t / DAY:.2f} -- orbit unstable, skipping policies",
              flush=True)
        return row

    row["impacted"] = False
    print(f"  truth N={truth_deg}: wall {truth_wall:7.1f}s "
          f"rhs {truth_rhs.n_calls}", flush=True)

    class _Truth:
        y = truth_y
    truth = _Truth()
    alt = np.linalg.norm(truth.y[:3], axis=0) - model.r_ref
    row["truth_altitude_range_km"] = [float(alt.min() / 1e3),
                                      float(alt.max() / 1e3)]

    row["policies"] = {}

    def run_policy(name, degfun):
        sol, rhs, wall = propagate(model, y0, dur, t_grid, degfun, args,
                                   L_RTOL, L_ATOL, max_step=L_MAX_STEP)
        st = err_stats(sol.y, truth.y)
        st.update({
            "n_rhs": rhs.n_calls, "wall_s": wall,
            "grav_s": rhs.grav_ns / 1e9,
            "mean_deg_sq": rhs.sum_deg_sq / max(rhs.n_calls, 1),
            "checkpoints": _checkpoint_stats(sol.y, truth.y, t_grid, dur),
        })
        row["policies"][name] = st
        print(f"  {name:11s}: rms {st['pos_rms_m']:11.3f} m  "
              f"final-I {st['ric_final_m']['in_track']:12.3f} m  "
              f"grav {st['grav_s']:7.1f}s", flush=True)
        return st

    st_emp = run_policy("sched_emp", alt_sched(tables[("emp", truth_deg)]))
    n_work = int(round(math.sqrt(st_emp["mean_deg_sq"])))
    row["n_work_matched"] = n_work
    run_policy("fixed_work", lambda t, h: n_work)
    run_policy("fixed_crit", lambda t, h: n_crit)
    run_policy("sched_down", alt_sched(tables[("down", truth_deg)]))
    run_policy("sched_up", alt_sched(tables[("up", truth_deg)]))

    e_crit = row["policies"]["fixed_crit"]["pos_rms_m"]
    e_work = row["policies"]["fixed_work"]["pos_rms_m"]
    g_crit = row["policies"]["fixed_crit"]["grav_s"]
    row["ratios"] = {}
    for s in ("sched_emp", "sched_down", "sched_up"):
        e_s = row["policies"][s]["pos_rms_m"]
        row["ratios"][s] = {
            "rho_crit": e_crit / e_s,
            "rho_work": e_work / e_s,
            "grav_time_saving_vs_crit":
                1.0 - row["policies"][s]["grav_s"] / g_crit,
        }
    # Track the actual best of the three schedules by full-arc position RMS.
    # Earlier archived output labelled sched_emp as the "best schedule" even
    # when sched_up supplied the best-schedule table entry.  Retain the
    # empirical diagnostic, but keep it separate from the selected policy.
    schedule_names = ("sched_emp", "sched_down", "sched_up")
    best_name = min(schedule_names,
                    key=lambda name: row["policies"][name]["pos_rms_m"])
    row["best_schedule_name"] = best_name
    cps = row["policies"][best_name]["checkpoints"]
    it = [abs(cps[f"d{d}"]["in_track_at_checkpoint_m"])
          for d in CHECKPOINT_DAYS if f"d{d}" in cps]
    row["best_schedule_in_track_growth"] = {
        "checkpoint_days": [d for d in CHECKPOINT_DAYS
                            if f"d{d}" in cps],
        "abs_in_track_m": it,
        "d28_over_d7": (it[-1] / it[0]) if len(it) >= 2 and it[0] > 0 else None,
    }
    emp_cps = row["policies"]["sched_emp"]["checkpoints"]
    emp_it = [abs(emp_cps[f"d{d}"]["in_track_at_checkpoint_m"])
              for d in CHECKPOINT_DAYS if f"d{d}" in emp_cps]
    row["sched_emp_in_track_growth"] = {
        "checkpoint_days": [d for d in CHECKPOINT_DAYS
                            if f"d{d}" in emp_cps],
        "abs_in_track_m": emp_it,
        "d28_over_d7": ((emp_it[-1] / emp_it[0])
                         if len(emp_it) >= 2 and emp_it[0] > 0 else None),
    }
    row["orbit_wall_s"] = time.perf_counter() - t0
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="two orbits x 6 h to validate the pipeline")
    a = ap.parse_args()

    global OUT_STEP
    if a.smoke:
        dur = 6.0 * 3600.0
        OUT_STEP = 120.0
        names = ("c2_50x300_polar", "c6_lro_30x216")
        out_name = "r7_doe_stage3_longarc_smoke.json"
    else:
        dur = 28.0 * DAY
        names = STAGE3_NAMES
        out_name = "r7_doe_stage3_longarc.json"

    orbits = [_orbit_dict(n) for n in names]
    need_degrees = sorted({600 if o["hp_km"] < 50.0 else 300 for o in orbits})
    models, powers, tables = {}, {}, {}
    for d in need_degrees:
        m = load_model(d)
        ar = kernel_args(m)
        warmup(m, ar)
        models[d] = (m, ar)
        powers[d] = degree_power(m)
        tables[("down", d)] = kaula_table(m, "down")
        tables[("up", d)] = kaula_table(m, "up")
        tables[("emp", d)] = emp_table(m, powers[d])
        print(f"[model] degree {d} loaded", flush=True)

    scenario = {
        "plan": "ROBUSTNESS_EXPERIMENT_PLAN.md Test 9 / Stage 3 long-arc",
        "duration_s": dur, "output_step_s": OUT_STEP,
        "checkpoint_days": list(CHECKPOINT_DAYS),
        "integrator": "DOP853", "rtol": L_RTOL, "atol": L_ATOL,
        "max_step_s": L_MAX_STEP,
        "tolerance_note": "long-arc family (Test 3 seven-day settings); "
                          "looser than rev7 stage-1/2 rtol=1e-12 to keep the "
                          "28-day adaptive step from collapsing",
        "surface_impact_handling": "terminal event at r = R_ref; impacting "
                                   "reference orbits are flagged and their "
                                   "policy comparison is omitted",
        "rotation": "uniform sidereal about polar axis; epoch degenerate "
                    "with body-fixed perilune longitude",
        "truth_rule": "N=600 if perilune < 50 km else N=300 "
                      "(LRO N=600->900 adequacy shown in Test 3)",
        "representative_orbits": {
            "c2_50x300_polar": "polar, scheduling loses badly",
            "c6_lro_30x216": "LRO-like, scheduling beats work-matched fixed",
            "c3_50x300_i60": "intermediate inclination",
            "c4_50x300_i5": "near-equatorial / non-polar",
        },
        "smoke": a.smoke,
    }

    rows = []
    for k, orb in enumerate(orbits):
        print(f"[{k + 1}/{len(orbits)}] {orb['name']} "
              f"hp={orb['hp_km']:.1f} ha={orb['ha_km']:.1f} "
              f"i={orb['incl_deg']:.1f}", flush=True)
        rows.append(run_orbit_long(orb, models, tables, powers, dur))
        dump(out_name, {"scenario": scenario, "rows": rows,
                        "complete": k + 1 == len(orbits)})

    print("done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
