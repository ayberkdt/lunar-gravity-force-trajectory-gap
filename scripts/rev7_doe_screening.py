"""R7 / Test-9 Stage 1: 24-orbit DOE screening of the scheduling penalty.

Robustness plan (revision/ROBUSTNESS_EXPERIMENT_PLAN.md, Test 9, Asama 1):
a stratified design of experiments over lunar orbit geometry instead of
phase-only variation of a single geometry. Six canonical orbits plus an
18-point maximin Latin-hypercube design over

    (perilune 30-120 km, apolune 180-500 km, inclination 0-100 deg,
     argument of periapsis 0-360 deg, RAAN 0-360 deg),

each propagated 24 h with the r1-r6 harness (uniformly rotating body-fixed
JGGRX field, DOP853, rtol 1e-12 / atol 1e-5). Under uniform rotation the
epoch is degenerate with the body-fixed perilune longitude, so epoch is not
a separate design dimension in this stage (deferred to the MOON_PA stage-2
matrix).

Pre-registered per-orbit comparators (no shared fixed-N comparator):
  fixed_crit  N = empirical tail-criterion N_min at the perilune altitude
  fixed_work  nearest integer to sqrt(mean N^2) of the empirical-lookup
              schedule run of the same orbit (rev4 cost-match convention)
  sched_down / sched_up   Kaula-criterion altitude schedules, q=10, cap 250
  sched_emp   empirical-lookup altitude schedule, q=10, cap 250

Truth: N=300 for perilune >= 50 km, N=600 below (plan section 9.5); each
policy run uses the same model object as its truth run. Results stream to
metrics/r7_doe_screening_stage1.json after every orbit.
"""

from __future__ import annotations

import argparse
import math
import sys
import time

import numpy as np

from rev3_common import (DAY, SEED, dump, err_stats, kernel_args, load_model,
                         propagate, warmup, degree_power)

# importing rev3_common above put the pinned Lunaris repo's src on sys.path
from lunaris.common.math_utils import (  # noqa: E402
    LUNAR_DEGREE_POWER_EXPONENT,
    coe_to_rv,
    recommended_sh_degree,
)

RTOL, ATOL = 1e-12, 1e-5
DUR = 1.0 * DAY
OUT_STEP = 60.0
EPS = 1.0e-3
FLOOR, CAP, Q = 60, 250, 10
H_GRID_KM = np.arange(20.0, 561.0, 10.0)

# ------------------------------------------------------------------ design
CANONICAL = [
    # name, hp_km, ha_km, incl_deg, argp_deg, raan_deg
    ("c1_circ100_polar", 100.0, 100.0, 90.0, 0.0, 0.0),
    ("c2_50x300_polar", 50.0, 300.0, 90.0, 0.0, 0.0),
    ("c3_50x300_i60", 50.0, 300.0, 60.0, 0.0, 0.0),
    ("c4_50x300_i5", 50.0, 300.0, 5.0, 0.0, 0.0),
    ("c5_100x300_polar", 100.0, 300.0, 90.0, 0.0, 0.0),
    ("c6_lro_30x216", 30.0, 216.0, 90.0, 270.0, 0.0),
]
LHS_BOUNDS = {  # plan section 9.1
    "hp_km": (30.0, 120.0),
    "ha_km": (180.0, 500.0),
    "incl_deg": (0.0, 100.0),
    "argp_deg": (0.0, 360.0),
    "raan_deg": (0.0, 360.0),
}
N_LHS = 18
N_LHS_CANDIDATES = 300


def lhs_maximin(n_pts: int, n_dim: int, seed: int) -> np.ndarray:
    """Best-of-N_LHS_CANDIDATES maximin Latin hypercube in the unit cube."""
    from scipy.stats import qmc
    rng = np.random.default_rng(seed)
    best, best_d = None, -1.0
    for _ in range(N_LHS_CANDIDATES):
        sampler = qmc.LatinHypercube(d=n_dim, seed=rng)
        pts = sampler.random(n_pts)
        d = np.min(
            np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
            + np.eye(n_pts) * 1e9
        )
        if d > best_d:
            best, best_d = pts, d
    return best


def build_design() -> list[dict]:
    orbits = [
        {"name": n, "family": "canonical", "hp_km": hp, "ha_km": ha,
         "incl_deg": i, "argp_deg": w, "raan_deg": om}
        for n, hp, ha, i, w, om in CANONICAL
    ]
    keys = list(LHS_BOUNDS)
    pts = lhs_maximin(N_LHS, len(keys), SEED)
    for j, row in enumerate(pts):
        o = {"name": f"lhs_{j:02d}", "family": "lhs"}
        for k, u in zip(keys, row):
            lo, hi = LHS_BOUNDS[k]
            o[k] = float(lo + u * (hi - lo))
        orbits.append(o)
    return orbits


def initial_state(model, orb: dict) -> np.ndarray:
    rp = model.r_ref + orb["hp_km"] * 1e3
    ra = model.r_ref + orb["ha_km"] * 1e3
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    r, v = coe_to_rv(a, e, math.radians(orb["incl_deg"]),
                     math.radians(orb["raan_deg"]),
                     math.radians(orb["argp_deg"]), 0.0, model.mu)
    return np.concatenate([np.asarray(r, float), np.asarray(v, float)])


def perilune_geography(y0: np.ndarray) -> dict:
    """Body-fixed perilune longitude/latitude at t=0 (frames aligned)."""
    x, y, z = y0[:3]
    r = math.sqrt(x * x + y * y + z * z)
    return {"perilune_lon_deg_bodyfixed_t0": math.degrees(math.atan2(y, x)),
            "perilune_lat_deg_bodyfixed_t0": math.degrees(math.asin(z / r))}


# --------------------------------------------------------------- schedules
def emp_nmin_exact(power: np.ndarray, r_ref: float, h_m: float,
                   eps: float = EPS) -> int:
    """Empirical tail-criterion N_min at an exact altitude (rev4-style,
    integer, unquantized)."""
    n_arr = np.arange(len(power), dtype=np.float64)
    r = r_ref + h_m
    ratio_n = np.exp(n_arr * math.log(r_ref / r))
    sig = np.sqrt((n_arr + 1.0) * (2.0 * n_arr + 1.0)) * ratio_n * np.sqrt(power)
    sq = sig ** 2
    total = float(np.sum(sq[2:]))
    budget = eps * eps * total
    tail = total
    nmin = len(sq) - 1
    for n in range(2, len(sq)):
        if tail <= budget:
            nmin = n - 1
            break
        tail -= sq[n]
    return max(FLOOR, nmin)


def kaula_table(model, policy: str) -> dict:
    table = {}
    for hk in H_GRID_KM:
        n = recommended_sh_degree(hk, model.r_ref,
                                  kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                  kaula_tail_fraction=EPS)
        n = max(FLOOR, min(CAP, n))
        if policy == "down":
            nq = (n // Q) * Q
        else:
            nq = ((n + Q - 1) // Q) * Q
        table[float(hk)] = max(FLOOR, min(CAP, nq))
    return table


def emp_table(model, power: np.ndarray) -> dict:
    table = {}
    for hk in H_GRID_KM:
        nmin = emp_nmin_exact(power, model.r_ref, hk * 1e3)
        table[float(hk)] = max(FLOOR, min(CAP, (nmin // Q) * Q))
    return table


def alt_sched(table: dict):
    hmax, hmin = max(table), min(table)
    def f(t, h_m):
        hb = min(hmax, max(hmin, 10.0 * math.floor(h_m / 1e3 / 10.0)))
        return table[hb]
    return f


# -------------------------------------------------------------------- main
def run_orbit(orb: dict, models: dict, tables: dict, powers: dict,
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
    truth, truth_rhs, truth_wall = propagate(
        model, y0, dur, t_grid, lambda t, h: truth_deg, args, RTOL, ATOL)
    print(f"  truth N={truth_deg}: wall {truth_wall:6.1f}s "
          f"rhs {truth_rhs.n_calls}", flush=True)

    alt = np.linalg.norm(truth.y[:3], axis=0) - model.r_ref
    row["truth_altitude_range_km"] = [float(alt.min() / 1e3),
                                      float(alt.max() / 1e3)]

    def run_policy(name, degfun):
        sol, rhs, wall = propagate(model, y0, dur, t_grid, degfun, args,
                                   RTOL, ATOL)
        st = err_stats(sol.y, truth.y)
        st.update({
            "n_rhs": rhs.n_calls, "wall_s": wall,
            "grav_s": rhs.grav_ns / 1e9,
            "mean_deg_sq": rhs.sum_deg_sq / max(rhs.n_calls, 1),
        })
        row["policies"][name] = st
        print(f"  {name:11s}: rms {st['pos_rms_m']:10.3f} m  "
              f"final-I {st['ric_final_m']['in_track']:10.3f} m  "
              f"grav {st['grav_s']:6.1f}s", flush=True)
        return st

    row["policies"] = {}
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
    row["orbit_wall_s"] = time.perf_counter() - t0
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="2 orbits x 2 h to validate the pipeline")
    ap.add_argument("--stage", type=int, default=1, choices=(1, 2),
                    help="1: 24 h screening (60 s grid); "
                         "2: 7-day main matrix (120 s grid, plan 9.4)")
    a = ap.parse_args()

    global OUT_STEP
    if a.smoke:
        dur = 2.0 * 3600.0
        out_name = "r7_doe_screening_stage1_smoke.json"
    elif a.stage == 2:
        dur = 7.0 * DAY
        OUT_STEP = 120.0
        out_name = "r7_doe_matrix_stage2.json"
    else:
        dur = DUR
        out_name = "r7_doe_screening_stage1.json"

    orbits = build_design()
    if a.smoke:
        orbits = [o for o in orbits
                  if o["name"] in ("c2_50x300_polar", "c6_lro_30x216")]

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
        "plan": f"ROBUSTNESS_EXPERIMENT_PLAN.md Test 9 / Stage {a.stage}"
                + (" screening" if a.stage == 1 else " main 7-day matrix"),
        "duration_s": dur, "output_step_s": OUT_STEP,
        "integrator": "DOP853", "rtol": RTOL, "atol": ATOL,
        "rotation": "uniform sidereal about polar axis; epoch degenerate "
                    "with body-fixed perilune longitude in this stage",
        "truth_rule": "N=600 if perilune < 50 km else N=300",
        "schedule_params": {"eps_tail_fraction": EPS, "floor": FLOOR,
                            "cap": CAP, "quantum": Q,
                            "alt_grid_km": [float(H_GRID_KM[0]),
                                            float(H_GRID_KM[-1]), 10.0]},
        "comparator_rule": {
            "fixed_crit": "empirical tail-criterion N_min at perilune "
                          "altitude (rev4 convention, unquantized)",
            "fixed_work": "nearest integer to sqrt(mean N^2) of the "
                          "sched_emp run of the same orbit",
        },
        "lhs": {"n_points": N_LHS, "n_candidates": N_LHS_CANDIDATES,
                "seed": SEED, "bounds": LHS_BOUNDS},
        "smoke": a.smoke,
    }

    rows = []
    for k, orb in enumerate(orbits):
        print(f"[{k + 1}/{len(orbits)}] {orb['name']} "
              f"hp={orb['hp_km']:.1f} ha={orb['ha_km']:.1f} "
              f"i={orb['incl_deg']:.1f} w={orb['argp_deg']:.1f} "
              f"W={orb['raan_deg']:.1f}", flush=True)
        rows.append(run_orbit(orb, models, tables, powers, dur))
        dump(out_name, {"scenario": scenario, "design": orbits,
                        "rows": rows, "complete": k + 1 == len(orbits)})

    # ------------------------------------------------------------ summary
    def frac_better(sched, comp):
        vals = [r["ratios"][sched][comp] for r in rows]
        return {"n_orbits": len(vals),
                "frac_schedule_better": float(np.mean([v > 1.0 for v in vals])),
                "median_rho": float(np.median(vals)),
                "p10_rho": float(np.percentile(vals, 10)),
                "p90_rho": float(np.percentile(vals, 90))}

    summary = {s: {"vs_fixed_crit": frac_better(s, "rho_crit"),
                   "vs_fixed_work": frac_better(s, "rho_work")}
               for s in ("sched_emp", "sched_down", "sched_up")}
    print("summary:", summary, flush=True)
    dump(out_name, {"scenario": scenario, "design": orbits, "rows": rows,
                    "summary": summary, "complete": True})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
