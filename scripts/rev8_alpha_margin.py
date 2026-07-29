"""R8: alpha-margin robustness control on the existing 24-orbit DOE.

Referee question: is the long-arc scheduling loss merely a consequence of the
nominal spectrum-tail budget selecting degrees that are too low?  This study
applies a uniform multiplicative safety margin alpha to the archived
empirical-lookup altitude schedule (the Stage-2 ``sched_emp`` policy) and
measures whether -- and at what compute cost -- the margin recovers
fixed-degree fidelity.

Design (specified before aggregate inspection; see revision notes):
  * Orbits: the existing 24-orbit stratified DOE (seed 20260719), unchanged.
  * PRIMARY family (``ladder``): N_alpha(h) = clamp(up_q10(ceil(alpha*N0(h))),
    floor 60, cap truth_degree-10) with N0(h) the archived empirical-lookup
    table value (eps=1e-3, floor 60, cap 250, q=10, downward quantization).
    alpha in {1.00, 1.10, 1.20, 1.30, 1.50}.  At alpha=1.00 the table is
    identical to the archived sched_emp table, so the run must reproduce the
    archived per-orbit error exactly (hard sanity gate).
  * CONTROL family (``exact``): no-ladder variant
    N_alpha(h) = clamp(ceil(alpha*min(250, N_emp_min_raw(h))), 60,
    truth_degree-10) for alpha in {1.00, 1.20}, separating the margin effect
    from the 10-degree ladder.
  * Comparators: archived per-orbit fixed_work / fixed_crit errors are reused
    (deterministic; not re-ranked).  fixed_crit is RE-RUN once per orbit in
    this session purely so gravity-time ratios compare like-for-like timings;
    its error doubles as a reproduction check.
  * Integrator/tolerances IDENTICAL to the archived Stage-2 matrix:
    DOP853, rtol 1e-12, scalar atol 1e-5, no max-step cap, 120 s output grid,
    7-day arcs.  (The review sheet's "vector atol / 60 s max step" spec would
    contradict its own alpha=1.00 reproduction gate; the archived settings
    win.  Deviations are recorded in the output scenario.)

Metrics per (orbit, alpha): 7-day Cartesian RMS error vs the orbit's own
truth, in-track fraction, rho_work = E_work/E_sched and
rho_crit = E_crit/E_sched (paper Table convention: rho > 1 means the schedule
wins), mean N^2 work proxy, RHS call count, measured gravity-kernel time and
its ratio to the re-run fixed_crit, and the per-orbit smallest sufficient
margin alpha* = min{alpha : E_sched(alpha) <= E_work}.

Stages:
  full run  : python rev8_alpha_margin.py            (~2 h, streams JSON)
  smoke     : python rev8_alpha_margin.py --smoke    (2 orbits x 2 h)
  floor pass: python rev8_alpha_margin.py --guard    (re-integrates flagged
              (orbit, alpha) pairs at rtol 1e-13 / atol 1e-6 and records the
              run-to-run envelope E_num; Table-9 discipline)
  resolution: python rev8_alpha_margin.py --resolution-v2 (fills the remaining
              schedule envelopes and applies the conservative truth-inclusive
              pairwise resolution criterion)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from rev3_common import (DAY, SEED, dump, err_stats, kernel_args, load_model,
                         propagate, warmup, degree_power)
from rev7_doe_screening import (CAP, EPS, FLOOR, H_GRID_KM, Q, alt_sched,
                                build_design, emp_nmin_exact, emp_table,
                                initial_state, kaula_table)

RTOL, ATOL = 1e-12, 1e-5          # identical to the archived Stage-2 matrix
OUT_STEP = 120.0
DUR = 7.0 * DAY
GUARD_RTOL, GUARD_ATOL = 1e-13, 1e-6

ALPHAS_LADDER = (1.00, 1.10, 1.20, 1.30, 1.50)
ALPHAS_EXACT = (1.00, 1.20)

METRICS = Path(__file__).resolve().parents[1] / "metrics"
STAGE2 = METRICS / "r7_doe_matrix_stage2.json"
OUT_NAME = "r8_alpha_margin.json"
GATE_REL_TOL = 1e-9

# Floor proxy for the guard stage: the largest 7-day truth envelope measured
# in robustness_numerical_floor_check.json is ~13 m; pairs with
# E_sched < 5 * proxy get their own envelope measured directly.
FLOOR_PROXY_M = 13.0


def gravity_file_sha256() -> dict:
    from lunaris.common.lunar_data import resolve_lunar_gravity_path
    p = Path(resolve_lunar_gravity_path(None))
    return {"gravity_file": str(p),
            "gravity_file_sha256":
                hashlib.sha256(p.read_bytes()).hexdigest()}


def up_q(n: int, q: int = Q) -> int:
    return ((n + q - 1) // q) * q


def ladder_table(base_emp: dict, alpha: float, truth_deg: int) -> dict:
    """Inflate the archived empirical-lookup table by alpha, up-quantize to
    the q=10 ladder, clamp to [FLOOR, truth_deg - 10]."""
    cap_t = truth_deg - Q
    return {h: max(FLOOR, min(cap_t, up_q(math.ceil(alpha * n0))))
            for h, n0 in base_emp.items()}


def exact_table(power: np.ndarray, r_ref: float, alpha: float,
                truth_deg: int) -> dict:
    """No-ladder control: inflate the raw (unquantized) empirical N_min."""
    cap_t = truth_deg - Q
    table = {}
    for hk in H_GRID_KM:
        n_raw = min(CAP, emp_nmin_exact(power, r_ref, hk * 1e3))
        table[float(hk)] = max(FLOOR, min(cap_t, math.ceil(alpha * n_raw)))
    return table


def load_stage2() -> dict:
    d = json.loads(STAGE2.read_text())
    return {r["name"]: r for r in d["rows"]}


def run_one(model, y0, dur, t_grid, degfun, args,
            rtol=RTOL, atol=ATOL) -> tuple:
    sol, rhs, wall = propagate(model, y0, dur, t_grid, degfun, args,
                               rtol, atol)
    return sol.y, rhs, wall


def policy_stats(sol_y, truth_y, rhs, wall) -> dict:
    st = err_stats(sol_y, truth_y)
    st.update({
        "n_rhs": rhs.n_calls, "wall_s": wall, "grav_s": rhs.grav_ns / 1e9,
        "mean_deg_sq": rhs.sum_deg_sq / max(rhs.n_calls, 1),
        "n_deg_changes": rhs.n_deg_changes,
        "in_track_fraction":
            st["ric_rms_m"]["in_track"] / max(st["pos_rms_m"], 1e-300),
    })
    return st


def run_orbit(orb: dict, arch: dict, models: dict, powers: dict,
              base_tables: dict, dur: float, is_smoke: bool,
              fam_name: str = "ladder", arch_key: str = "sched_emp",
              include_exact: bool = True) -> dict:
    truth_deg = 600 if orb["hp_km"] < 50.0 else 300
    model, args = models[truth_deg]
    y0 = initial_state(model, orb)
    t_grid = np.arange(0.0, dur + 1.0, OUT_STEP)
    base_emp = base_tables[truth_deg]

    n_crit = min(CAP, emp_nmin_exact(powers[truth_deg], model.r_ref,
                                     orb["hp_km"] * 1e3))
    row = {k: orb[k] for k in ("name", "family", "hp_km", "ha_km",
                               "incl_deg", "argp_deg", "raan_deg")}
    row["truth_degree"] = truth_deg
    row["n_critical_empirical"] = int(n_crit)

    a_row = arch[orb["name"]]
    if int(n_crit) != int(a_row["n_critical_empirical"]):
        raise RuntimeError(f"{orb['name']}: n_crit {n_crit} != archived "
                           f"{a_row['n_critical_empirical']}")
    e_work = a_row["policies"]["fixed_work"]["pos_rms_m"]
    e_crit = a_row["policies"]["fixed_crit"]["pos_rms_m"]
    row["archived"] = {
        f"{arch_key}_pos_rms_m": a_row["policies"][arch_key]["pos_rms_m"],
        "fixed_work_pos_rms_m": e_work,
        "fixed_crit_pos_rms_m": e_crit,
        "fixed_crit_grav_s": a_row["policies"]["fixed_crit"]["grav_s"],
        "n_work_matched": a_row["n_work_matched"],
    }

    t0 = time.perf_counter()
    truth_y, truth_rhs, truth_wall = run_one(
        model, y0, dur, t_grid, lambda t, h: truth_deg, args)
    print(f"  truth N={truth_deg}: wall {truth_wall:7.1f}s "
          f"rhs {truth_rhs.n_calls}", flush=True)

    # fixed_crit re-run: same-session gravity-time baseline + repro check
    crit_y, crit_rhs, crit_wall = run_one(
        model, y0, dur, t_grid, lambda t, h, n=int(n_crit): n, args)
    crit_st = policy_stats(crit_y, truth_y, crit_rhs, crit_wall)
    row["fixed_crit_rerun"] = crit_st
    if not is_smoke:
        rel = abs(crit_st["pos_rms_m"] - e_crit) / max(e_crit, 1e-300)
        row["fixed_crit_rerun"]["rel_diff_vs_archived"] = rel
        if rel > 1e-6:
            print(f"  WARNING fixed_crit repro rel diff {rel:.3e}",
                  flush=True)
    g_crit = crit_st["grav_s"]
    print(f"  fixed_crit N={int(n_crit)}: rms {crit_st['pos_rms_m']:11.3f} m"
          f"  grav {g_crit:7.1f}s", flush=True)

    row["runs"] = {}

    def run_family(fam: str, alpha: float, table: dict) -> dict:
        sol_y, rhs, wall = run_one(model, y0, dur, t_grid,
                                   alt_sched(table), args)
        st = policy_stats(sol_y, truth_y, rhs, wall)
        e_s = st["pos_rms_m"]
        st.update({
            "family": fam, "alpha": alpha,
            "rho_work": e_work / e_s, "rho_crit": e_crit / e_s,
            "grav_time_ratio_vs_crit_rerun": st["grav_s"] / g_crit,
            "grav_time_saving_vs_crit_rerun": 1.0 - st["grav_s"] / g_crit,
            "degree_table_min": min(table.values()),
            "degree_table_max": max(table.values()),
        })
        row["runs"][f"{fam}_a{alpha:.2f}"] = st
        print(f"  {fam:6s} a={alpha:.2f}: rms {e_s:11.3f} m  "
              f"rho_w {st['rho_work']:7.3f} rho_c {st['rho_crit']:7.3f}  "
              f"grav {st['grav_s']:7.1f}s (x{st['grav_time_ratio_vs_crit_rerun']:.2f})",
              flush=True)
        return st

    for alpha in ALPHAS_LADDER:
        table = ladder_table(base_emp, alpha, truth_deg)
        if alpha == 1.00 and table != base_emp:
            raise RuntimeError("alpha=1.00 ladder table != archived "
                               f"{arch_key} table")
        st = run_family(fam_name, alpha, table)
        # ----- hard sanity gate: alpha=1.00 == archived base schedule ------
        if alpha == 1.00 and not is_smoke:
            e_arch = a_row["policies"][arch_key]["pos_rms_m"]
            rel = abs(st["pos_rms_m"] - e_arch) / max(e_arch, 1e-300)
            st["gate_rel_diff_vs_archived"] = rel
            st["gate_passed"] = rel <= GATE_REL_TOL
            if not st["gate_passed"]:
                raise RuntimeError(
                    f"HARD GATE FAILED on {orb['name']}: alpha=1.00 gives "
                    f"{st['pos_rms_m']!r} m, archived {arch_key} is "
                    f"{e_arch!r} m (rel {rel:.3e} > {GATE_REL_TOL})")
            print(f"  gate OK: alpha=1.00 reproduces archived {arch_key} "
                  f"(rel {rel:.1e})", flush=True)

    if include_exact:
        for alpha in ALPHAS_EXACT:
            run_family("exact", alpha,
                       exact_table(powers[truth_deg], model.r_ref, alpha,
                                   truth_deg))

    # per-orbit smallest sufficient margin (primary family, vs comparators)
    alpha_star = None
    for alpha in ALPHAS_LADDER:
        if row["runs"][f"{fam_name}_a{alpha:.2f}"]["pos_rms_m"] <= e_work:
            alpha_star = alpha
            break
    row["alpha_star_vs_work"] = alpha_star
    alpha_star_c = None
    for alpha in ALPHAS_LADDER:
        if row["runs"][f"{fam_name}_a{alpha:.2f}"]["pos_rms_m"] <= e_crit:
            alpha_star_c = alpha
            break
    row["alpha_star_vs_crit"] = alpha_star_c
    row["orbit_wall_s"] = time.perf_counter() - t0
    return row


def summarize(rows: list[dict],
              families=(("ladder", ALPHAS_LADDER),
                        ("exact", ALPHAS_EXACT))) -> dict:
    out = {}
    for fam, alphas in families:
        for alpha in alphas:
            key = f"{fam}_a{alpha:.2f}"
            sel = [r for r in rows if key in r.get("runs", {})]
            if not sel:
                continue
            rw = [r["runs"][key]["rho_work"] for r in sel]
            rc = [r["runs"][key]["rho_crit"] for r in sel]
            gr = [r["runs"][key]["grav_time_ratio_vs_crit_rerun"]
                  for r in sel]
            itf = [r["runs"][key]["in_track_fraction"] for r in sel]
            out[key] = {
                "n_orbits": len(sel),
                "win_rate_vs_work": float(np.mean([v > 1.0 for v in rw])),
                "win_rate_vs_crit": float(np.mean([v > 1.0 for v in rc])),
                "median_rho_work": float(np.median(rw)),
                "p10_rho_work": float(np.percentile(rw, 10)),
                "p90_rho_work": float(np.percentile(rw, 90)),
                "median_rho_crit": float(np.median(rc)),
                "p10_rho_crit": float(np.percentile(rc, 10)),
                "p90_rho_crit": float(np.percentile(rc, 90)),
                "median_grav_time_ratio_vs_crit": float(np.median(gr)),
                "p10_grav_time_ratio_vs_crit": float(np.percentile(gr, 10)),
                "p90_grav_time_ratio_vs_crit": float(np.percentile(gr, 90)),
                "median_in_track_fraction": float(np.median(itf)),
            }
    stars = [r["alpha_star_vs_work"] for r in rows]
    out["alpha_star_vs_work"] = {
        "values": stars,
        "n_recovered": int(sum(1 for s in stars if s is not None)),
        "n_not_recovered_at_1p5":
            int(sum(1 for s in stars if s is None)),
    }
    out["alpha_star_vs_crit"] = {
        "values": [r["alpha_star_vs_crit"] for r in rows]}
    return out


def scenario_block(is_smoke: bool, dur: float) -> dict:
    sc = {
        "schema": "r8_alpha_margin_v1",
        "purpose": "alpha-margin robustness control: does a uniform safety "
                   "margin on the archived empirical-lookup altitude "
                   "schedule recover fixed-degree fidelity, and at what "
                   "gravity-time cost?",
        "duration_s": dur, "output_step_s": OUT_STEP,
        "integrator": "DOP853", "rtol": RTOL, "atol": ATOL,
        "max_step_s": None,
        "rotation": "uniform sidereal about polar axis; epoch degenerate "
                    "with body-fixed perilune longitude",
        "truth_rule": "N=600 if perilune < 50 km else N=300 (recomputed "
                      "deterministically; alpha=1.00 gate certifies "
                      "reproduction of the archived Stage-2 numbers)",
        "seed": SEED,
        "alphas_ladder": list(ALPHAS_LADDER),
        "alphas_exact": list(ALPHAS_EXACT),
        "margin_family": {
            "ladder": "N_a(h) = clamp(up_q10(ceil(a*N0(h))), 60, truth-10); "
                      "N0 = archived empirical-lookup table "
                      "(eps=1e-3, floor 60, cap 250, q10 down)",
            "exact": "N_a(h) = clamp(ceil(a*min(250, N_emp_min_raw(h))), "
                     "60, truth-10); no-ladder control",
        },
        "comparators": "archived Stage-2 per-orbit fixed_work / fixed_crit "
                       "errors reused; fixed_crit re-run in-session only "
                       "for like-for-like gravity timing",
        "rho_convention": "rho_work = E_work/E_sched, rho_crit = "
                          "E_crit/E_sched; rho > 1 means schedule wins "
                          "(paper Table 11 convention)",
        "deviations_from_review_sheet": {
            "eps": "review sheet said eps=1e-2; the archived Stage-2 "
                   "sched_emp (and hence the alpha=1.00 gate) uses "
                   "eps=1e-3, which is what Eq.(4)'s empirical lookup used "
                   "in the DOE. eps=1e-3 kept.",
            "tolerances": "review sheet asked for vector atol 1e-6/1e-9 + "
                          "60 s max step ('7-day tight matrix'); the "
                          "archived 24-orbit matrix ran at rtol 1e-12 / "
                          "scalar atol 1e-5 / no max step, and the "
                          "alpha=1.00 reproduction gate is only "
                          "satisfiable at those settings. Archived "
                          "settings kept; floors handled by the --guard "
                          "stage at rtol 1e-13 / atol 1e-6.",
        },
        "guard_stage": {"rtol": GUARD_RTOL, "atol": GUARD_ATOL,
                        "floor_proxy_m": FLOOR_PROXY_M,
                        "flag_rule": "E_sched < 5 * floor_proxy"},
        "smoke": is_smoke,
    }
    sc.update(gravity_file_sha256())
    import numba
    import scipy
    sc["versions"] = {"numpy": np.__version__, "scipy": scipy.__version__,
                      "numba": numba.__version__}
    return sc


def main_run(is_smoke: bool, secondary: bool = False) -> int:
    dur = 2.0 * 3600.0 if is_smoke else DUR
    if secondary:
        out_name = ("r8_alpha_margin_down_smoke.json" if is_smoke
                    else "r8_alpha_margin_down.json")
        fam_name, arch_key, include_exact = "kdown", "sched_down", False
        families = ((fam_name, ALPHAS_LADDER),)
    else:
        out_name = "r8_alpha_margin_smoke.json" if is_smoke else OUT_NAME
        fam_name, arch_key, include_exact = "ladder", "sched_emp", True
        families = (("ladder", ALPHAS_LADDER), ("exact", ALPHAS_EXACT))

    orbits = build_design()
    if is_smoke:
        orbits = [o for o in orbits
                  if o["name"] in ("c2_50x300_polar", "c6_lro_30x216")]
    arch = load_stage2()

    need = sorted({600 if o["hp_km"] < 50.0 else 300 for o in orbits})
    models, powers, base_tables = {}, {}, {}
    for d in need:
        m = load_model(d)
        ar = kernel_args(m)
        warmup(m, ar)
        models[d] = (m, ar)
        powers[d] = degree_power(m)
        base_tables[d] = (kaula_table(m, "down") if secondary
                          else emp_table(m, powers[d]))
        print(f"[model] degree {d} loaded", flush=True)

    scenario = scenario_block(is_smoke, dur)
    if secondary:
        scenario["margin_family"] = {
            "kdown": "SECONDARY worst-case baseline: N_a(h) = "
                     "clamp(up_q10(ceil(a*N0(h))), 60, truth-10); N0 = "
                     "archived Kaula-criterion downward-quantized table "
                     "(eps=1e-3, floor 60, cap 250, q10 down); alpha=1.00 "
                     "reproduces archived sched_down"}
    rows = []
    for k, orb in enumerate(orbits):
        print(f"[{k + 1}/{len(orbits)}] {orb['name']} "
              f"hp={orb['hp_km']:.1f} ha={orb['ha_km']:.1f} "
              f"i={orb['incl_deg']:.1f}", flush=True)
        rows.append(run_orbit(orb, arch, models, powers, base_tables, dur,
                              is_smoke, fam_name=fam_name, arch_key=arch_key,
                              include_exact=include_exact))
        dump(out_name, {"scenario": scenario, "rows": rows,
                        "complete": k + 1 == len(orbits)})

    summary = summarize(rows, families)
    dump(out_name, {"scenario": scenario, "rows": rows, "summary": summary,
                    "complete": True})
    print("summary:", json.dumps(summary, indent=1), flush=True)
    return 0


def main_guard(secondary: bool = False) -> int:
    """Numerical-floor pass: re-integrate flagged (orbit, alpha) pairs at
    one tighter tolerance step and record the run-to-run envelope."""
    src = "r8_alpha_margin_down.json" if secondary else OUT_NAME
    out = ("r8_alpha_margin_down_guard.json" if secondary
           else "r8_alpha_margin_guard.json")
    d = json.loads((METRICS / src).read_text())
    thr = 5.0 * FLOOR_PROXY_M
    flagged = []
    for r in d["rows"]:
        keys = [k for k, st in r["runs"].items()
                if st["pos_rms_m"] < thr]
        if keys:
            flagged.append((r, keys))
    print(f"[guard] {sum(len(k) for _, k in flagged)} flagged pairs on "
          f"{len(flagged)} orbits (E < {thr:.0f} m)", flush=True)

    orbits = {o["name"]: o for o in build_design()}
    need = sorted({r["truth_degree"] for r, _ in flagged} or {300})
    models, powers, base_tables = {}, {}, {}
    for deg in need:
        m = load_model(deg)
        ar = kernel_args(m)
        warmup(m, ar)
        models[deg] = (m, ar)
        powers[deg] = degree_power(m)
        base_tables[deg] = (kaula_table(m, "down") if secondary
                            else emp_table(m, powers[deg]))
        print(f"[model] degree {deg} loaded", flush=True)

    guard_rows = []
    for r, keys in flagged:
        orb = orbits[r["name"]]
        truth_deg = r["truth_degree"]
        model, args = models[truth_deg]
        y0 = initial_state(model, orb)
        t_grid = np.arange(0.0, DUR + 1.0, OUT_STEP)

        def rms(a, b):
            dd = np.linalg.norm(a[:3] - b[:3], axis=0)
            return float(np.sqrt(np.mean(dd * dd)))

        truth_b, _, _ = run_one(model, y0, DUR, t_grid,
                                lambda t, h: truth_deg, args)
        truth_t, _, _ = run_one(model, y0, DUR, t_grid,
                                lambda t, h: truth_deg, args,
                                GUARD_RTOL, GUARD_ATOL)
        g = {"name": r["name"], "truth_degree": truth_deg,
             "truth_envelope_rms_m": rms(truth_b, truth_t), "pairs": {}}
        print(f"[guard] {r['name']}: truth envelope "
              f"{g['truth_envelope_rms_m']:.3f} m", flush=True)
        for key in keys:
            st = r["runs"][key]
            fam, alpha = st["family"], st["alpha"]
            if fam in ("ladder", "kdown"):
                table = ladder_table(base_tables[truth_deg], alpha,
                                     truth_deg)
            else:
                table = exact_table(powers[truth_deg], model.r_ref, alpha,
                                    truth_deg)
            y_b, _, _ = run_one(model, y0, DUR, t_grid, alt_sched(table),
                                args)
            y_t, _, _ = run_one(model, y0, DUR, t_grid, alt_sched(table),
                                args, GUARD_RTOL, GUARD_ATOL)
            e_num = rms(y_b, y_t)
            e_pol = rms(y_b, truth_b)
            g["pairs"][key] = {
                "E_policy_rms_m": e_pol,
                "E_numerical_envelope_rms_m": e_num,
                "ratio_policy_over_envelope": e_pol / max(e_num, 1e-12),
                "resolved_above_floor": e_pol > 5.0 * e_num,
            }
            print(f"  {key}: E {e_pol:9.3f} m  E_num {e_num:9.3f} m  "
                  f"ratio {e_pol / max(e_num, 1e-12):7.1f}", flush=True)
        guard_rows.append(g)
        dump(out, {
            "schema": "r8_alpha_margin_guard_v1",
            "flag_rule": f"E_sched < {thr} m",
            "baseline": {"rtol": RTOL, "atol": ATOL},
            "tight": {"rtol": GUARD_RTOL, "atol": GUARD_ATOL},
            "rows": guard_rows,
            "complete": len(guard_rows) == len(flagged)})
    print("[guard] done", flush=True)
    return 0


def _uncapped_ladder_value(n0: int, alpha: float) -> int:
    return max(FLOOR, up_q(math.ceil(alpha * n0)))


def main_workmatch() -> int:
    """Per-(orbit, family, alpha) re-matched work comparators and pairwise
    resolution envelopes.

    Fixes the equal-work criticism: the archived fixed_work degree is matched
    to the *nominal* (alpha=1) schedule's mean N^2, so for alpha>1 the
    inflated schedule does more work than that comparator and
    rho_work = E_work/E_sched overstates the schedule.  Here the comparator
    is re-matched per run, N_work(alpha) = round(sqrt(<N_alpha^2>)) from the
    archived alpha-run's measured mean_deg_sq, then propagated at baseline
    AND tight tolerance so that rankings can use the pairwise criterion
    |E_A - E_B| > E_num,A + E_num,B.  fixed_crit is re-run at both
    tolerances for the same reason.  The truth run's altitude profile also
    gives, per (family, alpha), the fraction of arc time during which the
    truth-degree cap is binding (cap-artifact audit).

    Output: metrics/r8_alpha_margin_workmatch.json
    """
    srcs = {"ladder": OUT_NAME, "exact": OUT_NAME,
            "kdown": "r8_alpha_margin_down.json"}
    guard_files = {"ladder": "r8_alpha_margin_guard.json",
                   "exact": "r8_alpha_margin_guard.json",
                   "kdown": "r8_alpha_margin_down_guard.json"}
    data, genv = {}, {}
    for fam, f in srcs.items():
        p = METRICS / f
        if p.exists():
            d = json.loads(p.read_text())
            if d.get("complete"):
                data[fam] = {r["name"]: r for r in d["rows"]}
    for fam, f in set(guard_files.items()):
        p = METRICS / f
        if p.exists():
            for gr in json.loads(p.read_text())["rows"]:
                for k, pr in gr["pairs"].items():
                    genv[(gr["name"], k)] = \
                        pr["E_numerical_envelope_rms_m"]
    print(f"[workmatch] families: {sorted(data)}", flush=True)

    orbits = build_design()
    need = sorted({600 if o["hp_km"] < 50.0 else 300 for o in orbits})
    models, powers, tabs = {}, {}, {}
    for deg in need:
        m = load_model(deg)
        ar = kernel_args(m)
        warmup(m, ar)
        models[deg] = (m, ar)
        powers[deg] = degree_power(m)
        tabs[("emp", deg)] = emp_table(m, powers[deg])
        tabs[("kdown", deg)] = kaula_table(m, "down")
        print(f"[model] degree {deg} loaded", flush=True)

    def rms(a, b):
        dd = np.linalg.norm(a[:3] - b[:3], axis=0)
        return float(np.sqrt(np.mean(dd * dd)))

    out_rows = []
    for oi, orb in enumerate(orbits):
        truth_deg = 600 if orb["hp_km"] < 50.0 else 300
        model, args = models[truth_deg]
        y0 = initial_state(model, orb)
        t_grid = np.arange(0.0, DUR + 1.0, OUT_STEP)
        cap_t = truth_deg - Q
        print(f"[{oi + 1}/24] {orb['name']}", flush=True)

        truth_y, _, tw = run_one(model, y0, DUR, t_grid,
                                 lambda t, h: truth_deg, args)
        alt_km = (np.linalg.norm(truth_y[:3], axis=0) - model.r_ref) / 1e3
        print(f"  truth: {tw:.0f}s", flush=True)

        row = {"name": orb["name"], "truth_degree": truth_deg,
               "cap_degree": cap_t, "families": {}, "cap_audit": {},
               "fixed_runs": {}}

        # ---- cap audit from the truth altitude profile -------------------
        hbins = np.clip(10.0 * np.floor(alt_km / 10.0),
                        float(min(H_GRID_KM)), float(max(H_GRID_KM)))
        for fam in data:
            base = tabs[("kdown" if fam == "kdown" else "emp", truth_deg)]
            arow = data[fam][orb["name"]]
            alphas = sorted({arow["runs"][k]["alpha"]
                             for k in arow["runs"]
                             if arow["runs"][k]["family"] == fam})
            for alpha in alphas:
                if fam == "exact":
                    des = {h: math.ceil(alpha * min(CAP, emp_nmin_exact(
                        powers[truth_deg], model.r_ref, h * 1e3)))
                        for h in base}
                else:
                    des = {h: _uncapped_ladder_value(n0, alpha)
                           for h, n0 in base.items()}
                binding = np.array([des[float(h)] > cap_t for h in hbins])
                row["cap_audit"][f"{fam}_a{alpha:.2f}"] = {
                    "frac_time_cap_binding": float(np.mean(binding)),
                    "max_desired_in_visited_range":
                        int(max(des[float(h)] for h in set(hbins))),
                }

        # ---- comparator degrees to run (dedup) ---------------------------
        n_crit = min(CAP, emp_nmin_exact(powers[truth_deg], model.r_ref,
                                         orb["hp_km"] * 1e3))
        degs = {int(n_crit)}
        wanted = {}
        for fam in data:
            arow = data[fam][orb["name"]]
            for k, st in arow["runs"].items():
                if st["family"] != fam:
                    continue
                nw = int(round(math.sqrt(st["mean_deg_sq"])))
                wanted[(fam, k)] = nw
                degs.add(nw)

        for n in sorted(degs):
            yb, rb, wb = run_one(model, y0, DUR, t_grid,
                                 lambda t, h, nn=n: nn, args)
            yt, _, wt = run_one(model, y0, DUR, t_grid,
                                lambda t, h, nn=n: nn, args,
                                GUARD_RTOL, GUARD_ATOL)
            row["fixed_runs"][str(n)] = {
                "E_vs_truth_rms_m": rms(yb, truth_y),
                "E_num_envelope_rms_m": rms(yb, yt),
                "grav_s": rb.grav_ns / 1e9,
                "n_rhs": rb.n_calls,
            }
            print(f"  fixed N={n:3d}: E {row['fixed_runs'][str(n)]['E_vs_truth_rms_m']:9.2f} m  "
                  f"E_num {row['fixed_runs'][str(n)]['E_num_envelope_rms_m']:8.2f} m",
                  flush=True)

        # ---- pairwise resolution per (family, alpha) ---------------------
        for (fam, k), nw in wanted.items():
            arow = data[fam][orb["name"]]
            st = arow["runs"][k]
            e_s = st["pos_rms_m"]
            fr = row["fixed_runs"][str(nw)]
            e_w, en_w = fr["E_vs_truth_rms_m"], fr["E_num_envelope_rms_m"]
            fc = row["fixed_runs"][str(int(n_crit))]
            e_c, en_c = fc["E_vs_truth_rms_m"], fc["E_num_envelope_rms_m"]
            en_s = genv.get((orb["name"], k))
            rec = {
                "n_work_alpha": nw,
                "E_sched_rms_m": e_s,
                "E_work_alpha_rms_m": e_w,
                "E_num_work_rms_m": en_w,
                "E_crit_rms_m": e_c,
                "E_num_crit_rms_m": en_c,
                "E_num_sched_rms_m": en_s,
                "rho_work_alpha": e_w / e_s,
                "rho_crit": e_c / e_s,
                "grav_s_sched": st["grav_s"],
                "grav_s_work_alpha": fr["grav_s"],
            }
            # 100 m conservatively bounds every schedule envelope observed
            # in the guard passes (max 94 m); used only to decide whether an
            # unmeasured envelope must be measured / can be bounded away.
            ENV_BOUND = 100.0
            if en_s is None:
                # schedule envelope unmeasured (error was above the flag
                # threshold); pairwise calls need it only if contested
                contested_w = abs(e_s - e_w) < en_w + ENV_BOUND
                contested_c = abs(e_s - e_c) < en_c + ENV_BOUND
                if contested_w or contested_c:
                    if fam == "kdown":
                        table = ladder_table(tabs[("kdown", truth_deg)],
                                             st["alpha"], truth_deg)
                    elif fam == "ladder":
                        table = ladder_table(tabs[("emp", truth_deg)],
                                             st["alpha"], truth_deg)
                    else:
                        table = exact_table(powers[truth_deg], model.r_ref,
                                            st["alpha"], truth_deg)
                    ysb, _, _ = run_one(model, y0, DUR, t_grid,
                                        alt_sched(table), args)
                    yst, _, _ = run_one(model, y0, DUR, t_grid,
                                        alt_sched(table), args,
                                        GUARD_RTOL, GUARD_ATOL)
                    en_s = rms(ysb, yst)
                    rec["E_num_sched_rms_m"] = en_s
                    rec["sched_envelope_source"] = "measured_here"
                else:
                    rec["sched_envelope_source"] = "not_needed_gap_large"
            else:
                rec["sched_envelope_source"] = "guard"
            if en_s is not None:
                rec["work_rank_resolved"] = bool(
                    abs(e_s - e_w) > en_s + en_w)
                rec["crit_rank_resolved"] = bool(
                    abs(e_s - e_c) > en_s + en_c)
            else:
                rec["work_rank_resolved"] = bool(
                    abs(e_s - e_w) > en_w + ENV_BOUND)
                rec["crit_rank_resolved"] = bool(
                    abs(e_s - e_c) > en_c + ENV_BOUND)
            rec["sched_wins_work_alpha"] = bool(e_s < e_w)
            rec["sched_wins_crit"] = bool(e_s < e_c)
            row["families"][k] = rec

        out_rows.append(row)
        dump("r8_alpha_margin_workmatch.json", {
            "schema": "r8_alpha_margin_workmatch_v1",
            "purpose": "per-alpha re-matched work comparators, pairwise "
                       "rank resolution, cap-binding audit",
            "comparator_rule": "N_work(alpha) = round(sqrt(mean N^2 over "
                               "RHS calls of the archived alpha run)); "
                               "baseline+tight runs give E_num for the "
                               "pairwise criterion |E_A-E_B| > "
                               "E_num,A + E_num,B",
            "baseline": {"rtol": RTOL, "atol": ATOL},
            "tight": {"rtol": GUARD_RTOL, "atol": GUARD_ATOL},
            "rows": out_rows,
            "complete": len(out_rows) == len(orbits)})
    print("[workmatch] done", flush=True)
    return 0


def main_resolution_v2() -> int:
    """Complete all schedule envelopes and apply truth-inclusive resolution.

    The original workmatch pass measured every fixed comparator at both
    tolerances and measured schedule envelopes whenever a comparison was
    contested. This post-review pass fills the remaining large-gap schedule
    envelopes, imports the independently measured truth self-difference for
    every orbit from the guard artifacts, and writes a non-destructive v2
    artifact. For policy P, E_num,P = E_self,P + E_self,truth, so a pair is
    resolved only when

        |E_A-E_B| > E_self,A + E_self,B + 2 E_self,truth.
    """
    source_path = METRICS / "r8_alpha_margin_workmatch.json"
    out_name = "r8_alpha_margin_workmatch_v2.json"
    out_path = METRICS / out_name
    source = json.loads(source_path.read_text())
    if not source.get("complete"):
        raise SystemExit("incomplete workmatch source JSON")

    if out_path.exists():
        candidate = json.loads(out_path.read_text())
        rows = candidate["rows"] if candidate.get("schema") == \
            "r8_alpha_margin_workmatch_v2" else source["rows"]
    else:
        rows = source["rows"]
    by_name = {r["name"]: r for r in rows}

    truth_env = {}
    for filename in ("r8_alpha_margin_guard.json",
                     "r8_alpha_margin_down_guard.json"):
        guard = json.loads((METRICS / filename).read_text())
        if not guard.get("complete"):
            raise SystemExit(f"incomplete guard artifact: {filename}")
        for row in guard["rows"]:
            value = float(row["truth_envelope_rms_m"])
            previous = truth_env.get(row["name"])
            if previous is not None and not math.isclose(
                    previous, value, rel_tol=1e-12, abs_tol=1e-9):
                raise SystemExit(
                    f"inconsistent truth envelope for {row['name']}")
            truth_env[row["name"]] = value

    design = {o["name"]: o for o in build_design()}
    if set(truth_env) != set(design):
        raise SystemExit("truth-envelope coverage does not match the DOE")

    srcs = {"ladder": OUT_NAME, "exact": OUT_NAME,
            "kdown": "r8_alpha_margin_down.json"}
    data = {}
    for fam, filename in srcs.items():
        artifact = json.loads((METRICS / filename).read_text())
        if not artifact.get("complete"):
            raise SystemExit(f"incomplete alpha artifact: {filename}")
        data[fam] = {r["name"]: r for r in artifact["rows"]}

    def save(complete: bool) -> None:
        dump(out_name, {
            "schema": "r8_alpha_margin_workmatch_v2",
            "purpose": "complete two-tolerance schedule/comparator envelopes "
                       "and truth-inclusive pairwise rank resolution",
            "source_artifact": source_path.name,
            "truth_envelope_sources": ["r8_alpha_margin_guard.json",
                                       "r8_alpha_margin_down_guard.json"],
            "comparator_rule": "N_work(alpha) = round(sqrt(mean N^2 over "
                               "RHS calls of the archived alpha run))",
            "resolution_rule": "E_num,P = E_self,P + E_self,truth; resolved "
                               "iff |E_A-E_B| > E_num,A + E_num,B",
            "baseline": {"rtol": RTOL, "atol": ATOL},
            "tight": {"rtol": GUARD_RTOL, "atol": GUARD_ATOL},
            "rows": rows,
            "complete": complete,
        })

    pending = []
    for name, row in by_name.items():
        for key, rec in row["families"].items():
            if rec.get("E_num_sched_rms_m") is None:
                fam, _ = key.rsplit("_a", 1)
                pending.append((name, fam, key))
    print(f"[resolution-v2] {len(pending)} schedule envelopes to fill",
          flush=True)

    models = {}
    powers = {}
    tabs = {}
    needed_degrees = sorted({int(by_name[name]["truth_degree"])
                             for name, _, _ in pending})
    for degree in needed_degrees:
        model = load_model(degree)
        args = kernel_args(model)
        warmup(model, args)
        models[degree] = (model, args)
        powers[degree] = degree_power(model)
        tabs[("emp", degree)] = emp_table(model, powers[degree])
        tabs[("kdown", degree)] = kaula_table(model, "down")

    def rms(a, b):
        dd = np.linalg.norm(a[:3] - b[:3], axis=0)
        return float(np.sqrt(np.mean(dd * dd)))

    for index, (name, fam, key) in enumerate(pending, 1):
        row = by_name[name]
        degree = int(row["truth_degree"])
        model, args = models[degree]
        orbit = design[name]
        y0 = initial_state(model, orbit)
        t_grid = np.arange(0.0, DUR + 1.0, OUT_STEP)
        st = data[fam][name]["runs"][key]
        alpha = float(st["alpha"])
        if fam == "kdown":
            table = ladder_table(tabs[("kdown", degree)], alpha, degree)
        elif fam == "ladder":
            table = ladder_table(tabs[("emp", degree)], alpha, degree)
        else:
            table = exact_table(powers[degree], model.r_ref, alpha, degree)
        y_base, _, _ = run_one(model, y0, DUR, t_grid, alt_sched(table), args)
        y_tight, _, _ = run_one(model, y0, DUR, t_grid, alt_sched(table),
                                args, GUARD_RTOL, GUARD_ATOL)
        envelope = rms(y_base, y_tight)
        rec = row["families"][key]
        rec["E_num_sched_rms_m"] = envelope
        rec["sched_envelope_source"] = "resolution_v2_backfill"
        print(f"[{index}/{len(pending)}] {name} {key}: "
              f"E_self={envelope:.3f} m", flush=True)
        save(False)

    for row in rows:
        truth = truth_env[row["name"]]
        row["truth_envelope_rms_m"] = truth
        for rec in row["families"].values():
            sched = rec.get("E_num_sched_rms_m")
            if sched is None:
                raise SystemExit("schedule-envelope backfill is incomplete")
            work = float(rec["E_num_work_rms_m"])
            crit = float(rec["E_num_crit_rms_m"])
            rec["work_rank_resolved_without_truth"] = \
                bool(rec["work_rank_resolved"])
            rec["crit_rank_resolved_without_truth"] = \
                bool(rec["crit_rank_resolved"])
            rec["E_num_truth_rms_m"] = truth
            rec["E_num_sched_including_truth_rms_m"] = sched + truth
            rec["E_num_work_including_truth_rms_m"] = work + truth
            rec["E_num_crit_including_truth_rms_m"] = crit + truth
            rec["work_rank_resolved"] = bool(
                abs(rec["E_sched_rms_m"] - rec["E_work_alpha_rms_m"])
                > sched + work + 2.0 * truth)
            rec["crit_rank_resolved"] = bool(
                abs(rec["E_sched_rms_m"] - rec["E_crit_rms_m"])
                > sched + crit + 2.0 * truth)
    save(True)
    print(f"[resolution-v2] complete: metrics/{out_name}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="2 orbits x 2 h pipeline shakedown")
    ap.add_argument("--guard", action="store_true",
                    help="numerical-floor pass on flagged pairs")
    ap.add_argument("--secondary", action="store_true",
                    help="worst-case baseline: Kaula downward-quantized "
                         "table x alpha (archived sched_down at alpha=1)")
    ap.add_argument("--workmatch", action="store_true",
                    help="per-alpha re-matched work comparators, pairwise "
                         "resolution envelopes, cap audit")
    ap.add_argument("--resolution-v2", action="store_true",
                    help="fill all schedule envelopes and apply the "
                         "truth-inclusive conservative resolution rule")
    a = ap.parse_args()
    if a.resolution_v2:
        return main_resolution_v2()
    if a.workmatch:
        return main_workmatch()
    if a.guard:
        return main_guard(a.secondary)
    return main_run(a.smoke, a.secondary)


if __name__ == "__main__":
    raise SystemExit(main())
