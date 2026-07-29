"""P0-5: multi-geometry long-arc scheduling matrix, incl. a MOON_PA case.

Seven-day arcs at the tight tolerance (rtol 1e-12, atol 1e-5), each case
with its own N=300 reference and six degree policies:

  M_phaseB  : 50 x 300 km polar, apolune start (phase control)
  M_inc60   : 50 x 300 km, 60-deg inclination, perilune start
  M_100x300 : 100 x 300 km polar, perilune start (altitude-range control)
  M_moonpa  : 50 x 300 km polar, perilune start, full DE440 MOON_PA
              orientation via SPICE (instead of uniform rotation)

Policies: fixed 138 (reference cap degree), fixed 106 (comparable measured
cost), downward dwell-aware schedule, upward-quantized schedule, 600 s
minimum-dwell schedule, empirical-lookup schedule. Direct accepted/rejected
counts from the instrumented solver.
"""

from __future__ import annotations

import math
import time

import numpy as np

from rev3_common import (DAY, OMEGA_MOON, Rhs, dump, err_stats, kernel_args,
                         load_model, degree_power, make_p_table, make_emp_table,
                         alt_sched, mindwell_from_profile, eccentric_state,
                         propagate, propagate_instr, warmup)

import sys
from pathlib import Path

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba  # noqa: E402

RTOL, ATOL = 1e-12, 1e-5
KERNEL_DIR = Path(r"C:\Users\ayber\Desktop\lunaris\data\ephemeris_models")
SPICE_KERNELS = ["naif0012.tls", "pck00011.tpc", "gm_de440.tpc",
                 "de440s.bsp", "moon_de440_250416.tf",
                 "moon_pa_de440_200625.bpc"]
ET0_TDB_S = 788961600.0  # 2025-01-01 00:00:00 TDB, matches the Tudat contract


class RhsMoonPA:
    """Inertial RHS with full DE440 MOON_PA orientation via SPICE pxform."""

    def __init__(self, model, degree_of, args):
        import spiceypy as sp
        self.sp = sp
        self.model = model
        self.degree_of = degree_of
        self.args = args
        self.n_calls = 0
        self.sum_deg_sq = 0.0
        self.grav_ns = 0
        self.deg_counts: dict[int, int] = {}

    def __call__(self, t, y):
        self.n_calls += 1
        x, yy, z, vx, vy, vz = y
        M = self.sp.pxform("J2000", "MOON_PA", ET0_TDB_S + t)
        rb = M @ np.array([x, yy, z])
        r = math.sqrt(x * x + yy * yy + z * z)
        n = self.degree_of(t, r - self.model.r_ref)
        self.sum_deg_sq += float(n) * float(n)
        self.deg_counts[n] = self.deg_counts.get(n, 0) + 1
        t0 = time.perf_counter_ns()
        axb, ayb, azb = sh_accel_fixed_numba(rb[0], rb[1], rb[2], n, *self.args)
        self.grav_ns += time.perf_counter_ns() - t0
        ai = M.T @ np.array([axb, ayb, azb])
        return (vx, vy, vz, ai[0], ai[1], ai[2])


def growth_slope(sol_y, truth_y, t_grid):
    dp = np.linalg.norm(sol_y[:3] - truth_y[:3], axis=0)
    days = np.arange(1, 8)
    rms = []
    for k in days:
        m = (t_grid > (k - 1) * DAY) & (t_grid <= k * DAY)
        rms.append(float(np.sqrt(np.mean(dp[m] ** 2))))
    lr = np.polyfit(np.log(days), np.log(np.maximum(rms, 1e-12)), 1)
    return float(lr[0]), rms


def count_switches(model, t_grid, Y, degfun):
    tg = np.arange(0.0, t_grid[-1], 10.0)
    rg = np.interp(tg, t_grid, np.linalg.norm(Y[:3], axis=0))
    degs = np.array([degfun(float(t), float(r) - model.r_ref)
                     for t, r in zip(tg, rg)])
    return int(np.sum(degs[1:] != degs[:-1]))


def run_case(model, args, y0, dur, t_grid, policies, rhs_cls=Rhs):
    # reference
    if rhs_cls is Rhs:
        Yt, rhs_t, info_t = propagate_instr(model, y0, dur, t_grid,
                                            lambda t, h: 300, args, RTOL, ATOL)
    else:
        from scipy.integrate import solve_ivp
        rhs_t = rhs_cls(model, lambda t, h: 300, args)
        t0 = time.perf_counter()
        sol = solve_ivp(rhs_t, (0.0, dur), y0, method="DOP853",
                        rtol=RTOL, atol=ATOL, t_eval=t_grid)
        info_t = {"n_rhs": rhs_t.n_calls, "wall_s": time.perf_counter() - t0,
                  "grav_s": rhs_t.grav_ns / 1e9}
        Yt = sol.y
    rows = []
    for pname, degfun in policies.items():
        if rhs_cls is Rhs:
            Y, rhs, info = propagate_instr(model, y0, dur, t_grid, degfun,
                                           args, RTOL, ATOL)
        else:
            from scipy.integrate import solve_ivp
            rhs = rhs_cls(model, degfun, args)
            t0 = time.perf_counter()
            sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853",
                            rtol=RTOL, atol=ATOL, t_eval=t_grid)
            info = {"n_rhs": rhs.n_calls,
                    "wall_s": time.perf_counter() - t0,
                    "grav_s": rhs.grav_ns / 1e9}
            Y = sol.y
        st = err_stats(Y, Yt)
        slope, daily = growth_slope(Y, Yt, t_grid)
        st.update({
            "run": pname, **info,
            "mean_deg_sq": rhs.sum_deg_sq / rhs.n_calls,
            "proxy_rel_fixed138": (rhs.sum_deg_sq / rhs.n_calls) / 138.0 ** 2,
            "n_switches_10s_grid": count_switches(model, t_grid, Y, degfun),
            "growth_slope_loglog": slope, "daily_rms_m": daily,
        })
        rows.append(st)
        print(f"  {pname}: RMS {st['pos_rms_m']:9.2f} m  final "
              f"{st['pos_final_m']:9.2f} m  slope {slope:.2f}  "
              f"wall {info['wall_s']:.1f} s")
    return {"truth_info": info_t, "rows": rows}


def build_policies(model, power, truth_for_mindwell, dur):
    tab_down = make_p_table(model, 1e-3, 60, policy="down")
    tab_up = make_p_table(model, 1e-3, 60, policy="up")
    tab_emp = make_emp_table(model, power, 1e-3, 60)
    mind, seg_t, seg_n = mindwell_from_profile(
        truth_for_mindwell, model, tab_down, dur)
    return {
        "fixed_138": lambda t, h: 138,
        "fixed_106": lambda t, h: 106,
        "sched_down": alt_sched(tab_down),
        "sched_up": alt_sched(tab_up),
        "sched_mindwell600": mind,
        "sched_emp": alt_sched(tab_emp),
    }, {"down": tab_down, "up": tab_up, "emp": tab_emp,
        "mindwell_segments": {"t_s": seg_t, "degree": seg_n}}


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    power = degree_power(load_model(1800))
    dur = 7.0 * DAY
    t_grid = np.arange(0.0, dur + 1.0, 120.0)

    cases = {
        "M_phaseB": dict(y0=eccentric_state(model, 50.0, 300.0,
                                            at_apolune=True), rhs=Rhs),
        "M_inc60": dict(y0=eccentric_state(model, 50.0, 300.0,
                                           incl_deg=60.0), rhs=Rhs),
        "M_100x300": dict(y0=eccentric_state(model, 100.0, 300.0), rhs=Rhs),
        "M_moonpa": dict(y0=eccentric_state(model, 50.0, 300.0),
                         rhs=RhsMoonPA),
    }

    out = {}
    sched_meta = {}
    for cname, spec in cases.items():
        print(f"== case {cname} ==")
        if spec["rhs"] is RhsMoonPA:
            import spiceypy as sp
            for k in SPICE_KERNELS:
                sp.furnsh(str(KERNEL_DIR / k))
        # a quick uniform-rotation truth just to precompute the min-dwell
        # time schedule from the altitude profile (open-loop in time)
        truth_sched, _, _ = propagate(model, spec["y0"], dur, t_grid,
                                      lambda t, h: 60, args, 1e-9, 1e-2)
        policies, meta = build_policies(model, power, truth_sched, dur)
        sched_meta[cname] = meta
        out[cname] = run_case(model, args, spec["y0"], dur, t_grid,
                              policies, rhs_cls=spec["rhs"])

    dump("r3_longarc_matrix.json", {
        "scenario": {"duration_s": dur, "integrator": "DOP853",
                     "rtol": RTOL, "atol": ATOL, "output_step_s": 120.0,
                     "truth_degree": 300,
                     "cases": {
                         "M_phaseB": "50x300 polar, apolune start, uniform rotation",
                         "M_inc60": "50x300, i=60 deg, perilune start, uniform rotation",
                         "M_100x300": "100x300 polar, perilune start, uniform rotation",
                         "M_moonpa": "50x300 polar, perilune start, DE440 MOON_PA via SPICE",
                     },
                     "moonpa_epoch_tdb_j2000_s": ET0_TDB_S,
                     "spice_kernels": SPICE_KERNELS},
        "schedules": {k: {kk: ({f"{a:.0f}": b for a, b in vv.items()}
                               if isinstance(vv, dict) and kk != "mindwell_segments"
                               else vv)
                          for kk, vv in v.items()}
                      for k, v in sched_meta.items()},
        "cases": out,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
