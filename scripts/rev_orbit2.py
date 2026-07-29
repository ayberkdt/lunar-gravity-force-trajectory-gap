"""Revision-1 follow-up: tight-tolerance 7-day long arc.

Motivated by the convergence control in rev_orbit.py: at the baseline
tolerance (rtol 1e-11, atol 1e-4) the integrator floor on the 24 h circular
case is 4.4 m RMS, which contaminates both the high-degree end of the
truncation mapping and the 7-day comparisons. This script reruns the 7-day
eccentric arc at rtol 1e-12, atol 1e-5, reports the integration floor on
the same arc, and adds two switch-rate-mitigation variants:

  * coarse3    : three-level altitude schedule (130 below 80 km, 100 in
                 80-150 km, 60 above 150 km), ~4 switches per orbit;
  * mindwell600: the dwell-aware schedule with a 600 s minimum dwell,
                 implemented as a time-based schedule precomputed from the
                 truth altitude profile (open-loop in time, as a mission
                 would schedule from the predicted orbit).
"""

from __future__ import annotations

import bisect
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.common.lunar_data import resolve_lunar_gravity_path  # noqa: E402
from lunaris.common.math_utils import (  # noqa: E402
    LUNAR_DEGREE_POWER_EXPONENT,
    recommended_sh_degree,
)
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel,
    sh_accel_fixed_numba,
)

OUT = Path(__file__).resolve().parents[1] / "metrics"
OMEGA_MOON = 2.0 * math.pi / (27.321661 * 86400.0)
DAY = 86400.0
RTOL, ATOL = 1e-12, 1e-5


def commit_sha() -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return r.stdout.strip()


class Rhs:
    def __init__(self, model, degree_of, args):
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
        th = OMEGA_MOON * t
        c, s = math.cos(th), math.sin(th)
        xb = c * x + s * yy
        yb = -s * x + c * yy
        r = math.sqrt(x * x + yy * yy + z * z)
        n = self.degree_of(t, r - self.model.r_ref)
        self.sum_deg_sq += float(n) * float(n)
        self.deg_counts[n] = self.deg_counts.get(n, 0) + 1
        t0 = time.perf_counter_ns()
        axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *self.args)
        self.grav_ns += time.perf_counter_ns() - t0
        return (vx, vy, vz, c * axb - s * ayb, s * axb + c * ayb, azb)


def kernel_args(model):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


def propagate(model, y0, dur, t_eval, degree_of, args, rtol=RTOL, atol=ATOL):
    rhs = Rhs(model, degree_of, args)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853", rtol=rtol, atol=atol,
                    t_eval=t_eval)
    wall = time.perf_counter() - t0
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs, wall


def ric_errors(sol, truth):
    d = (sol.y[:3] - truth.y[:3]).T
    r = truth.y[:3].T
    v = truth.y[3:].T
    R = r / np.linalg.norm(r, axis=1, keepdims=True)
    C = np.cross(r, v)
    C = C / np.linalg.norm(C, axis=1, keepdims=True)
    I = np.cross(C, R)
    return (np.sum(d * R, axis=1), np.sum(d * I, axis=1), np.sum(d * C, axis=1))


def err_stats(sol, truth):
    dp = np.linalg.norm(sol.y[:3] - truth.y[:3], axis=0)
    dv = np.linalg.norm(sol.y[3:] - truth.y[3:], axis=0)
    er, ei, ec = ric_errors(sol, truth)
    return {
        "pos_rms_m": float(np.sqrt(np.mean(dp ** 2))),
        "pos_max_m": float(np.max(dp)), "pos_final_m": float(dp[-1]),
        "vel_rms_m_s": float(np.sqrt(np.mean(dv ** 2))),
        "vel_max_m_s": float(np.max(dv)), "vel_final_m_s": float(dv[-1]),
        "ric_rms_m": {"radial": float(np.sqrt(np.mean(er ** 2))),
                      "in_track": float(np.sqrt(np.mean(ei ** 2))),
                      "cross_track": float(np.sqrt(np.mean(ec ** 2)))},
    }


def make_p_table(model, eps, floor, cap=138, q=10):
    table = {}
    for hk in np.arange(40.0, 321.0, 10.0):
        n = recommended_sh_degree(hk, model.r_ref,
                                  kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                  kaula_tail_fraction=eps)
        n = max(floor, min(cap, n))
        table[hk] = max(floor, min(cap, (n // q) * q))
    return table


def alt_sched(table):
    def f(t, h_m):
        hb = min(320.0, max(40.0, 10.0 * math.floor(h_m / 1e3 / 10.0)))
        return table[hb]
    return f


def coarse3_sched(t, h_m):
    hk = h_m / 1e3
    if hk < 80.0:
        return 130
    if hk < 150.0:
        return 100
    return 60


def main() -> int:
    commit = commit_sha()
    print("commit:", commit)
    gf = resolve_lunar_gravity_path(None)
    model = GravityModel.from_file(str(gf), requested_degree=300)
    args = kernel_args(model)
    sh_accel_fixed_numba(model.r_ref + 100e3, 0.0, 0.0, 60, *args)

    dur = 7.0 * DAY
    rp, ra = model.r_ref + 50e3, model.r_ref + 300e3
    a_sma = 0.5 * (rp + ra)
    vp = math.sqrt(model.mu * (2.0 / rp - 1.0 / a_sma))
    y0 = np.array([rp, 0.0, 0.0, 0.0, 0.0, vp])
    t_eval = np.arange(0.0, dur + 1.0, 120.0)

    print("truth tight...")
    truth, rhs_t, wall_t = propagate(model, y0, dur, t_eval, lambda t, h: 300, args)
    print(f"  {rhs_t.n_calls} RHS, {wall_t:.1f} s")
    print("truth baseline tolerance (integration floor)...")
    truth_base, _, _ = propagate(model, y0, dur, t_eval, lambda t, h: 300, args,
                                 rtol=1e-11, atol=1e-4)
    floor_stats = err_stats(truth_base, truth)
    print(f"  7-day integration floor: RMS {floor_stats['pos_rms_m']:.2f} m, "
          f"final {floor_stats['pos_final_m']:.2f} m")

    tab_dwell = make_p_table(model, 1e-3, 60)
    tab_naive = make_p_table(model, 1e-2, 37)

    # time-based min-dwell schedule from the truth altitude profile
    t10 = np.arange(0.0, dur + 10.0, 10.0)
    rr = np.interp(t10, truth.t, np.linalg.norm(truth.y[:3], axis=0))
    raw = np.array([alt_sched(tab_dwell)(0.0, float(r) - model.r_ref) for r in rr])
    seg_t, seg_n = [0.0], [int(raw[0])]
    last_switch = 0.0
    for i in range(1, len(t10)):
        if raw[i] != seg_n[-1] and (t10[i] - last_switch) >= 600.0:
            seg_t.append(float(t10[i]))
            seg_n.append(int(raw[i]))
            last_switch = float(t10[i])

    def mindwell(t, h_m):
        i = bisect.bisect_right(seg_t, t) - 1
        return seg_n[max(i, 0)]

    runs = {
        "fixed_138": lambda t, h: 138,
        "fixed_106": lambda t, h: 106,
        "sched_dwell_alt": alt_sched(tab_dwell),
        "sched_naive_alt": alt_sched(tab_naive),
        "sched_dwell_up": alt_sched(make_p_table(model, 1e-3, 60, q=10)),  # placeholder replaced below
        "sched_coarse3": coarse3_sched,
        "sched_mindwell600": mindwell,
    }
    # true up-quantized table
    tab_up = {}
    for hk in np.arange(40.0, 321.0, 10.0):
        n = recommended_sh_degree(hk, model.r_ref,
                                  kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                  kaula_tail_fraction=1e-3)
        n = max(60, min(138, n))
        tab_up[hk] = max(60, min(138, ((n + 9) // 10) * 10))
    runs["sched_dwell_up"] = alt_sched(tab_up)

    rows = []
    series = {}
    for name, degfun in runs.items():
        sol, rhs, wall = propagate(model, y0, dur, t_eval, degfun, args)
        st = err_stats(sol, truth)
        tg = np.arange(0.0, dur, 10.0)
        rg = np.interp(tg, sol.t, np.linalg.norm(sol.y[:3], axis=0))
        degs = np.array([degfun(float(t), float(r) - model.r_ref)
                         for t, r in zip(tg, rg)])
        st.update({
            "run": name, "n_rhs": rhs.n_calls, "wall_s": wall,
            "grav_s": rhs.grav_ns / 1e9,
            "mean_deg_sq": rhs.sum_deg_sq / rhs.n_calls,
            "proxy_rel_fixed138": (rhs.sum_deg_sq / rhs.n_calls) / 138.0 ** 2,
            "n_switches_10s_grid": int(np.sum(degs[1:] != degs[:-1])),
            "dwell_rhs_fraction": {str(k): v / rhs.n_calls
                                   for k, v in sorted(rhs.deg_counts.items())},
        })
        rows.append(st)
        er, ei, ec = ric_errors(sol, truth)
        series[name] = {"t_s": t_eval[::5].tolist(),
                        "radial_m": er[::5].tolist(),
                        "in_track_m": ei[::5].tolist(),
                        "cross_track_m": ec[::5].tolist()}
        print(f"  {name}: RMS {st['pos_rms_m']:8.2f} m  final {st['pos_final_m']:8.2f} m"
              f"  proxy {st['proxy_rel_fixed138']:.3f}"
              f"  switches {st['n_switches_10s_grid']}"
              f"  wall {wall:.1f} s (grav {st['grav_s']:.1f} s)")

    payload = {
        "repo_commit_sha": commit,
        "scenario": {"type": "eccentric_polar", "perilune_km": 50.0,
                     "apolune_km": 300.0, "duration_s": dur,
                     "truth_degree": 300, "integrator": "DOP853",
                     "rtol": RTOL, "atol": ATOL, "output_step_s": 120.0,
                     "rotation": "uniform sidereal about polar axis"},
        "integration_floor_7day_baseline_vs_tight": floor_stats,
        "schedules": {"dwell": {f"{k:.0f}": v for k, v in tab_dwell.items()},
                      "naive": {f"{k:.0f}": v for k, v in tab_naive.items()},
                      "up": {f"{k:.0f}": v for k, v in tab_up.items()},
                      "mindwell_segments": {"t_s": seg_t, "degree": seg_n}},
        "truth": {"n_rhs": rhs_t.n_calls, "wall_s": wall_t,
                  "grav_s": rhs_t.grav_ns / 1e9},
        "rows": rows,
        "error_series": series,
    }
    (OUT / "r1_longarc_tight.json").write_text(json.dumps(payload, indent=2),
                                               encoding="utf-8")
    print("[written] r1_longarc_tight.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
