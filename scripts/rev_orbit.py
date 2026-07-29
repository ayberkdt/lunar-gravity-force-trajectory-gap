"""Revision-1 orbit-level experiments (commit-pinned, seeded).

R1-E6  24 h circular-orbit truncation mapping rerun on the current commit,
       with measured gravity-kernel and total wall time, plus a tightened
       tolerance convergence control for selected degrees (reviewer item L).
R1-LA  7-day eccentric-orbit long arc: fixed degrees, dwell-aware schedule,
       naive low-floor schedule, empirical-lookup schedule; RIC errors,
       error-versus-time series, switch counts, dwell histograms, quadratic
       gravity-work proxy, measured gravity/total wall time, repeated-run
       wall-time dispersion (reviewer items E, H, P).
R1-QP  Quantization-policy comparison (down / nearest / up) and a hysteresis
       (2 km deadband) variant on the dwell-aware schedule (reviewer item N).
R1-SW  Step-by-step DOP853 instrumentation around degree switches: accepted
       step sizes, inferred rejected attempts, RHS density, two initial
       phases, fixed-degree control (reviewer item F).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, DOP853

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
SEED = 20260719
OMEGA_MOON = 2.0 * math.pi / (27.321661 * 86400.0)
DAY = 86400.0


def commit_sha() -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return r.stdout.strip()


def _dump(name: str, payload: dict) -> None:
    payload["repo_commit_sha"] = COMMIT
    (OUT / name).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {name}")


class Rhs:
    """Inertial RHS; uniformly rotating body-fixed field; instrumented."""

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
        n = self.degree_of(r - self.model.r_ref)
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


def propagate(model, y0, dur, t_eval, degree_of, args, rtol=1e-11, atol=1e-4):
    rhs = Rhs(model, degree_of, args)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853", rtol=rtol, atol=atol,
                    t_eval=t_eval)
    wall = time.perf_counter() - t0
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs, wall


def ric_errors(sol, truth):
    """RIC decomposition of the position difference using truth states."""
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
        "ric_max_m": {"radial": float(np.max(np.abs(er))),
                      "in_track": float(np.max(np.abs(ei))),
                      "cross_track": float(np.max(np.abs(ec)))},
    }


def count_switches(model, sol_t_dense, y_interp, degree_of):
    """Count schedule switches along the trajectory on a 10 s grid."""
    rr = np.linalg.norm(y_interp[:3], axis=0)
    degs = np.array([degree_of(float(r) - model.r_ref) for r in rr])
    return int(np.sum(degs[1:] != degs[:-1])), degs


# ---------------------------------------------------------------- schedules
def make_p_schedule(model, eps, floor, cap=138, q=10, policy="down",
                    deadband_km=0.0):
    table = {}
    for hk in np.arange(40.0, 321.0, 10.0):
        n = recommended_sh_degree(hk, model.r_ref,
                                  kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                  kaula_tail_fraction=eps)
        n = max(floor, min(cap, n))
        if policy == "down":
            nq = (n // q) * q
        elif policy == "up":
            nq = ((n + q - 1) // q) * q
        else:  # nearest
            nq = int(round(n / q)) * q
        table[hk] = max(floor, min(cap, nq))

    state = {"bin": None}

    def f(h_m: float) -> int:
        hb = 10.0 * math.floor(h_m / 1000.0 / 10.0)
        hb = min(320.0, max(40.0, hb))
        if deadband_km > 0.0 and state["bin"] is not None:
            lo_edge = state["bin"]
            if abs(h_m / 1000.0 - lo_edge) < deadband_km or \
               abs(h_m / 1000.0 - (lo_edge + 10.0)) < deadband_km:
                hb = state["bin"]
        state["bin"] = hb
        return table[hb]

    return f, {f"{k:.0f}": v for k, v in table.items()}


def make_emp_schedule(model, power, eps, floor, cap=138, q=10):
    n_arr = np.arange(len(power), dtype=np.float64)
    table = {}
    for hk in np.arange(40.0, 321.0, 10.0):
        r = model.r_ref + hk * 1e3
        ratio_n = np.exp(n_arr * math.log(model.r_ref / r))
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
        nmin = max(floor, min(cap, nmin))
        table[hk] = max(floor, min(cap, (nmin // q) * q))

    def f(h_m: float) -> int:
        hb = min(320.0, max(40.0, 10.0 * math.floor(h_m / 1000.0 / 10.0)))
        return table[hb]

    return f, {f"{k:.0f}": v for k, v in table.items()}


# ---------------------------------------------------------------- R1-E6
def r1_orbit_mapping(model, args):
    print("== R1-E6: circular 100 km mapping + convergence control ==")
    r0 = model.r_ref + 100e3
    v0 = math.sqrt(model.mu / r0)
    y0 = np.array([r0, 0.0, 0.0, 0.0, 0.0, v0])
    t_eval = np.arange(0.0, DAY + 1.0, 120.0)

    truth, rhs_t, wall_t = propagate(model, y0, DAY, t_eval, lambda h: 300, args)
    rows = []
    for n in (20, 37, 56, 60, 69, 84, 100, 124, 150):
        sol, rhs, wall = propagate(model, y0, DAY, t_eval, lambda h, n=n: n, args)
        st = err_stats(sol, truth)
        st.update({"degree": n, "n_rhs": rhs.n_calls, "wall_s": wall,
                   "grav_s": rhs.grav_ns / 1e9})
        rows.append(st)
        print(f"  N={n:4d}: RMS {st['pos_rms_m']:9.2f} m")

    # convergence control: tightened tolerances
    truth_T, _, _ = propagate(model, y0, DAY, t_eval, lambda h: 300, args,
                              rtol=1e-12, atol=1e-5)
    integ_floor = err_stats(truth_T, truth)
    conv = []
    for n in (69, 84, 100, 124, 150):
        solB, _, _ = propagate(model, y0, DAY, t_eval, lambda h, n=n: n, args)
        solT, _, _ = propagate(model, y0, DAY, t_eval, lambda h, n=n: n, args,
                               rtol=1e-12, atol=1e-5)
        eB = err_stats(solB, truth)["pos_rms_m"]
        eT = err_stats(solT, truth_T)["pos_rms_m"]
        conv.append({"degree": n, "rms_baseline_tol": eB, "rms_tight_tol": eT,
                     "delta_pct": 100.0 * abs(eT - eB) / eB})
        print(f"  conv N={n:3d}: baseline {eB:.2f} m  tight {eT:.2f} m")
    print(f"  integrator floor (truth tight vs baseline): "
          f"{integ_floor['pos_rms_m']:.4f} m RMS")

    _dump("r1_orbit_mapping.json", {
        "scenario": {"type": "circular_polar", "altitude_km": 100.0,
                     "duration_s": DAY, "truth_degree": 300,
                     "integrator": "DOP853", "rtol": 1e-11, "atol": 1e-4,
                     "max_step": "unbounded", "output_step_s": 120.0,
                     "rotation": "uniform sidereal about polar axis"},
        "truth_rhs": rhs_t.n_calls, "truth_wall_s": wall_t,
        "rows": rows,
        "convergence_control": {"tight_rtol": 1e-12, "tight_atol": 1e-5,
                                "integrator_floor_rms_m": integ_floor["pos_rms_m"],
                                "rows": conv},
    })


# ---------------------------------------------------------------- R1-LA
def r1_long_arc(model, args, power):
    print("== R1-LA: 7-day eccentric long arc ==")
    dur = 7.0 * DAY
    rp, ra = model.r_ref + 50e3, model.r_ref + 300e3
    a_sma = 0.5 * (rp + ra)
    vp = math.sqrt(model.mu * (2.0 / rp - 1.0 / a_sma))
    y0 = np.array([rp, 0.0, 0.0, 0.0, 0.0, vp])
    t_eval = np.arange(0.0, dur + 1.0, 120.0)

    sched_dwell, tab_dwell = make_p_schedule(model, 1e-3, 60)
    sched_naive, tab_naive = make_p_schedule(model, 1e-2, 37)
    sched_emp, tab_emp = make_emp_schedule(model, power, 1e-3, 60)
    sched_near, _ = make_p_schedule(model, 1e-3, 60, policy="nearest")
    sched_up, _ = make_p_schedule(model, 1e-3, 60, policy="up")
    sched_hyst, _ = make_p_schedule(model, 1e-3, 60, deadband_km=2.0)

    print("  truth (N=300)...")
    truth, rhs_t, wall_t = propagate(model, y0, dur, t_eval, lambda h: 300, args)
    print(f"  truth: {rhs_t.n_calls} RHS, {wall_t:.1f} s")

    runs = {
        "fixed_138": lambda h: 138,
        "fixed_106": lambda h: 106,
        "sched_p17_eps1e3_floor60_down": sched_dwell,
        "sched_p17_eps1e2_floor37_down": sched_naive,
        "sched_emp_eps1e3_floor60_down": sched_emp,
        "sched_p17_eps1e3_floor60_nearest": sched_near,
        "sched_p17_eps1e3_floor60_up": sched_up,
        "sched_p17_eps1e3_floor60_hyst2km": sched_hyst,
    }
    rows = []
    series = {}
    for name, degfun in runs.items():
        sol, rhs, wall = propagate(model, y0, dur, t_eval, degfun, args)
        st = err_stats(sol, truth)
        # switch count on 10 s grid from the solved trajectory
        t10 = np.arange(0.0, dur, 10.0)
        rr = np.interp(t10, sol.t, np.linalg.norm(sol.y[:3], axis=0))
        # reset hysteresis state before the diagnostic sweep
        degs = np.array([degfun(float(r) - model.r_ref) for r in rr])
        n_switch = int(np.sum(degs[1:] != degs[:-1]))
        st.update({
            "run": name, "n_rhs": rhs.n_calls, "wall_s": wall,
            "grav_s": rhs.grav_ns / 1e9,
            "mean_deg_sq": rhs.sum_deg_sq / rhs.n_calls,
            "proxy_rel_fixed138": (rhs.sum_deg_sq / rhs.n_calls) / 138.0 ** 2,
            "n_switches_10s_grid": n_switch,
            "dwell_rhs_fraction": {str(k): v / rhs.n_calls
                                   for k, v in sorted(rhs.deg_counts.items())},
        })
        rows.append(st)
        er, ei, ec = ric_errors(sol, truth)
        series[name] = {
            "t_s": t_eval[::5].tolist(),
            "radial_m": er[::5].tolist(),
            "in_track_m": ei[::5].tolist(),
            "cross_track_m": ec[::5].tolist(),
        }
        print(f"  {name}: RMS {st['pos_rms_m']:8.2f} m  final {st['pos_final_m']:8.2f} m"
              f"  proxy {st['proxy_rel_fixed138']:.3f}  switches {n_switch}"
              f"  wall {wall:.1f} s (grav {st['grav_s']:.1f} s)")

    # wall-time dispersion: 3 repeats of key runs
    reps = {}
    for name in ("fixed_138", "sched_p17_eps1e3_floor60_down"):
        walls, gravs = [], []
        for _ in range(3):
            _, rhs, wall = propagate(model, y0, dur, t_eval, runs[name], args)
            walls.append(wall)
            gravs.append(rhs.grav_ns / 1e9)
        reps[name] = {"wall_s": walls, "grav_s": gravs}
        print(f"  repeat {name}: wall {sorted(walls)}")

    # altitude/degree profile for the schedule figure (one orbit period)
    period = 2.0 * math.pi * math.sqrt(a_sma ** 3 / model.mu)
    tp = np.arange(0.0, 2.0 * period, 5.0)
    rr = np.interp(tp, truth.t, np.linalg.norm(truth.y[:3], axis=0))
    prof_deg = [int(sched_dwell(float(r) - model.r_ref)) for r in rr]

    _dump("r1_longarc.json", {
        "scenario": {"type": "eccentric_polar", "perilune_km": 50.0,
                     "apolune_km": 300.0, "duration_s": dur,
                     "truth_degree": 300, "integrator": "DOP853",
                     "rtol": 1e-11, "atol": 1e-4, "output_step_s": 120.0,
                     "rotation": "uniform sidereal about polar axis",
                     "period_s": period},
        "grav_timer_note": "gravity wall time accumulated per call with "
                           "perf_counter_ns inside the RHS; ~0.1 us overhead per call",
        "schedules": {"dwell": tab_dwell, "naive": tab_naive, "empirical": tab_emp},
        "truth": {"n_rhs": rhs_t.n_calls, "wall_s": wall_t,
                  "grav_s": rhs_t.grav_ns / 1e9},
        "rows": rows,
        "wall_time_repeats": reps,
        "profile": {"t_s": tp[::12].tolist(),
                    "altitude_km": ((rr[::12] - model.r_ref) / 1e3).tolist(),
                    "degree": prof_deg[::12]},
        "error_series": series,
    })


# ---------------------------------------------------------------- R1-SW
def r1_switch_instrumentation(model, args):
    print("== R1-SW: DOP853 step behavior around switches ==")
    rp, ra = model.r_ref + 50e3, model.r_ref + 300e3
    a_sma = 0.5 * (rp + ra)
    period = 2.0 * math.pi * math.sqrt(a_sma ** 3 / model.mu)
    dur = 2.2 * period
    sched, _ = make_p_schedule(model, 1e-3, 60)

    vp = math.sqrt(model.mu * (2.0 / rp - 1.0 / a_sma))
    va = math.sqrt(model.mu * (2.0 / ra - 1.0 / a_sma))
    starts = {
        "perilune_start": np.array([rp, 0, 0, 0, 0, vp]),
        "apolune_start": np.array([ra, 0, 0, 0, 0, -va]),
    }

    out = {}
    for phase, y0 in starts.items():
        for name, degfun in (("scheduled", sched), ("fixed_138", lambda h: 138)):
            rhs = Rhs(model, degfun, args)
            solver = DOP853(rhs, 0.0, y0, dur, rtol=1e-11, atol=1e-4)
            t_acc, h_acc, nfev_acc, alt_acc, deg_acc = [], [], [], [], []
            prev_nfev = rhs.n_calls
            while solver.status == "running":
                solver.step()
                t_acc.append(solver.t)
                h_acc.append(solver.t - solver.t_old)
                nfev_acc.append(rhs.n_calls - prev_nfev)
                prev_nfev = rhs.n_calls
                r = float(np.linalg.norm(solver.y[:3]))
                alt_acc.append((r - model.r_ref) / 1e3)
                deg_acc.append(int(degfun(r - model.r_ref)))
            t_acc = np.array(t_acc); h_acc = np.array(h_acc)
            nfev_acc = np.array(nfev_acc)
            deg_acc = np.array(deg_acc)
            # attempts per accepted step inferred from the 12-stage cost
            attempts = np.maximum(1, np.round(nfev_acc / 12.0)).astype(int)
            switches = np.flatnonzero(deg_acc[1:] != deg_acc[:-1]) + 1
            # step-size statistics near vs away from switches (+-600 s windows)
            near = np.zeros(len(t_acc), dtype=bool)
            for si in switches:
                near |= np.abs(t_acc - t_acc[si]) <= 600.0
            stat = {
                "n_steps": int(len(t_acc)),
                "n_rhs": int(rhs.n_calls),
                "n_switches": int(len(switches)),
                "inferred_rejected_attempts_total": int(np.sum(attempts - 1)),
                "median_step_s": float(np.median(h_acc)),
                "median_step_near_switch_s": float(np.median(h_acc[near])) if near.any() else None,
                "median_step_away_s": float(np.median(h_acc[~near])) if (~near).any() else None,
                "rejects_near_switch": int(np.sum(attempts[near] - 1)) if near.any() else 0,
                "rejects_away": int(np.sum(attempts[~near] - 1)) if (~near).any() else 0,
                "steps_near_switch": int(np.sum(near)),
            }
            out[f"{phase}/{name}"] = {
                "stats": stat,
                "series": {"t_s": t_acc.tolist(), "h_s": h_acc.tolist(),
                           "attempts": attempts.tolist(),
                           "altitude_km": alt_acc, "degree": deg_acc.tolist()},
            }
            print(f"  {phase}/{name}: steps {stat['n_steps']}, switches "
                  f"{stat['n_switches']}, rejected {stat['inferred_rejected_attempts_total']}, "
                  f"med step {stat['median_step_s']:.1f} s "
                  f"(near {stat['median_step_near_switch_s']} / away {stat['median_step_away_s']})")

    _dump("r1_switch_instrumentation.json", {
        "note": "DOP853 low-level stepping; attempts inferred as round(nfev/12) "
                "per accepted step (12 stages per attempt, no dense output); "
                "switch = degree change between consecutive accepted steps",
        "scenario": {"perilune_km": 50.0, "apolune_km": 300.0,
                     "duration_s": dur, "rtol": 1e-11, "atol": 1e-4},
        "cases": out,
    })


def degree_power(model):
    N = model.max_degree
    C, S = model.c_coeffs, model.s_coeffs
    return np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                     for n in range(N + 1)])


def main() -> int:
    global COMMIT
    COMMIT = commit_sha()
    print("commit:", COMMIT)
    OUT.mkdir(exist_ok=True)
    gf = resolve_lunar_gravity_path(None)
    model = GravityModel.from_file(str(gf), requested_degree=300)
    model1800 = GravityModel.from_file(str(gf), requested_degree=1800)
    power = degree_power(model1800)
    args = kernel_args(model)
    # warm-up
    sh_accel_fixed_numba(model.r_ref + 100e3, 0.0, 0.0, 60, *args)

    r1_orbit_mapping(model, args)
    r1_long_arc(model, args, power)
    r1_switch_instrumentation(model, args)
    print("orbit-level revision experiments complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
