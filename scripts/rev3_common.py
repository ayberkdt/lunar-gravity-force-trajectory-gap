"""Shared harness for the rev3 (publication-readiness) experiment set.

Everything is pinned to the Lunaris repository at D:\\Masaustu\\LUNAR_SIMULATION;
the commit SHA is recorded in every output payload. The RHS, error metrics,
and schedule builders reproduce the r1/r2 harness (rev_orbit.py/rev_orbit2.py)
so that new results are directly comparable with the archived ones.
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
    r = subprocess.run(["git", "-c", f"safe.directory={REPO.as_posix()}",
                        "rev-parse", "HEAD"], cwd=REPO,
                       capture_output=True, text=True, check=True)
    return r.stdout.strip()


def working_tree_clean() -> bool:
    r = subprocess.run(["git", "-c", f"safe.directory={REPO.as_posix()}",
                        "status", "--porcelain"], cwd=REPO,
                       capture_output=True, text=True, check=True)
    return r.stdout.strip() == ""


def dump(name: str, payload: dict) -> None:
    payload.setdefault("repo_commit_sha", commit_sha())
    payload.setdefault("repo_working_tree_clean", working_tree_clean())
    (OUT / name).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {name}")


def load_model(requested_degree: int = 300):
    gf = resolve_lunar_gravity_path(None)
    return GravityModel.from_file(str(gf), requested_degree=requested_degree)


def kernel_args(model):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


def warmup(model, args) -> None:
    sh_accel_fixed_numba(model.r_ref + 100e3, 0.0, 0.0, 60, *args)


def degree_power(model) -> np.ndarray:
    N = model.max_degree
    C, S = model.c_coeffs, model.s_coeffs
    return np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                     for n in range(N + 1)])


class Rhs:
    """Inertial RHS with uniformly rotating body-fixed field; instrumented."""

    def __init__(self, model, degree_of, args):
        self.model = model
        self.degree_of = degree_of
        self.args = args
        self.n_calls = 0
        self.sum_deg_sq = 0.0
        self.grav_ns = 0
        self.deg_counts: dict[int, int] = {}
        self.prev_deg: int | None = None
        self.n_deg_changes = 0

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
        if self.prev_deg is not None and n != self.prev_deg:
            self.n_deg_changes += 1
        self.prev_deg = n
        t0 = time.perf_counter_ns()
        axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *self.args)
        self.grav_ns += time.perf_counter_ns() - t0
        return (vx, vy, vz, c * axb - s * ayb, s * axb + c * ayb, azb)


def propagate(model, y0, dur, t_eval, degree_of, args, rtol, atol,
              max_step=np.inf):
    rhs = Rhs(model, degree_of, args)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853", rtol=rtol, atol=atol,
                    max_step=max_step, t_eval=t_eval)
    wall = time.perf_counter() - t0
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs, wall


def ric_errors(sol_y, truth_y):
    d = (sol_y[:3] - truth_y[:3]).T
    r = truth_y[:3].T
    v = truth_y[3:].T
    R = r / np.linalg.norm(r, axis=1, keepdims=True)
    C = np.cross(r, v)
    C = C / np.linalg.norm(C, axis=1, keepdims=True)
    I = np.cross(C, R)
    return (np.sum(d * R, axis=1), np.sum(d * I, axis=1), np.sum(d * C, axis=1))


def err_stats(sol_y, truth_y) -> dict:
    dp = np.linalg.norm(sol_y[:3] - truth_y[:3], axis=0)
    dv = np.linalg.norm(sol_y[3:] - truth_y[3:], axis=0)
    er, ei, ec = ric_errors(sol_y, truth_y)
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
        "ric_final_m": {"radial": float(er[-1]), "in_track": float(ei[-1]),
                        "cross_track": float(ec[-1])},
    }


# ------------------------------------------------------------------ orbits
def eccentric_state(model, hp_km: float, ha_km: float, incl_deg: float = 90.0,
                    at_apolune: bool = False) -> np.ndarray:
    """Initial state at perilune (or apolune) on the +x axis; the orbit
    normal makes the requested inclination with the +z rotation axis
    (cos i = v_y-direction cosine of the velocity at the node)."""
    rp = model.r_ref + hp_km * 1e3
    ra = model.r_ref + ha_km * 1e3
    a_sma = 0.5 * (rp + ra)
    ci = math.cos(math.radians(incl_deg))
    si = math.sin(math.radians(incl_deg))
    if at_apolune:
        va = math.sqrt(model.mu * (2.0 / ra - 1.0 / a_sma))
        return np.array([ra, 0.0, 0.0, 0.0, -ci * va, -si * va])
    vp = math.sqrt(model.mu * (2.0 / rp - 1.0 / a_sma))
    return np.array([rp, 0.0, 0.0, 0.0, ci * vp, si * vp])


def orbit_period(model, hp_km: float, ha_km: float) -> float:
    a_sma = model.r_ref + 0.5 * (hp_km + ha_km) * 1e3
    return 2.0 * math.pi * math.sqrt(a_sma ** 3 / model.mu)


# --------------------------------------------------------------- schedules
def make_p_table(model, eps, floor, cap=138, q=10, policy="down"):
    table = {}
    for hk in np.arange(40.0, 521.0, 10.0):
        n = recommended_sh_degree(hk, model.r_ref,
                                  kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                  kaula_tail_fraction=eps)
        n = max(floor, min(cap, n))
        if policy == "down":
            nq = (n // q) * q
        elif policy == "up":
            nq = ((n + q - 1) // q) * q
        else:
            nq = int(round(n / q)) * q
        table[hk] = max(floor, min(cap, nq))
    return table


def alt_sched(table):
    hmax = max(table)
    def f(t, h_m):
        hb = min(hmax, max(40.0, 10.0 * math.floor(h_m / 1e3 / 10.0)))
        return table[hb]
    return f


def make_emp_table(model, power, eps, floor, cap=138, q=10):
    n_arr = np.arange(len(power), dtype=np.float64)
    table = {}
    for hk in np.arange(40.0, 521.0, 10.0):
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
    return table


def mindwell_from_profile(truth, model, table, dur, min_dwell_s=600.0):
    """Time-based minimum-dwell schedule precomputed from the truth
    altitude profile (open-loop in time)."""
    import bisect
    t10 = np.arange(0.0, dur + 10.0, 10.0)
    rr = np.interp(t10, truth.t, np.linalg.norm(truth.y[:3], axis=0))
    f = alt_sched(table)
    raw = np.array([f(0.0, float(r) - model.r_ref) for r in rr])
    seg_t, seg_n = [0.0], [int(raw[0])]
    last = 0.0
    for i in range(1, len(t10)):
        if raw[i] != seg_n[-1] and (t10[i] - last) >= min_dwell_s:
            seg_t.append(float(t10[i]))
            seg_n.append(int(raw[i]))
            last = float(t10[i])

    def g(t, h_m):
        i = bisect.bisect_right(seg_t, t) - 1
        return seg_n[max(i, 0)]

    return g, seg_t, seg_n


# ------------------------------------------- direct DOP853 instrumentation
class InstrumentedDOP853(DOP853):
    """DOP853 with a direct attempted/rejected trial-step counter.

    Every iteration of the acceptance loop in RungeKutta._step_impl calls
    _estimate_error_norm exactly once, so wrapping it counts attempted
    trial steps directly; rejected trials are attempts that returned an
    error norm >= 1. This does not rely on the 12-stage RHS-count
    inference used in the r1 instrumentation.
    """

    def __init__(self, *a, **k):
        self.attempt_t: list[float] = []
        self.attempt_h: list[float] = []
        self.attempt_err: list[float] = []
        super().__init__(*a, **k)

    def _estimate_error_norm(self, K, h, scale):
        en = super()._estimate_error_norm(K, h, scale)
        self.attempt_t.append(float(self.t))
        self.attempt_h.append(float(h))
        self.attempt_err.append(float(en))
        return en

    @property
    def n_attempts(self) -> int:
        return len(self.attempt_err)

    @property
    def n_rejected(self) -> int:
        return int(sum(1 for e in self.attempt_err if not (e < 1.0)))


def propagate_instr(model, y0, dur, t_grid, degree_of, args, rtol, atol,
                    keep_steps: bool = False, max_step=np.inf, rhs_obj=None):
    """Low-level instrumented propagation: direct accepted/rejected trial
    counts plus states interpolated onto t_grid via DOP853 dense output.

    Returns (Y[6, len(t_grid)], rhs, info dict). A pre-built RHS object
    (e.g. a MOON_PA SPICE RHS) can be supplied via rhs_obj."""
    rhs = rhs_obj if rhs_obj is not None else Rhs(model, degree_of, args)
    t0w = time.perf_counter()
    solver = InstrumentedDOP853(rhs, 0.0, np.asarray(y0, float), dur,
                                rtol=rtol, atol=atol, max_step=max_step)
    Y = np.empty((6, len(t_grid)))
    filled = 0
    if t_grid[0] == 0.0:
        Y[:, 0] = y0
        filled = 1
    step_h: list[float] = []
    step_t: list[float] = []
    step_rej: list[int] = []
    prev_attempts = solver.n_attempts
    while solver.status == "running":
        solver.step()
        if solver.status == "failed":
            raise RuntimeError("DOP853 failed")
        step_t.append(float(solver.t))
        step_h.append(float(solver.t - solver.t_old))
        step_rej.append(solver.n_attempts - prev_attempts - 1)
        prev_attempts = solver.n_attempts
        hi = solver.t
        n_new = 0
        while filled + n_new < len(t_grid) and t_grid[filled + n_new] <= hi + 1e-9:
            n_new += 1
        if n_new:
            dense = solver.dense_output()
            for k in range(n_new):
                Y[:, filled + k] = dense(t_grid[filled + k])
            filled += n_new
    if filled < len(t_grid):
        dense = solver.dense_output()
        for k in range(filled, len(t_grid)):
            Y[:, k] = dense(min(t_grid[k], solver.t))
    wall = time.perf_counter() - t0w
    info = {
        "n_accepted_steps": len(step_t),
        "n_attempted_steps": solver.n_attempts,
        "n_rejected_trials": solver.n_rejected,
        "n_rhs": rhs.n_calls,
        "wall_s": wall,
        "grav_s": rhs.grav_ns / 1e9,
    }
    if keep_steps:
        info["step_t_s"] = step_t
        info["step_h_s"] = step_h
        info["step_rejected_trials"] = step_rej
    return Y, rhs, info
