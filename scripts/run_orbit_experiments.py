"""E6/E7: orbit-level mapping of truncation error and adaptive switching.

E6: circular polar orbit at 100 km altitude, uniformly rotating Moon,
propagated 24 h with DOP853 at several truncation degrees; position error
against an N=300 truth run.

E7: eccentric polar orbit (50 x 300 km), 24 h: fixed-degree baselines vs
altitude-scheduled discrete degree switching (p=1.7 criterion, downward
quantization q=10); error and cost (RHS-weighted mean N^2) comparison.

All gravity evaluations use the production Lunaris serial Kahan kernel.
"""

from __future__ import annotations

import json
import math
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
SEED = 20260719
OMEGA_MOON = 2.0 * math.pi / (27.321661 * 86400.0)  # rad/s, sidereal
DAY_S = 86400.0


def kernel_args(model: GravityModel):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


class Rhs:
    """Inertial-frame RHS with a uniformly rotating body-fixed gravity field."""

    def __init__(self, model: GravityModel, degree_of, args):
        self.model = model
        self.degree_of = degree_of  # callable: altitude_m -> int
        self.args = args
        self.n_calls = 0
        self.sum_deg_sq = 0.0

    def __call__(self, t, y):
        self.n_calls += 1
        x, yy, z, vx, vy, vz = y
        th = OMEGA_MOON * t
        c, s = math.cos(th), math.sin(th)
        # inertial -> body (R_z(-theta))
        xb = c * x + s * yy
        yb = -s * x + c * yy
        r = math.sqrt(x * x + yy * yy + z * z)
        n = self.degree_of(r - self.model.r_ref)
        self.sum_deg_sq += float(n) * float(n)
        axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *self.args)
        # body -> inertial (R_z(theta))
        ax = c * axb - s * ayb
        ay = s * axb + c * ayb
        return (vx, vy, vz, ax, ay, azb)


def propagate(model, y0, t_span, t_eval, degree_of, args):
    rhs = Rhs(model, degree_of, args)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, t_span, y0, method="DOP853",
                    rtol=1e-11, atol=1e-4, t_eval=t_eval, dense_output=False)
    wall = time.perf_counter() - t0
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs, wall


def pos_err_stats(sol, sol_truth):
    d = sol.y[:3] - sol_truth.y[:3]
    err = np.linalg.norm(d, axis=0)
    return {"rms_m": float(np.sqrt(np.mean(err ** 2))),
            "final_m": float(err[-1]),
            "max_m": float(np.max(err))}


def main() -> int:
    OUT.mkdir(exist_ok=True)
    gravity_file = resolve_lunar_gravity_path(None)
    model = GravityModel.from_file(str(gravity_file), requested_degree=300)
    args = kernel_args(model)
    R, mu = model.r_ref, model.mu

    t_eval = np.arange(0.0, DAY_S + 1.0, 120.0)

    # ---------------- E6: circular polar orbit, 100 km ----------------
    r0 = R + 100e3
    v0 = math.sqrt(mu / r0)
    y0 = np.array([r0, 0.0, 0.0, 0.0, 0.0, v0])

    print("E6: truth run (N=300)...")
    truth, rhs_t, wall_t = propagate(model, y0, (0.0, DAY_S), t_eval,
                                     lambda h: 300, args)
    print(f"  truth: {rhs_t.n_calls} evals, {wall_t:.1f}s")

    degrees = [20, 37, 56, 60, 69, 84, 100, 124, 150]
    rows = []
    for n in degrees:
        sol, rhs, wall = propagate(model, y0, (0.0, DAY_S), t_eval,
                                   lambda h, n=n: n, args)
        st = pos_err_stats(sol, truth)
        st.update({"degree": n, "n_rhs": rhs.n_calls, "wall_s": wall})
        rows.append(st)
        print(f"  N={n:4d}: RMS {st['rms_m']:10.2f} m  final {st['final_m']:10.2f} m")

    e6 = {
        "seed": SEED,
        "scenario": {"type": "circular_polar", "altitude_km": 100.0,
                     "duration_s": DAY_S, "truth_degree": 300,
                     "integrator": "DOP853 rtol=1e-11 atol=1e-4",
                     "rotation_rad_s": OMEGA_MOON},
        "criteria_at_100km": {"empirical": 69, "kaula_p1_7": 70,
                              "kaula_p2_0": 56, "attenuation_only": 124},
        "truth_rhs_calls": rhs_t.n_calls,
        "rows": rows,
    }
    (OUT / "e6_orbit_mapping.json").write_text(json.dumps(e6, indent=2),
                                               encoding="utf-8")
    print("[written] e6_orbit_mapping.json")

    # ---------------- E7: eccentric polar orbit, 50 x 300 km ----------------
    rp, ra = R + 50e3, R + 300e3
    a_sma = 0.5 * (rp + ra)
    vp = math.sqrt(mu * (2.0 / rp - 1.0 / a_sma))
    y0e = np.array([rp, 0.0, 0.0, 0.0, 0.0, vp])

    print("E7: truth run (N=300)...")
    truth_e, rhs_te, wall_te = propagate(model, y0e, (0.0, DAY_S), t_eval,
                                         lambda h: 300, args)
    print(f"  truth: {rhs_te.n_calls} evals, {wall_te:.1f}s")

    # altitude-scheduled degree table: p=1.7 criterion, quantized down to q=10
    q = 10
    bins_km = np.arange(40.0, 321.0, 10.0)
    table = {}
    for hk in bins_km:
        n_rec = recommended_sh_degree(hk, R,
                                      kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                      kaula_tail_fraction=1e-2)
        n_rec = max(37, min(138, n_rec))
        table[hk] = max(37, (n_rec // q) * q)

    def sched(h_m: float) -> int:
        hk = min(320.0, max(40.0, 10.0 * math.floor(h_m / 1000.0 / 10.0)))
        return table[hk]

    runs = {
        "fixed_138": (lambda h: 138),
        "fixed_106": (lambda h: 106),
        "adaptive_p17_q10": sched,
    }
    rows_e = []
    for name, degfun in runs.items():
        sol, rhs, wall = propagate(model, y0e, (0.0, DAY_S), t_eval, degfun, args)
        st = pos_err_stats(sol, truth_e)
        mean_deg_sq = rhs.sum_deg_sq / rhs.n_calls
        st.update({"run": name, "n_rhs": rhs.n_calls, "wall_s": wall,
                   "mean_degree_sq": mean_deg_sq,
                   "cost_rel_fixed138": mean_deg_sq / (138.0 ** 2)})
        rows_e.append(st)
        print(f"  {name}: RMS {st['rms_m']:9.2f} m  final {st['final_m']:9.2f} m  "
              f"cost/f138 {st['cost_rel_fixed138']:.3f}  wall {wall:.1f}s")

    e7 = {
        "seed": SEED,
        "scenario": {"type": "eccentric_polar", "perilune_km": 50.0,
                     "apolune_km": 300.0, "duration_s": DAY_S,
                     "truth_degree": 300, "quantization_step": q,
                     "degree_bounds": [37, 138],
                     "integrator": "DOP853 rtol=1e-11 atol=1e-4"},
        "degree_table_km": {f"{k:.0f}": v for k, v in table.items()},
        "rows": rows_e,
    }
    (OUT / "e7_adaptive_orbit.json").write_text(json.dumps(e7, indent=2),
                                                encoding="utf-8")
    print("[written] e7_adaptive_orbit.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
