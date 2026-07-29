"""Robustness Test 4: variational (STM) mechanism of the in-track error.

The scheduling penalty is dominated by a coherent in-track drift.  This test
shows that the drift is the linear variational response of the reference
trajectory to the *measured* truncation force defect, i.e. that

    delta_x(T) ~= integral_0^T Phi(T, tau) B Delta_a(tau) d tau,

where Phi is the reference state-transition matrix, B maps acceleration into
the velocity block, and Delta_a(tau) = a_policy(x_ref) - a_ref(x_ref) is the
policy-minus-reference acceleration evaluated *along the reference*.

Rather than form Phi explicitly, we integrate the forced variational equation

    d/dt delta_r = delta_v
    d/dt delta_v = G(tau) delta_r + Delta_a(tau)

alongside the reference, with G(tau) = d a_ref / d r the reference gravity
gradient (central finite differences at N=300, 1 m step).  All policies share
one reference and one gradient, so their linear predictions are obtained in a
single augmented integration.  The linear prediction delta_x_linear is then
compared with the true nonlinear difference delta_x_actual = x_policy - x_ref.

Same uniformly rotating N=300 field as the truncation experiments and Test 2
(so the polar downward-schedule case reproduces the 568.8 m headline entry).

Full run:  ``.venv\\Scripts\\python.exe robustness_test4_variational_mechanism.py``
Smoke run: add ``--smoke`` (four hours; writes a separate artifact).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from rev3_common import (DAY, OMEGA_MOON, OUT, alt_sched, degree_power,
                         eccentric_state, kernel_args, load_model,
                         make_emp_table, make_p_table, propagate, warmup,
                         commit_sha, working_tree_clean)
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

FD_STEP_M = 1.0
OUTPUT_STEP_S = 120.0
RTOL = 1.0e-10
ATOL = 1.0e-8
MAX_STEP_S = 60.0


def accel_inertial(r_vec, t, n, args, r_ref):
    """Inertial acceleration of the uniformly rotating body-fixed N-field."""
    x, y, z = r_vec
    th = OMEGA_MOON * t
    c, s = math.cos(th), math.sin(th)
    xb = c * x + s * y
    yb = -s * x + c * y
    axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *args)
    return np.array([c * axb - s * ayb, s * axb + c * ayb, azb])


def gravity_gradient(r_vec, t, n, args, r_ref, h=FD_STEP_M):
    """3x3 d a_inertial / d r by central finite differences at degree n."""
    g = np.empty((3, 3))
    for j in range(3):
        rp = r_vec.copy(); rp[j] += h
        rm = r_vec.copy(); rm[j] -= h
        g[:, j] = (accel_inertial(rp, t, n, args, r_ref)
                   - accel_inertial(rm, t, n, args, r_ref)) / (2.0 * h)
    return g


def ric_axes(r, v):
    radial = r / np.linalg.norm(r)
    cross = np.cross(r, v)
    cross /= np.linalg.norm(cross)
    intrack = np.cross(cross, radial)
    return radial, intrack, cross


def project_ric(ref_states, delta):
    """Project a (3, N) position difference onto the reference RIC frame."""
    out = np.empty((delta.shape[1], 3))
    for k in range(delta.shape[1]):
        radial, intrack, cross = ric_axes(ref_states[:3, k], ref_states[3:, k])
        d = delta[:, k]
        out[k] = (radial @ d, intrack @ d, cross @ d)
    return out  # (N, 3): radial, in_track, cross_track


def compare_series(ref_states, actual_pos, linear_pos):
    """Compare actual vs linear position differences (both (3, N))."""
    a_ric = project_ric(ref_states, actual_pos)
    l_ric = project_ric(ref_states, linear_pos)
    res_ric = a_ric - l_ric
    labels = ("radial", "in_track", "cross_track")

    def rms(col):
        return float(np.sqrt(np.mean(col * col)))

    ai, li, ri = a_ric[:, 1], l_ric[:, 1], res_ric[:, 1]
    denom = math.sqrt(float(np.sum(ai * ai)) * float(np.sum(li * li)))
    corr_it = float(np.sum(ai * li) / denom) if denom > 0 else float("nan")
    actual_norm = np.linalg.norm(actual_pos, axis=0)
    linear_norm = np.linalg.norm(linear_pos, axis=0)
    resid_norm = np.linalg.norm(actual_pos - linear_pos, axis=0)
    return {
        "actual_ric_rms_m": {l: rms(a_ric[:, i]) for i, l in enumerate(labels)},
        "linear_ric_rms_m": {l: rms(l_ric[:, i]) for i, l in enumerate(labels)},
        "residual_ric_rms_m": {l: rms(res_ric[:, i]) for i, l in enumerate(labels)},
        "actual_ric_final_m": {l: float(a_ric[-1, i]) for i, l in enumerate(labels)},
        "linear_ric_final_m": {l: float(l_ric[-1, i]) for i, l in enumerate(labels)},
        "actual_pos_rms_m": rms(actual_norm),
        "linear_pos_rms_m": rms(linear_norm),
        "residual_pos_rms_m": rms(resid_norm),
        "in_track_explained_fraction": float(1.0 - rms(ri) / rms(ai)) if rms(ai) > 0 else float("nan"),
        "in_track_correlation": corr_it,
        "in_track_final_actual_m": float(a_ric[-1, 1]),
        "in_track_final_linear_m": float(l_ric[-1, 1]),
        "pos_rms_explained_fraction": float(1.0 - rms(resid_norm) / rms(actual_norm)) if rms(actual_norm) > 0 else float("nan"),
    }


def run_geometry(name, y0, model, args, policies, degree_funcs, duration, grid):
    """One augmented integration -> linear delta_x for every policy."""
    r_ref = model.r_ref
    policy_names = list(policies)
    npol = len(policy_names)

    def augmented(t, Y):
        r = Y[0:3]
        v = Y[3:6]
        a_ref = accel_inertial(r, t, 300, args, r_ref)
        g = gravity_gradient(r, t, 300, args, r_ref)
        alt = float(np.linalg.norm(r)) - r_ref
        dY = np.empty_like(Y)
        dY[0:3] = v
        dY[3:6] = a_ref
        for p, pname in enumerate(policy_names):
            n_p = int(degree_funcs[pname](t, alt))
            a_p = accel_inertial(r, t, n_p, args, r_ref) if n_p != 300 else a_ref
            delta_a = a_p - a_ref
            base = 6 + 6 * p
            dr = Y[base:base + 3]
            dv = Y[base + 3:base + 6]
            dY[base:base + 3] = dv
            dY[base + 3:base + 6] = g @ dr + delta_a
        return dY

    Y0 = np.zeros(6 + 6 * npol)
    Y0[0:6] = y0
    print(f"  [{name}] augmented variational integration ({npol} policies)",
          flush=True)
    sol = solve_ivp(augmented, (0.0, duration), Y0, method="DOP853",
                    t_eval=grid, rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S)
    if not sol.success:
        raise RuntimeError(sol.message)
    ref_states = sol.y[0:6]              # (6, N)
    linear_pos = {pname: sol.y[6 + 6 * p:6 + 6 * p + 3]
                  for p, pname in enumerate(policy_names)}

    # Actual nonlinear differences: propagate each policy, subtract reference.
    results = {}
    for pname in policy_names:
        print(f"  [{name}] nonlinear propagation: {pname}", flush=True)
        sol_p, _rhs, _wall = propagate(model, y0, duration, grid,
                                       degree_funcs[pname], args,
                                       rtol=1.0e-11, atol=1.0e-4,
                                       max_step=MAX_STEP_S)
        actual_pos = sol_p.y[0:3] - ref_states[0:3]
        results[pname] = compare_series(ref_states, actual_pos,
                                        linear_pos[pname])
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="run four hours and write a _smoke artifact")
    cli = parser.parse_args()
    duration = 4.0 * 3600.0 if cli.smoke else 7.0 * DAY
    grid = np.arange(0.0, duration + 0.1, OUTPUT_STEP_S)

    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)

    down = make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="down")
    up = make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="up")
    empirical = make_emp_table(model, degree_power(model), 1.0e-3,
                               60, cap=138, q=10)

    geometries = {
        "50x300_polar": {
            "y0": eccentric_state(model, 50.0, 300.0, 90.0),
            "fixed_critical_degree": 138,
        },
        "lro_like_30x216": {
            "y0": eccentric_state(model, 30.0, 216.0, 90.0),
            "fixed_critical_degree": 194,
        },
    }

    payload = {
        "schema": "robustness_test4_variational_mechanism_v1",
        "formal_run": not cli.smoke,
        "repo_commit_sha": commit_sha(),
        "repo_working_tree_clean": working_tree_clean(),
        "method": {
            "field": "uniformly rotating N=300 body-fixed field (rev3 harness)",
            "gravity_gradient": f"central finite differences, {FD_STEP_M} m step, N=300",
            "forcing": "Delta_a = a_policy(x_ref) - a_ref(x_ref) along reference",
            "integrator": "SciPy DOP853",
            "rtol": RTOL, "atol": ATOL, "max_step_s": MAX_STEP_S,
            "output_step_s": OUTPUT_STEP_S,
        },
        "scenario": {"duration_s": duration},
        "geometries": {},
    }

    for name, geo in geometries.items():
        print(f"== Test 4 {name} ==", flush=True)
        nc = int(geo["fixed_critical_degree"])
        degree_funcs = {
            "fixed_critical": (lambda _t, _h, n=nc: n),
            "schedule_down": alt_sched(down),
            "schedule_up": alt_sched(up),
            "schedule_emp": alt_sched(empirical),
        }
        policies = list(degree_funcs)
        results = run_geometry(name, np.asarray(geo["y0"], float), model, args,
                               policies, degree_funcs, duration, grid)
        payload["geometries"][name] = {
            "fixed_critical_degree": nc,
            "policies": results,
        }

    suffix = "_smoke" if cli.smoke else ""
    out = OUT / f"robustness_test4_variational_mechanism{suffix}.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {out.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
