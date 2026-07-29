"""Geometric (symplectic) verification of the propagation contract (R11).

Reviewer objection
------------------
"DOP853 is not enough; the results sit under the noise floor, so
geometry-preserving symplectic runs are required."

What this script establishes
----------------------------
The production dynamics are an inertial integration of a uniformly rotating
body-fixed field.  In the co-rotating frame this is an *autonomous* Hamiltonian

    H(q, p) = |p|^2/2 - w . (q x p) - U(q)                              (1)

with ``q`` the body-fixed position, ``p`` the conjugate momentum (the inertial
velocity resolved on body axes), ``w`` the lunar rotation vector, and ``U`` the
geodesy potential (``a = +grad U``).  Its value is the Jacobi integral, and it
is exactly the quantity a symplectic method is supposed to protect.

Equation (1) splits into three *exactly solvable* pieces,

    A: |p|^2/2      drift    q += h p
    B: -U(q)        kick     p += h a_body(q, N)
    C: -w.(q x p)   rotation q, p  ->  Rz(-|w| h) q,  Rz(-|w| h) p

so a symmetric composition of A, B, C is an explicit, genuinely symplectic
integrator *for the production dynamics* -- not a simplified surrogate.  We
build order 2 (Strang), order 4 and order 6 (Yoshida triple jump) from it.

Four blocks:

  G1  Jacobi conservation, fixed degree.  The smooth-Hamiltonian case.  Report
      max |dH/H| over seven days at three step sizes and confirm it is bounded
      and non-secular -- i.e. the integrator really is geometry preserving.

  G2  Cross-integrator agreement.  Refine the symplectic step and show the
      geometric trajectory converges to the DOP853 vector-tolerance trajectory.
      The residual bounds how much of the reported policy differences could be
      DOP853 artifacts.  Also report DOP853's own Jacobi drift.

  G3  The decisive control.  Apply the SAME symplectic integrator to a
      degree-*switching* schedule.  Because U changes at a switch, H is no
      longer a single smooth Hamiltonian, the backward-error theory behind
      bounded energy drift does not apply, and the Jacobi error is expected to
      jump and to NOT improve under step refinement.  Symplecticity therefore
      cannot supply the protection the objection asks for, for exactly the
      policies the study is about.

  G4  Work-normalized accuracy.  Force evaluations needed to reach a given
      position accuracy, symplectic vs DOP853, so the choice of DOP853 is
      justified on cost as well as on correctness.

Usage
-----
    python rev11_geometric_verification.py run
    python rev11_geometric_verification.py smoke
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev3_common import DAY, OMEGA_MOON, load_model, kernel_args, warmup
from rev7_doe_screening import CANONICAL, initial_state, emp_table, alt_sched
from rev3_common import degree_power
from lunaris.physics.spherical_harmonics import (
    sh_accel_fixed_numba, sh_potential_accel_fixed,
)

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r11_geometric_verification.json"
SMOKE_OUTPUT = METRICS / "r11_geometric_verification_smoke.json"

OUTPUT_STEP = 120.0
TRUTH_RTOL = 1.0e-13
TRUTH_ATOL = np.array([1.0e-6] * 3 + [1.0e-9] * 3)

# Geometries: circular polar, eccentric polar, low-perilune LRO-like.
CASES = {
    "c1_circ100_polar": {"hp_km": 100.0, "ha_km": 100.0, "incl_deg": 90.0,
                         "argp_deg": 0.0, "raan_deg": 0.0, "degree": 120},
    "c2_50x300_polar": {"hp_km": 50.0, "ha_km": 300.0, "incl_deg": 90.0,
                        "argp_deg": 0.0, "raan_deg": 0.0, "degree": 140},
    "c6_lro_30x216": {"hp_km": 30.0, "ha_km": 216.0, "incl_deg": 90.0,
                      "argp_deg": 270.0, "raan_deg": 0.0, "degree": 140},
}
STEPS_S = (20.0, 10.0, 5.0, 2.5)
SWITCH_STEPS_S = (20.0, 10.0, 5.0, 2.5)

def _composition(order: int) -> np.ndarray:
    """Sub-step weights of a symmetric composition of the given order.

    Yoshida triple jump.  The coefficients depend on the order of the method
    being composed: lifting a symmetric method of order ``p`` to ``p+2`` uses
    ``z1 = 1/(2 - 2^(1/(p+1)))``.  Reusing the ``p=2`` value for every jump
    silently caps the composition at order 4.
    """
    if order < 2 or order % 2:
        raise ValueError("order must be even and >= 2")
    weights = np.array([1.0])
    p = 2
    while p < order:
        z1 = 1.0 / (2.0 - 2.0 ** (1.0 / (p + 1.0)))
        z0 = 1.0 - 2.0 * z1
        weights = np.concatenate([z1 * weights, z0 * weights, z1 * weights])
        p += 2
    return weights


# ------------------------------------------------------------------ integrator
def symplectic_propagate(model, args, y0_inertial, duration, t_grid,
                         degree_of, step, order):
    """Explicit symplectic A/B/C splitting of Eq. (1) in the rotating frame.

    ``y0_inertial`` is the usual inertial state; at t=0 the body and inertial
    frames coincide, so q0 = r0 and p0 = v0.  Output is returned in the
    inertial frame on ``t_grid`` so it is directly comparable with DOP853.
    """
    weights = _composition(order)
    r_ref = model.r_ref
    q = np.asarray(y0_inertial[:3], dtype=float).copy()
    p = np.asarray(y0_inertial[3:], dtype=float).copy()

    n_steps = int(math.ceil(duration / step))
    out_q = np.empty((3, len(t_grid)))
    out_p = np.empty((3, len(t_grid)))
    out_deg = np.empty(len(t_grid), dtype=int)
    filled = 0
    n_force = 0
    deg_changes = 0
    prev_deg = None
    sum_deg_sq = 0.0

    def kick(qv, pv, dt, t):
        nonlocal n_force, prev_deg, deg_changes, sum_deg_sq
        rad = math.sqrt(qv[0] * qv[0] + qv[1] * qv[1] + qv[2] * qv[2])
        n = degree_of(t, rad - r_ref)
        if prev_deg is not None and n != prev_deg:
            deg_changes += 1
        prev_deg = n
        sum_deg_sq += float(n) * float(n)
        n_force += 1
        ax, ay, az = sh_accel_fixed_numba(qv[0], qv[1], qv[2], n, *args)
        pv[0] += dt * ax
        pv[1] += dt * ay
        pv[2] += dt * az
        return n

    def rotate(qv, pv, dt):
        ang = -OMEGA_MOON * dt
        c, s = math.cos(ang), math.sin(ang)
        for vec in (qv, pv):
            x, y = vec[0], vec[1]
            vec[0] = c * x - s * y
            vec[1] = s * x + c * y

    t = 0.0
    current_deg = degree_of(0.0, float(np.linalg.norm(q)) - r_ref)
    wall0 = time.perf_counter_ns()
    while filled < len(t_grid) and t_grid[filled] <= 1e-9:
        out_q[:, filled] = q
        out_p[:, filled] = p
        out_deg[filled] = current_deg
        filled += 1
    for _ in range(n_steps):
        h = min(step, duration - t)
        if h <= 0.0:
            break
        for w in weights:
            hw = w * h
            rotate(q, p, 0.5 * hw)
            current_deg = kick(q, p, 0.5 * hw, t)
            q += hw * p
            current_deg = kick(q, p, 0.5 * hw, t + hw)
            rotate(q, p, 0.5 * hw)
        t += h
        # sample onto the output grid (grid is a multiple of the step)
        while filled < len(t_grid) and t_grid[filled] <= t + 1e-9:
            ang = OMEGA_MOON * t
            c, s = math.cos(ang), math.sin(ang)
            out_q[0, filled] = c * q[0] - s * q[1]
            out_q[1, filled] = s * q[0] + c * q[1]
            out_q[2, filled] = q[2]
            out_p[0, filled] = c * p[0] - s * p[1]
            out_p[1, filled] = s * p[0] + c * p[1]
            out_p[2, filled] = p[2]
            out_deg[filled] = current_deg
            filled += 1
    wall_ns = time.perf_counter_ns() - wall0
    telemetry = {"n_force_evaluations": n_force, "n_steps": n_steps,
                 "order": order, "step_s": step,
                 "switch_count_at_force_samples": deg_changes,
                 "mean_degree_sq": sum_deg_sq / max(n_force, 1),
                 "total_wall_ns": int(wall_ns)}
    return out_q[:, :filled], out_p[:, :filled], out_deg[:filled], telemetry


# ------------------------------------------------------------------- diagnostic
def jacobi_series(model, q_inertial, p_inertial, t_grid, degrees):
    """Jacobi integral of Eq. (1) sampled on the output grid.

    Positions arrive in the inertial frame; rotate them back to body axes,
    evaluate U in one batch per distinct degree, and assemble H.
    """
    ang = OMEGA_MOON * np.asarray(t_grid, dtype=float)
    c, s = np.cos(ang), np.sin(ang)
    qb = np.empty_like(q_inertial)
    qb[0] = c * q_inertial[0] + s * q_inertial[1]
    qb[1] = -s * q_inertial[0] + c * q_inertial[1]
    qb[2] = q_inertial[2]
    pb = np.empty_like(p_inertial)
    pb[0] = c * p_inertial[0] + s * p_inertial[1]
    pb[1] = -s * p_inertial[0] + c * p_inertial[1]
    pb[2] = p_inertial[2]

    potential = np.empty(qb.shape[1])
    for degree in np.unique(degrees):
        mask = degrees == degree
        V, _ = sh_potential_accel_fixed(
            np.ascontiguousarray(qb[:, mask].T), model.c_coeffs,
            model.s_coeffs, model.mu, model.r_ref, int(degree))
        potential[mask] = V
    kinetic = 0.5 * np.sum(pb * pb, axis=0)
    # w . (q x p) with w along +z reduces to |w| (qx py - qy px)
    angular = OMEGA_MOON * (qb[0] * pb[1] - qb[1] * pb[0])
    return kinetic - angular - potential


def jacobi_report(h_series):
    h0 = h_series[0]
    rel = (h_series - h0) / abs(h0)
    return {"H0_m2_s2": float(h0),
            "max_abs_relative_drift": float(np.max(np.abs(rel))),
            "final_relative_drift": float(rel[-1]),
            "std_relative": float(np.std(rel)),
            # a secular (drifting) error has |mean| comparable to |max|;
            # a bounded oscillation has |mean| much smaller
            "mean_relative": float(np.mean(rel)),
            "secular_fraction": float(abs(np.mean(rel)) /
                                      max(np.max(np.abs(rel)), 1e-300))}


def pos_rms(a, b):
    n = min(a.shape[1], b.shape[1])
    d = np.linalg.norm(a[:, :n] - b[:, :n], axis=0)
    return float(np.sqrt(np.mean(d * d))), float(np.max(d))


# ------------------------------------------------------------------------ main
def run(smoke: bool) -> int:
    duration = 0.5 * DAY if smoke else 7.0 * DAY
    steps = (20.0, 10.0) if smoke else STEPS_S
    switch_steps = (20.0, 10.0) if smoke else SWITCH_STEPS_S
    cases = dict(list(CASES.items())[:1]) if smoke else CASES
    t_grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)

    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    schedule_table = emp_table(model, degree_power(model))
    schedule = alt_sched(schedule_table)

    payload = {
        "schema": "r11_geometric_verification_v1",
        "created_utc": base.utc_now(),
        "smoke": smoke,
        "hamiltonian": ("H = |p|^2/2 - w.(q x p) - U(q); rotating-frame "
                        "Jacobi integral of the production dynamics"),
        "splitting": "A drift / B kick / C exact rotation, symmetric composition",
        "composition": {"order_4_weights": _composition(4).tolist(),
                        "order_6_weights": _composition(6).tolist()},
        "omega_moon_rad_s": OMEGA_MOON,
        "duration_s": duration, "output_step_s": OUTPUT_STEP,
        "reference_integrator": {"method": "DOP853", "rtol": TRUTH_RTOL,
                                 "atol_position_m": 1e-6,
                                 "atol_velocity_m_s": 1e-9,
                                 "max_step_s": 60.0},
        "cases": {},
        "source": base.provenance(),
    }

    for name, spec in cases.items():
        degree = int(spec["degree"])
        y0 = initial_state(model, {**spec, "name": name})
        print(f"[{name}] N={degree} duration={duration/DAY:.2f} d", flush=True)
        entry = {"geometry": {k: spec[k] for k in
                              ("hp_km", "ha_km", "incl_deg", "argp_deg")},
                 "fixed_degree": degree}

        # --- DOP853 reference at the tighter vector tolerance ---------------
        t_ref, y_ref, status, _, failure, tel_ref = \
            base.propagate_event_instrumented(
                model, y0, duration, t_grid, lambda t, h: degree, args,
                TRUTH_RTOL, TRUTH_ATOL, max_step=60.0)
        if status == "numerical_failure":
            entry["reference_failure"] = failure
            payload["cases"][name] = entry
            continue
        ref_q, ref_p = y_ref[:3], y_ref[3:]
        ref_deg = np.full(ref_q.shape[1], degree, dtype=int)
        ref_h = jacobi_series(model, ref_q, ref_p, t_ref, ref_deg)
        entry["dop853_reference"] = {
            "telemetry": tel_ref, "status": status,
            "jacobi": jacobi_report(ref_h)}
        print(f"   DOP853  nrhs={tel_ref['n_rhs']:7d} "
              f"|dH/H|max={entry['dop853_reference']['jacobi']['max_abs_relative_drift']:.3e}",
              flush=True)

        # --- G1/G2: fixed degree, symplectic step refinement ----------------
        entry["symplectic_fixed_degree"] = []
        for order in (4, 6):
            for step in steps:
                q, p, deg, tel = symplectic_propagate(
                    model, args, y0, duration, t_grid,
                    lambda t, h, n=degree: n, step, order)
                hs = jacobi_series(model, q, p, t_grid[:q.shape[1]], deg)
                rms, mx = pos_rms(q, ref_q)
                record = {"order": order, "step_s": step,
                          "telemetry": tel, "jacobi": jacobi_report(hs),
                          "position_rms_vs_dop853_m": rms,
                          "position_max_vs_dop853_m": mx}
                entry["symplectic_fixed_degree"].append(record)
                print(f"   sympl o{order} h={step:5.1f}s  "
                      f"|dH/H|max={record['jacobi']['max_abs_relative_drift']:.3e}  "
                      f"secfrac={record['jacobi']['secular_fraction']:.2f}  "
                      f"dPos_rms={rms:.4e} m  nforce={tel['n_force_evaluations']}",
                      flush=True)

        # --- G3: the same integrator on a switching schedule ----------------
        entry["symplectic_switching_schedule"] = []
        t_sw, y_sw, sw_status, _, sw_fail, tel_sw = \
            base.propagate_event_instrumented(
                model, y0, duration, t_grid, schedule, args,
                TRUTH_RTOL, TRUTH_ATOL, max_step=60.0)
        if sw_status != "numerical_failure":
            # The degree actually used at each output epoch, so the piecewise
            # Jacobi integral is formed with the policy's own potential.
            sw_deg = np.array([schedule(float(tt), float(np.linalg.norm(rr)) -
                                        model.r_ref)
                               for tt, rr in zip(t_sw, y_sw[:3].T)], dtype=int)
            sw_h = jacobi_series(model, y_sw[:3], y_sw[3:], t_sw, sw_deg)
            entry["dop853_switching"] = {
                "telemetry": tel_sw, "status": sw_status,
                "jacobi": jacobi_report(sw_h),
                "distinct_degrees": [int(x) for x in np.unique(sw_deg)]}
            print(f"   DOP853 switching |dH/H|max="
                  f"{entry['dop853_switching']['jacobi']['max_abs_relative_drift']:.3e}",
                  flush=True)
        for step in switch_steps:
            q, p, deg, tel = symplectic_propagate(
                model, args, y0, duration, t_grid, schedule, step, 4)
            hs = jacobi_series(model, q, p, t_grid[:q.shape[1]], deg)
            record = {"order": 4, "step_s": step, "telemetry": tel,
                      "jacobi": jacobi_report(hs),
                      "distinct_degrees": [int(x) for x in np.unique(deg)]}
            entry["symplectic_switching_schedule"].append(record)
            print(f"   switch o4 h={step:5.1f}s  "
                  f"|dH/H|max={record['jacobi']['max_abs_relative_drift']:.3e}  "
                  f"switches={tel['switch_count_at_force_samples']}", flush=True)

        payload["cases"][name] = entry

    # ---------------- G4 aggregate: refinement behavior -------------------
    def refinement(records, key="max_abs_relative_drift"):
        by_step = sorted(records, key=lambda r: -r["step_s"])
        return [{"step_s": r["step_s"], key: r["jacobi"][key]} for r in by_step]

    summary = {}
    for name, entry in payload["cases"].items():
        if "symplectic_fixed_degree" not in entry:
            continue
        fixed4 = [r for r in entry["symplectic_fixed_degree"] if r["order"] == 4]
        fixed6 = [r for r in entry["symplectic_fixed_degree"] if r["order"] == 6]
        switch = entry.get("symplectic_switching_schedule", [])
        finest6 = min(fixed6, key=lambda r: r["step_s"]) if fixed6 else None
        summary[name] = {
            "fixed_degree_order4_drift_vs_step": refinement(fixed4),
            "fixed_degree_order6_drift_vs_step": refinement(fixed6),
            "switching_order4_drift_vs_step": refinement(switch),
            "dop853_jacobi_max_abs_relative_drift":
                entry["dop853_reference"]["jacobi"]["max_abs_relative_drift"],
            "finest_symplectic_vs_dop853_position_rms_m":
                finest6["position_rms_vs_dop853_m"] if finest6 else None,
            "force_evaluations_finest_symplectic":
                finest6["telemetry"]["n_force_evaluations"] if finest6 else None,
            "force_evaluations_dop853":
                entry["dop853_reference"]["telemetry"]["n_rhs"],
        }
        if switch:
            drifts = [r["jacobi"]["max_abs_relative_drift"] for r in
                      sorted(switch, key=lambda r: -r["step_s"])]
            summary[name]["switching_drift_reduction_coarse_to_fine"] = (
                drifts[0] / drifts[-1] if drifts[-1] > 0 else None)
        if fixed4:
            d4 = [r["jacobi"]["max_abs_relative_drift"] for r in
                  sorted(fixed4, key=lambda r: -r["step_s"])]
            summary[name]["fixed_drift_reduction_coarse_to_fine"] = (
                d4[0] / d4[-1] if d4[-1] > 0 else None)
    payload["summary"] = summary

    target = SMOKE_OUTPUT if smoke else OUTPUT
    base.atomic_json(target, payload)
    print(f"[written] {target.name}", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("run", "smoke"))
    args = parser.parse_args()
    return run(args.command == "smoke")


if __name__ == "__main__":
    raise SystemExit(main())
