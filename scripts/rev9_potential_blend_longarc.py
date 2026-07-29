"""Long-arc accuracy test of the corrected potential-level blend.

The test uses the same rotating lunar field, tolerances, output cadence, and
two representative geometries as the Stage-3 long-arc audit.  It compares a
fixed N=120 field, a discrete N=30/120 switch, the naive acceleration blend,
and the conservative potential blend against the geometry's fixed high-degree
truth trajectory.  The 28-day run includes 1, 7, 14, 21, and 28-day
checkpoints, so it also supplies the requested 24-hour and seven-day evidence.

Full run:
    .venv\Scripts\python.exe rev9_potential_blend_longarc.py

Smoke run (two geometries, six hours):
    .venv\Scripts\python.exe rev9_potential_blend_longarc.py --smoke
"""

from __future__ import annotations

import argparse
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from rev3_common import (
    DAY,
    OMEGA_MOON,
    dump,
    err_stats,
    kernel_args,
    load_model,
    propagate,
    warmup,
)
from rev7_doe_screening import CANONICAL, initial_state

from lunaris.physics.spherical_harmonics import (
    sh_accel_fixed_numba,
    sh_potential_accel_fixed,
)


N_LO, N_HI = 30, 120
ALT_NEAR, ALT_FAR = 50e3, 200e3
RTOL, ATOL, MAX_STEP = 1.0e-11, 1.0e-6, 120.0
OUTPUT_STEP = 300.0
FULL_DAYS = 28.0
ORBIT_NAMES = ("c2_50x300_polar", "c6_lro_30x216")
CHECKPOINT_DAYS = (1, 7, 14, 21, 28)


def orbit_dict(name: str) -> dict:
    for n, hp, ha, incl, argp, raan in CANONICAL:
        if n == name:
            return {
                "name": n,
                "family": "canonical",
                "hp_km": hp,
                "ha_km": ha,
                "incl_deg": incl,
                "argp_deg": argp,
                "raan_deg": raan,
            }
    raise KeyError(name)


def weight_and_derivative(radius: float, reference_radius: float) -> tuple[float, float]:
    """Return high-degree weight and its radial derivative."""
    altitude = radius - reference_radius
    if altitude <= ALT_NEAR:
        return 1.0, 0.0
    if altitude >= ALT_FAR:
        return 0.0, 0.0
    t = (ALT_FAR - altitude) / (ALT_FAR - ALT_NEAR)
    weight = t * t * (3.0 - 2.0 * t)
    derivative = 6.0 * t * (1.0 - t) * (-1.0 / (ALT_FAR - ALT_NEAR))
    return weight, derivative


class BlendRhs:
    """Rotating-field RHS for the four N=30/120 comparison policies."""

    def __init__(self, model, args, policy: str):
        self.model = model
        self.args = args
        self.policy = policy
        self.n_calls = 0
        self.grav_ns = 0

    def _potential_acceleration(self, xyz: np.ndarray, degree: int):
        potential, acceleration = sh_potential_accel_fixed(
            xyz.reshape(1, 3),
            self.model.c_coeffs,
            self.model.s_coeffs,
            self.model.mu,
            self.model.r_ref,
            degree,
            -1,
        )
        return float(potential[0]), acceleration[0]

    def __call__(self, t, y):
        self.n_calls += 1
        x, y_inertial, z, vx, vy, vz = y
        angle = OMEGA_MOON * t
        cosine, sine = math.cos(angle), math.sin(angle)
        xb = cosine * x + sine * y_inertial
        yb = -sine * x + cosine * y_inertial
        body_position = np.array([xb, yb, z], dtype=float)
        radius = float(np.linalg.norm(body_position))
        weight, derivative = weight_and_derivative(radius, self.model.r_ref)

        start = time.perf_counter_ns()
        if self.policy == "fixed_N120":
            acceleration = np.asarray(
                sh_accel_fixed_numba(xb, yb, z, N_HI, *self.args), dtype=float
            )
        elif self.policy == "switch_N30_N120":
            degree = N_HI if weight >= 0.5 else N_LO
            acceleration = np.asarray(
                sh_accel_fixed_numba(xb, yb, z, degree, *self.args), dtype=float
            )
        elif self.policy == "blend_acceleration":
            a_lo = np.asarray(
                sh_accel_fixed_numba(xb, yb, z, N_LO, *self.args), dtype=float
            )
            a_hi = np.asarray(
                sh_accel_fixed_numba(xb, yb, z, N_HI, *self.args), dtype=float
            )
            acceleration = (1.0 - weight) * a_lo + weight * a_hi
        elif self.policy == "blend_potential_corrected":
            u_lo, a_lo = self._potential_acceleration(body_position, N_LO)
            u_hi, a_hi = self._potential_acceleration(body_position, N_HI)
            acceleration = (1.0 - weight) * a_lo + weight * a_hi
            if derivative != 0.0:
                acceleration = acceleration + (
                    (u_hi - u_lo) * derivative * body_position / radius
                )
        else:
            raise ValueError(self.policy)
        self.grav_ns += time.perf_counter_ns() - start

        axb, ayb, azb = acceleration
        return (
            vx,
            vy,
            vz,
            cosine * axb - sine * ayb,
            sine * axb + cosine * ayb,
            azb,
        )


def propagate_policy(model, args, policy, y0, duration, time_grid):
    rhs = BlendRhs(model, args, policy)
    start = time.perf_counter()
    solution = solve_ivp(
        rhs,
        (0.0, duration),
        y0,
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        max_step=MAX_STEP,
        t_eval=time_grid,
    )
    wall = time.perf_counter() - start
    if not solution.success:
        raise RuntimeError(solution.message)
    if float(np.min(np.linalg.norm(solution.y[:3], axis=0))) <= model.r_ref:
        raise RuntimeError(f"{policy} reached the reference surface")
    return solution, rhs, wall


def checkpoint_stats(solution_y, truth_y, time_grid, duration):
    result = {}
    for day in CHECKPOINT_DAYS:
        checkpoint = day * DAY
        if checkpoint > duration + 1.0:
            continue
        index = int(np.searchsorted(time_grid, checkpoint - 1e-6))
        index = min(index, time_grid.size - 1)
        result[f"d{day}"] = err_stats(
            solution_y[:, : index + 1], truth_y[:, : index + 1]
        )
    return result


def run_orbit(orbit, duration, output_name, rows, scenario):
    truth_degree = 600 if orbit["hp_km"] < 50.0 else 300
    model = load_model(truth_degree)
    args = kernel_args(model)
    warmup(model, args)
    y0 = initial_state(model, orbit)
    time_grid = np.arange(0.0, duration + 1.0, OUTPUT_STEP)

    row = dict(orbit)
    row["truth_degree"] = truth_degree
    row["policies"] = {}
    orbit_start = time.perf_counter()

    print(
        f"[{len(rows) + 1}/{len(ORBIT_NAMES)}] {orbit['name']} "
        f"truth N={truth_degree}",
        flush=True,
    )
    truth, truth_rhs, truth_wall = propagate(
        model,
        y0,
        duration,
        time_grid,
        lambda t, h: truth_degree,
        args,
        RTOL,
        ATOL,
        max_step=MAX_STEP,
    )
    row["truth"] = {
        "wall_s": truth_wall,
        "n_rhs": truth_rhs.n_calls,
        "grav_s": truth_rhs.grav_ns / 1e9,
    }
    print(
        f"  truth: wall {truth_wall:.1f}s rhs {truth_rhs.n_calls}", flush=True
    )

    policies = (
        "fixed_N120",
        "switch_N30_N120",
        "blend_acceleration",
        "blend_potential_corrected",
    )
    for policy in policies:
        solution, rhs, wall = propagate_policy(
            model, args, policy, y0, duration, time_grid
        )
        statistics = err_stats(solution.y, truth.y)
        statistics.update(
            {
                "wall_s": wall,
                "n_rhs": rhs.n_calls,
                "grav_s": rhs.grav_ns / 1e9,
                "checkpoints": checkpoint_stats(
                    solution.y, truth.y, time_grid, duration
                ),
            }
        )
        row["policies"][policy] = statistics
        print(
            f"  {policy:26s} rms {statistics['pos_rms_m']:11.3f} m "
            f"wall {wall:8.1f}s rhs {rhs.n_calls}",
            flush=True,
        )
        dump(
            output_name,
            {
                "scenario": scenario,
                "rows": rows + [row],
                "complete": False,
                "active_orbit": orbit["name"],
                "last_completed_policy": policy,
            },
        )

    corrected = row["policies"]["blend_potential_corrected"]
    naive = row["policies"]["blend_acceleration"]
    row["comparison"] = {
        "corrected_over_naive_pos_rms": (
            corrected["pos_rms_m"] / naive["pos_rms_m"]
            if naive["pos_rms_m"] > 0.0
            else None
        ),
        "corrected_minus_naive_pos_rms_m": (
            corrected["pos_rms_m"] - naive["pos_rms_m"]
        ),
    }
    row["orbit_wall_s"] = time.perf_counter() - orbit_start
    rows.append(row)
    dump(
        output_name,
        {"scenario": scenario, "rows": rows, "complete": False},
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    duration = 6.0 * 3600.0 if args.smoke else FULL_DAYS * DAY
    output_name = (
        "r9_potential_blend_longarc_smoke.json"
        if args.smoke
        else "r9_potential_blend_longarc.json"
    )
    scenario = {
        "purpose": "long-arc accuracy of corrected potential-level blend",
        "duration_s": duration,
        "output_step_s": OUTPUT_STEP,
        "integrator": "DOP853",
        "rtol": RTOL,
        "atol": ATOL,
        "max_step_s": MAX_STEP,
        "rotation": "uniform lunar sidereal rotation about polar axis",
        "degrees": {"low": N_LO, "high": N_HI},
        "transition_altitude_m": {"near": ALT_NEAR, "far": ALT_FAR},
        "weight": "smoothstep; high-degree weight is 1 below near and 0 above far",
        "corrected_term": "(U_hi-U_lo) * dw/dr * r_hat",
        "truth_rule": "N=600 for perilune below 50 km, otherwise N=300",
        "orbits": list(ORBIT_NAMES),
        "checkpoint_days": [
            day for day in CHECKPOINT_DAYS if day * DAY <= duration + 1.0
        ],
        "potential_path": (
            "production sh_potential_accel_fixed called separately at N=30 "
            "and N=120 on every corrected-blend RHS evaluation"
        ),
        "smoke": args.smoke,
    }

    rows = []
    for name in ORBIT_NAMES:
        run_orbit(orbit_dict(name), duration, output_name, rows, scenario)

    dump(output_name, {"scenario": scenario, "rows": rows, "complete": True})
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
