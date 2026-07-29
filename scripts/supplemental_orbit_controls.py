"""Orbit controls requested after the publication-readiness experiment set.

The script measures: (1) the seven-day N=300 versus N=600 reference error,
(2) vector-atol tight/tighter convergence for two narrow-gap geometries and
geometry-specific fixed degrees, and (3) max-step sensitivity of a switching
case. All output comparisons use the common 120-s grid.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from rev3_common import (
    DAY,
    Rhs,
    alt_sched,
    dump,
    eccentric_state,
    err_stats,
    kernel_args,
    load_model,
    make_p_table,
    propagate_instr,
    warmup,
)
from rev3_longarc_matrix import ET0_TDB_S, KERNEL_DIR, RhsMoonPA, SPICE_KERNELS


DURATION = 7.0 * DAY
GRID = np.arange(0.0, DURATION + 1.0, 120.0)
TIGHT = {"rtol": 1e-12, "atol": np.array([1e-5] * 3 + [1e-8] * 3)}
TIGHTER = {"rtol": 1e-13, "atol": np.array([1e-6] * 3 + [1e-9] * 3)}


def run_uniform(model, args, y0, degree_of, tol, max_step):
    Y, rhs, info = propagate_instr(
        model, y0, DURATION, GRID, degree_of, args,
        tol["rtol"], tol["atol"], max_step=max_step,
    )
    return Y, {**info, "degree_counts": rhs.deg_counts}


def run_moonpa(model, args, y0, degree, tol, max_step):
    rhs = RhsMoonPA(model, lambda t, h: degree, args)
    start = time.perf_counter()
    sol = solve_ivp(
        rhs, (0.0, DURATION), y0, method="DOP853", t_eval=GRID,
        rtol=tol["rtol"], atol=tol["atol"], max_step=max_step,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y, {
        "n_rhs": rhs.n_calls,
        "wall_s": time.perf_counter() - start,
        "grav_s": rhs.grav_ns / 1e9,
    }


def tolerance_case(model, args, y0, policies):
    states = {}
    runs = {}
    for label, tolerance in (("tight", TIGHT), ("tighter", TIGHTER)):
        for policy, degree_of in policies.items():
            print(f"{policy} {label}")
            states[(policy, label)], runs[(policy, label)] = run_uniform(
                model, args, y0, degree_of, tolerance, 60.0,
            )
    out = {"runs": {}, "self_convergence": {}, "error_vs_same_tolerance_ref300": {}}
    for (policy, label), info in runs.items():
        out["runs"].setdefault(policy, {})[label] = info
    for policy in policies:
        out["self_convergence"][policy] = err_stats(
            states[(policy, "tight")], states[(policy, "tighter")]
        )
        if policy != "ref300":
            out["error_vs_same_tolerance_ref300"][policy] = {
                label: err_stats(states[(policy, label)], states[("ref300", label)])
                for label in ("tight", "tighter")
            }
    return out


def main() -> None:
    model600 = load_model(600)
    args600 = kernel_args(model600)
    warmup(model600, args600)
    y0_50 = eccentric_state(model600, 50.0, 300.0)

    print("N=300 versus N=600, uniform rotation")
    y300, i300 = run_uniform(model600, args600, y0_50, lambda t, h: 300, TIGHTER, 60.0)
    y600, i600 = run_uniform(model600, args600, y0_50, lambda t, h: 600, TIGHTER, 60.0)
    reference_degree = {
        "uniform_rotation": {
            "N300": i300, "N600": i600,
            "difference_N300_minus_N600": err_stats(y300, y600),
        }
    }

    print("N=300 versus N=600, DE440 MOON_PA")
    import spiceypy as sp
    for kernel in SPICE_KERNELS:
        sp.furnsh(str(Path(KERNEL_DIR) / kernel))
    ym300, im300 = run_moonpa(model600, args600, y0_50, 300, TIGHTER, 60.0)
    ym600, im600 = run_moonpa(model600, args600, y0_50, 600, TIGHTER, 60.0)
    reference_degree["moon_pa"] = {
        "epoch_tdb_j2000_s": ET0_TDB_S,
        "N300": im300, "N600": im600,
        "difference_N300_minus_N600": err_stats(ym300, ym600),
    }

    model300 = load_model(300)
    args300 = kernel_args(model300)
    warmup(model300, args300)
    down = alt_sched(make_p_table(model300, 1e-3, 60, policy="down"))
    up = alt_sched(make_p_table(model300, 1e-3, 60, policy="up"))

    print("narrow-gap: 60-degree inclination")
    inc60 = tolerance_case(model300, args300,
        eccentric_state(model300, 50.0, 300.0, incl_deg=60.0), {
            "ref300": lambda t, h: 300,
            "fixed130_calibrated": lambda t, h: 130,
            "fixed138_conservative": lambda t, h: 138,
            "sched_up": up,
        })
    print("narrow-gap: 100x300 km")
    alt100 = tolerance_case(model300, args300,
        eccentric_state(model300, 100.0, 300.0), {
            "ref300": lambda t, h: 300,
            "fixed70_calibrated": lambda t, h: 70,
            "fixed80_percentile": lambda t, h: 80,
            "fixed106": lambda t, h: 106,
            "sched_down": down,
        })

    print("max-step sensitivity")
    step_states = {}
    step_runs = {}
    for max_step in (30.0, 60.0, 120.0, math.inf):
        key = "unbounded" if math.isinf(max_step) else f"{max_step:.0f}s"
        for policy, degree_of in (("ref300", lambda t, h: 300), ("sched_down", down)):
            print(key, policy)
            step_states[(key, policy)], step_runs[(key, policy)] = run_uniform(
                model300, args300, eccentric_state(model300, 50.0, 300.0),
                degree_of, TIGHT, max_step,
            )
    max_step_out = {"runs": {}, "schedule_error_vs_same_step_ref300": {},
                    "state_difference_vs_30s": {}}
    for (key, policy), info in step_runs.items():
        max_step_out["runs"].setdefault(key, {})[policy] = info
    for key in ("30s", "60s", "120s", "unbounded"):
        max_step_out["schedule_error_vs_same_step_ref300"][key] = err_stats(
            step_states[(key, "sched_down")], step_states[(key, "ref300")]
        )
        if key != "30s":
            max_step_out["state_difference_vs_30s"][key] = {
                policy: err_stats(step_states[(key, policy)], step_states[("30s", policy)])
                for policy in ("ref300", "sched_down")
            }

    dump("supplemental_orbit_controls.json", {
        "scenario": {
            "duration_s": DURATION,
            "output_step_s": 120.0,
            "integrator": "SciPy DOP853",
            "tight": {"rtol": TIGHT["rtol"], "position_atol_m": 1e-5,
                      "velocity_atol_m_s": 1e-8},
            "tighter": {"rtol": TIGHTER["rtol"], "position_atol_m": 1e-6,
                        "velocity_atol_m_s": 1e-9},
        },
        "reference_degree_control": reference_degree,
        "narrow_gap_convergence": {"inclination_60deg": inc60,
                                   "100x300km": alt100},
        "max_step_sensitivity": max_step_out,
    })


if __name__ == "__main__":
    main()
