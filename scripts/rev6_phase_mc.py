"""R6-C: initial-phase Monte Carlo of the 7-day scheduling penalty.

The core 50x300 km polar scheduling result rests on specific initial
phases. This experiment repeats the 7-day comparison from eight initial
true anomalies (0..315 deg in 45-deg steps), each with its own N=300
reference, for fixed 138/106 and the downward and upward schedules, at
the scalar tight tolerance of the screening matrix (rtol 1e-12,
atol 1e-5). Reported: per-phase errors and the distribution of the
schedule-to-fixed penalty ratio.
"""

from __future__ import annotations

import math

import numpy as np

from rev3_common import (DAY, dump, err_stats, kernel_args, load_model,
                         make_p_table, alt_sched, propagate, warmup)

RTOL, ATOL = 1e-12, 1e-5
PHASES_DEG = list(range(0, 360, 45))


def state_at_anomaly(model, hp_km, ha_km, nu_deg):
    """Polar-plane (x-z) state at true anomaly nu on the ellipse."""
    rp = model.r_ref + hp_km * 1e3
    ra = model.r_ref + ha_km * 1e3
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    p = a * (1.0 - e * e)
    nu = math.radians(nu_deg)
    r = p / (1.0 + e * math.cos(nu))
    vr = math.sqrt(model.mu / p) * e * math.sin(nu)
    vt = math.sqrt(model.mu / p) * (1.0 + e * math.cos(nu))
    rhat = np.array([math.cos(nu), 0.0, math.sin(nu)])
    that = np.array([-math.sin(nu), 0.0, math.cos(nu)])
    pos = r * rhat
    vel = vr * rhat + vt * that
    return np.concatenate([pos, vel])


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    dur = 7.0 * DAY
    t_grid = np.arange(0.0, dur + 1.0, 120.0)

    policies = {
        "fixed_138": lambda t, h: 138,
        "fixed_106": lambda t, h: 106,
        "sched_down": alt_sched(make_p_table(model, 1e-3, 60, policy="down")),
        "sched_up": alt_sched(make_p_table(model, 1e-3, 60, policy="up")),
    }

    rows = []
    for nu in PHASES_DEG:
        y0 = state_at_anomaly(model, 50.0, 300.0, nu)
        truth, _, _ = propagate(model, y0, dur, t_grid, lambda t, h: 300,
                                args, RTOL, ATOL)
        phase_row = {"nu_deg": nu, "policies": {}}
        for pname, degfun in policies.items():
            sol, rhs, wall = propagate(model, y0, dur, t_grid, degfun, args,
                                       RTOL, ATOL)
            st = err_stats(sol.y, truth.y)
            st.update({"n_rhs": rhs.n_calls, "wall_s": wall,
                       "grav_s": rhs.grav_ns / 1e9})
            phase_row["policies"][pname] = st
        f138 = phase_row["policies"]["fixed_138"]["pos_rms_m"]
        for s in ("sched_down", "sched_up"):
            phase_row["policies"][s]["penalty_vs_fixed138"] = \
                phase_row["policies"][s]["pos_rms_m"] / f138
        rows.append(phase_row)
        print(f"nu={nu:3d}: f138 {f138:7.2f}  f106 "
              f"{phase_row['policies']['fixed_106']['pos_rms_m']:7.2f}  "
              f"down {phase_row['policies']['sched_down']['pos_rms_m']:8.2f} "
              f"(x{phase_row['policies']['sched_down']['penalty_vs_fixed138']:.1f})  "
              f"up {phase_row['policies']['sched_up']['pos_rms_m']:8.2f} "
              f"(x{phase_row['policies']['sched_up']['penalty_vs_fixed138']:.1f})")

    summary = {}
    for pname in policies:
        vals = np.array([r["policies"][pname]["pos_rms_m"] for r in rows])
        summary[pname] = {"min_m": float(vals.min()),
                          "median_m": float(np.median(vals)),
                          "max_m": float(vals.max())}
    for s in ("sched_down", "sched_up"):
        pen = np.array([r["policies"][s]["penalty_vs_fixed138"] for r in rows])
        summary[s]["penalty_min"] = float(pen.min())
        summary[s]["penalty_median"] = float(np.median(pen))
        summary[s]["penalty_max"] = float(pen.max())
    print("summary:", summary)

    dump("r6_phase_mc.json", {
        "scenario": {"type": "eccentric_polar_50x300", "duration_s": dur,
                     "truth_degree": 300, "integrator": "DOP853",
                     "rtol": RTOL, "atol": ATOL, "output_step_s": 120.0,
                     "phases_deg": PHASES_DEG,
                     "rotation": "uniform sidereal about polar axis"},
        "rows": rows,
        "summary": summary,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
