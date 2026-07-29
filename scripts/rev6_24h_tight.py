"""R6-B: 24-hour eccentric scheduling comparison at vector tolerances.

The archived 24-hour scheduling comparison ran at the scalar baseline
tolerance, whose directly measured 24-hour floor (4.4 m RMS on the
circular control) is comparable to the meter-level results it reported.
This experiment re-runs the 24-hour 50x300 km polar comparison with
component-scaled vector tolerances and measures each configuration's own
floor by tolerance tightening, so the short-arc conclusion stands on
numbers above their own numerical noise.

Settings: DOP853, rtol 1e-12 with atol (1e-5 m, 1e-8 m/s), max step 60 s
(tight); rtol 1e-13 with atol (1e-6 m, 1e-9 m/s) (tighter control).
"""

from __future__ import annotations

import numpy as np

from rev3_common import (DAY, dump, err_stats, kernel_args, load_model,
                         make_p_table, alt_sched, eccentric_state,
                         propagate_instr, warmup)

TIGHT = dict(rtol=1e-12, atol=np.array([1e-5] * 3 + [1e-8] * 3))
TIGHTER = dict(rtol=1e-13, atol=np.array([1e-6] * 3 + [1e-9] * 3))
MAX_STEP = 60.0


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    dur = 1.0 * DAY
    y0 = eccentric_state(model, 50.0, 300.0)
    t_grid = np.arange(0.0, dur + 1.0, 120.0)

    policies = {
        "ref_300": lambda t, h: 300,
        "fixed_138": lambda t, h: 138,
        "fixed_106": lambda t, h: 106,
        "sched_down": alt_sched(make_p_table(model, 1e-3, 60, policy="down")),
        "sched_up": alt_sched(make_p_table(model, 1e-3, 60, policy="up")),
    }

    states = {}
    info_all = {}
    for pname, degfun in policies.items():
        for lvl, cfg in (("tight", TIGHT), ("tighter", TIGHTER)):
            Y, rhs, info = propagate_instr(model, y0, dur, t_grid, degfun,
                                           args, cfg["rtol"], cfg["atol"],
                                           max_step=MAX_STEP)
            states[(pname, lvl)] = Y
            info["mean_deg_sq"] = rhs.sum_deg_sq / rhs.n_calls
            info["proxy_rel_fixed138"] = (rhs.sum_deg_sq / rhs.n_calls) / 138.0**2
            info_all[f"{pname}/{lvl}"] = info
            print(f"{pname}/{lvl}: rhs {info['n_rhs']} acc "
                  f"{info['n_accepted_steps']} rej {info['n_rejected_trials']}")

    rows = []
    for pname in policies:
        floor = err_stats(states[(pname, "tight")], states[(pname, "tighter")])
        row = {"policy": pname,
               "floor_tight_vs_tighter_rms_m": floor["pos_rms_m"],
               "floor_max_m": floor["pos_max_m"],
               **{k: info_all[f"{pname}/tight"][k]
                  for k in ("n_rhs", "n_accepted_steps", "n_rejected_trials",
                            "wall_s", "grav_s", "proxy_rel_fixed138")}}
        if pname != "ref_300":
            for lvl in ("tight", "tighter"):
                e = err_stats(states[(pname, lvl)], states[("ref_300", lvl)])
                row[f"err_{lvl}"] = e
            print(f"{pname}: RMS {row['err_tight']['pos_rms_m']:.3f} m "
                  f"(tighter {row['err_tighter']['pos_rms_m']:.3f}), "
                  f"floor {row['floor_tight_vs_tighter_rms_m']:.3f} m")
        rows.append(row)

    dump("r6_24h_tight.json", {
        "scenario": {"type": "eccentric_polar_50x300", "duration_s": dur,
                     "truth_degree": 300, "integrator": "DOP853",
                     "tight": {"rtol": 1e-12, "atol_pos_m": 1e-5,
                               "atol_vel_m_s": 1e-8},
                     "tighter": {"rtol": 1e-13, "atol_pos_m": 1e-6,
                                 "atol_vel_m_s": 1e-9},
                     "max_step_s": MAX_STEP, "output_step_s": 120.0,
                     "rotation": "uniform sidereal about polar axis"},
        "rows": rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
