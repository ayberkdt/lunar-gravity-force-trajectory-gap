"""R6-D: direct accepted/rejected telemetry for the MOON_PA scheduling case.

The multi-geometry screening matrix reported '-' for the MOON_PA
accept/reject columns because the low-level counter was not wired to the
ephemeris-driven RHS. This experiment re-runs the full DE440/MOON_PA
seven-day case with the instrumented solver and the SPICE RHS, filling
those entries and cross-checking the error statistics against the
archived matrix run (same scalar tight tolerance).
"""

from __future__ import annotations

import numpy as np

from rev3_common import (DAY, dump, err_stats, kernel_args, load_model,
                         degree_power, make_p_table, make_emp_table,
                         alt_sched, mindwell_from_profile, eccentric_state,
                         propagate, propagate_instr, warmup)
from rev3_longarc_matrix import RhsMoonPA, KERNEL_DIR, SPICE_KERNELS

RTOL, ATOL = 1e-12, 1e-5


def main() -> int:
    import spiceypy as sp
    for k in SPICE_KERNELS:
        sp.furnsh(str(KERNEL_DIR / k))

    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    power = degree_power(load_model(1800))
    dur = 7.0 * DAY
    t_grid = np.arange(0.0, dur + 1.0, 120.0)
    y0 = eccentric_state(model, 50.0, 300.0)

    # min-dwell time schedule from a quick uniform-rotation altitude profile
    truth_sched, _, _ = propagate(model, y0, dur, t_grid, lambda t, h: 60,
                                  args, 1e-9, 1e-2)
    tab_down = make_p_table(model, 1e-3, 60, policy="down")
    mind, _, _ = mindwell_from_profile(truth_sched, model, tab_down, dur)
    policies = {
        "fixed_138": lambda t, h: 138,
        "fixed_106": lambda t, h: 106,
        "sched_down": alt_sched(tab_down),
        "sched_up": alt_sched(make_p_table(model, 1e-3, 60, policy="up")),
        "sched_mindwell600": mind,
        "sched_emp": alt_sched(make_emp_table(model, power, 1e-3, 60)),
    }

    rhs_t = RhsMoonPA(model, lambda t, h: 300, args)
    Yt, _, info_t = propagate_instr(model, y0, dur, t_grid, None, args,
                                    RTOL, ATOL, rhs_obj=rhs_t)
    print(f"ref_300: rhs {info_t['n_rhs']} acc {info_t['n_accepted_steps']} "
          f"rej {info_t['n_rejected_trials']}")

    rows = []
    for pname, degfun in policies.items():
        rhs = RhsMoonPA(model, degfun, args)
        Y, _, info = propagate_instr(model, y0, dur, t_grid, None, args,
                                     RTOL, ATOL, rhs_obj=rhs)
        st = err_stats(Y, Yt)
        st.update({"policy": pname, **info,
                   "mean_deg_sq": rhs.sum_deg_sq / rhs.n_calls})
        rows.append(st)
        print(f"{pname}: RMS {st['pos_rms_m']:8.2f} m  acc "
              f"{info['n_accepted_steps']} rej {info['n_rejected_trials']} "
              f"rhs {info['n_rhs']}")

    dump("r6_moonpa_telemetry.json", {
        "scenario": {"type": "eccentric_polar_50x300_moonpa",
                     "duration_s": dur, "truth_degree": 300,
                     "integrator": "DOP853 (instrumented)",
                     "rtol": RTOL, "atol": ATOL, "output_step_s": 120.0,
                     "rotation": "DE440 MOON_PA via SPICE pxform per RHS"},
        "reference_info": info_t,
        "rows": rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
