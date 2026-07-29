"""P0-3: direct 7-day integration-floor convergence.

Five force configurations (N=300 reference, fixed 106, fixed 138, best
upward-quantized schedule, dwell-aware downward schedule) are each run at
three DOP853 settings (rtol/atol 1e-11/1e-4, 1e-12/1e-5, 1e-13/1e-6) on the
50 x 300 km polar arc. Pairwise differences between successive tolerance
levels give the direct per-configuration integration floor; direct
accepted/rejected trial-step counts come from the instrumented solver.
"""

from __future__ import annotations

import numpy as np

from rev3_common import (DAY, dump, err_stats, kernel_args, load_model,
                         make_p_table, alt_sched, eccentric_state,
                         propagate_instr, warmup)

TOLS = [(1e-11, 1e-4), (1e-12, 1e-5), (1e-13, 1e-6)]


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    dur = 7.0 * DAY
    y0 = eccentric_state(model, 50.0, 300.0)
    t_grid = np.arange(0.0, dur + 1.0, 120.0)

    tab_down = make_p_table(model, 1e-3, 60, policy="down")
    tab_up = make_p_table(model, 1e-3, 60, policy="up")
    configs = {
        "ref_300": lambda t, h: 300,
        "fixed_106": lambda t, h: 106,
        "fixed_138": lambda t, h: 138,
        "sched_up": alt_sched(tab_up),
        "sched_down": alt_sched(tab_down),
    }

    results = {}
    states = {}
    for cname, degfun in configs.items():
        results[cname] = {}
        for rtol, atol in TOLS:
            key = f"rtol{rtol:.0e}"
            Y, rhs, info = propagate_instr(model, y0, dur, t_grid, degfun,
                                           args, rtol, atol)
            states[(cname, key)] = Y
            results[cname][key] = {"rtol": rtol, "atol": atol, **info}
            print(f"{cname} {key}: rhs {info['n_rhs']}, steps "
                  f"{info['n_accepted_steps']}, rejected "
                  f"{info['n_rejected_trials']}, wall {info['wall_s']:.1f} s")

    # pairwise tolerance-convergence differences per configuration
    for cname in configs:
        keys = [f"rtol{r:.0e}" for r, _ in TOLS]
        for a, b in zip(keys[:-1], keys[1:]):
            d = err_stats(states[(cname, a)], states[(cname, b)])
            results[cname][f"diff_{a}_vs_{b}"] = d
            print(f"{cname} {a} vs {b}: RMS {d['pos_rms_m']:.3e} m, "
                  f"max {d['pos_max_m']:.3e} m, final {d['pos_final_m']:.3e} m")

    # truncation errors against the tightest common reference
    trunc = {}
    ref = states[("ref_300", "rtol1e-13")]
    for cname in ("fixed_106", "fixed_138", "sched_up", "sched_down"):
        trunc[cname] = {}
        for key in [f"rtol{r:.0e}" for r, _ in TOLS]:
            trunc[cname][key] = err_stats(states[(cname, key)], ref)

    dump("r3_floor7day.json", {
        "scenario": {"type": "eccentric_polar", "perilune_km": 50.0,
                     "apolune_km": 300.0, "duration_s": dur,
                     "integrator": "DOP853 (instrumented)",
                     "tolerance_levels": TOLS, "output_step_s": 120.0,
                     "max_step": "unbounded",
                     "dense_output": "per accepted step, grid sampling",
                     "rotation": "uniform sidereal about polar axis"},
        "schedules": {"down": {f"{k:.0f}": v for k, v in tab_down.items()},
                      "up": {f"{k:.0f}": v for k, v in tab_up.items()}},
        "results": results,
        "truncation_error_vs_ref300_rtol1e-13": trunc,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
