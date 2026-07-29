"""P0-6: direct DOP853 rejection statistics around degree switches.

Replaces the r1 inference (attempts = round(nfev/12)) with direct counts
from an instrumented solver whose error-norm evaluations are intercepted
(one per attempted trial step). The counter is validated on a controlled
ODE where the RHS-call bookkeeping nfev = 2 + 12 x attempts holds exactly.

Outputs an event-aligned aggregate over ALL switches (not a single
representative switch): median accepted step size and IQR versus time from
the switch, rejection probability per bin, upward and downward switches
separated, for both initial phases.
"""

from __future__ import annotations

import math

import numpy as np

from rev3_common import (DAY, InstrumentedDOP853, Rhs, dump, kernel_args,
                         load_model, make_p_table, alt_sched,
                         eccentric_state, orbit_period, warmup)

RTOL, ATOL = 1e-11, 1e-4
WINDOW_S = 1800.0
BIN_S = 120.0


def counter_validation() -> dict:
    calls = {"n": 0}

    def f(t, y):
        calls["n"] += 1
        return [y[1], -y[0] * (1.0 + 5.0 * np.tanh(50.0 * np.sin(3.0 * t)))]

    s = InstrumentedDOP853(f, 0.0, [1.0, 0.0], 60.0, rtol=1e-10, atol=1e-10)
    acc = 0
    while s.status == "running":
        s.step()
        acc += 1
    return {
        "ode": "y'' = -y (1 + 5 tanh(50 sin 3t)), rtol=atol=1e-10, t in [0,60]",
        "accepted": acc, "attempts": s.n_attempts, "rejected": s.n_rejected,
        "identity_attempts_eq_accepted_plus_rejected":
            bool(s.n_attempts == acc + s.n_rejected),
        "nfev": calls["n"], "nfev_predicted_2_plus_12x_attempts":
            2 + 12 * s.n_attempts,
        "nfev_match": bool(calls["n"] == 2 + 12 * s.n_attempts),
    }


def run_case(model, args, y0, dur, degfun):
    rhs = Rhs(model, degfun, args)
    solver = InstrumentedDOP853(rhs, 0.0, np.asarray(y0, float), dur,
                                rtol=RTOL, atol=ATOL)
    t_acc, h_acc, rej_acc, deg_acc = [], [], [], []
    prev = solver.n_attempts
    while solver.status == "running":
        solver.step()
        t_acc.append(float(solver.t))
        h_acc.append(float(solver.t - solver.t_old))
        rej_acc.append(solver.n_attempts - prev - 1)
        prev = solver.n_attempts
        r = float(np.linalg.norm(solver.y[:3]))
        deg_acc.append(int(degfun(float(solver.t), r - model.r_ref)))
    return (np.array(t_acc), np.array(h_acc), np.array(rej_acc),
            np.array(deg_acc), rhs, solver)


def event_aligned(t_acc, h_acc, rej_acc, deg_acc):
    """Stack accepted steps around every switch; per-bin medians/IQR and
    rejection probability, split by switch direction."""
    switches = np.flatnonzero(deg_acc[1:] != deg_acc[:-1]) + 1
    directions = np.sign(deg_acc[switches] - deg_acc[switches - 1])
    edges = np.arange(-WINDOW_S, WINDOW_S + BIN_S, BIN_S)
    centers = 0.5 * (edges[:-1] + edges[1:])
    out = {}
    for label, mask in (("down", directions < 0), ("up", directions > 0)):
        sw_times = t_acc[switches[mask]]
        h_bins = [[] for _ in centers]
        rej_bins = [[] for _ in centers]
        for ts in sw_times:
            rel = t_acc - ts
            sel = np.abs(rel) <= WINDOW_S
            idx = np.clip(((rel[sel] + WINDOW_S) // BIN_S).astype(int),
                          0, len(centers) - 1)
            for i, hh, rr in zip(idx, h_acc[sel], rej_acc[sel]):
                h_bins[i].append(hh)
                rej_bins[i].append(rr)
        med = [float(np.median(b)) if b else None for b in h_bins]
        q1 = [float(np.percentile(b, 25)) if b else None for b in h_bins]
        q3 = [float(np.percentile(b, 75)) if b else None for b in h_bins]
        rejp = [float(np.mean(np.array(b) > 0)) if b else None
                for b in rej_bins]
        out[label] = {"n_switches": int(mask.sum()),
                      "bin_center_s": centers.tolist(),
                      "median_step_s": med, "q1_step_s": q1,
                      "q3_step_s": q3, "rejection_probability": rejp}
    return out, int(len(switches))


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    period = orbit_period(model, 50.0, 300.0)
    dur = 2.2 * period
    sched = alt_sched(make_p_table(model, 1e-3, 60, policy="down"))

    cases = {
        "perilune_start": eccentric_state(model, 50.0, 300.0),
        "apolune_start": eccentric_state(model, 50.0, 300.0,
                                         at_apolune=True),
    }
    payload = {"counter_validation": counter_validation(),
               "scenario": {"perilune_km": 50.0, "apolune_km": 300.0,
                            "duration_s": dur, "rtol": RTOL, "atol": ATOL,
                            "window_s": WINDOW_S, "bin_s": BIN_S,
                            "rotation": "uniform sidereal about polar axis"},
               "cases": {}}
    print("counter validation:", payload["counter_validation"])

    for phase, y0 in cases.items():
        for name, degfun in (("scheduled", sched),
                             ("fixed_138", lambda t, h: 138)):
            t_acc, h_acc, rej_acc, deg_acc, rhs, solver = run_case(
                model, args, y0, dur, degfun)
            aligned, n_sw = event_aligned(t_acc, h_acc, rej_acc, deg_acc)
            near = np.zeros(len(t_acc), dtype=bool)
            sw_idx = np.flatnonzero(deg_acc[1:] != deg_acc[:-1]) + 1
            for si in sw_idx:
                near |= np.abs(t_acc - t_acc[si]) <= 600.0
            stats = {
                "n_accepted": int(len(t_acc)),
                "n_attempts": int(solver.n_attempts),
                "n_rejected_direct": int(solver.n_rejected),
                "n_rhs": int(rhs.n_calls),
                "n_switches": n_sw,
                "median_step_s": float(np.median(h_acc)),
                "median_step_near_switch_s":
                    float(np.median(h_acc[near])) if near.any() else None,
                "median_step_away_s":
                    float(np.median(h_acc[~near])) if (~near).any() else None,
                "rejected_near_switch": int(rej_acc[near].sum()),
                "rejected_away": int(rej_acc[~near].sum()),
                "steps_near_switch": int(near.sum()),
            }
            payload["cases"][f"{phase}/{name}"] = {
                "stats": stats, "event_aligned": aligned,
                "series": {"t_s": t_acc.tolist(), "h_s": h_acc.tolist(),
                           "rejected_trials": rej_acc.tolist(),
                           "degree": deg_acc.tolist()},
            }
            print(f"{phase}/{name}: accepted {stats['n_accepted']}, "
                  f"rejected {stats['n_rejected_direct']} (direct), "
                  f"switches {stats['n_switches']}")

    dump("r3_switch_direct.json", payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
