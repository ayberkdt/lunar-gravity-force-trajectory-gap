"""Robustness Test 2: exact event-resolved degree switching.

The existing production-like control selects the active degree at every RHS
evaluation, including Runge--Kutta stages.  This runner implements a second,
independent switching path: each constant-degree arc is integrated separately,
the altitude crossing is located by ``solve_ivp`` root finding, and DOP853 is
restarted from the continuous event state.

Four variants are compared on the 50 x 300 km polar and LRO-like 30 x 216 km
geometries:

* stage-evaluated switching (the existing harness behavior),
* exact event-resolved switching,
* exact switching with a 2 km hysteresis deadband, and
* exact switching with a 600 s minimum dwell.

Full run:  ``.venv\\Scripts\\python.exe robustness_test2_event_switching.py``
Smoke run: add ``--smoke`` (four hours; never overwrites the full artifact).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from rev3_common import (DAY, Rhs, alt_sched, eccentric_state, err_stats,
                         kernel_args, load_model, make_p_table, warmup)


BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
RTOL = 1.0e-11
ATOL = np.array([1.0e-4] * 3 + [1.0e-7] * 3)
MAX_STEP_S = 60.0
OUT_STEP_S = 60.0
EVENT_TIME_GUARD_S = 1.0e-7
EVENT_VALUE_GUARD_M = 1.0e-4


@dataclass(frozen=True)
class Boundary:
    altitude_m: float
    below_degree: int
    above_degree: int


def schedule_degree(table: dict[float, int], altitude_m: float) -> int:
    """Mirror ``rev3_common.alt_sched`` without closing over hidden state."""
    hmin, hmax = min(table), max(table)
    key = min(hmax, max(hmin, 10.0 * math.floor(altitude_m / 1.0e3 / 10.0)))
    return int(table[key])


def schedule_boundaries(table: dict[float, int]) -> list[Boundary]:
    """Return only bin edges where the active degree actually changes."""
    keys = sorted(table)
    out: list[Boundary] = []
    for lower, upper in zip(keys[:-1], keys[1:]):
        n0, n1 = int(table[lower]), int(table[upper])
        if n0 != n1:
            out.append(Boundary(float(upper * 1.0e3), n0, n1))
    return out


def _fixed_rhs(model, args, degree: int):
    return Rhs(model, lambda _t, _h, n=int(degree): n, args)


def propagate_event_resolved(model, args, y0: np.ndarray, duration_s: float,
                             grid: np.ndarray, table: dict[float, int], *,
                             deadband_m: float = 0.0,
                             min_dwell_s: float = 0.0) -> tuple[np.ndarray, dict]:
    """Integrate constant-degree arcs and restart at exact altitude events.

    Hysteresis is a Schmitt trigger: an upward transition occurs at
    ``boundary + deadband`` and its reverse at ``boundary - deadband``.
    Minimum dwell suppresses events until the dwell expires, then immediately
    reconciles the active degree with the current altitude.
    """
    boundaries = schedule_boundaries(table)
    state = np.asarray(y0, dtype=float).copy()
    t_now = 0.0
    current_degree = schedule_degree(table, np.linalg.norm(state[:3]) - model.r_ref)
    eligible_t = 0.0
    out = np.empty((6, grid.size), dtype=float)
    out[:, 0] = state
    next_grid = 1
    events_log: list[dict] = []
    segments = 0
    rhs_calls = 0
    gravity_ns = 0
    started = time.perf_counter()

    while t_now < duration_s - 1.0e-9:
        # Enforce a dwell as a fixed-degree time arc.  At its end, reconcile
        # with the instantaneous altitude before enabling crossings again.
        dwell_only = t_now + 1.0e-9 < eligible_t
        segment_end = min(duration_s, eligible_t) if dwell_only else duration_s
        rhs = _fixed_rhs(model, args, current_degree)
        event_specs: list[tuple[Boundary, int, float]] = []
        event_functions = []

        if not dwell_only:
            start_t = t_now
            start_h = float(np.linalg.norm(state[:3]) - model.r_ref)
            for boundary in boundaries:
                direction = 0
                threshold = boundary.altitude_m
                if current_degree == boundary.below_degree:
                    direction = +1
                    threshold += deadband_m
                elif current_degree == boundary.above_degree:
                    direction = -1
                    threshold -= deadband_m
                else:
                    continue

                # Avoid a spurious zero at a freshly restarted event.  The
                # returned sign is the post-event side and is used only for a
                # sub-microsecond guard interval.
                post_sign = 1.0 if direction > 0 else -1.0

                def event(t, y, *, threshold=threshold, start_t=start_t,
                          start_h=start_h, post_sign=post_sign):
                    value = float(np.linalg.norm(y[:3]) - model.r_ref - threshold)
                    if (t - start_t <= EVENT_TIME_GUARD_S
                            and abs(start_h - threshold) <= 1.0e-3):
                        return post_sign * EVENT_VALUE_GUARD_M
                    return value

                event.terminal = True
                event.direction = direction
                event_functions.append(event)
                event_specs.append((boundary, direction, threshold))

        sol = solve_ivp(
            rhs, (t_now, segment_end), state, method="DOP853",
            rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S, dense_output=True,
            events=event_functions if event_functions else None,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        segments += 1
        rhs_calls += rhs.n_calls
        gravity_ns += rhs.grav_ns
        t_stop = float(sol.t[-1])
        while next_grid < grid.size and grid[next_grid] <= t_stop + 1.0e-8:
            out[:, next_grid] = sol.sol(float(grid[next_grid]))
            next_grid += 1
        state = np.asarray(sol.y[:, -1], dtype=float)
        previous_t = t_now
        t_now = t_stop

        hit_index = None
        if event_functions:
            for index, times in enumerate(sol.t_events):
                if len(times):
                    hit_index = index
                    state = np.asarray(sol.y_events[index][-1], dtype=float)
                    t_now = float(times[-1])
                    break

        if hit_index is not None:
            boundary, direction, threshold = event_specs[hit_index]
            old_degree = current_degree
            current_degree = (boundary.above_degree if direction > 0
                              else boundary.below_degree)
            altitude = float(np.linalg.norm(state[:3]) - model.r_ref)
            events_log.append({
                "t_s": t_now,
                "altitude_m": altitude,
                "target_altitude_m": threshold,
                "root_residual_m": altitude - threshold,
                "direction": "ascending" if direction > 0 else "descending",
                "from_degree": int(old_degree),
                "to_degree": int(current_degree),
            })
            eligible_t = t_now + min_dwell_s
            continue

        if dwell_only and t_now >= eligible_t - 1.0e-8:
            altitude = float(np.linalg.norm(state[:3]) - model.r_ref)
            desired = schedule_degree(table, altitude)
            if desired != current_degree:
                events_log.append({
                    "t_s": t_now, "altitude_m": altitude,
                    "target_altitude_m": None, "root_residual_m": None,
                    "direction": "dwell_expiry",
                    "from_degree": int(current_degree),
                    "to_degree": int(desired),
                })
                current_degree = desired
                eligible_t = t_now + min_dwell_s
            continue

        if t_now <= previous_t + 1.0e-10:
            raise RuntimeError("event-resolved propagation made no time progress")
        if t_now >= duration_s - 1.0e-9:
            break

    if next_grid != grid.size:
        raise RuntimeError(f"filled {next_grid}/{grid.size} output samples")
    residuals = [abs(e["root_residual_m"]) for e in events_log
                 if e["root_residual_m"] is not None]
    info = {
        "segments": segments,
        "events": events_log,
        "event_count": len(events_log),
        "max_root_residual_m": float(max(residuals, default=0.0)),
        "n_rhs": rhs_calls,
        "gravity_wall_s": gravity_ns / 1.0e9,
        "wall_s": time.perf_counter() - started,
    }
    return out, info


def run_geometry(name: str, y0: np.ndarray, model, args, duration_s: float,
                 table: dict[float, int]) -> dict:
    grid = np.arange(0.0, duration_s + 0.1, OUT_STEP_S)
    truth_rhs = _fixed_rhs(model, args, 300)
    truth = solve_ivp(
        truth_rhs, (0.0, duration_s), y0, method="DOP853", t_eval=grid,
        rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S,
    )
    if not truth.success:
        raise RuntimeError(truth.message)

    stage_rhs = Rhs(model, alt_sched(table), args)
    stage = solve_ivp(
        stage_rhs, (0.0, duration_s), y0, method="DOP853", t_eval=grid,
        rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S,
    )
    if not stage.success:
        raise RuntimeError(stage.message)

    variants = {
        "exact_event": {"deadband_m": 0.0, "min_dwell_s": 0.0},
        "exact_hysteresis_2km": {"deadband_m": 2.0e3, "min_dwell_s": 0.0},
        "exact_min_dwell_600s": {"deadband_m": 0.0, "min_dwell_s": 600.0},
    }
    states = {"stage_evaluated": stage.y}
    runs = {
        "stage_evaluated": {
            "n_rhs": int(stage_rhs.n_calls),
            "gravity_wall_s": stage_rhs.grav_ns / 1.0e9,
            "switch_count_at_rhs_samples": int(stage_rhs.n_deg_changes),
            "distinct_degrees": int(len(stage_rhs.deg_counts)),
        }
    }
    for variant, settings in variants.items():
        print(f"  {name}: {variant}", flush=True)
        state, info = propagate_event_resolved(
            model, args, y0, duration_s, grid, table, **settings,
        )
        states[variant] = state
        runs[variant] = info

    errors = {key: err_stats(value, truth.y) for key, value in states.items()}
    stage_rms = errors["stage_evaluated"]["pos_rms_m"]
    exact_rms = errors["exact_event"]["pos_rms_m"]
    exact_ric = errors["exact_event"]["ric_rms_m"]
    ric_norm = math.sqrt(sum(value * value for value in exact_ric.values()))
    return {
        "name": name,
        "truth_run": {"degree": 300, "n_rhs": int(truth_rhs.n_calls)},
        "runs": runs,
        "errors_vs_N300": errors,
        "exact_minus_stage": err_stats(states["exact_event"], states["stage_evaluated"]),
        "diagnostics": {
            "exact_to_stage_rms_ratio": float(exact_rms / max(stage_rms, 1.0e-30)),
            "exact_in_track_rms_fraction": float(abs(exact_ric["in_track"]) / max(ric_norm, 1.0e-30)),
            "policy_bias_persists_after_exact_switching": bool(
                exact_rms >= 0.5 * stage_rms and exact_rms >= 10.0
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="run four hours and write a _smoke artifact")
    args_cli = parser.parse_args()
    duration = 4.0 * 3600.0 if args_cli.smoke else 7.0 * DAY
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    table = make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="down")
    geometries = [
        ("50x300_polar", eccentric_state(model, 50.0, 300.0, 90.0)),
        ("lro_like_30x216", eccentric_state(model, 30.0, 216.0, 90.0)),
    ]
    payload = {
        "schema": "robustness_test2_event_switching_v1",
        "formal_run": not args_cli.smoke,
        "scenario": {
            "duration_s": duration, "output_step_s": OUT_STEP_S,
            "integrator": "SciPy DOP853",
            "rtol": RTOL, "atol_position_m": float(ATOL[0]),
            "atol_velocity_m_s": float(ATOL[3]), "max_step_s": MAX_STEP_S,
            "rotation": "uniform lunar sidereal rotation",
            "schedule": "Kaula p-tail, eps=1e-3, floor=60, cap=138, q=10 down",
            "event_method": "terminal altitude roots; continuous state restart with fixed degree per segment",
        },
        "geometries": [],
    }
    for name, y0 in geometries:
        print(f"== {name} ==", flush=True)
        payload["geometries"].append(
            run_geometry(name, y0, model, args, duration, table)
        )
    suffix = "_smoke" if args_cli.smoke else ""
    out = METRICS / f"robustness_test2_event_switching{suffix}.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {out.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
