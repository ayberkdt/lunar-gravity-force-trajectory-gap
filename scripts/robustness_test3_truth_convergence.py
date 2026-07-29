"""Robustness Test 3: LRO-like reference-degree convergence.

Runs the same MOON_PA LRO-like initial condition at N=300, 600, and 900 for
both gravity-only and expanded Earth/Sun/SRP/eclipsed force models.  Scheduled
and fixed critical-degree policies are also evaluated against N=900 so the
reference self-convergence can be expressed as a fraction of the interpreted
policy discrepancy.

Full run:  ``.venv\\Scripts\\python.exe robustness_test3_truth_convergence.py``
Smoke run: add ``--smoke`` (four hours; writes a separate artifact).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from rev3_common import (DAY, alt_sched, err_stats, kernel_args, load_model,
                         make_p_table, warmup)
from rev4_robustness_controls import (
    _run_expanded,
    build_ephemeris,
    lro_like_state,
)


BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
REFERENCE_DEGREES = (300, 600, 900)
POLICY_FRACTION_LIMIT = 0.05
ABSOLUTE_GAP_LIMIT_M = 5.0


def _run_force_set(name: str, model, args, ephem, y0: np.ndarray,
                   duration_s: float, *, expanded: bool) -> dict:
    use_common_forces = bool(expanded)
    states: dict[str, np.ndarray] = {}
    runs: dict[str, dict] = {}
    policies = {
        **{f"reference_N{degree}": (lambda _t, _h, n=degree: n)
           for degree in REFERENCE_DEGREES},
        "fixed_N194_critical": lambda _t, _h: 194,
        "schedule_down": alt_sched(
            make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="down")
        ),
        "schedule_up": alt_sched(
            make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="up")
        ),
    }
    for policy, degree_of in policies.items():
        print(f"  {name}: {policy}", flush=True)
        sol, info = _run_expanded(
            model, args, ephem, y0, degree_of, duration_s,
            use_third_body=use_common_forces, use_srp=use_common_forces,
            track_components=(policy == "reference_N900"),
        )
        states[policy] = sol.y
        runs[policy] = info

    n900 = states["reference_N900"]
    convergence = {
        "N300_minus_N600": err_stats(states["reference_N300"],
                                      states["reference_N600"]),
        "N600_minus_N900": err_stats(states["reference_N600"], n900),
        "N300_minus_N900": err_stats(states["reference_N300"], n900),
    }
    policy_errors = {
        policy: err_stats(states[policy], n900)
        for policy in ("fixed_N194_critical", "schedule_down", "schedule_up")
    }
    gap_rms = convergence["N600_minus_N900"]["pos_rms_m"]
    policy_rms = {key: value["pos_rms_m"] for key, value in policy_errors.items()}
    smallest_policy_rms = min(policy_rms.values())
    allowed = min(ABSOLUTE_GAP_LIMIT_M,
                  POLICY_FRACTION_LIMIT * smallest_policy_rms)
    gap_ric = convergence["N600_minus_N900"]["ric_rms_m"]
    gap_ric_norm = math.sqrt(sum(value * value for value in gap_ric.values()))
    return {
        "force_set": name,
        "runs": runs,
        "reference_self_convergence": convergence,
        "policy_errors_vs_N900": policy_errors,
        "acceptance": {
            "absolute_limit_m": ABSOLUTE_GAP_LIMIT_M,
            "policy_fraction_limit": POLICY_FRACTION_LIMIT,
            "smallest_policy_rms_m": smallest_policy_rms,
            "effective_rms_limit_m": allowed,
            "N600_to_N900_rms_m": gap_rms,
            "N600_to_N900_fraction_of_each_policy": {
                key: gap_rms / max(value, 1.0e-30)
                for key, value in policy_rms.items()
            },
            "N600_to_N900_in_track_fraction": (
                abs(gap_ric["in_track"]) / max(gap_ric_norm, 1.0e-30)
            ),
            "passes": bool(gap_rms <= allowed),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="run four hours and write a _smoke artifact")
    args_cli = parser.parse_args()
    duration = 4.0 * 3600.0 if args_cli.smoke else 7.0 * DAY

    # One N=900 model is intentional: every runtime truncation sees the same
    # parsed source coefficients and kernel implementation, isolating N only.
    model = load_model(900)
    args = kernel_args(model)
    warmup(model, args)
    ephem = build_ephemeris(duration)
    y0, orbit_metadata = lro_like_state(model, ephem)

    payload = {
        "schema": "robustness_test3_truth_convergence_v1",
        "formal_run": not args_cli.smoke,
        "scenario": {
            "duration_s": duration,
            "reference_degrees": list(REFERENCE_DEGREES),
            "orientation": "DE440/MOON_PA, 60 s quaternion table",
            "orbit": orbit_metadata,
            "integrator": "SciPy DOP853",
            "rtol": 1.0e-11,
            "atol_position_m": 1.0e-4,
            "atol_velocity_m_s": 1.0e-7,
            "max_step_s": 60.0,
            "coefficient_control": "single N=900 model parsed once; runtime truncation only",
        },
        "force_sets": [],
    }
    payload["force_sets"].append(
        _run_force_set("moon_sh_only", model, args, ephem, y0, duration,
                       expanded=False)
    )
    payload["force_sets"].append(
        _run_force_set("moon_sh_plus_earth_sun_srp_eclipse", model, args,
                       ephem, y0, duration, expanded=True)
    )
    suffix = "_smoke" if args_cli.smoke else ""
    out = METRICS / f"robustness_test3_truth_convergence{suffix}.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {out.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
