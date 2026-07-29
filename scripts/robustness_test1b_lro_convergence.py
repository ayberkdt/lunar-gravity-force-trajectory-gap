"""Robustness Test 1b: LRO-like Tudat step-convergence with Richardson.

The LRO-like refinement envelopes in the main Test 1 comparison (Tudat
30-to-60 s step difference) are an appreciable fraction of the policy errors,
so the bare statement ``magnitudes agree within the refinement envelope'' is
not, by itself, a high-precision validation at that geometry. This script
resolves it: the independent TudatPy segmented propagation is run for the
LRO-like orbit at three maximum steps (60, 30, 15 s) for the N=300 reference
and each policy, the formal observed order

    p_obs = log2( ||y_h - y_{h/2}|| / ||y_{h/2} - y_{h/4}|| )

is computed per trajectory, and a Richardson-extrapolated trajectory

    y_R = y_{h/4} + (y_{h/4} - y_{h/2}) / (2^{p} - 1)

is formed. Policy errors are then reported against the Richardson-extrapolated
reference, and the residual numerical envelope ||y_{h/4} - y_R|| is shown to be
a small fraction of the policy error, so the LRO comparison is high precision.

Reuses the frozen Test 1 contract and the same independent Tudat kernel/parser
as robustness_test1_tudat_dynamic.py.  Launch under micromamba (see runbook):
  micromamba run --root-prefix ... --prefix ...\\.tudat-env python \
      robustness_test1b_lro_convergence.py
"""

from __future__ import annotations

import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np

from robustness_test1_tudat_dynamic import (
    BASE, DATA_ROOT, METRICS, REPO, TUDAT_RUNNER, OUTPUT_STEP_S,
    _table_from_json, _degree_at, _boundaries, state_metrics, sha256)

STEPS_S = (60.0, 30.0, 15.0)     # h, h/2, h/4
ORBIT_NAME = "lro_like_30x216"
DURATION_S = 7.0 * 86400.0


def build_propagator(contract):
    """Return a propagate_policy(y0, degree_policy, step_s) closure that
    replicates the independent Tudat segmented, event-resolved propagation."""
    sys.path.insert(0, str(TUDAT_RUNNER))
    from scipy.interpolate import CubicSpline
    from tudatpy.dynamics import environment_setup, propagation_setup, simulator
    from tudatpy.interface import spice
    from tudatpy.math import root_finders
    from run_tudat import load_jggrx_independently
    from validation_common import load_scenario, sha256_file

    scenario = load_scenario()
    gravity_path = Path(scenario["gravity_file"])
    if sha256_file(gravity_path) != scenario["expected_gravity_sha256"]:
        raise ValueError("gravity hash mismatch")
    max_degree = int(contract["gravity_degree_loaded"])
    mu, radius, cosine, sine, parser_meta = load_jggrx_independently(
        gravity_path, max_degree)
    spice.clear_kernels()
    for kernel in scenario["kernel_files"]:
        spice.load_kernel(str(kernel))
    start_epoch = float(spice.convert_date_string_to_ephemeris_time(
        contract["start_utc"]))
    duration = float(contract["duration_s"])

    def create_bodies():
        settings = environment_setup.get_default_body_settings(
            ["Moon"], "Moon", "J2000")
        moon = settings.get("Moon")
        moon.gravity_field_settings = environment_setup.gravity_field.spherical_harmonic(
            mu, radius, cosine, sine, "MOON_PA")
        moon.rotation_model_settings = environment_setup.rotation_model.spice(
            "J2000", "MOON_PA", "MOON_PA")
        moon.shape_deformation_settings = []
        moon.shape_settings = environment_setup.shape.spherical(radius)
        bodies = environment_setup.create_system_of_bodies(settings)
        bodies.create_empty_body("Vehicle")
        bodies.get("Vehicle").mass = 100.0
        return bodies

    def propagate_policy(y0, degree_policy, step_s):
        bodies = create_bodies()
        is_fixed = isinstance(degree_policy, int)
        table = None if is_fixed else degree_policy
        state = y0.copy()
        t_now = start_epoch
        t_end = start_epoch + duration
        altitude = float(np.linalg.norm(state[:3]) - radius)
        degree = int(degree_policy) if is_fixed else _degree_at(table, altitude)
        boundary_rows = [] if is_fixed else _boundaries(table)
        history = {t_now: state.copy()}
        event_rows = []
        segments = 0
        started = time.perf_counter()
        root = root_finders.bisection(
            maximum_iteration=40,
            maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result)
        while t_now < t_end - 1.0e-7:
            acceleration = propagation_setup.acceleration.spherical_harmonic_gravity(
                degree, degree)
            models = propagation_setup.create_acceleration_models(
                bodies, {"Vehicle": {"Moon": [acceleration]}},
                ["Vehicle"], ["Moon"])
            control = (propagation_setup.integrator
                       .step_size_control_elementwise_scalar_tolerance(1.0e-12, 1.0e-6))
            validation = propagation_setup.integrator.step_size_validation(
                1.0e-3, float(step_s))
            integrator = propagation_setup.integrator.runge_kutta_variable_step(
                initial_time_step=min(10.0, float(step_s)),
                coefficient_set=propagation_setup.integrator.CoefficientSets.rkf_78,
                step_size_control_settings=control,
                step_size_validation_settings=validation)
            terminations = [propagation_setup.propagator.time_termination(
                t_end, terminate_exactly_on_final_condition=True)]
            altitude_now = float(np.linalg.norm(state[:3]) - radius)
            if table is not None:
                altitude_variable = propagation_setup.dependent_variable.altitude(
                    "Vehicle", "Moon")
                for boundary in boundary_rows:
                    b = float(boundary["altitude_m"])
                    if (degree == int(boundary["below_degree"])
                            and b + 1.0e-3 > altitude_now):
                        terminations.append(
                            propagation_setup.propagator.dependent_variable_termination(
                                altitude_variable, b + 1.0e-3, use_as_lower_limit=False,
                                terminate_exactly_on_final_condition=True,
                                termination_root_finder_settings=root))
                    elif (degree == int(boundary["above_degree"])
                          and b - 1.0e-3 < altitude_now):
                        terminations.append(
                            propagation_setup.propagator.dependent_variable_termination(
                                altitude_variable, b - 1.0e-3, use_as_lower_limit=True,
                                terminate_exactly_on_final_condition=True,
                                termination_root_finder_settings=root))
            termination = propagation_setup.propagator.hybrid_termination(
                terminations, fulfill_single_condition=True
            ) if len(terminations) > 1 else terminations[0]
            propagator = propagation_setup.propagator.translational(
                central_bodies=["Moon"], acceleration_models=models,
                bodies_to_integrate=["Vehicle"], initial_states=state,
                initial_time=t_now, integrator_settings=integrator,
                termination_settings=termination)
            propagator.processing_settings.results_save_frequency_in_steps = 1
            simulation = simulator.create_dynamics_simulator(bodies, propagator)
            segment_history = simulation.propagation_results.state_history
            for epoch, value in segment_history.items():
                history[float(epoch)] = np.asarray(value, dtype=float).copy()
            epochs = sorted(segment_history)
            new_t = float(epochs[-1])
            state = np.asarray(segment_history[epochs[-1]], dtype=float).copy()
            segments += 1
            if new_t <= t_now + 1.0e-8:
                raise RuntimeError("Tudat event segment made no progress")
            t_now = new_t
            if t_now >= t_end - 1.0e-7 or table is None:
                continue
            h = float(np.linalg.norm(state[:3]) - radius)
            radial_rate = float(state[:3] @ state[3:] / np.linalg.norm(state[:3]))
            probe_h = h + math.copysign(1.0, radial_rate if radial_rate else 1.0)
            new_degree = _degree_at(table, probe_h)
            if new_degree == degree:
                raise RuntimeError(f"Tudat stopped at h={h} m without transition")
            nearest = min(boundary_rows,
                          key=lambda row: abs(float(row["altitude_m"]) - h))
            event_rows.append({
                "epoch_rel_s": t_now - start_epoch, "altitude_m": h,
                "from_degree": degree, "to_degree": int(new_degree)})
            degree = int(new_degree)

        epochs = np.asarray(sorted(history), dtype=float)
        states = np.vstack([history[float(epoch)] for epoch in epochs])
        requested = start_epoch + np.arange(0.0, duration + 0.1, OUTPUT_STEP_S)
        interpolated = CubicSpline(epochs, states, axis=0)(requested)
        return interpolated, {
            "maximum_step_s": step_s, "segments": segments,
            "event_count": len(event_rows),
            "runtime_s": time.perf_counter() - started}

    return propagate_policy, gravity_path, parser_meta


def _rms(a):
    return float(np.sqrt(np.mean(np.sum(a * a, axis=1))))


def convergence_stats(y_h, y_h2, y_h4):
    """p_obs and Richardson extrapolation from three step sizes (h, h/2, h/4).
    States are (N, 6); position block used for the order estimate."""
    d1 = _rms(y_h[:, :3] - y_h2[:, :3])
    d2 = _rms(y_h2[:, :3] - y_h4[:, :3])
    p_obs = math.log2(d1 / d2) if d2 > 0 and d1 > 0 else float("nan")
    p_eff = p_obs if (math.isfinite(p_obs) and p_obs > 0.5) else 7.0
    y_rich = y_h4 + (y_h4 - y_h2) / (2.0 ** p_eff - 1.0)
    return {
        "successive_diff_h_to_h2_rms_m": d1,
        "successive_diff_h2_to_h4_rms_m": d2,
        "p_obs_position": p_obs,
        "p_used_for_richardson": p_eff,
        "richardson_minus_finest_rms_m": _rms(y_rich[:, :3] - y_h4[:, :3]),
    }, y_rich


def main() -> int:
    contract_path = DATA_ROOT / "contract.json"
    if not contract_path.exists():
        raise FileNotFoundError("run robustness_test1 --backend lunaris first")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    orbit = contract["orbits"][ORBIT_NAME]
    y0 = np.asarray(orbit["initial_state_j2000_m_m_s"], dtype=float)

    tables = {key: _table_from_json(contract["schedule_definition"][short])
              for key, short in (("schedule_down", "down"),
                                 ("schedule_up", "up"),
                                 ("schedule_emp", "emp"))}
    policies = {"reference_N300": 300,
                "fixed_critical": int(orbit["fixed_critical_degree"]),
                **tables}

    propagate_policy, gravity_path, parser_meta = build_propagator(contract)

    print(f"== Test 1b LRO convergence: steps {STEPS_S} s ==", flush=True)
    runs = {name: {} for name in policies}
    infos = {name: {} for name in policies}
    for name, degree_policy in policies.items():
        for step_s in STEPS_S:
            print(f"  {name} / max {step_s:g} s", flush=True)
            state, info = propagate_policy(y0, degree_policy, step_s)
            runs[name][step_s] = state
            infos[name][step_s] = info

    # Richardson per trajectory
    rich = {}
    conv = {}
    for name in policies:
        c, y_r = convergence_stats(runs[name][STEPS_S[0]],
                                   runs[name][STEPS_S[1]],
                                   runs[name][STEPS_S[2]])
        conv[name] = c
        rich[name] = y_r

    ref_rich = rich["reference_N300"]
    ref_finest = runs["reference_N300"][STEPS_S[2]]
    results = {}
    for name in policies:
        if name == "reference_N300":
            continue
        finest = runs[name][STEPS_S[2]]
        results[name] = {
            "convergence": conv[name],
            "error_vs_richardson_reference": state_metrics(ref_rich, rich[name]),
            "error_vs_finest_reference_15s": state_metrics(ref_finest, finest),
            "residual_numerical_envelope_rms_m":
                conv[name]["richardson_minus_finest_rms_m"],
            "event_counts": {str(s): infos[name][s]["event_count"] for s in STEPS_S},
        }

    payload = {
        "schema": "robustness_test1b_lro_convergence_v1",
        "orbit": ORBIT_NAME, "steps_s": list(STEPS_S),
        "python": sys.version, "platform": platform.platform(),
        "gravity_sha256": sha256(gravity_path),
        "reference_convergence": conv["reference_N300"],
        "policies": results,
    }
    out = METRICS / "robustness_test1b_lro_convergence.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"[written] {out.name}", flush=True)
    for name, r in results.items():
        rr = r["error_vs_richardson_reference"]["pos_rms_m"]
        env = r["residual_numerical_envelope_rms_m"]
        print(f"  {name:15s} err_vs_richardson={rr:8.2f} m  "
              f"num_envelope={env:6.3f} m  p_obs="
              f"{r['convergence']['p_obs_position']:.2f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
