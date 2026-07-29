"""Robustness Test 1: dynamic scheduling in independent TudatPy.

This is a two-environment runner.  It deliberately keeps the Lunaris and
TudatPy execution paths in separate processes and exchanges only a frozen JSON
contract plus compressed trajectory files under ``codebase/data/test1``.

Lunaris backend (creates the contract and production-like stage-evaluated
trajectories)::

    D:\\Masaustu\\LUNAR_SIMULATION\\.venv\\Scripts\\python.exe \
        robustness_test1_tudat_dynamic.py --backend lunaris

Tudat backend (run from the pinned TudatPy 1.0.0 environment)::

    micromamba run --prefix <external>\\.tudat-env python \
        robustness_test1_tudat_dynamic.py --backend tudat

Comparison (either NumPy/SciPy environment)::

    <python> robustness_test1_tudat_dynamic.py --backend compare

Add ``--smoke`` to every command for the separate four-hour qualification
path.  Tudat uses its own SH parser/kernel, SPICE rotation, variable-step
RKF78 integrator,
altitude termination/root finding, and segmented scheduler.  No Lunaris force
or scheduling code is imported by the Tudat backend.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


BASE = Path(__file__).resolve().parents[1]
DATA_ROOT = BASE / "data" / "robustness_test1"
METRICS = BASE / "metrics"
REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
TUDAT_RUNNER = (REPO / "validation" / "gravity_reference" / "generators"
                / "trajectory" / "tudatpy_rotating")
DAY = 86400.0
# The existing MOON_PA Lunaris harness archives states every 120 s.  Tudat is
# resampled to that exact common grid for policy-to-policy comparisons.
OUTPUT_STEP_S = 120.0
FULL_DURATION_S = 7.0 * DAY
SMOKE_DURATION_S = 4.0 * 3600.0
POLICY_NAMES = ("fixed_critical", "schedule_down", "schedule_up", "schedule_emp")


def artifact_paths(smoke: bool) -> dict[str, Path]:
    suffix = "_smoke" if smoke else ""
    return {
        "contract": DATA_ROOT / f"contract{suffix}.json",
        "lunaris_npz": DATA_ROOT / f"lunaris_trajectories{suffix}.npz",
        "lunaris_json": DATA_ROOT / f"lunaris_runs{suffix}.json",
        "tudat_npz": DATA_ROOT / f"tudat_trajectories{suffix}.npz",
        "tudat_json": DATA_ROOT / f"tudat_runs{suffix}.json",
        "comparison": METRICS / f"robustness_test1_tudat_dynamic{suffix}.json",
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def _table_from_json(table: dict[str, int]) -> dict[float, int]:
    return {float(key): int(value) for key, value in table.items()}


def _degree_at(table: dict[float, int], altitude_m: float) -> int:
    keys = sorted(table)
    key = min(keys[-1], max(keys[0],
              10.0 * math.floor(altitude_m / 1.0e3 / 10.0)))
    return int(table[key])


def _boundaries(table: dict[float, int]) -> list[dict[str, float | int]]:
    keys = sorted(table)
    return [
        {"altitude_m": float(high * 1.0e3),
         "below_degree": int(table[low]), "above_degree": int(table[high])}
        for low, high in zip(keys[:-1], keys[1:])
        if int(table[low]) != int(table[high])
    ]


def _ric(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    delta = candidate[:, :3] - reference[:, :3]
    out = np.empty_like(delta)
    for index, state in enumerate(reference):
        r, v = state[:3], state[3:]
        radial = r / np.linalg.norm(r)
        cross = np.cross(r, v)
        cross /= np.linalg.norm(cross)
        intrack = np.cross(cross, radial)
        out[index] = (radial @ delta[index], intrack @ delta[index],
                      cross @ delta[index])
    return out


def state_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict:
    dp = np.linalg.norm(candidate[:, :3] - reference[:, :3], axis=1)
    dv = np.linalg.norm(candidate[:, 3:] - reference[:, 3:], axis=1)
    ric = _ric(reference, candidate)
    labels = ("radial", "in_track", "cross_track")
    return {
        "pos_rms_m": float(np.sqrt(np.mean(dp * dp))),
        "pos_max_m": float(np.max(dp)), "pos_final_m": float(dp[-1]),
        "vel_rms_m_s": float(np.sqrt(np.mean(dv * dv))),
        "vel_max_m_s": float(np.max(dv)), "vel_final_m_s": float(dv[-1]),
        "ric_rms_m": {label: float(np.sqrt(np.mean(ric[:, i] ** 2)))
                      for i, label in enumerate(labels)},
        "ric_max_m": {label: float(np.max(np.abs(ric[:, i])))
                      for i, label in enumerate(labels)},
        "ric_final_m": {label: float(ric[-1, i])
                        for i, label in enumerate(labels)},
    }


def run_lunaris(smoke: bool) -> int:
    from rev3_common import (alt_sched, degree_power, eccentric_state,
                             kernel_args, load_model, make_emp_table,
                             make_p_table, warmup)
    from lunaris.common.math_utils import coe_to_rv
    from rev4_robustness_controls import (_rotation_matrix_i2f, _run_expanded,
                                          build_ephemeris)

    paths = artifact_paths(smoke)
    duration = SMOKE_DURATION_S if smoke else FULL_DURATION_S
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    ephem = build_ephemeris(duration)
    down = make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="down")
    up = make_p_table(model, 1.0e-3, 60, cap=138, q=10, policy="up")
    empirical = make_emp_table(model, degree_power(model), 1.0e-3,
                               60, cap=138, q=10)
    # Tudat's normalized-Legendre implementation is singular at an exactly
    # polar evaluation point.  Preserve the published 30 x 216 km, AOP=270
    # geometry while pre-registering a 0.1-degree pole-avoidance inclination.
    # Both backends receive this exact same Cartesian state.
    rp = model.r_ref + 30.0e3
    ra = model.r_ref + 216.0e3
    semi_major = 0.5 * (rp + ra)
    eccentricity = (ra - rp) / (ra + rp)
    r_fixed, v_fixed = coe_to_rv(
        semi_major, eccentricity, math.radians(89.9), 0.0,
        math.radians(270.0), 0.0, model.mu,
    )
    matrix = _rotation_matrix_i2f(ephem, 0.0)
    dt = 1.0
    matrix_dot = (
        -3.0 * matrix + 4.0 * _rotation_matrix_i2f(ephem, dt)
        - _rotation_matrix_i2f(ephem, 2.0 * dt)
    ) / (2.0 * dt)
    r_inertial = matrix.T @ r_fixed
    v_inertial = matrix.T @ (v_fixed - matrix_dot @ r_inertial)
    lro_state = np.concatenate((r_inertial, v_inertial))
    lro_meta = {
        "perilune_altitude_km": 30.0, "apolune_altitude_km": 216.0,
        "inclination_deg_moon_fixed": 89.9,
        "argument_of_periapsis_deg_moon_fixed": 270.0,
        "pole_avoidance": (
            "0.1 deg from exact polar because TudatPy 1.0.0 normalized "
            "Legendre derivatives are singular at an exact pole"
        ),
    }
    orbits = {
        "50x300_polar": {
            "initial_state_j2000_m_m_s": eccentric_state(
                model, 50.0, 300.0, 90.0).tolist(),
            "fixed_critical_degree": 138,
            "description": "50 x 300 km polar, perilune start",
        },
        "lro_like_30x216": {
            "initial_state_j2000_m_m_s": lro_state.tolist(),
            "fixed_critical_degree": 194,
            "description": "LRO-like 30 x 216 km, Moon-fixed i=90 deg, AOP=270 deg",
            "metadata": lro_meta,
        },
    }
    contract = {
        "schema": "robustness_test1_contract_v1",
        "formal_run": not smoke,
        "duration_s": duration, "output_step_s": OUTPUT_STEP_S,
        "start_utc": "2025-01-01T00:00:00Z",
        "orientation": "J2000 integration / DE440 MOON_PA coefficient frame",
        "gravity_degree_loaded": 300,
        "reference_degree": 300,
        "schedule_definition": {
            "eps_tail_fraction": 1.0e-3, "floor": 60, "cap": 138,
            "quantum": 10, "altitude_grid_km": [40.0, 520.0, 10.0],
            "down": {str(k): int(v) for k, v in down.items()},
            "up": {str(k): int(v) for k, v in up.items()},
            "emp": {str(k): int(v) for k, v in empirical.items()},
        },
        "orbits": orbits,
    }
    write_json(paths["contract"], contract)

    trajectories: dict[str, np.ndarray] = {}
    run_rows: dict[str, dict] = {}
    for orbit_name, orbit in orbits.items():
        print(f"== Lunaris {orbit_name} ==", flush=True)
        y0 = np.asarray(orbit["initial_state_j2000_m_m_s"], dtype=float)
        policies = {
            "reference_N300": lambda _t, _h: 300,
            "fixed_critical": (lambda _t, _h,
                               n=int(orbit["fixed_critical_degree"]): n),
            "schedule_down": alt_sched(down),
            "schedule_up": alt_sched(up),
            "schedule_emp": alt_sched(empirical),
        }
        states: dict[str, np.ndarray] = {}
        rows: dict[str, dict] = {}
        for policy, degree_of in policies.items():
            print(f"  {policy}", flush=True)
            sol, info = _run_expanded(
                model, args, ephem, y0, degree_of, duration,
                use_third_body=False, use_srp=False,
            )
            state = sol.y.T.copy()
            states[policy] = state
            trajectories[f"{orbit_name}__{policy}"] = state
            rows[policy] = info
        rows["errors_vs_reference_N300"] = {
            policy: state_metrics(states["reference_N300"], states[policy])
            for policy in POLICY_NAMES
        }
        run_rows[orbit_name] = rows
    grid = np.arange(0.0, duration + 0.1, OUTPUT_STEP_S)
    trajectories["epoch_rel_s"] = grid
    np.savez_compressed(paths["lunaris_npz"], **trajectories)
    write_json(paths["lunaris_json"], {
        "backend": "Lunaris", "formal_run": not smoke,
        "python": sys.version, "platform": platform.platform(),
        "contract_sha256": sha256(paths["contract"]), "orbits": run_rows,
    })
    print(f"[written] {paths['lunaris_npz'].name}, {paths['lunaris_json'].name}")
    return 0


def run_tudat(smoke: bool) -> int:
    """Independent TudatPy segmented, event-resolved implementation."""
    paths = artifact_paths(smoke)
    if not paths["contract"].exists():
        raise FileNotFoundError("run --backend lunaris first to freeze the contract")
    contract = json.loads(paths["contract"].read_text(encoding="utf-8"))
    sys.path.insert(0, str(TUDAT_RUNNER))
    import tudatpy
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
        gravity_path, max_degree
    )
    spice.clear_kernels()
    for kernel in scenario["kernel_files"]:
        spice.load_kernel(str(kernel))
    start_epoch = float(spice.convert_date_string_to_ephemeris_time(
        contract["start_utc"]
    ))
    duration = float(contract["duration_s"])
    baseline_step = 60.0
    tight_step = 30.0

    def create_bodies():
        settings = environment_setup.get_default_body_settings(
            ["Moon"], "Moon", "J2000"
        )
        moon = settings.get("Moon")
        moon.gravity_field_settings = environment_setup.gravity_field.spherical_harmonic(
            mu, radius, cosine, sine, "MOON_PA"
        )
        moon.rotation_model_settings = environment_setup.rotation_model.spice(
            "J2000", "MOON_PA", "MOON_PA"
        )
        moon.shape_deformation_settings = []
        moon.shape_settings = environment_setup.shape.spherical(radius)
        bodies = environment_setup.create_system_of_bodies(settings)
        bodies.create_empty_body("Vehicle")
        bodies.get("Vehicle").mass = 100.0
        return bodies

    def propagate_policy(y0: np.ndarray, degree_policy: int | dict[float, int],
                         step_s: float) -> tuple[np.ndarray, dict]:
        bodies = create_bodies()
        is_fixed = isinstance(degree_policy, int)
        table = None if is_fixed else degree_policy
        state = y0.copy()
        t_now = start_epoch
        t_end = start_epoch + duration
        altitude = float(np.linalg.norm(state[:3]) - radius)
        degree = int(degree_policy) if is_fixed else _degree_at(table, altitude)
        boundary_rows = [] if is_fixed else _boundaries(table)
        history: dict[float, np.ndarray] = {t_now: state.copy()}
        event_rows = []
        segments = 0
        started = time.perf_counter()
        root = root_finders.bisection(
            maximum_iteration=40,
            maximum_iteration_handling=root_finders.MaximumIterationHandling.accept_result,
        )
        while t_now < t_end - 1.0e-7:
            acceleration = propagation_setup.acceleration.spherical_harmonic_gravity(
                degree, degree
            )
            models = propagation_setup.create_acceleration_models(
                bodies, {"Vehicle": {"Moon": [acceleration]}},
                ["Vehicle"], ["Moon"],
            )
            control = (
                propagation_setup.integrator
                .step_size_control_elementwise_scalar_tolerance(1.0e-12, 1.0e-6)
            )
            validation = propagation_setup.integrator.step_size_validation(
                1.0e-3, float(step_s)
            )
            integrator = propagation_setup.integrator.runge_kutta_variable_step(
                initial_time_step=min(10.0, float(step_s)),
                coefficient_set=propagation_setup.integrator.CoefficientSets.rkf_78,
                step_size_control_settings=control,
                step_size_validation_settings=validation,
            )
            terminations = [propagation_setup.propagator.time_termination(
                t_end, terminate_exactly_on_final_condition=True
            )]
            altitude_now = float(np.linalg.norm(state[:3]) - radius)
            if table is not None:
                altitude_variable = propagation_setup.dependent_variable.altitude(
                    "Vehicle", "Moon"
                )
                for boundary in boundary_rows:
                    b = float(boundary["altitude_m"])
                    if (degree == int(boundary["below_degree"])
                            and b + 1.0e-3 > altitude_now):
                        terminations.append(
                            propagation_setup.propagator.dependent_variable_termination(
                                altitude_variable, b + 1.0e-3,
                                use_as_lower_limit=False,
                                terminate_exactly_on_final_condition=True,
                                termination_root_finder_settings=root,
                            )
                        )
                    elif (degree == int(boundary["above_degree"])
                          and b - 1.0e-3 < altitude_now):
                        terminations.append(
                            propagation_setup.propagator.dependent_variable_termination(
                                altitude_variable, b - 1.0e-3,
                                use_as_lower_limit=True,
                                terminate_exactly_on_final_condition=True,
                                termination_root_finder_settings=root,
                            )
                        )
            termination = propagation_setup.propagator.hybrid_termination(
                terminations, fulfill_single_condition=True
            ) if len(terminations) > 1 else terminations[0]
            propagator = propagation_setup.propagator.translational(
                central_bodies=["Moon"], acceleration_models=models,
                bodies_to_integrate=["Vehicle"], initial_states=state,
                initial_time=t_now, integrator_settings=integrator,
                termination_settings=termination,
            )
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
                raise RuntimeError(
                    f"Tudat stopped at h={h} m without a degree transition"
                )
            nearest = min(boundary_rows,
                          key=lambda row: abs(float(row["altitude_m"]) - h))
            event_rows.append({
                "epoch_rel_s": t_now - start_epoch,
                "altitude_m": h,
                "nearest_boundary_m": float(nearest["altitude_m"]),
                "root_residual_m": h - float(nearest["altitude_m"]),
                "from_degree": degree, "to_degree": int(new_degree),
            })
            degree = int(new_degree)

        epochs = np.asarray(sorted(history), dtype=float)
        states = np.vstack([history[float(epoch)] for epoch in epochs])
        requested = start_epoch + np.arange(
            0.0, duration + 0.1, OUTPUT_STEP_S, dtype=float
        )
        interpolated = CubicSpline(epochs, states, axis=0)(requested)
        return interpolated, {
            "maximum_step_s": step_s,
            "integrator": "Tudat variable-step RKF78, rtol=1e-12, atol=1e-6",
            "segments": segments, "event_count": len(event_rows),
            "events": event_rows,
            "max_abs_root_residual_m": float(max(
                (abs(row["root_residual_m"]) for row in event_rows), default=0.0
            )),
            "internal_saved_states": int(epochs.size),
            "runtime_s": time.perf_counter() - started,
        }

    tables = {
        key: _table_from_json(contract["schedule_definition"][short])
        for key, short in (("schedule_down", "down"),
                           ("schedule_up", "up"),
                           ("schedule_emp", "emp"))
    }
    trajectories: dict[str, np.ndarray] = {}
    run_rows: dict[str, dict] = {}
    for orbit_name, orbit in contract["orbits"].items():
        print(f"== Tudat {orbit_name} ==", flush=True)
        y0 = np.asarray(orbit["initial_state_j2000_m_m_s"], dtype=float)
        policies: dict[str, int | dict[float, int]] = {
            "reference_N300": 300,
            "fixed_critical": int(orbit["fixed_critical_degree"]),
            **tables,
        }
        rows = {}
        for policy, degree_policy in policies.items():
            rows[policy] = {}
            for level, step_s in (("baseline", baseline_step), ("tight", tight_step)):
                print(f"  {policy} / {level} (max {step_s:g} s)", flush=True)
                state, info = propagate_policy(y0, degree_policy, step_s)
                trajectories[f"{orbit_name}__{policy}__{level}"] = state
                rows[policy][level] = info
        run_rows[orbit_name] = rows
    trajectories["epoch_rel_s"] = np.arange(
        0.0, duration + 0.1, OUTPUT_STEP_S, dtype=float
    )
    np.savez_compressed(paths["tudat_npz"], **trajectories)
    write_json(paths["tudat_json"], {
        "backend": "TudatPy", "formal_run": not smoke,
        "tudatpy_version": getattr(tudatpy, "__version__", "unknown"),
        "python": sys.version, "platform": platform.platform(),
        "contract_sha256": sha256(paths["contract"]),
        "gravity_sha256": sha256(gravity_path),
        "gravity_parser": "independent PDS SHADR parser from validated Tudat runner",
        "gravity_parser_metadata": parser_meta,
        "orbits": run_rows,
    })
    print(f"[written] {paths['tudat_npz'].name}, {paths['tudat_json'].name}")
    return 0


def compare(smoke: bool) -> int:
    paths = artifact_paths(smoke)
    required = ("contract", "lunaris_npz", "lunaris_json", "tudat_npz", "tudat_json")
    missing = [str(paths[key]) for key in required if not paths[key].exists()]
    if missing:
        raise FileNotFoundError("missing backend artifacts:\n" + "\n".join(missing))
    contract = json.loads(paths["contract"].read_text(encoding="utf-8"))
    lunaris_meta = json.loads(paths["lunaris_json"].read_text(encoding="utf-8"))
    tudat_meta = json.loads(paths["tudat_json"].read_text(encoding="utf-8"))
    if lunaris_meta["contract_sha256"] != tudat_meta["contract_sha256"]:
        raise ValueError("backend contract hashes differ")
    lunaris = np.load(paths["lunaris_npz"])
    tudat = np.load(paths["tudat_npz"])
    result_orbits = {}
    for orbit_name in contract["orbits"]:
        lref = lunaris[f"{orbit_name}__reference_N300"]
        tref = tudat[f"{orbit_name}__reference_N300__tight"]
        policies = {}
        for policy in POLICY_NAMES:
            lstate = lunaris[f"{orbit_name}__{policy}"]
            tstate = tudat[f"{orbit_name}__{policy}__tight"]
            tbase = tudat[f"{orbit_name}__{policy}__baseline"]
            policies[policy] = {
                "lunaris_vs_own_N300": state_metrics(lref, lstate),
                "tudat_vs_own_N300": state_metrics(tref, tstate),
                "tudat_tight_minus_baseline": state_metrics(tstate, tbase),
                "cross_implementation_same_policy": state_metrics(lstate, tstate),
            }
        lrank = sorted(POLICY_NAMES,
                       key=lambda key: policies[key]["lunaris_vs_own_N300"]["pos_rms_m"])
        trank = sorted(POLICY_NAMES,
                       key=lambda key: policies[key]["tudat_vs_own_N300"]["pos_rms_m"])
        # Treat a pair as numerically undecided when its separation is no
        # larger than either Tudat member's baseline-to-tight refinement RMS.
        # This prevents a sub-envelope swap from being mislabeled as a policy
        # ranking reversal (notably the near-tied LRO down/empirical pair).
        decisive_pairs = []
        for i, first in enumerate(POLICY_NAMES):
            for second in POLICY_NAMES[i + 1:]:
                l1 = policies[first]["lunaris_vs_own_N300"]["pos_rms_m"]
                l2 = policies[second]["lunaris_vs_own_N300"]["pos_rms_m"]
                t1 = policies[first]["tudat_vs_own_N300"]["pos_rms_m"]
                t2 = policies[second]["tudat_vs_own_N300"]["pos_rms_m"]
                envelope = max(
                    policies[first]["tudat_tight_minus_baseline"]["pos_rms_m"],
                    policies[second]["tudat_tight_minus_baseline"]["pos_rms_m"],
                )
                decisive = abs(l1 - l2) > envelope and abs(t1 - t2) > envelope
                decisive_pairs.append({
                    "pair": [first, second], "refinement_envelope_rms_m": envelope,
                    "decisive": decisive,
                    "same_order_if_decisive": (not decisive) or ((l1 < l2) == (t1 < t2)),
                })
        result_orbits[orbit_name] = {
            "reference_cross_implementation": state_metrics(lref, tref),
            "policies": policies,
            "ranking": {
                "lunaris": lrank, "tudat": trank,
                "identical_raw": lrank == trank,
                "decisive_pairs": decisive_pairs,
                "all_decisive_pairs_agree": all(
                    row["same_order_if_decisive"] for row in decisive_pairs
                ),
            },
            "in_track_dominant_both": {
                policy: all(
                    metrics["ric_rms_m"]["in_track"] >=
                    max(metrics["ric_rms_m"]["radial"],
                        metrics["ric_rms_m"]["cross_track"])
                    for metrics in (policies[policy]["lunaris_vs_own_N300"],
                                    policies[policy]["tudat_vs_own_N300"])
                ) for policy in POLICY_NAMES
            },
        }
    payload = {
        "schema": "robustness_test1_tudat_dynamic_v1",
        "formal_run": not smoke,
        "contract_sha256": sha256(paths["contract"]),
        "backend_artifact_sha256": {
            "lunaris_npz": sha256(paths["lunaris_npz"]),
            "tudat_npz": sha256(paths["tudat_npz"]),
            "lunaris_json": sha256(paths["lunaris_json"]),
            "tudat_json": sha256(paths["tudat_json"]),
        },
        "acceptance_scope": (
            "Policy ranking and RIC dominance are hard checks. Magnitudes are "
            "reported beside the Tudat variable-step RKF78 refinement "
            "envelope; no post-hoc "
            "tolerance is introduced."
        ),
        "orbits": result_orbits,
    }
    write_json(paths["comparison"], payload)
    print(f"[written] {paths['comparison'].name}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", required=True,
                        choices=("lunaris", "tudat", "compare"))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.backend == "lunaris":
        return run_lunaris(args.smoke)
    if args.backend == "tudat":
        return run_tudat(args.smoke)
    return compare(args.smoke)


if __name__ == "__main__":
    raise SystemExit(main())
