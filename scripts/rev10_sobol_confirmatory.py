"""R10 independent scrambled-Sobol confirmatory population.

This driver freezes the pre-analysis protocol and two 64-point designs, runs
short smoke checks, and executes the resumable seven-day Sobol A baseline.
Sobol B is design-only in this revision.

Commands
--------
    python rev10_sobol_confirmatory.py design
    python rev10_sobol_confirmatory.py smoke
    python rev10_sobol_confirmatory.py baseline
    python rev10_sobol_confirmatory.py status
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import qmc

from rev3_common import (
    DAY,
    REPO,
    Rhs,
    InstrumentedDOP853,
    commit_sha,
    degree_power,
    err_stats,
    kernel_args,
    load_model,
    warmup,
    working_tree_clean,
)
from rev7_doe_screening import (
    ATOL,
    CAP,
    EPS,
    FLOOR,
    H_GRID_KM,
    OUT_STEP,
    Q,
    RTOL,
    alt_sched,
    emp_nmin_exact,
    emp_table,
    initial_state,
    kaula_table,
)

from lunaris.common.lunar_data import resolve_lunar_gravity_path


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PROTOCOL_DIR = ROOT / "experiments" / "protocols"
CASE_ROOT = METRICS / "r10_cases"
RAW_ROOT = METRICS / "r10_raw"

PROTOCOL_PATH = PROTOCOL_DIR / "sobol_confirmatory_protocol.json"
DESIGN_A_PATH = METRICS / "r10_sobolA_design.json"
DESIGN_B_PATH = METRICS / "r10_sobolB_design_frozen.json"
BASELINE_PATH = METRICS / "r10_sobolA_baseline.json"
ACTIVE_PATH = METRICS / "r10_sobolA_baseline_active.json"
SMOKE_PATH = METRICS / "r10_sobolA_smoke.json"

SOBOL_SEED_A = 20260723
SOBOL_SEED_B = 20260724
N_POINTS = 64
DURATION_S = 7.0 * DAY
OUTPUT_STEP_S = 120.0
SMOKE_DURATION_S = 2.0 * 3600.0
SURFACE_ROOT_XTOL_S = 1.0e-8
ROUNDTRIP_LIMIT_DEG = 1.0e-10
DESIGN_DUPLICATE_TOL = 1.0e-14


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def object_hash(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def wrapped_delta_deg(actual: float, requested: float) -> float:
    return (actual - requested + 180.0) % 360.0 - 180.0


def experiment_script_sha() -> str:
    return file_hash(Path(__file__).resolve())


def provenance() -> dict:
    gravity_path = Path(resolve_lunar_gravity_path(None)).resolve()
    kernel_path = (REPO / "src" / "lunaris" / "physics" /
                   "spherical_harmonics.py").resolve()
    return {
        "lunaris_repo": str(REPO.resolve()),
        "lunaris_commit": commit_sha(),
        "lunaris_working_tree_clean": working_tree_clean(),
        "release_tag": "paper-truncation-v1.0",
        "kernel_path": str(kernel_path),
        "kernel_sha256": file_hash(kernel_path),
        "gravity_path": str(gravity_path),
        "gravity_sha256": file_hash(gravity_path),
        "experiment_script": str(Path(__file__).resolve()),
        "experiment_script_sha256": experiment_script_sha(),
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
    }


def build_protocol() -> dict:
    payload = {
        "schema": "r10_sobol_confirmatory_protocol_v1",
        "created_utc": utc_now(),
        "status": "frozen_before_formal_trajectory_inspection",
        "purpose": (
            "independent expanded-domain coverage confirmation of the "
            "24-orbit exploratory population"
        ),
        "sobol": {
            "dimension": 5,
            "sample_count": N_POINTS,
            "generator": "scipy.stats.qmc.Sobol(scramble=True).random_base2(m=6)",
            "seed_A": SOBOL_SEED_A,
            "seed_B": SOBOL_SEED_B,
            "seed_B_role": "frozen_unpropagated_future_replicate",
            "no_rejection": True,
            "no_replacement": True,
            "no_candidate_selection": True,
            "duplicate_tolerance_u_maxnorm": DESIGN_DUPLICATE_TOL,
        },
        "transformations": {
            "hp_km": "30 + 120*u0",
            "ha_min_km": "max(180, hp_km + 30)",
            "ha_km": "ha_min_km + u1*(600-ha_min_km)",
            "incl_deg": "180*u2",
            "argp_deg": "360*u3",
            "perilune_lon_deg_bodyfixed_t0": "360*u4",
            "raan_deg": (
                "lambda_p - atan2(cos(i)*sin(omega), cos(omega)) mod 360"
            ),
            "nu0_deg": 0.0,
            "roundtrip_limit_deg": ROUNDTRIP_LIMIT_DEG,
        },
        "force_contract": {
            "gravity": "JGGRX_1800F",
            "dynamics": "gravity-only",
            "orientation": "uniform lunar sidereal rotation about polar axis",
            "epoch": "t=0 body-fixed and inertial frames aligned",
            "truth_rule": "N=600 when hp<50 km, otherwise N=300",
        },
        "population_baseline": {
            "duration_s": DURATION_S,
            "output_step_s": OUTPUT_STEP_S,
            "integrator": "DOP853 production kernel with direct trial-step instrumentation",
            "rtol": RTOL,
            "atol": ATOL,
            "atol_kind": "scalar",
            "max_step_s": None,
            "surface_event": "norm(r)-r_ref, terminal downward crossing",
        },
        "policies": {
            "primary": "empirical lookup schedule",
            "primary_comparators": ["fixed_work", "fixed_critical"],
            "secondary": ["schedule_up", "schedule_down"],
            "empirical_lookup": {
                "altitude_grid_km": [float(H_GRID_KM[0]),
                                     float(H_GRID_KM[-1]), 10.0],
                "tail_fraction": EPS,
                "floor": FLOOR,
                "cap": CAP,
                "degree_quantum": Q,
                "degree_quantization": "downward",
                "altitude_bin_convention": "floor",
                "evaluation": "at every gravity RHS call including stages",
            },
            "fixed_critical": (
                "unquantized empirical Nmin at own perilune, capped at 250"
            ),
            "fixed_work": (
                "round(sqrt(mean N_call^2)) from empirical schedule actual RHS calls"
            ),
            "schedule_up_down": (
                "existing Kaula p=2 upward/downward quantized sensitivity families"
            ),
        },
        "primary_metric": "full-common-grid Cartesian position RMS against truth",
        "primary_ratios": {
            "rho_work": "E_fixed_work/E_empirical",
            "rho_crit": "E_fixed_critical/E_empirical",
            "interpretation": "rho>1 favors the empirical schedule",
        },
        "failure_contract": {
            "truth_impact": "exclude from full-seven-day accuracy ratios; retain in survival counts",
            "policy_impact": "noncompetitive; report impact and last-common-arc error",
            "numerical_failure": "separate from impact; preserve attempt; never replace design point",
        },
        "planned_convergence": {
            "vector_tight": {
                "rtol": 1.0e-12,
                "atol_position_m": 1.0e-5,
                "atol_velocity_m_s": 1.0e-8,
                "max_step_s": 60.0,
            },
            "vector_tighter": {
                "rtol": 1.0e-13,
                "atol_position_m": 1.0e-6,
                "atol_velocity_m_s": 1.0e-9,
                "max_step_s": 60.0,
            },
            "selection": (
                "frozen union A-F from the implementation plan; minimum 12; "
                "no truncation of raw-win or close cases"
            ),
            "resolution": (
                "abs(EA-EB) > (Eself,A+Eself,truth) + "
                "(Eself,B+Eself,truth)"
            ),
        },
        "planned_truth_audit": {
            "cases": (
                "four lowest-perilune truth survivors plus sub-50-km raw wins"
            ),
            "degrees": [600, 900],
            "acceptance": (
                "E600_900 < min(5 m, 0.05*smallest interpreted policy RMS vs N900)"
            ),
        },
        "planned_aggregates": [
            "survival and impact counts",
            "rho_work and rho_crit median, p10, p90",
            "raw and resolved outcome counts",
            "median gravity-time saving versus fixed critical",
            "median in-track fraction",
            "declared descriptive regime summaries",
        ],
        "timing": {
            "wall_and_kernel_clock": "time.perf_counter_ns",
            "cpu_clock": "time.process_time_ns",
            "formal_baseline_requires_no_other_python_process": True,
            "parallel_trajectory_timing_allowed_in_manuscript": False,
        },
        "provenance": provenance(),
        "doi_created": False,
        "archival_deposit_created": False,
    }
    payload["protocol_sha256"] = object_hash(payload)
    return payload


def protocol_payload() -> dict:
    if not PROTOCOL_PATH.exists():
        raise FileNotFoundError("run the design command first")
    payload = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    claimed = payload.pop("protocol_sha256")
    actual = object_hash(payload)
    payload["protocol_sha256"] = claimed
    if actual != claimed:
        raise RuntimeError("protocol hash mismatch")
    if payload["provenance"]["experiment_script_sha256"] != experiment_script_sha():
        raise RuntimeError("driver changed after protocol freeze")
    return payload


def orbit_from_u(index: int, row: np.ndarray, family: str, model) -> dict:
    hp_km = 30.0 + 120.0 * float(row[0])
    ha_min_km = max(180.0, hp_km + 30.0)
    ha_km = ha_min_km + float(row[1]) * (600.0 - ha_min_km)
    incl_deg = 180.0 * float(row[2])
    argp_deg = 360.0 * float(row[3])
    requested_lon = 360.0 * float(row[4])
    inclination = math.radians(incl_deg)
    argument = math.radians(argp_deg)
    longitude_offset = math.degrees(
        math.atan2(math.cos(inclination) * math.sin(argument),
                   math.cos(argument))
    )
    raan_deg = (requested_lon - longitude_offset) % 360.0
    orbit = {
        "name": f"{family}_{index:03d}",
        "family": family,
        "sobol_index": index,
        "u": [float(x) for x in row],
        "hp_km": hp_km,
        "ha_km": ha_km,
        "incl_deg": incl_deg,
        "argp_deg": argp_deg,
        "requested_perilune_lon_deg_bodyfixed_t0": requested_lon,
        "raan_deg": raan_deg,
        "nu0_deg": 0.0,
    }
    rp = model.r_ref + hp_km * 1000.0
    ra = model.r_ref + ha_km * 1000.0
    orbit["semimajor_axis_m"] = 0.5 * (rp + ra)
    orbit["eccentricity"] = (ra - rp) / (ra + rp)
    state = initial_state(model, orbit)
    actual_lon = math.degrees(math.atan2(state[1], state[0])) % 360.0
    actual_lat = math.degrees(math.asin(state[2] / np.linalg.norm(state[:3])))
    orbit["initial_state_si"] = [float(x) for x in state]
    orbit["reconstructed_perilune_lon_deg_bodyfixed_t0"] = actual_lon
    orbit["reconstructed_perilune_lat_deg_bodyfixed_t0"] = actual_lat
    orbit["longitude_roundtrip_error_deg"] = wrapped_delta_deg(
        actual_lon, requested_lon
    )
    orbit["truth_degree"] = 600 if hp_km < 50.0 else 300
    return orbit


def inclination_regime(value: float) -> str:
    if value < 60.0:
        return "prograde"
    if value <= 120.0:
        return "high_inclination"
    return "retrograde"


def make_design(seed: int, family: str, role: str, model,
                protocol_sha: str) -> dict:
    points = qmc.Sobol(d=5, scramble=True, seed=seed).random_base2(m=6)
    orbits = [orbit_from_u(index, row, family, model)
              for index, row in enumerate(points)]
    counts = {name: 0 for name in
              ("prograde", "high_inclination", "retrograde")}
    for orbit in orbits:
        counts[inclination_regime(orbit["incl_deg"])] += 1
    payload = {
        "schema": "r10_sobol_design_v1",
        "family": family,
        "role": role,
        "seed": seed,
        "sample_count": len(orbits),
        "dimension": 5,
        "generator": "scrambled Sobol random_base2(m=6)",
        "protocol_sha256": protocol_sha,
        "n_rejected": 0,
        "n_replaced": 0,
        "optimized": False,
        "propagation_status": (
            "planned_confirmatory_baseline" if family == "sobolA"
            else "frozen_unpropagated"
        ),
        "realized_inclination_range_deg": [
            min(o["incl_deg"] for o in orbits),
            max(o["incl_deg"] for o in orbits),
        ],
        "inclination_regime_counts": counts,
        "orbits": orbits,
    }
    payload["design_sha256"] = object_hash(payload)
    return payload


def validate_designs(design_a: dict, design_b: dict) -> dict:
    failures: list[str] = []
    for design in (design_a, design_b):
        orbits = design["orbits"]
        if len(orbits) != N_POINTS:
            failures.append(f"{design['family']}: expected 64 points")
        matrix = np.asarray([o["u"] for o in orbits], dtype=float)
        if len(np.unique(matrix, axis=0)) != N_POINTS:
            failures.append(f"{design['family']}: duplicate point")
        for orbit in orbits:
            if not (30.0 <= orbit["hp_km"] <= 150.0):
                failures.append(f"{orbit['name']}: hp out of domain")
            if not (180.0 <= orbit["ha_km"] <= 600.0):
                failures.append(f"{orbit['name']}: ha out of domain")
            if orbit["ha_km"] - orbit["hp_km"] < 30.0 - 1.0e-12:
                failures.append(f"{orbit['name']}: insufficient altitude gap")
            if abs(orbit["longitude_roundtrip_error_deg"]) > ROUNDTRIP_LIMIT_DEG:
                failures.append(f"{orbit['name']}: longitude roundtrip failed")
            if not np.all(np.isfinite(orbit["initial_state_si"])):
                failures.append(f"{orbit['name']}: nonfinite state")
        if min(design["inclination_regime_counts"].values()) < 2:
            failures.append(f"{design['family']}: inclination regime underfilled")
    a = np.asarray([o["u"] for o in design_a["orbits"]])
    b = np.asarray([o["u"] for o in design_b["orbits"]])
    minimum_cross_distance = float(np.min(np.max(np.abs(a[:, None, :] -
                                                        b[None, :, :]), axis=2)))
    if minimum_cross_distance <= DESIGN_DUPLICATE_TOL:
        failures.append("A-B duplicate within tolerance")
    return {
        "passes": not failures,
        "failures": failures,
        "minimum_A_B_u_maxnorm": minimum_cross_distance,
        "duplicate_tolerance": DESIGN_DUPLICATE_TOL,
    }


def command_design(force: bool) -> int:
    if (CASE_ROOT.exists() or BASELINE_PATH.exists()) and force:
        raise RuntimeError("refusing to refreeze a design after baseline artifacts exist")
    if PROTOCOL_PATH.exists() and not force:
        protocol = protocol_payload()
        print(f"protocol already frozen: {protocol['protocol_sha256']}")
    else:
        protocol = build_protocol()
        atomic_json(PROTOCOL_PATH, protocol)
        print(f"[written] {PROTOCOL_PATH}")
    model = load_model(300)
    design_a = make_design(SOBOL_SEED_A, "sobolA", "confirmatory", model,
                           protocol["protocol_sha256"])
    design_b = make_design(SOBOL_SEED_B, "sobolB", "frozen_future_replicate",
                           model, protocol["protocol_sha256"])
    validation = validate_designs(design_a, design_b)
    design_a["cross_design_validation"] = validation
    design_b["cross_design_validation"] = validation
    if not validation["passes"]:
        raise RuntimeError("design validation failed: " +
                           "; ".join(validation["failures"]))
    for payload in (design_a, design_b):
        payload.pop("design_sha256", None)
        payload["design_sha256"] = object_hash(payload)
    atomic_json(DESIGN_A_PATH, design_a)
    atomic_json(DESIGN_B_PATH, design_b)
    print(f"[written] {DESIGN_A_PATH}")
    print(f"[written] {DESIGN_B_PATH}")
    print(json.dumps({
        "A_inclination_range_deg": design_a["realized_inclination_range_deg"],
        "A_regimes": design_a["inclination_regime_counts"],
        "B_inclination_range_deg": design_b["realized_inclination_range_deg"],
        "B_regimes": design_b["inclination_regime_counts"],
        "minimum_A_B_u_maxnorm": validation["minimum_A_B_u_maxnorm"],
    }, indent=2))
    return 0


def load_design_a() -> dict:
    protocol = protocol_payload()
    payload = json.loads(DESIGN_A_PATH.read_text(encoding="utf-8"))
    if payload["protocol_sha256"] != protocol["protocol_sha256"]:
        raise RuntimeError("design/protocol hash mismatch")
    if not payload.get("cross_design_validation", {}).get("passes"):
        raise RuntimeError("design validation is not passing")
    claimed = payload.pop("design_sha256")
    actual = object_hash(payload)
    payload["design_sha256"] = claimed
    if claimed != actual:
        raise RuntimeError("design A hash mismatch")
    return payload


def other_python_processes() -> list[dict]:
    try:
        import psutil
    except ImportError:
        return []
    found = []
    own = os.getpid()
    related = {own}
    try:
        related.update(parent.pid for parent in psutil.Process(own).parents())
    except (psutil.AccessDenied, psutil.NoSuchProcess):
        pass
    for proc in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            if proc.info["pid"] in related:
                continue
            name = (proc.info.get("name") or "").lower()
            if name in ("python.exe", "pythonw.exe", "python"):
                found.append({
                    "pid": proc.info["pid"],
                    "cmdline": proc.info.get("cmdline") or [],
                })
        except (psutil.AccessDenied, psutil.NoSuchProcess):
            continue
    return found


def propagate_event_instrumented(model, y0: np.ndarray, duration: float,
                                 t_grid: np.ndarray,
                                 degree_of: Callable[[float, float], int],
                                 args: tuple, rtol: float, atol,
                                 max_step: float = np.inf) -> tuple:
    rhs = Rhs(model, degree_of, args)
    wall_start = time.perf_counter_ns()
    cpu_start = time.process_time_ns()
    solver = InstrumentedDOP853(
        rhs, 0.0, np.asarray(y0, float), duration,
        rtol=rtol, atol=atol, max_step=max_step,
    )
    output = np.empty((6, len(t_grid)), dtype=float)
    output[:, 0] = y0
    filled = 1
    accepted_steps = 0
    impact_t = None
    impact_state = None
    status = "complete"
    failure_message = None
    previous_g = float(np.linalg.norm(y0[:3]) - model.r_ref)
    if previous_g <= 0.0:
        raise ValueError("initial state is on or below the reference surface")

    try:
        while solver.status == "running":
            old_t = float(solver.t)
            solver.step()
            if solver.status == "failed":
                raise RuntimeError("DOP853 failed")
            accepted_steps += 1
            new_t = float(solver.t)
            dense = solver.dense_output()
            new_g = float(np.linalg.norm(solver.y[:3]) - model.r_ref)
            interval_end = new_t
            if previous_g > 0.0 and new_g <= 0.0:
                impact_t = float(brentq(
                    lambda t: float(np.linalg.norm(dense(t)[:3]) - model.r_ref),
                    old_t,
                    new_t,
                    xtol=SURFACE_ROOT_XTOL_S,
                    rtol=4.0 * np.finfo(float).eps,
                ))
                impact_state = np.asarray(dense(impact_t), dtype=float)
                interval_end = impact_t
                status = "surface_impact"
            while filled < len(t_grid) and t_grid[filled] <= interval_end + 1e-9:
                output[:, filled] = dense(float(t_grid[filled]))
                filled += 1
            if impact_t is not None:
                break
            previous_g = new_g
        if status == "complete" and filled != len(t_grid):
            raise RuntimeError(
                f"completed integration filled {filled}/{len(t_grid)} epochs"
            )
    except Exception as exc:
        status = "numerical_failure"
        failure_message = f"{type(exc).__name__}: {exc}"

    cpu_ns = time.process_time_ns() - cpu_start
    wall_ns = time.perf_counter_ns() - wall_start
    times = np.asarray(t_grid[:filled], dtype=float)
    states = np.asarray(output[:, :filled], dtype=float)
    counts = {int(k): int(v) for k, v in rhs.deg_counts.items()}
    n_rhs = max(rhs.n_calls, 1)
    telemetry = {
        "n_rhs": int(rhs.n_calls),
        "n_accepted_steps": int(accepted_steps),
        "n_attempted_steps": int(solver.n_attempts),
        "n_rejected_trials": int(solver.n_rejected),
        "switch_count_at_rhs_samples": int(rhs.n_deg_changes),
        "degree_counts": counts,
        "degree_range": ([min(counts), max(counts)] if counts else None),
        "mean_degree": (
            float(sum(k * v for k, v in counts.items()) / n_rhs)
            if counts else None
        ),
        "mean_degree_sq": (
            float(rhs.sum_deg_sq / n_rhs) if rhs.n_calls else None
        ),
        "gravity_kernel_ns": int(rhs.grav_ns),
        "process_cpu_ns": int(cpu_ns),
        "total_wall_ns": int(wall_ns),
    }
    event = None
    if impact_t is not None:
        event = {
            "type": "reference_surface_downward_crossing",
            "epoch_s": impact_t,
            "state_si": [float(x) for x in impact_state],
            "root_residual_m": float(np.linalg.norm(impact_state[:3]) -
                                     model.r_ref),
        }
    return times, states, status, event, failure_message, telemetry


def trajectory_paths(run_kind: str, orbit_index: int, policy: str) -> tuple:
    case_dir = CASE_ROOT / run_kind / f"sobolA_{orbit_index:03d}"
    raw_dir = RAW_ROOT / run_kind / f"sobolA_{orbit_index:03d}"
    return case_dir / f"{policy}_baseline.json", raw_dir / f"{policy}_baseline.npz"


def valid_cached(sidecar_path: Path, raw_path: Path,
                 expected_config_sha: str, expected_duration: float) -> bool:
    try:
        meta = json.loads(sidecar_path.read_text(encoding="utf-8"))
        if meta["config_sha256"] != expected_config_sha:
            return False
        if meta["status"] not in ("complete", "surface_impact"):
            return False
        if meta["raw_sha256"] != file_hash(raw_path):
            return False
        with np.load(raw_path) as data:
            times = data["t_s"]
            states = data["state_si"]
            if states.shape != (6, len(times)) or len(times) == 0:
                return False
            if not np.all(np.isfinite(times)) or not np.all(np.isfinite(states)):
                return False
            if np.any(np.diff(times) <= 0.0):
                return False
            if meta["status"] == "complete" and abs(times[-1] - expected_duration) > 1e-6:
                return False
            if meta["status"] == "surface_impact" and not meta.get("event"):
                return False
        return True
    except (OSError, KeyError, ValueError, json.JSONDecodeError):
        return False


def preserve_invalid(path: Path) -> None:
    if path.exists():
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        target = path.with_name(path.name + f".invalid.{stamp}")
        os.replace(path, target)


def load_raw(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path) as data:
        return np.asarray(data["t_s"], dtype=float), np.asarray(
            data["state_si"], dtype=float
        )


def run_trajectory(run_kind: str, orbit: dict, policy: str, model, args,
                   degree_of: Callable[[float, float], int], policy_spec: dict,
                   duration: float, output_step: float,
                   timing_comparable: bool) -> tuple[dict, np.ndarray, np.ndarray]:
    protocol = protocol_payload()
    config = {
        "schema": "r10_trajectory_config_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "design_sha256": load_design_a()["design_sha256"],
        "sobol_seed": SOBOL_SEED_A,
        "sobol_index": orbit["sobol_index"],
        "original_sobol_coordinates": orbit["u"],
        "initial_state_si": orbit["initial_state_si"],
        "truth_degree_rule_result": orbit["truth_degree"],
        "policy": policy,
        "policy_spec": policy_spec,
        "duration_s": duration,
        "output_step_s": output_step,
        "integrator": "InstrumentedDOP853",
        "rtol": RTOL,
        "atol": ATOL,
        "atol_kind": "scalar",
        "max_step_s": None,
        "surface_event": "norm(r)-r_ref downward terminal",
        "timing_comparable": timing_comparable,
        "source": protocol["provenance"],
    }
    config_sha = object_hash(config)
    sidecar_path, raw_path = trajectory_paths(run_kind, orbit["sobol_index"], policy)
    if sidecar_path.exists() and raw_path.exists() and valid_cached(
            sidecar_path, raw_path, config_sha, duration):
        meta = json.loads(sidecar_path.read_text(encoding="utf-8"))
        times, states = load_raw(raw_path)
        print(f"  {policy:16s} cached {meta['status']}", flush=True)
        return meta, times, states
    if sidecar_path.exists() or raw_path.exists():
        preserve_invalid(sidecar_path)
        preserve_invalid(raw_path)

    active = {
        "updated_utc": utc_now(),
        "run_kind": run_kind,
        "sobol_index": orbit["sobol_index"],
        "orbit": orbit["name"],
        "policy": policy,
        "config_sha256": config_sha,
    }
    atomic_json(ACTIVE_PATH, active)
    grid = np.arange(0.0, duration + 0.5 * output_step, output_step)
    times, states, status, event, failure, telemetry = propagate_event_instrumented(
        model,
        np.asarray(orbit["initial_state_si"], dtype=float),
        duration,
        grid,
        degree_of,
        args,
        RTOL,
        ATOL,
        max_step=np.inf,
    )
    if status == "numerical_failure":
        failure_path = sidecar_path.with_name(
            sidecar_path.stem + f".failure.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
        )
        atomic_json(failure_path, {
            "schema": "r10_trajectory_failure_v1",
            "created_utc": utc_now(),
            "config": config,
            "config_sha256": config_sha,
            "status": status,
            "failure_message": failure,
            "telemetry": telemetry,
        })
        raise RuntimeError(f"{orbit['name']} {policy}: {failure}")

    arrays = {"t_s": times, "state_si": states}
    if event is not None:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"], dtype=float)
    atomic_npz(raw_path, **arrays)
    meta = {
        "schema": "r10_trajectory_result_v1",
        "created_utc": utc_now(),
        "config": config,
        "config_sha256": config_sha,
        "status": status,
        "event": event,
        "failure_message": None,
        "telemetry": telemetry,
        "raw_path": str(raw_path.relative_to(ROOT)),
        "raw_sha256": file_hash(raw_path),
        "n_output_epochs": int(len(times)),
        "last_output_epoch_s": float(times[-1]),
    }
    atomic_json(sidecar_path, meta)
    print(
        f"  {policy:16s} {status:14s} wall "
        f"{telemetry['total_wall_ns']/1e9:8.1f}s rhs {telemetry['n_rhs']}",
        flush=True,
    )
    return meta, times, states


def common_error(policy_t: np.ndarray, policy_y: np.ndarray,
                 truth_t: np.ndarray, truth_y: np.ndarray) -> dict:
    n = min(len(policy_t), len(truth_t))
    if n == 0 or not np.allclose(policy_t[:n], truth_t[:n], rtol=0.0, atol=1e-9):
        raise RuntimeError("policy/truth output grids do not share a common prefix")
    stats = err_stats(policy_y[:, :n], truth_y[:, :n])
    stats["common_epoch_count"] = n
    stats["common_arc_end_s"] = float(policy_t[n - 1])
    stats["in_track_fraction"] = (
        stats["ric_rms_m"]["in_track"] / max(stats["pos_rms_m"], 1e-300)
    )
    return stats


def add_error_to_meta(meta: dict, sidecar_path: Path, error: dict,
                      is_full_arc: bool) -> dict:
    meta["error_against_truth"] = error
    meta["error_scope"] = "full_seven_day" if is_full_arc else "last_common_arc"
    atomic_json(sidecar_path, meta)
    return meta


def run_orbit(run_kind: str, orbit: dict, models: dict, tables: dict,
              powers: dict, duration: float, output_step: float,
              timing_comparable: bool) -> dict:
    truth_degree = int(orbit["truth_degree"])
    model, args = models[truth_degree]
    truth_fn = lambda t, h, n=truth_degree: n
    truth_meta, truth_t, truth_y = run_trajectory(
        run_kind, orbit, "truth", model, args, truth_fn,
        {"kind": "fixed_truth", "degree": truth_degree},
        duration, output_step, timing_comparable,
    )

    empirical_fn = alt_sched(tables[("emp", truth_degree)])
    emp_meta, emp_t, emp_y = run_trajectory(
        run_kind, orbit, "schedule_empirical", model, args, empirical_fn,
        {"kind": "empirical_lookup", "floor": FLOOR, "cap": CAP,
         "quantum": Q, "degree_quantization": "downward",
         "altitude_bin": "10-km floor"},
        duration, output_step, timing_comparable,
    )
    mean_n2 = emp_meta["telemetry"]["mean_degree_sq"]
    if mean_n2 is None:
        raise RuntimeError(f"{orbit['name']}: empirical schedule has no RHS calls")
    n_work = int(round(math.sqrt(mean_n2)))
    n_critical = int(min(
        CAP,
        emp_nmin_exact(powers[truth_degree], model.r_ref,
                       orbit["hp_km"] * 1000.0),
    ))
    policy_defs = [
        ("fixed_work", lambda t, h, n=n_work: n,
         {"kind": "fixed_work", "degree": n_work,
          "derived_from": "empirical actual baseline RHS-call mean N^2"}),
        ("fixed_critical", lambda t, h, n=n_critical: n,
         {"kind": "fixed_critical", "degree": n_critical,
          "definition": "unquantized empirical Nmin at own perilune, cap 250"}),
        ("schedule_up", alt_sched(tables[("up", truth_degree)]),
         {"kind": "Kaula_p2_upward_quantized", "floor": FLOOR,
          "cap": CAP, "quantum": Q}),
        ("schedule_down", alt_sched(tables[("down", truth_degree)]),
         {"kind": "Kaula_p2_downward_quantized", "floor": FLOOR,
          "cap": CAP, "quantum": Q}),
    ]
    results: dict[str, dict] = {}
    emp_sidecar, _ = trajectory_paths(run_kind, orbit["sobol_index"],
                                      "schedule_empirical")
    empirical_error = common_error(emp_t, emp_y, truth_t, truth_y)
    emp_meta = add_error_to_meta(
        emp_meta, emp_sidecar, empirical_error,
        emp_meta["status"] == truth_meta["status"] == "complete",
    )
    results["schedule_empirical"] = emp_meta

    for name, function, spec in policy_defs:
        meta, times, states = run_trajectory(
            run_kind, orbit, name, model, args, function, spec,
            duration, output_step, timing_comparable,
        )
        error = common_error(times, states, truth_t, truth_y)
        sidecar, _ = trajectory_paths(run_kind, orbit["sobol_index"], name)
        meta = add_error_to_meta(
            meta, sidecar, error,
            meta["status"] == truth_meta["status"] == "complete",
        )
        results[name] = meta

    eligible = truth_meta["status"] == "complete"
    ratios = None
    if eligible:
        empirical_complete = results["schedule_empirical"]["status"] == "complete"
        work_complete = results["fixed_work"]["status"] == "complete"
        critical_complete = results["fixed_critical"]["status"] == "complete"
        if empirical_complete and work_complete and critical_complete:
            empirical_error_value = empirical_error["pos_rms_m"]
            rho_work = (results["fixed_work"]["error_against_truth"]["pos_rms_m"] /
                        empirical_error_value)
            rho_crit = (results["fixed_critical"]["error_against_truth"]["pos_rms_m"] /
                        empirical_error_value)
            ratios = {
                "rho_work": float(rho_work),
                "rho_crit": float(rho_crit),
                "empirical_raw_win_vs_work": bool(rho_work > 1.0),
                "empirical_raw_win_vs_critical": bool(rho_crit > 1.0),
                "gravity_time_saving_vs_critical": float(
                    1.0 - results["schedule_empirical"]["telemetry"]["gravity_kernel_ns"] /
                    max(results["fixed_critical"]["telemetry"]["gravity_kernel_ns"], 1)
                ),
            }
    trajectory_meta = {"truth": truth_meta, **results}
    return {
        "name": orbit["name"],
        "sobol_index": orbit["sobol_index"],
        "design_point": orbit,
        "truth_degree": truth_degree,
        "truth_survives_full_arc": truth_meta["status"] == "complete",
        "n_work": n_work,
        "n_critical": n_critical,
        "trajectory_status": {k: v["status"] for k, v in trajectory_meta.items()},
        "policies": {
            k: {
                "status": v["status"],
                "error_scope": v.get("error_scope"),
                "error_against_truth": v.get("error_against_truth"),
                "telemetry": v["telemetry"],
                "event": v.get("event"),
                "config_sha256": v["config_sha256"],
                "raw_path": v["raw_path"],
                "raw_sha256": v["raw_sha256"],
            }
            for k, v in trajectory_meta.items()
        },
        "primary_ratios": ratios,
        "orbit_total_trajectory_wall_s": float(sum(
            value["telemetry"]["total_wall_ns"] for value in trajectory_meta.values()
        ) / 1e9),
    }


def summarize_rows(rows: list[dict]) -> dict:
    eligible = [r for r in rows if r.get("primary_ratios") is not None]
    summary = {
        "completed_orbit_records": len(rows),
        "eligible_primary_comparisons": len(eligible),
        "truth_surviving": sum(r.get("truth_survives_full_arc") is True
                               for r in rows),
        "truth_impacts": sum(r.get("truth_survives_full_arc") is False
                             for r in rows),
        "orbit_numerical_failures": sum(r.get("batch_status") ==
                                          "numerical_failure" for r in rows),
        "policy_impact_counts": {},
    }
    for policy in ("schedule_empirical", "fixed_work", "fixed_critical",
                   "schedule_up", "schedule_down"):
        summary["policy_impact_counts"][policy] = sum(
            r.get("trajectory_status", {}).get(policy) == "surface_impact"
            for r in rows
        )
    if eligible:
        for key in ("rho_work", "rho_crit"):
            values = np.asarray([r["primary_ratios"][key] for r in eligible])
            summary[key] = {
                "median": float(np.median(values)),
                "p10": float(np.percentile(values, 10)),
                "p90": float(np.percentile(values, 90)),
                "raw_empirical_win_count": int(np.sum(values > 1.0)),
            }
    return summary


def load_existing_rows(run_kind: str) -> list[dict]:
    target = SMOKE_PATH if run_kind == "smoke" else BASELINE_PATH
    if not target.exists():
        return []
    payload = json.loads(target.read_text(encoding="utf-8"))
    return payload.get("rows", [])


def write_batch_index(run_kind: str, rows: list[dict], design: dict,
                      started_utc: str, timing_comparable: bool,
                      complete: bool, concurrency_events: list[dict],
                      session_wall_ns: int, session_cpu_ns: int) -> None:
    target = SMOKE_PATH if run_kind == "smoke" else BASELINE_PATH
    payload = {
        "schema": "r10_sobolA_baseline_index_v1",
        "run_kind": run_kind,
        "protocol_sha256": protocol_payload()["protocol_sha256"],
        "design_sha256": design["design_sha256"],
        "started_utc": started_utc,
        "updated_utc": utc_now(),
        "complete": complete,
        "timing_comparable": timing_comparable,
        "concurrency_events": concurrency_events,
        "current_session_wall_ns": int(session_wall_ns),
        "current_session_cpu_ns": int(session_cpu_ns),
        "planned_orbit_count": len(design["orbits"]) if run_kind == "baseline" else 2,
        "rows": rows,
        "summary": summarize_rows(rows),
    }
    if complete:
        payload["ended_utc"] = utc_now()
        payload["actual_total_trajectory_wall_s"] = float(sum(
            r["orbit_total_trajectory_wall_s"] for r in rows
        ))
    atomic_json(target, payload)


def command_run(run_kind: str, allow_concurrent: bool) -> int:
    design = load_design_a()
    initial_other = other_python_processes()
    timing_comparable = not initial_other
    if initial_other and not allow_concurrent:
        raise RuntimeError(
            "other Python processes are active; refusing a timing-comparable run: " +
            json.dumps(initial_other)
        )
    concurrency_events: list[dict] = []
    if initial_other:
        concurrency_events.append({"utc": utc_now(), "processes": initial_other})
    orbits = design["orbits"]
    duration = DURATION_S
    output_step = OUTPUT_STEP_S
    if run_kind == "smoke":
        lowest = min(orbits, key=lambda o: o["hp_km"])
        middle = min(orbits, key=lambda o: abs(o["hp_km"] - 90.0))
        orbits = [lowest, middle]
        duration = SMOKE_DURATION_S
    rows = load_existing_rows(run_kind)
    existing = {r["sobol_index"]: r for r in rows}
    started_utc = utc_now()
    session_wall_start = time.perf_counter_ns()
    session_cpu_start = time.process_time_ns()

    needed_degrees = sorted({int(o["truth_degree"]) for o in orbits})
    models: dict[int, tuple] = {}
    powers: dict[int, np.ndarray] = {}
    tables: dict[tuple[str, int], dict] = {}
    for degree in needed_degrees:
        model = load_model(degree)
        args = kernel_args(model)
        warmup(model, args)
        models[degree] = (model, args)
        powers[degree] = degree_power(model)
        tables[("emp", degree)] = emp_table(model, powers[degree])
        tables[("up", degree)] = kaula_table(model, "up")
        tables[("down", degree)] = kaula_table(model, "down")
        print(f"[model] N={degree} loaded and warmed", flush=True)

    print(
        f"[{run_kind}] {len(orbits)} orbits, timing_comparable={timing_comparable}",
        flush=True,
    )
    for position, orbit in enumerate(orbits, start=1):
        current_other = other_python_processes()
        if current_other:
            timing_comparable = False
            event = {"utc": utc_now(), "before_orbit": orbit["sobol_index"],
                     "processes": current_other}
            if not concurrency_events or concurrency_events[-1]["processes"] != current_other:
                concurrency_events.append(event)
            if not allow_concurrent:
                write_batch_index(run_kind, list(existing.values()), design,
                                  started_utc, False, False, concurrency_events,
                                  time.perf_counter_ns() - session_wall_start,
                                  time.process_time_ns() - session_cpu_start)
                raise RuntimeError("another Python process appeared during the formal run")
        print(
            f"[{position}/{len(orbits)}] {orbit['name']} hp={orbit['hp_km']:.3f} "
            f"ha={orbit['ha_km']:.3f} i={orbit['incl_deg']:.3f}",
            flush=True,
        )
        try:
            row = run_orbit(
                run_kind, orbit, models, tables, powers, duration, output_step,
                timing_comparable,
            )
            row["batch_status"] = "completed"
        except Exception as exc:
            row = {
                "name": orbit["name"],
                "sobol_index": orbit["sobol_index"],
                "design_point": orbit,
                "batch_status": "numerical_failure",
                "failure_message": f"{type(exc).__name__}: {exc}",
                "truth_survives_full_arc": None,
                "trajectory_status": {},
                "primary_ratios": None,
                "orbit_total_trajectory_wall_s": 0.0,
            }
            print(f"  [failure preserved] {row['failure_message']}", flush=True)
        existing[orbit["sobol_index"]] = row
        ordered = [existing[key] for key in sorted(existing)]
        write_batch_index(run_kind, ordered, design, started_utc,
                          timing_comparable, False, concurrency_events,
                          time.perf_counter_ns() - session_wall_start,
                          time.process_time_ns() - session_cpu_start)
    ordered = [existing[key] for key in sorted(existing)]
    write_batch_index(run_kind, ordered, design, started_utc,
                      timing_comparable, True, concurrency_events,
                      time.perf_counter_ns() - session_wall_start,
                      time.process_time_ns() - session_cpu_start)
    if ACTIVE_PATH.exists():
        ACTIVE_PATH.unlink()
    print(f"[{run_kind}] complete", flush=True)
    return 0


def command_status() -> int:
    if BASELINE_PATH.exists():
        payload = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
        print(json.dumps({
            "complete": payload.get("complete"),
            "completed_orbit_records": len(payload.get("rows", [])),
            "planned_orbit_count": payload.get("planned_orbit_count"),
            "timing_comparable": payload.get("timing_comparable"),
            "updated_utc": payload.get("updated_utc"),
            "summary": payload.get("summary"),
        }, indent=2))
    else:
        print("no formal baseline index yet")
    if ACTIVE_PATH.exists():
        print("active:")
        print(ACTIVE_PATH.read_text(encoding="utf-8"))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    design_parser = sub.add_parser("design")
    design_parser.add_argument("--force", action="store_true")
    smoke_parser = sub.add_parser("smoke")
    smoke_parser.add_argument("--allow-concurrent", action="store_true")
    baseline_parser = sub.add_parser("baseline")
    baseline_parser.add_argument("--allow-concurrent", action="store_true")
    sub.add_parser("status")
    arguments = parser.parse_args()
    if arguments.command == "design":
        return command_design(arguments.force)
    if arguments.command == "smoke":
        return command_run("smoke", arguments.allow_concurrent)
    if arguments.command == "baseline":
        return command_run("baseline", arguments.allow_concurrent)
    return command_status()


if __name__ == "__main__":
    raise SystemExit(main())
