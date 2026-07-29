"""28-day LRO-like corrected-potential-blend convergence batch."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import rev10_sobol_confirmatory as base
from rev7_doe_screening import CANONICAL, initial_state
from rev9_potential_blend_longarc import BlendRhs


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r10_blend_lro_convergence.json"
SMOKE_OUTPUT = METRICS / "r10_blend_lro_convergence_smoke.json"
ACTIVE = METRICS / "r10_blend_lro_convergence_active.json"
CASE_ROOT = METRICS / "r10_cases" / "blend_lro_convergence"
RAW_ROOT = METRICS / "r10_raw" / "blend_lro_convergence"

DURATION = 28.0 * base.DAY
OUTPUT_STEP = 300.0
MAX_STEP = 120.0
LEVELS = {
    "baseline": {"rtol": 1.0e-11, "atol": 1.0e-6},
    "tighter": {"rtol": 1.0e-12, "atol": 1.0e-7},
}
# Heavy-first. Estimates are conservative start guards, not reported timings.
RUN_ORDER = (
    ("corrected_blend", "tighter", 15000.0),
    ("corrected_blend", "baseline", 10000.0),
    ("truth_N600", "tighter", 3600.0),
    ("truth_N600", "baseline", 2200.0),
    ("fixed_N120", "tighter", 300.0),
    ("fixed_N120", "baseline", 300.0),
)


def orbit_dict() -> dict:
    for name, hp, ha, inc, argp, raan in CANONICAL:
        if name == "c6_lro_30x216":
            return {
                "name": name, "hp_km": hp, "ha_km": ha,
                "incl_deg": inc, "argp_deg": argp, "raan_deg": raan,
            }
    raise RuntimeError("LRO canonical orbit not found")


def parse_deadline(value: str | None) -> datetime | None:
    if value is None:
        return None
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        raise ValueError("deadline must contain an explicit UTC offset")
    return result.astimezone(timezone.utc)


def remaining_seconds(deadline: datetime | None) -> float:
    if deadline is None:
        return math.inf
    return (deadline - datetime.now(timezone.utc)).total_seconds()


def generic_propagate(rhs, model, y0, duration, grid, rtol, atol):
    wall0 = time.perf_counter_ns()
    cpu0 = time.process_time_ns()
    solver = base.InstrumentedDOP853(
        rhs, 0.0, np.asarray(y0, float), duration,
        rtol=rtol, atol=atol, max_step=MAX_STEP,
    )
    output = np.empty((6, len(grid)), dtype=float)
    output[:, 0] = y0
    filled = 1
    accepted = 0
    impact_t = None
    impact_state = None
    status = "complete"
    failure = None
    previous_g = float(np.linalg.norm(y0[:3]) - model.r_ref)
    try:
        while solver.status == "running":
            old_t = float(solver.t)
            solver.step()
            if solver.status == "failed":
                raise RuntimeError("DOP853 failed")
            accepted += 1
            new_t = float(solver.t)
            dense = solver.dense_output()
            new_g = float(np.linalg.norm(solver.y[:3]) - model.r_ref)
            end = new_t
            if previous_g > 0.0 and new_g <= 0.0:
                impact_t = float(brentq(
                    lambda t: float(np.linalg.norm(dense(t)[:3]) - model.r_ref),
                    old_t, new_t, xtol=1.0e-8,
                    rtol=4.0 * np.finfo(float).eps,
                ))
                impact_state = np.asarray(dense(impact_t), dtype=float)
                end = impact_t
                status = "surface_impact"
            while filled < len(grid) and grid[filled] <= end + 1.0e-9:
                output[:, filled] = dense(float(grid[filled]))
                filled += 1
            if impact_t is not None:
                break
            previous_g = new_g
        if status == "complete" and filled != len(grid):
            raise RuntimeError(f"filled {filled}/{len(grid)} output epochs")
    except Exception as exc:
        status = "numerical_failure"
        failure = f"{type(exc).__name__}: {exc}"
    telemetry = {
        "n_rhs": int(rhs.n_calls),
        "n_accepted_steps": int(accepted),
        "n_attempted_steps": int(solver.n_attempts),
        "n_rejected_trials": int(solver.n_rejected),
        "gravity_kernel_ns": int(rhs.grav_ns),
        "process_cpu_ns": int(time.process_time_ns() - cpu0),
        "total_wall_ns": int(time.perf_counter_ns() - wall0),
    }
    event = None
    if impact_t is not None:
        event = {
            "type": "reference_surface_downward_crossing",
            "epoch_s": impact_t,
            "state_si": [float(x) for x in impact_state],
            "root_residual_m": float(np.linalg.norm(impact_state[:3]) - model.r_ref),
        }
    return (np.asarray(grid[:filled]), np.asarray(output[:, :filled]),
            status, event, failure, telemetry)


def paths(policy: str, level: str, smoke: bool) -> tuple[Path, Path]:
    suffix = "_smoke" if smoke else ""
    return (
        CASE_ROOT / f"{policy}_{level}{suffix}.json",
        RAW_ROOT / f"{policy}_{level}{suffix}.npz",
    )


def run_case(policy: str, level: str, model, args, y0, duration: float,
             smoke: bool) -> tuple[dict, np.ndarray, np.ndarray]:
    protocol = base.protocol_payload()
    tol = LEVELS[level]
    config = {
        "schema": "r10_blend_lro_convergence_config_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "orbit": orbit_dict(),
        "initial_state_si": [float(x) for x in y0],
        "policy": policy,
        "level": level,
        "duration_s": duration,
        "output_step_s": OUTPUT_STEP,
        "integrator": "InstrumentedDOP853",
        "rtol": tol["rtol"],
        "atol": tol["atol"],
        "atol_kind": "scalar",
        "max_step_s": MAX_STEP,
        "truth_degree": 600,
        "fixed_degree": 120,
        "blend_degrees": [30, 120],
        "transition_altitude_m": [50000.0, 200000.0],
        "corrected_term": "(U_hi-U_lo)*dw/dr*r_hat",
        "orientation": "uniform lunar sidereal rotation",
        "timing_comparable": True,
        "source": protocol["provenance"],
    }
    config_sha = base.object_hash(config)
    sidecar, raw = paths(policy, level, smoke)
    if sidecar.exists() and raw.exists() and base.valid_cached(
            sidecar, raw, config_sha, duration):
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        t, y = base.load_raw(raw)
        print(f"  {policy}/{level} cached {meta['status']}", flush=True)
        return meta, t, y
    if sidecar.exists() or raw.exists():
        base.preserve_invalid(sidecar)
        base.preserve_invalid(raw)
    if policy == "corrected_blend":
        rhs = BlendRhs(model, args, "blend_potential_corrected")
    elif policy == "truth_N600":
        rhs = base.Rhs(model, lambda t, h: 600, args)
    elif policy == "fixed_N120":
        rhs = base.Rhs(model, lambda t, h: 120, args)
    else:
        raise ValueError(policy)
    base.atomic_json(ACTIVE, {
        "updated_utc": base.utc_now(), "policy": policy, "level": level,
        "config_sha256": config_sha,
    })
    grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    t, y, status, event, failure, telemetry = generic_propagate(
        rhs, model, y0, duration, grid, tol["rtol"], tol["atol"]
    )
    if status == "numerical_failure":
        fail = sidecar.with_name(sidecar.stem + f".failure.{int(time.time())}.json")
        base.atomic_json(fail, {
            "config": config, "config_sha256": config_sha,
            "status": status, "failure_message": failure,
            "telemetry": telemetry,
        })
        raise RuntimeError(f"{policy}/{level}: {failure}")
    arrays = {"t_s": t, "state_si": y}
    if event:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"])
    base.atomic_npz(raw, **arrays)
    meta = {
        "schema": "r10_blend_lro_convergence_result_v1",
        "created_utc": base.utc_now(), "config": config,
        "config_sha256": config_sha, "status": status, "event": event,
        "failure_message": None, "telemetry": telemetry,
        "raw_path": str(raw.relative_to(ROOT)), "raw_sha256": base.file_hash(raw),
        "n_output_epochs": int(len(t)), "last_output_epoch_s": float(t[-1]),
    }
    base.atomic_json(sidecar, meta)
    print(f"  {policy}/{level} {status} wall={telemetry['total_wall_ns']/1e9:.1f}s", flush=True)
    return meta, t, y


def compute_summary(records: dict) -> dict | None:
    required = {(p, l) for p in ("truth_N600", "fixed_N120", "corrected_blend")
                for l in LEVELS}
    if not required.issubset(records):
        return None
    data = {}
    for key, record in records.items():
        sidecar, raw = paths(key[0], key[1], False)
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        t, y = base.load_raw(raw)
        data[key] = (meta, t, y)
    result = {"policies": {}}
    truth_self = base.common_error(
        data[("truth_N600", "baseline")][1], data[("truth_N600", "baseline")][2],
        data[("truth_N600", "tighter")][1], data[("truth_N600", "tighter")][2],
    )["pos_rms_m"]
    result["truth_self_difference_rms_m"] = truth_self
    for policy in ("fixed_N120", "corrected_blend"):
        errors = {}
        for level in LEVELS:
            errors[level] = base.common_error(
                data[(policy, level)][1], data[(policy, level)][2],
                data[("truth_N600", level)][1], data[("truth_N600", level)][2],
            )
        self_diff = base.common_error(
            data[(policy, "baseline")][1], data[(policy, "baseline")][2],
            data[(policy, "tighter")][1], data[(policy, "tighter")][2],
        )["pos_rms_m"]
        result["policies"][policy] = {
            "error_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self,
        }
    fixed = result["policies"]["fixed_N120"]
    blend = result["policies"]["corrected_blend"]
    difference = abs(
        blend["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
        - fixed["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
    )
    threshold = (blend["truth_inclusive_envelope_m"] +
                 fixed["truth_inclusive_envelope_m"])
    result["comparison"] = {
        "absolute_baseline_error_difference_m": difference,
        "resolution_threshold_m": threshold,
        "resolved": bool(difference > threshold),
        "winner_if_resolved": (
            "corrected_blend"
            if blend["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
            < fixed["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
            else "fixed_N120"
        ) if difference > threshold else None,
    }
    return result


def write_index(records: dict, complete: bool, skipped: list, started: str,
                wall0: int, deadline: datetime | None, smoke: bool) -> None:
    target = SMOKE_OUTPUT if smoke else OUTPUT
    payload = {
        "schema": "r10_blend_lro_convergence_index_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "started_utc": started, "updated_utc": base.utc_now(),
        "deadline_utc": deadline.isoformat() if deadline else None,
        "complete": complete, "timing_comparable": True,
        "records": [records[key] for key in sorted(records)],
        "skipped_for_deadline": skipped,
        "summary": None if smoke else compute_summary(records),
        "session_wall_ns": time.perf_counter_ns() - wall0,
    }
    if complete:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(target, payload)


def run(smoke: bool, deadline: datetime | None) -> int:
    if base.other_python_processes():
        raise RuntimeError("another Python process is active")
    protocol = base.protocol_payload()
    model = base.load_model(600)
    args = base.kernel_args(model)
    base.warmup(model, args)
    orbit = orbit_dict()
    y0 = initial_state(model, orbit)
    duration = 2.0 * 3600.0 if smoke else DURATION
    order = (("corrected_blend", "tighter", 600.0),) if smoke else RUN_ORDER
    records = {}
    skipped = []
    target = SMOKE_OUTPUT if smoke else OUTPUT
    if target.exists():
        old = json.loads(target.read_text(encoding="utf-8"))
        records = {(r["policy"], r["level"]): r for r in old.get("records", [])}
    started = base.utc_now()
    wall0 = time.perf_counter_ns()
    print(f"[blend] cases={len(order)} deadline={deadline} source={protocol['protocol_sha256']}", flush=True)
    for policy, level, estimate in order:
        if base.other_python_processes():
            raise RuntimeError("another Python process appeared during blend batch")
        remaining = remaining_seconds(deadline)
        if remaining < estimate + 600.0:
            skipped.append({"policy": policy, "level": level,
                            "remaining_s": remaining, "guard_estimate_s": estimate})
            print(f"  skip {policy}/{level}: deadline guard", flush=True)
            continue
        meta, _, _ = run_case(policy, level, model, args, y0, duration, smoke)
        records[(policy, level)] = {
            "policy": policy, "level": level, "status": meta["status"],
            "config_sha256": meta["config_sha256"],
            "raw_path": meta["raw_path"], "raw_sha256": meta["raw_sha256"],
            "telemetry": meta["telemetry"],
        }
        write_index(records, False, skipped, started, wall0, deadline, smoke)
    required_count = len(order)
    complete = len(records) == required_count
    write_index(records, complete, skipped, started, wall0, deadline, smoke)
    if ACTIVE.exists():
        ACTIVE.unlink()
    print(f"[blend] finished complete={complete} records={len(records)}/{required_count}", flush=True)
    return 0 if complete else 3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("smoke", "run", "status"))
    parser.add_argument("--deadline")
    args = parser.parse_args()
    if args.command == "status":
        if OUTPUT.exists():
            data = json.loads(OUTPUT.read_text(encoding="utf-8"))
            print(json.dumps({"complete": data["complete"],
                              "records": len(data["records"]),
                              "skipped": data["skipped_for_deadline"],
                              "summary": data["summary"]}, indent=2))
        else:
            print("no blend convergence output")
        return 0
    return run(args.command == "smoke", parse_deadline(args.deadline))


if __name__ == "__main__":
    raise SystemExit(main())
