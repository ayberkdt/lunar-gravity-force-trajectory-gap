"""Selective N=600/N=900 truth-degree audit for Sobol A.

The six audit cases are the frozen union of the four lowest-perilune truth
survivors and every sub-50-km raw empirical-schedule win in the completed
R10 baseline. Existing validated N=600 and policy trajectories are reused;
only the six N=900 truth trajectories are newly propagated.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
BASELINE_PATH = METRICS / "r10_sobolA_baseline.json"
OUTPUT_PATH = METRICS / "r10_sobolA_truth_audit.json"
SMOKE_PATH = METRICS / "r10_sobolA_truth_audit_smoke.json"
ACTIVE_PATH = METRICS / "r10_sobolA_truth_audit_active.json"
CASE_ROOT = METRICS / "r10_cases"
RAW_ROOT = METRICS / "r10_raw"

AUDIT_INDICES = (4, 20, 27, 28, 36, 59)
POLICIES = (
    "schedule_empirical",
    "fixed_work",
    "fixed_critical",
    "schedule_up",
    "schedule_down",
)


def load_baseline() -> dict:
    payload = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    if not payload.get("complete") or len(payload.get("rows", [])) != 64:
        raise RuntimeError("the 64-orbit baseline is not complete")
    if payload.get("summary", {}).get("orbit_numerical_failures") != 0:
        raise RuntimeError("the baseline contains numerical failures")
    if not payload.get("timing_comparable"):
        raise RuntimeError("the baseline timing contract is not comparable")
    protocol = base.protocol_payload()
    design = base.load_design_a()
    if payload["protocol_sha256"] != protocol["protocol_sha256"]:
        raise RuntimeError("baseline/protocol hash mismatch")
    if payload["design_sha256"] != design["design_sha256"]:
        raise RuntimeError("baseline/design hash mismatch")
    return payload


def derive_audit_indices(baseline: dict) -> tuple[int, ...]:
    rows = [row for row in baseline["rows"]
            if row.get("truth_survives_full_arc")]
    lowest = {
        row["sobol_index"]
        for row in sorted(rows, key=lambda item: item["design_point"]["hp_km"])[:4]
    }
    sub50_raw_wins = {
        row["sobol_index"]
        for row in rows
        if row["design_point"]["hp_km"] < 50.0
        and row.get("primary_ratios") is not None
        and (
            row["primary_ratios"]["rho_work"] > 1.0
            or row["primary_ratios"]["rho_crit"] > 1.0
        )
    }
    derived = tuple(sorted(lowest | sub50_raw_wins))
    if derived != AUDIT_INDICES:
        raise RuntimeError(
            f"frozen audit indices {AUDIT_INDICES} differ from derived {derived}"
        )
    return derived


def audit_paths(run_kind: str, index: int) -> tuple[Path, Path]:
    case_dir = CASE_ROOT / run_kind / f"sobolA_{index:03d}"
    raw_dir = RAW_ROOT / run_kind / f"sobolA_{index:03d}"
    return case_dir / "truth_N900_baseline.json", raw_dir / "truth_N900_baseline.npz"


def load_validated_baseline_trajectory(index: int, policy: str) -> tuple:
    sidecar, raw = base.trajectory_paths("baseline", index, policy)
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    if meta["status"] not in ("complete", "surface_impact"):
        raise RuntimeError(f"baseline {index}/{policy} is not terminal")
    if meta["raw_sha256"] != base.file_hash(raw):
        raise RuntimeError(f"baseline {index}/{policy} raw hash mismatch")
    times, states = base.load_raw(raw)
    return meta, times, states


def run_n900(run_kind: str, orbit: dict, model, args, duration: float,
             timing_comparable: bool) -> tuple[dict, np.ndarray, np.ndarray]:
    protocol = base.protocol_payload()
    design = base.load_design_a()
    config = {
        "schema": "r10_truth_audit_trajectory_config_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "design_sha256": design["design_sha256"],
        "baseline_sha256": base.file_hash(BASELINE_PATH),
        "audit_script": str(Path(__file__).resolve()),
        "audit_script_sha256": base.file_hash(Path(__file__).resolve()),
        "sobol_seed": base.SOBOL_SEED_A,
        "sobol_index": orbit["sobol_index"],
        "original_sobol_coordinates": orbit["u"],
        "initial_state_si": orbit["initial_state_si"],
        "policy": "fixed_truth_N900",
        "degree": 900,
        "duration_s": duration,
        "output_step_s": base.OUTPUT_STEP_S,
        "integrator": "InstrumentedDOP853",
        "rtol": base.RTOL,
        "atol": base.ATOL,
        "atol_kind": "scalar",
        "max_step_s": None,
        "surface_event": "norm(r)-r_ref downward terminal",
        "timing_comparable": timing_comparable,
        "source": protocol["provenance"],
    }
    config_sha = base.object_hash(config)
    sidecar_path, raw_path = audit_paths(run_kind, orbit["sobol_index"])
    if sidecar_path.exists() and raw_path.exists() and base.valid_cached(
            sidecar_path, raw_path, config_sha, duration):
        meta = json.loads(sidecar_path.read_text(encoding="utf-8"))
        times, states = base.load_raw(raw_path)
        print(f"  N=900 cached {meta['status']}", flush=True)
        return meta, times, states
    if sidecar_path.exists() or raw_path.exists():
        base.preserve_invalid(sidecar_path)
        base.preserve_invalid(raw_path)

    base.atomic_json(ACTIVE_PATH, {
        "updated_utc": base.utc_now(),
        "run_kind": run_kind,
        "sobol_index": orbit["sobol_index"],
        "orbit": orbit["name"],
        "policy": "truth_N900",
        "config_sha256": config_sha,
    })
    grid = np.arange(0.0, duration + 0.5 * base.OUTPUT_STEP_S,
                     base.OUTPUT_STEP_S)
    degree = lambda t, h: 900
    times, states, status, event, failure, telemetry = (
        base.propagate_event_instrumented(
            model,
            np.asarray(orbit["initial_state_si"], dtype=float),
            duration,
            grid,
            degree,
            args,
            base.RTOL,
            base.ATOL,
            max_step=np.inf,
        )
    )
    if status == "numerical_failure":
        failure_path = sidecar_path.with_name(
            sidecar_path.stem + f".failure.{int(time.time())}.json"
        )
        base.atomic_json(failure_path, {
            "schema": "r10_truth_audit_failure_v1",
            "created_utc": base.utc_now(),
            "config": config,
            "config_sha256": config_sha,
            "status": status,
            "failure_message": failure,
            "telemetry": telemetry,
        })
        raise RuntimeError(f"{orbit['name']} N=900: {failure}")

    arrays = {"t_s": times, "state_si": states}
    if event is not None:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"], dtype=float)
    base.atomic_npz(raw_path, **arrays)
    meta = {
        "schema": "r10_truth_audit_trajectory_result_v1",
        "created_utc": base.utc_now(),
        "config": config,
        "config_sha256": config_sha,
        "status": status,
        "event": event,
        "failure_message": None,
        "telemetry": telemetry,
        "raw_path": str(raw_path.relative_to(ROOT)),
        "raw_sha256": base.file_hash(raw_path),
        "n_output_epochs": int(len(times)),
        "last_output_epoch_s": float(times[-1]),
    }
    base.atomic_json(sidecar_path, meta)
    print(
        f"  N=900 {status:14s} wall {telemetry['total_wall_ns']/1e9:8.1f}s "
        f"rhs {telemetry['n_rhs']}",
        flush=True,
    )
    return meta, times, states


def audit_orbit(run_kind: str, orbit: dict, model, args, duration: float,
                timing_comparable: bool) -> dict:
    n600_meta, n600_t, n600_y = load_validated_baseline_trajectory(
        orbit["sobol_index"], "truth"
    )
    n900_meta, n900_t, n900_y = run_n900(
        run_kind, orbit, model, args, duration, timing_comparable
    )
    truth_difference = base.common_error(n600_t, n600_y, n900_t, n900_y)
    policy_errors_n900 = {}
    policy_sources = {}
    for policy in POLICIES:
        meta, times, states = load_validated_baseline_trajectory(
            orbit["sobol_index"], policy
        )
        policy_errors_n900[policy] = base.common_error(
            times, states, n900_t, n900_y
        )
        policy_sources[policy] = {
            "config_sha256": meta["config_sha256"],
            "raw_sha256": meta["raw_sha256"],
        }
    smallest_policy = min(
        value["pos_rms_m"] for value in policy_errors_n900.values()
    )
    threshold = min(5.0, 0.05 * smallest_policy)
    full_arc = n600_meta["status"] == n900_meta["status"] == "complete"
    passes = bool(
        full_arc and truth_difference["pos_rms_m"] < threshold
    )
    empirical = policy_errors_n900["schedule_empirical"]["pos_rms_m"]
    ratios_n900 = {
        "rho_work": (
            policy_errors_n900["fixed_work"]["pos_rms_m"] / empirical
        ),
        "rho_crit": (
            policy_errors_n900["fixed_critical"]["pos_rms_m"] / empirical
        ),
    }
    return {
        "name": orbit["name"],
        "sobol_index": orbit["sobol_index"],
        "hp_km": orbit["hp_km"],
        "ha_km": orbit["ha_km"],
        "incl_deg": orbit["incl_deg"],
        "N600_source": {
            "config_sha256": n600_meta["config_sha256"],
            "raw_sha256": n600_meta["raw_sha256"],
            "status": n600_meta["status"],
        },
        "N900_source": {
            "config_sha256": n900_meta["config_sha256"],
            "raw_sha256": n900_meta["raw_sha256"],
            "status": n900_meta["status"],
            "telemetry": n900_meta["telemetry"],
        },
        "N600_to_N900": truth_difference,
        "policy_errors_against_N900": policy_errors_n900,
        "policy_sources": policy_sources,
        "smallest_interpreted_policy_rms_m": smallest_policy,
        "acceptance_threshold_m": threshold,
        "acceptance_rule": (
            "E600_900 < min(5 m, 0.05*smallest interpreted policy RMS vs N900)"
        ),
        "passes": passes,
        "adopted_truth_degree": 600 if passes else 900,
        "primary_ratios_against_N900": ratios_n900,
    }


def write_index(path: Path, rows: list[dict], complete: bool,
                timing_comparable: bool, started_utc: str,
                wall_start: int, cpu_start: int) -> None:
    payload = {
        "schema": "r10_sobolA_truth_audit_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "design_sha256": base.load_design_a()["design_sha256"],
        "baseline_sha256": base.file_hash(BASELINE_PATH),
        "audit_script_sha256": base.file_hash(Path(__file__).resolve()),
        "audit_indices": list(AUDIT_INDICES),
        "started_utc": started_utc,
        "updated_utc": base.utc_now(),
        "complete": complete,
        "timing_comparable": timing_comparable,
        "rows": rows,
        "summary": {
            "completed_cases": len(rows),
            "passes": sum(row.get("passes") is True for row in rows),
            "failures": sum(row.get("passes") is False for row in rows),
            "numerical_failures": sum(
                row.get("batch_status") == "numerical_failure" for row in rows
            ),
            "adopted_N900_indices": [
                row["sobol_index"] for row in rows
                if row.get("adopted_truth_degree") == 900
            ],
        },
        "current_session_wall_ns": time.perf_counter_ns() - wall_start,
        "current_session_cpu_ns": time.process_time_ns() - cpu_start,
    }
    if complete:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(path, payload)


def run(smoke: bool) -> int:
    baseline = load_baseline()
    indices = derive_audit_indices(baseline)
    other = base.other_python_processes()
    if other:
        raise RuntimeError(
            "other Python processes are active; refusing timing-comparable audit: "
            + json.dumps(other)
        )
    timing_comparable = True
    output_path = SMOKE_PATH if smoke else OUTPUT_PATH
    run_kind = "truth_audit_smoke" if smoke else "truth_audit"
    duration = 2.0 * 3600.0 if smoke else base.DURATION_S
    run_indices = indices[:1] if smoke else indices
    existing = {}
    if output_path.exists():
        previous = json.loads(output_path.read_text(encoding="utf-8"))
        existing = {row["sobol_index"]: row for row in previous.get("rows", [])}
    design = base.load_design_a()
    orbit_by_index = {orbit["sobol_index"]: orbit for orbit in design["orbits"]}
    model = base.load_model(900)
    args = base.kernel_args(model)
    base.warmup(model, args)
    started_utc = base.utc_now()
    wall_start = time.perf_counter_ns()
    cpu_start = time.process_time_ns()
    print(
        f"[{run_kind}] {len(run_indices)} cases, timing_comparable=True",
        flush=True,
    )
    for position, index in enumerate(run_indices, start=1):
        current_other = base.other_python_processes()
        if current_other:
            ordered = [existing[key] for key in sorted(existing)]
            write_index(output_path, ordered, False, False, started_utc,
                        wall_start, cpu_start)
            raise RuntimeError("another Python process appeared during audit")
        orbit = orbit_by_index[index]
        print(
            f"[{position}/{len(run_indices)}] {orbit['name']} "
            f"hp={orbit['hp_km']:.3f} km",
            flush=True,
        )
        try:
            row = audit_orbit(
                run_kind, orbit, model, args, duration, timing_comparable
            )
            row["batch_status"] = "completed"
            print(
                f"  E600-900={row['N600_to_N900']['pos_rms_m']:.6f} m, "
                f"limit={row['acceptance_threshold_m']:.6f} m, "
                f"passes={row['passes']}",
                flush=True,
            )
        except Exception as exc:
            row = {
                "name": orbit["name"],
                "sobol_index": index,
                "batch_status": "numerical_failure",
                "failure_message": f"{type(exc).__name__}: {exc}",
                "passes": None,
                "adopted_truth_degree": None,
            }
            print(f"  [failure preserved] {row['failure_message']}", flush=True)
        existing[index] = row
        ordered = [existing[key] for key in sorted(existing)]
        write_index(output_path, ordered, False, timing_comparable,
                    started_utc, wall_start, cpu_start)
    ordered = [existing[key] for key in sorted(existing)]
    write_index(output_path, ordered, True, timing_comparable,
                started_utc, wall_start, cpu_start)
    if ACTIVE_PATH.exists():
        ACTIVE_PATH.unlink()
    print(f"[{run_kind}] complete", flush=True)
    return 0


def status() -> int:
    if not OUTPUT_PATH.exists():
        print("no formal truth-audit index yet")
    else:
        payload = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
        print(json.dumps({
            "complete": payload.get("complete"),
            "completed_cases": len(payload.get("rows", [])),
            "planned_cases": len(payload.get("audit_indices", [])),
            "timing_comparable": payload.get("timing_comparable"),
            "summary": payload.get("summary"),
            "updated_utc": payload.get("updated_utc"),
        }, indent=2))
    if ACTIVE_PATH.exists():
        print("active:")
        print(ACTIVE_PATH.read_text(encoding="utf-8"))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("smoke", "run", "status"))
    args = parser.parse_args()
    if args.command == "smoke":
        return run(True)
    if args.command == "run":
        return run(False)
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
