"""28-day LRO-like corrected-blend convergence at VECTOR tolerance (R11).

Why this run exists
-------------------
``revision/REVIEW_FIX_REPORT_2026-07-22.md`` lists "a dedicated tighter-tolerance
rerun of the 28-day LRO-like corrected blend" as deliberately deferred, and the
abstract consequently qualifies that comparison as *nominal*.  The R10 batch
used scalar tolerances (``rtol 1e-11/atol 1e-6`` and ``rtol 1e-12/atol 1e-7``).
Because SciPy builds the error scale componentwise as ``atol + rtol*|y|``, and
the velocity components carry ``|v| ~ 1.6e3 m/s``, the scalar ``atol`` governs
the velocity error control; the resulting envelope produced a 154.5 m
resolution threshold against a 146.6 m error gap, leaving the comparison
unresolved.

This rerun repeats the identical experiment with position/velocity-split vector
tolerances and a tighter rtol ladder, which is what actually lowers the floor.
Geometry, force models, blend definition, duration, output grid, error metric,
and the resolution rule are unchanged and are imported from
``rev10_blend_lro_convergence`` so the two campaigns differ *only* in the
tolerance contract.

Usage
-----
    python rev11_blend_lro_vector.py run --deadline 2026-07-24T15:00:00+03:00
    python rev11_blend_lro_vector.py smoke
    python rev11_blend_lro_vector.py status
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev10_blend_lro_convergence as blend
from rev7_doe_screening import initial_state
from rev9_potential_blend_longarc import BlendRhs

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r11_blend_lro_vector.json"
SMOKE_OUTPUT = METRICS / "r11_blend_lro_vector_smoke.json"
ACTIVE = METRICS / "r11_blend_lro_vector_active.json"
CASE_ROOT = METRICS / "r11_cases" / "blend_lro_vector"
RAW_ROOT = METRICS / "r11_raw" / "blend_lro_vector"

DURATION = blend.DURATION
OUTPUT_STEP = blend.OUTPUT_STEP
MAX_STEP = blend.MAX_STEP

# The only substantive change against R10: split position/velocity tolerances.
LEVELS = {
    "baseline": {"rtol": 1.0e-12, "atol_position_m": 1.0e-5,
                 "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol_position_m": 1.0e-6,
                "atol_velocity_m_s": 1.0e-9},
}
# Heavy-first; guards scaled from the measured R10 walls (the vector ladder
# costs roughly 3-4x the scalar one per trajectory).
RUN_ORDER = (
    ("corrected_blend", "tighter", 42000.0),
    ("corrected_blend", "baseline", 30000.0),
    ("truth_N600", "tighter", 12000.0),
    ("truth_N600", "baseline", 8000.0),
    ("fixed_N120", "tighter", 1500.0),
    ("fixed_N120", "baseline", 1000.0),
)


def atol_vector(level: str) -> np.ndarray:
    tol = LEVELS[level]
    return np.array([tol["atol_position_m"]] * 3 +
                    [tol["atol_velocity_m_s"]] * 3)


def paths(policy: str, level: str, smoke: bool) -> tuple[Path, Path]:
    suffix = "_smoke" if smoke else ""
    return (CASE_ROOT / f"{policy}_{level}{suffix}.json",
            RAW_ROOT / f"{policy}_{level}{suffix}.npz")


def build_config(policy: str, level: str, y0, duration: float) -> dict:
    tol = LEVELS[level]
    protocol = base.protocol_payload()
    return {
        "schema": "r11_blend_lro_vector_config_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "reused_module_sha256": base.file_hash(
            Path(blend.__file__).resolve()),
        "orbit": blend.orbit_dict(),
        "initial_state_si": [float(x) for x in y0],
        "policy": policy, "level": level,
        "duration_s": duration, "output_step_s": OUTPUT_STEP,
        "integrator": "InstrumentedDOP853",
        "rtol": tol["rtol"],
        "atol_kind": "vector",
        "atol_position_m": tol["atol_position_m"],
        "atol_velocity_m_s": tol["atol_velocity_m_s"],
        "max_step_s": MAX_STEP,
        "truth_degree": 600, "fixed_degree": 120,
        "blend_degrees": [30, 120],
        "transition_altitude_m": [50000.0, 200000.0],
        "corrected_term": "(U_hi-U_lo)*dw/dr*r_hat",
        "orientation": "uniform lunar sidereal rotation",
        "timing_comparable": False,
        "timing_note": ("accuracy rerun executed alongside other R11 "
                        "campaigns; cost claims use the serial R10 baseline"),
        "supersedes": "r10_blend_lro_convergence (scalar atol)",
        "source": protocol["provenance"],
    }


def run_case(policy: str, level: str, model, args, y0, duration: float,
             smoke: bool):
    config = build_config(policy, level, y0, duration)
    config_sha = base.object_hash(config)
    sidecar, raw = paths(policy, level, smoke)
    if sidecar.exists() and raw.exists() and base.valid_cached(
            sidecar, raw, config_sha, duration):
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        print(f"  {policy}/{level} cached {meta['status']}", flush=True)
        return meta
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
    base.atomic_json(ACTIVE, {"updated_utc": base.utc_now(), "policy": policy,
                              "level": level, "config_sha256": config_sha})
    grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    t, y, status, event, failure, telemetry = blend.generic_propagate(
        rhs, model, y0, duration, grid, LEVELS[level]["rtol"],
        atol_vector(level))
    if status == "numerical_failure":
        base.atomic_json(
            sidecar.with_name(sidecar.stem + f".failure.{int(time.time())}.json"),
            {"config": config, "config_sha256": config_sha, "status": status,
             "failure_message": failure, "telemetry": telemetry})
        raise RuntimeError(f"{policy}/{level}: {failure}")
    arrays = {"t_s": t, "state_si": y}
    if event:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"])
    base.atomic_npz(raw, **arrays)
    meta = {
        "schema": "r11_blend_lro_vector_result_v1",
        "created_utc": base.utc_now(), "config": config,
        "config_sha256": config_sha, "status": status, "event": event,
        "failure_message": None, "telemetry": telemetry,
        "raw_path": str(raw.relative_to(ROOT)),
        "raw_sha256": base.file_hash(raw),
        "n_output_epochs": int(len(t)), "last_output_epoch_s": float(t[-1])}
    base.atomic_json(sidecar, meta)
    print(f"  {policy}/{level} {status} "
          f"wall={telemetry['total_wall_ns']/1e9:.1f}s", flush=True)
    return meta


def compute_summary(smoke: bool) -> dict | None:
    keys = [(p, l) for p in ("truth_N600", "fixed_N120", "corrected_blend")
            for l in LEVELS]
    data = {}
    for policy, level in keys:
        sidecar, raw = paths(policy, level, smoke)
        if not (sidecar.exists() and raw.exists()):
            return None
        t, y = base.load_raw(raw)
        data[(policy, level)] = (t, y)
    truth_self = base.common_error(
        *data[("truth_N600", "baseline")], *data[("truth_N600", "tighter")]
    )["pos_rms_m"]
    result = {"truth_self_difference_rms_m": truth_self, "policies": {}}
    for policy in ("fixed_N120", "corrected_blend"):
        errors = {level: base.common_error(*data[(policy, level)],
                                           *data[("truth_N600", level)])
                  for level in LEVELS}
        self_diff = base.common_error(*data[(policy, "baseline")],
                                      *data[(policy, "tighter")])["pos_rms_m"]
        result["policies"][policy] = {
            "error_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self}
    fixed = result["policies"]["fixed_N120"]
    blended = result["policies"]["corrected_blend"]
    e_blend = blended["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
    e_fixed = fixed["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]
    difference = abs(e_blend - e_fixed)
    threshold = (blended["truth_inclusive_envelope_m"] +
                 fixed["truth_inclusive_envelope_m"])
    result["comparison"] = {
        "corrected_blend_error_m": e_blend, "fixed_N120_error_m": e_fixed,
        "absolute_baseline_error_difference_m": difference,
        "resolution_threshold_m": threshold,
        "resolved": bool(difference > threshold),
        "winner_if_resolved": (("corrected_blend" if e_blend < e_fixed
                                else "fixed_N120")
                               if difference > threshold else None)}
    result["r10_scalar_comparison_for_reference"] = {
        "absolute_baseline_error_difference_m": 146.6,
        "resolution_threshold_m": 154.5, "resolved": False}
    return result


def write_index(records, complete, skipped, started, wall0, deadline, smoke):
    target = SMOKE_OUTPUT if smoke else OUTPUT
    payload = {
        "schema": "r11_blend_lro_vector_index_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "started_utc": started, "updated_utc": base.utc_now(),
        "deadline_utc": deadline.isoformat() if deadline else None,
        "complete": complete, "timing_comparable": False,
        "levels": LEVELS,
        "records": [records[key] for key in sorted(records)],
        "skipped_for_deadline": skipped,
        "summary": None if smoke else compute_summary(smoke),
        "session_wall_ns": time.perf_counter_ns() - wall0}
    if complete:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(target, payload)


def run(smoke: bool, deadline) -> int:
    model = base.load_model(600)
    args = base.kernel_args(model)
    base.warmup(model, args)
    y0 = initial_state(model, blend.orbit_dict())
    duration = 2.0 * 3600.0 if smoke else DURATION
    order = (("corrected_blend", "tighter", 600.0),
             ("corrected_blend", "baseline", 600.0),
             ("truth_N600", "tighter", 600.0),
             ("truth_N600", "baseline", 600.0),
             ("fixed_N120", "tighter", 600.0),
             ("fixed_N120", "baseline", 600.0)) if smoke else RUN_ORDER
    records, skipped = {}, []
    target = SMOKE_OUTPUT if smoke else OUTPUT
    if target.exists():
        old = json.loads(target.read_text(encoding="utf-8"))
        records = {(r["policy"], r["level"]): r for r in old.get("records", [])}
    started = base.utc_now()
    wall0 = time.perf_counter_ns()
    print(f"[r11-blend] cases={len(order)} deadline={deadline}", flush=True)
    for policy, level, estimate in order:
        left = blend.remaining_seconds(deadline)
        if left < estimate + 600.0:
            skipped.append({"policy": policy, "level": level,
                            "remaining_s": left, "guard_estimate_s": estimate})
            print(f"  skip {policy}/{level}: deadline guard", flush=True)
            continue
        meta = run_case(policy, level, model, args, y0, duration, smoke)
        records[(policy, level)] = {
            "policy": policy, "level": level, "status": meta["status"],
            "config_sha256": meta["config_sha256"],
            "raw_path": meta["raw_path"], "raw_sha256": meta["raw_sha256"],
            "telemetry": meta["telemetry"]}
        write_index(records, False, skipped, started, wall0, deadline, smoke)
    complete = len(records) == len(order)
    write_index(records, complete, skipped, started, wall0, deadline, smoke)
    if ACTIVE.exists():
        ACTIVE.unlink()
    print(f"[r11-blend] finished complete={complete} "
          f"records={len(records)}/{len(order)}", flush=True)
    return 0 if complete else 3


def run_one(policy: str, level: str, smoke: bool) -> int:
    """Run a single (policy, level) case.

    The six cases are independent, so the campaign is driven as six concurrent
    single-case processes; ``corrected_blend`` dominates the wall clock and the
    other five finish early and release their cores.
    """
    model = base.load_model(600)
    args = base.kernel_args(model)
    base.warmup(model, args)
    y0 = initial_state(model, blend.orbit_dict())
    duration = 2.0 * 3600.0 if smoke else DURATION
    print(f"[r11-blend-one] {policy}/{level}", flush=True)
    run_case(policy, level, model, args, y0, duration, smoke)
    return 0


def collect(smoke: bool) -> int:
    """Assemble the index from whatever single-case artifacts exist."""
    records = {}
    for policy, level, _ in RUN_ORDER:
        sidecar, _ = paths(policy, level, smoke)
        if not sidecar.exists():
            continue
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        records[(policy, level)] = {
            "policy": policy, "level": level, "status": meta["status"],
            "config_sha256": meta["config_sha256"],
            "raw_path": meta["raw_path"], "raw_sha256": meta["raw_sha256"],
            "telemetry": meta["telemetry"]}
    complete = len(records) == len(RUN_ORDER)
    write_index(records, complete, [], base.utc_now(),
                time.perf_counter_ns(), None, smoke)
    print(f"[r11-blend] collected {len(records)}/{len(RUN_ORDER)} "
          f"complete={complete}", flush=True)
    return 0 if complete else 3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command",
                        choices=("run", "smoke", "status", "one", "collect"))
    parser.add_argument("--deadline")
    parser.add_argument("--policy")
    parser.add_argument("--level")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.command == "one":
        return run_one(args.policy, args.level, args.smoke)
    if args.command == "collect":
        return collect(args.smoke)
    if args.command == "status":
        if OUTPUT.exists():
            data = json.loads(OUTPUT.read_text(encoding="utf-8"))
            print(json.dumps({"complete": data["complete"],
                              "records": len(data["records"]),
                              "summary": data["summary"]}, indent=2))
        else:
            print("no r11 blend output")
        return 0
    return run(args.command == "smoke", blend.parse_deadline(args.deadline))


if __name__ == "__main__":
    raise SystemExit(main())
