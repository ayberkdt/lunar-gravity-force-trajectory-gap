"""Initial-phase sweep of the 7-day scheduling penalty, at vector tolerance (R11).

Why this run exists
-------------------
``rev6_phase_mc.py`` sampled eight initial true anomalies at the scalar
tolerance ``rtol 1e-12 / atol 1e-5`` and reported policy errors of 5-18 m.
Under that setting the seven-day integration envelope is itself several meters
(see ``metrics/robustness_numerical_floor_check.json``), so most of the
reported spread was not separable from integration noise, and eight samples is
thin for a distributional claim.  ``REVIEW_FIX_REPORT_2026-07-22`` lists the
repeat "with the final vector-tolerance setup" as deferred.

This run replaces it: 24 initial phases instead of 8, position/velocity-split
vector tolerances at two levels, and a per-phase numerical envelope so every
phase carries its own resolved/unresolved verdict rather than a bare ratio.

Contract
--------
Same geometry (50x300 km polar), same truth degree (300), same policies
(fixed 138, fixed 106, downward and upward schedules) and the same seven-day
horizon as R6-C, so the numbers are directly comparable with the archived ones.
The tolerance ladder matches ``rev11_full_convergence`` exactly.

Usage
-----
    python rev11_phase_sweep.py run --workers 5
    python rev11_phase_sweep.py smoke
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev3_common import DAY, alt_sched, make_p_table
from rev6_phase_mc import state_at_anomaly

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r11_phase_sweep.json"
SMOKE_OUTPUT = METRICS / "r11_phase_sweep_smoke.json"
CASE_ROOT = METRICS / "r11_cases" / "phase_sweep"
RAW_ROOT = METRICS / "r11_raw" / "phase_sweep"

N_PHASES = 24
PHASES_DEG = [round(360.0 * i / N_PHASES, 4) for i in range(N_PHASES)]
HP_KM, HA_KM = 50.0, 300.0
TRUTH_DEGREE = 300
DURATION = 7.0 * DAY
OUTPUT_STEP = 120.0
MAX_STEP = 60.0
LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol_position_m": 1.0e-5,
              "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol_position_m": 1.0e-6,
                "atol_velocity_m_s": 1.0e-9},
}
POLICIES = ("truth", "fixed_138", "fixed_106", "sched_down", "sched_up")
COMPARED = ("fixed_138", "fixed_106", "sched_down", "sched_up")

_MODEL: dict = {}


def model_args():
    if "m" not in _MODEL:
        model = base.load_model(TRUTH_DEGREE)
        args = base.kernel_args(model)
        base.warmup(model, args)
        _MODEL["m"] = (model, args)
        _MODEL["down"] = alt_sched(make_p_table(model, 1e-3, 60, policy="down"))
        _MODEL["up"] = alt_sched(make_p_table(model, 1e-3, 60, policy="up"))
    return _MODEL["m"]


def degree_function(policy: str):
    model, _ = model_args()
    if policy == "truth":
        return (lambda t, h: TRUTH_DEGREE), {"kind": "fixed_truth",
                                             "degree": TRUTH_DEGREE}
    if policy == "fixed_138":
        return (lambda t, h: 138), {"kind": "fixed", "degree": 138}
    if policy == "fixed_106":
        return (lambda t, h: 106), {"kind": "fixed", "degree": 106}
    if policy == "sched_down":
        return _MODEL["down"], {"kind": "kaula_downward_quantized",
                                "epsilon": 1e-3, "floor": 60}
    if policy == "sched_up":
        return _MODEL["up"], {"kind": "kaula_upward_quantized",
                              "epsilon": 1e-3, "floor": 60}
    raise ValueError(policy)


def atol_vector(level: str) -> np.ndarray:
    tol = LEVELS[level]
    return np.array([tol["atol_position_m"]] * 3 +
                    [tol["atol_velocity_m_s"]] * 3)


def paths(phase_index: int, policy: str, level: str, smoke: bool):
    suffix = "_smoke" if smoke else ""
    return (CASE_ROOT / f"phase_{phase_index:02d}" / f"{policy}_{level}{suffix}.json",
            RAW_ROOT / f"phase_{phase_index:02d}" / f"{policy}_{level}{suffix}.npz")


def worker(task: dict) -> dict:
    index, policy, level = task["index"], task["policy"], task["level"]
    nu, duration, smoke = task["nu_deg"], task["duration"], task["smoke"]
    sidecar, raw = paths(index, policy, level, smoke)
    try:
        model, args = model_args()
        degree_of, spec = degree_function(policy)
        y0 = state_at_anomaly(model, HP_KM, HA_KM, nu)
        tol = LEVELS[level]
        config = {
            "schema": "r11_phase_sweep_config_v1",
            "protocol_sha256": base.protocol_payload()["protocol_sha256"],
            "script_sha256": task["script_sha"],
            "phase_index": index, "nu_deg": nu,
            "geometry": {"hp_km": HP_KM, "ha_km": HA_KM,
                         "plane": "polar (x-z)"},
            "initial_state_si": [float(x) for x in y0],
            "truth_degree": TRUTH_DEGREE,
            "policy": policy, "policy_spec": spec, "level": level,
            "duration_s": duration, "output_step_s": OUTPUT_STEP,
            "integrator": "InstrumentedDOP853",
            "rtol": tol["rtol"], "atol_kind": "vector",
            "atol_position_m": tol["atol_position_m"],
            "atol_velocity_m_s": tol["atol_velocity_m_s"],
            "max_step_s": MAX_STEP,
            "timing_comparable": False,
            "supersedes": "r6_phase_mc (8 phases, scalar atol)",
            "source": base.provenance(),
        }
        config_sha = base.object_hash(config)
        if sidecar.exists() and raw.exists() and base.valid_cached(
                sidecar, raw, config_sha, duration):
            return {"index": index, "policy": policy, "level": level,
                    "status": "cached", "wall_s": 0.0}
        if sidecar.exists() or raw.exists():
            base.preserve_invalid(sidecar)
            base.preserve_invalid(raw)
        grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        t, y, status, event, failure, telemetry = \
            base.propagate_event_instrumented(
                model, y0, duration, grid, degree_of, args, tol["rtol"],
                atol_vector(level), max_step=MAX_STEP)
        if status == "numerical_failure":
            return {"index": index, "policy": policy, "level": level,
                    "status": status, "message": failure, "wall_s": 0.0}
        arrays = {"t_s": t, "state_si": y}
        if event:
            arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
            arrays["impact_state_si"] = np.asarray(event["state_si"])
        base.atomic_npz(raw, **arrays)
        base.atomic_json(sidecar, {
            "schema": "r11_phase_sweep_result_v1",
            "created_utc": base.utc_now(), "config": config,
            "config_sha256": config_sha, "status": status, "event": event,
            "failure_message": None, "telemetry": telemetry,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
        return {"index": index, "policy": policy, "level": level,
                "status": status, "wall_s": telemetry["total_wall_ns"] / 1e9}
    except Exception as exc:
        return {"index": index, "policy": policy, "level": level,
                "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(), "wall_s": 0.0}


def phase_summary(index: int, nu: float, smoke: bool) -> dict:
    data = {}
    for policy in POLICIES:
        for level in LEVELS:
            sidecar, raw = paths(index, policy, level, smoke)
            t, y = base.load_raw(raw)
            data[(policy, level)] = (t, y)
    truth_self = base.common_error(*data[("truth", "tight")],
                                   *data[("truth", "tighter")])["pos_rms_m"]
    policies = {}
    for policy in COMPARED:
        errors = {level: base.common_error(*data[(policy, level)],
                                           *data[("truth", level)])
                  for level in LEVELS}
        self_diff = base.common_error(*data[(policy, "tight")],
                                      *data[(policy, "tighter")])["pos_rms_m"]
        policies[policy] = {
            "errors_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self}
    reference = policies["fixed_138"]
    e_ref = reference["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
    penalties = {}
    for policy in ("sched_down", "sched_up", "fixed_106"):
        entry = policies[policy]
        e = entry["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
        difference = abs(e - e_ref)
        threshold = (entry["truth_inclusive_envelope_m"] +
                     reference["truth_inclusive_envelope_m"])
        penalties[policy] = {
            "error_m": e, "reference_error_m": e_ref,
            "penalty_vs_fixed_138": e / e_ref if e_ref > 0 else None,
            "absolute_error_difference_m": difference,
            "resolution_threshold_m": threshold,
            "resolved": bool(difference > threshold),
            "winner_if_resolved": ((policy if e < e_ref else "fixed_138")
                                   if difference > threshold else None)}
    return {"phase_index": index, "nu_deg": nu,
            "truth_self_difference_rms_m": truth_self,
            "policies": policies, "penalties": penalties}


def summarize(rows: list[dict]) -> dict:
    def stat(values):
        arr = np.asarray([v for v in values if v is not None and np.isfinite(v)])
        if arr.size == 0:
            return None
        return {"n": int(arr.size), "median": float(np.median(arr)),
                "p10": float(np.percentile(arr, 10)),
                "p90": float(np.percentile(arr, 90)),
                "min": float(arr.min()), "max": float(arr.max())}

    out = {"phases": len(rows),
           "truth_self_difference_rms_m": stat(
               [r["truth_self_difference_rms_m"] for r in rows])}
    out["error_by_policy_m"] = {
        p: stat([r["policies"][p]["errors_against_same_tolerance_truth"]
                 ["tight"]["pos_rms_m"] for r in rows]) for p in COMPARED}
    out["self_difference_by_policy_m"] = {
        p: stat([r["policies"][p]["self_difference_rms_m"] for r in rows])
        for p in COMPARED}
    out["penalty_vs_fixed_138"] = {}
    for policy in ("sched_down", "sched_up", "fixed_106"):
        cs = [r["penalties"][policy] for r in rows]
        out["penalty_vs_fixed_138"][policy] = {
            "ratio": stat([c["penalty_vs_fixed_138"] for c in cs]),
            "resolved": sum(c["resolved"] for c in cs),
            "unresolved": sum(not c["resolved"] for c in cs),
            "resolved_fixed138_wins": sum(
                c["winner_if_resolved"] == "fixed_138" for c in cs),
            "resolved_policy_wins": sum(
                c["winner_if_resolved"] == policy for c in cs)}
    return out


def run(smoke: bool, workers: int) -> int:
    phases = PHASES_DEG[:2] if smoke else PHASES_DEG
    duration = 2.0 * 3600.0 if smoke else DURATION
    script_sha = base.file_hash(Path(__file__).resolve())
    for index in range(len(phases)):
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / f"phase_{index:02d}").mkdir(parents=True, exist_ok=True)
    tasks = [{"index": i, "nu_deg": nu, "policy": p, "level": l,
              "duration": duration, "smoke": smoke, "script_sha": script_sha}
             for i, nu in enumerate(phases) for p in POLICIES for l in LEVELS]
    started = base.utc_now()
    wall0 = time.perf_counter_ns()
    print(f"[r11-phase] phases={len(phases)} trajectories={len(tasks)} "
          f"workers={workers}", flush=True)
    failures = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(worker, t) for t in tasks]
        for n, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            if record["status"] in ("numerical_failure", "worker_error"):
                failures.append(record)
                print(f"  !! phase {record['index']} {record['policy']}/"
                      f"{record['level']}: {record.get('message')}", flush=True)
            if n % 10 == 0 or n == len(tasks):
                elapsed = time.time() - t0
                print(f"  [{n:3d}/{len(tasks)}] elapsed={elapsed/60:5.1f} min "
                      f"eta={(len(tasks)-n)*elapsed/n/60:5.1f} min", flush=True)
    rows = []
    for index, nu in enumerate(phases):
        try:
            rows.append(phase_summary(index, nu, smoke))
        except Exception as exc:
            failures.append({"index": index, "policy": "summary",
                             "status": "summary_incomplete",
                             "message": f"{type(exc).__name__}: {exc}"})
    payload = {
        "schema": "r11_phase_sweep_index_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "script_sha256": script_sha,
        "started_utc": started, "ended_utc": base.utc_now(),
        "complete": len(rows) == len(phases) and not failures,
        "timing_comparable": False,
        "scenario": {"type": "eccentric_polar_50x300",
                     "hp_km": HP_KM, "ha_km": HA_KM,
                     "duration_s": duration, "truth_degree": TRUTH_DEGREE,
                     "integrator": "DOP853", "levels": LEVELS,
                     "output_step_s": OUTPUT_STEP, "max_step_s": MAX_STEP,
                     "phases_deg": phases,
                     "rotation": "uniform sidereal about polar axis"},
        "rows": rows, "failures": failures, "summary": summarize(rows),
        "session_wall_ns": time.perf_counter_ns() - wall0}
    base.atomic_json(SMOKE_OUTPUT if smoke else OUTPUT, payload)
    print(f"[r11-phase] finished phases={len(rows)}/{len(phases)} "
          f"failures={len(failures)}", flush=True)
    return 0 if payload["complete"] else 3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("run", "smoke", "status"))
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()
    if args.command == "status":
        if OUTPUT.exists():
            data = json.loads(OUTPUT.read_text(encoding="utf-8"))
            print(json.dumps({"complete": data["complete"],
                              "phases": len(data["rows"]),
                              "summary": data["summary"]}, indent=2))
        else:
            print("no r11 phase sweep output")
        return 0
    return run(args.command == "smoke", args.workers)


if __name__ == "__main__":
    raise SystemExit(main())
