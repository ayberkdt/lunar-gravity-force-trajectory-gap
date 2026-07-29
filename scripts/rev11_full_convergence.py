"""Full-population vector-tolerance convergence for the Sobol A design (R11).

Motivation
----------
The R10 confirmatory campaign integrated the 64-orbit population with a
*scalar* tolerance ``atol = 1e-5``.  Because SciPy forms the error scale as
``atol + rtol*|y|`` componentwise, and the velocity components carry
``|v| ~ 1.6e3 m/s``, the velocity error control was effectively
``atol = 1e-5 m/s`` (``rtol*|v| = 1.6e-9`` is negligible against it).  A
persistent along-track velocity error of that size integrates to meters of
position drift per week, which is the same order as the policy differences the
study reports.  That, and not the absence of geometric structure in DOP853, is
the origin of the reported numerical floor: a two-day control on
``c1_circ100_polar`` gives a 0.854 m self-difference under the scalar setting
against 0.0075 m under a vector setting, a 114x reduction for 1.9x the cost.

R10 already demonstrated the fix on a selected 17-orbit subset
(``rev10_sobol_convergence.py``).  This script applies it to the *entire*
64-orbit population and to all five reported policies, so that every ratio in
the manuscript carries a per-case resolved verdict instead of a subset audit.

Contract
--------
* Two tolerance levels, identical to the frozen R10 convergence contract:
  ``tight``   rtol 1e-12, atol (1e-5 m, 1e-8 m/s), max_step 60 s
  ``tighter`` rtol 1e-13, atol (1e-6 m, 1e-9 m/s), max_step 60 s
* Policies: truth, schedule_empirical, schedule_up, schedule_down,
  fixed_critical, fixed_work.
* 64 orbits x 6 policies x 2 levels = 768 seven-day trajectories.
* Errors are always taken against the truth at the *same* tolerance level.
* The per-policy self-difference between levels is the numerical envelope; a
  comparison is resolved only when the error gap exceeds the sum of the two
  truth-inclusive envelopes.  This is the R10 rule, unchanged.

Timing
------
Trajectories run in a process pool, so wall clocks here are NOT comparable and
``timing_comparable`` is recorded false in every artifact.  All cost claims in
the manuscript continue to come from the serial R10 baseline; this campaign is
an accuracy campaign only.  Kernel-work proxies (``mean_degree_sq``,
``gravity_kernel_ns``) are still recorded per trajectory for cross-checking.

Usage
-----
    python rev11_full_convergence.py run --workers 5 --deadline 2026-07-24T13:00:00+03:00
    python rev11_full_convergence.py smoke --workers 3
    python rev11_full_convergence.py status
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev7_doe_screening import alt_sched, degree_power, emp_table, kaula_table


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

# Path configuration is read from the environment at IMPORT time so that it is
# honoured inside ProcessPoolExecutor workers, which (on Windows spawn) re-import
# this module fresh and do NOT inherit any monkey-patched module globals from the
# parent.  Environment variables ARE inherited by spawned children, so a driver
# that runs a different population (e.g. design B) sets these before creating the
# pool.  With no environment set, the design-A defaults below apply unchanged.
_TREE = os.environ.get("R11_TREE", "convergence")
CORRECTED = Path(os.environ.get(
    "R11_CORRECTED", METRICS / "r10_sobolA_baseline_truth_corrected.json"))
OUTPUT = Path(os.environ.get("R11_OUTPUT", METRICS / "r11_full_convergence.json"))
SMOKE_OUTPUT = Path(os.environ.get(
    "R11_SMOKE_OUTPUT", METRICS / "r11_full_convergence_smoke.json"))
CASE_ROOT = METRICS / "r11_cases" / _TREE
RAW_ROOT = METRICS / "r11_raw" / _TREE

LEVELS = {
    "tight": {
        "rtol": 1.0e-12,
        "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3),
        "atol_position_m": 1.0e-5,
        "atol_velocity_m_s": 1.0e-8,
    },
    "tighter": {
        "rtol": 1.0e-13,
        "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3),
        "atol_position_m": 1.0e-6,
        "atol_velocity_m_s": 1.0e-9,
    },
}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
POLICIES = ("truth", "schedule_empirical", "schedule_up", "schedule_down",
            "fixed_critical", "fixed_work")
COMPARED = ("schedule_empirical", "schedule_up", "schedule_down",
            "fixed_critical", "fixed_work")

# ----------------------------------------------------------------- worker side
_MODELS: dict[int, tuple] = {}
_TABLES: dict[tuple[str, int], dict] = {}


def _model(degree: int):
    if degree not in _MODELS:
        model = base.load_model(degree)
        args = base.kernel_args(model)
        base.warmup(model, args)
        _MODELS[degree] = (model, args)
    return _MODELS[degree]


def _table(kind: str, degree: int):
    key = (kind, degree)
    if key not in _TABLES:
        model, _ = _model(degree)
        if kind == "emp":
            _TABLES[key] = emp_table(model, degree_power(model))
        else:
            _TABLES[key] = kaula_table(model, kind)
    return _TABLES[key]


def paths(index: int, policy: str, level: str, smoke: bool,
          case_root=None, raw_root=None):
    # Roots may be passed explicitly (carried inside the pickled task dict) so a
    # spawned worker writes to the correct population's tree regardless of the
    # module globals it re-imported.  Falls back to the module globals, which is
    # what the design-A parent-side calls (read_meta) use.
    croot = Path(case_root) if case_root is not None else CASE_ROOT
    rroot = Path(raw_root) if raw_root is not None else RAW_ROOT
    suffix = "_smoke" if smoke else ""
    return (croot / f"sobolA_{index:03d}" / f"{policy}_{level}{suffix}.json",
            rroot / f"sobolA_{index:03d}" / f"{policy}_{level}{suffix}.npz")


def degree_function(row: dict, policy: str):
    """Return (degree_of, policy_spec) for one policy on one orbit."""
    original = int(row["original_truth_degree"])
    if policy == "truth":
        degree = int(row["adopted_truth_degree"])
        return (lambda t, h, n=degree: n), {"kind": "fixed_truth", "degree": degree}
    if policy == "fixed_critical":
        degree = int(row["n_critical"])
        return (lambda t, h, n=degree: n), {
            "kind": "fixed_critical", "degree": degree,
            "definition": "unquantized empirical Nmin at own perilune, cap 250"}
    if policy == "fixed_work":
        degree = int(row["n_work"])
        return (lambda t, h, n=degree: n), {
            "kind": "fixed_work", "degree": degree,
            "source": "scalar-baseline empirical RHS history (frozen from R10)"}
    kind = {"schedule_empirical": "emp", "schedule_up": "up",
            "schedule_down": "down"}[policy]
    table = _table(kind, original)
    return alt_sched(table), {"kind": f"altitude_lookup_{kind}",
                              "lookup_source_degree": original}


def build_config(row: dict, policy: str, level: str, spec: dict,
                 duration: float, script_sha: str, protocol_sha: str,
                 corrected_sha: str) -> dict:
    tol = LEVELS[level]
    return {
        "schema": "r11_full_convergence_config_v1",
        "protocol_sha256": protocol_sha,
        "corrected_baseline_sha256": corrected_sha,
        "script_sha256": script_sha,
        "sobol_index": row["sobol_index"],
        "initial_state_si": row["design_point"]["initial_state_si"],
        "adopted_truth_degree": row["adopted_truth_degree"],
        "original_truth_degree": row["original_truth_degree"],
        "policy": policy, "policy_spec": spec, "level": level,
        "duration_s": duration, "output_step_s": OUTPUT_STEP,
        "integrator": "InstrumentedDOP853",
        "rtol": tol["rtol"],
        "atol_kind": "vector",
        "atol_position_m": tol["atol_position_m"],
        "atol_velocity_m_s": tol["atol_velocity_m_s"],
        "max_step_s": MAX_STEP,
        "timing_comparable": False,
        "execution": "parallel_process_pool",
        "timing_note": ("wall clocks are not comparable under parallel "
                        "execution; manuscript cost claims use the serial "
                        "R10 baseline"),
        "surface_event": "norm(r)-r_ref downward terminal",
        "source": base.provenance(),
    }


def worker(task: dict) -> dict:
    """Run (or reuse) one trajectory. Returns a small status record."""
    row, policy, level = task["row"], task["policy"], task["level"]
    duration, smoke = task["duration"], task["smoke"]
    index = row["sobol_index"]
    sidecar, raw = paths(index, policy, level, smoke,
                         task.get("case_root"), task.get("raw_root"))
    try:
        degree_of, spec = degree_function(row, policy)
        config = build_config(row, policy, level, spec, duration,
                              task["script_sha"], task["protocol_sha"],
                              task["corrected_sha"])
        config_sha = base.object_hash(config)
        if sidecar.exists() and raw.exists() and base.valid_cached(
                sidecar, raw, config_sha, duration):
            return {"index": index, "policy": policy, "level": level,
                    "status": "cached", "wall_s": 0.0}
        if sidecar.exists() or raw.exists():
            base.preserve_invalid(sidecar)
            base.preserve_invalid(raw)
        model, args = _model(int(row["adopted_truth_degree"]))
        tol = LEVELS[level]
        grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        t, y, status, event, failure, telemetry = \
            base.propagate_event_instrumented(
                model, np.asarray(row["design_point"]["initial_state_si"]),
                duration, grid, degree_of, args, tol["rtol"], tol["atol"],
                max_step=MAX_STEP)
        if status == "numerical_failure":
            base.atomic_json(
                sidecar.with_name(sidecar.stem + f".failure.{int(time.time())}.json"),
                {"config": config, "config_sha256": config_sha,
                 "status": status, "failure_message": failure,
                 "telemetry": telemetry})
            return {"index": index, "policy": policy, "level": level,
                    "status": "numerical_failure", "message": failure,
                    "wall_s": telemetry["total_wall_ns"] / 1e9}
        arrays = {"t_s": t, "state_si": y}
        if event:
            arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
            arrays["impact_state_si"] = np.asarray(event["state_si"])
        base.atomic_npz(raw, **arrays)
        base.atomic_json(sidecar, {
            "schema": "r11_full_convergence_result_v1",
            "created_utc": base.utc_now(), "config": config,
            "config_sha256": config_sha, "status": status, "event": event,
            "failure_message": None, "telemetry": telemetry,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
        return {"index": index, "policy": policy, "level": level,
                "status": status, "wall_s": telemetry["total_wall_ns"] / 1e9}
    except Exception as exc:  # never let one trajectory kill the campaign
        return {"index": index, "policy": policy, "level": level,
                "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(), "wall_s": 0.0}


# ----------------------------------------------------------------- parent side
def load_rows() -> list[dict]:
    payload = json.loads(CORRECTED.read_text(encoding="utf-8"))
    rows = payload["rows"]
    if len(rows) != 64:
        raise RuntimeError(f"expected 64 corrected rows, found {len(rows)}")
    return rows


def cost_estimate(row: dict) -> float:
    """Rough serial seconds for one orbit's 12 trajectories (ordering only)."""
    degree = int(row["adopted_truth_degree"])
    truth = {300: 30.0, 600: 120.0, 900: 270.0}.get(degree, 30.0)
    policy = sum(
        (max(int(row["n_critical"]), int(row["n_work"]), 140) / 300.0) ** 2 * 30.0
        for _ in range(5))
    return 3.93 * 2.0 * (truth + policy)


def parse_deadline(value: str | None):
    if not value:
        return None
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        raise ValueError("deadline needs a UTC offset")
    return result.astimezone(timezone.utc)


def remaining(deadline) -> float:
    if deadline is None:
        return math.inf
    return (deadline - datetime.now(timezone.utc)).total_seconds()


def read_meta(index: int, policy: str, level: str, smoke: bool):
    sidecar, raw = paths(index, policy, level, smoke)
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    t, y = base.load_raw(raw)
    return meta, t, y


def orbit_summary(row: dict, smoke: bool) -> dict:
    index = row["sobol_index"]
    data = {}
    for policy in POLICIES:
        for level in LEVELS:
            data[(policy, level)] = read_meta(index, policy, level, smoke)
    truth_self = base.common_error(
        data[("truth", "tight")][1], data[("truth", "tight")][2],
        data[("truth", "tighter")][1], data[("truth", "tighter")][2],
    )["pos_rms_m"]
    policies = {}
    for policy in COMPARED:
        errors = {}
        for level in LEVELS:
            errors[level] = base.common_error(
                data[(policy, level)][1], data[(policy, level)][2],
                data[("truth", level)][1], data[("truth", level)][2])
        self_diff = base.common_error(
            data[(policy, "tight")][1], data[(policy, "tight")][2],
            data[(policy, "tighter")][1], data[(policy, "tighter")][2],
        )["pos_rms_m"]
        policies[policy] = {
            "errors_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self,
            "status": data[(policy, "tight")][0]["status"],
        }
    comparisons = {}
    for schedule_name in ("schedule_empirical", "schedule_up", "schedule_down"):
        schedule = policies[schedule_name]
        es = schedule["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
        for comparator in ("fixed_work", "fixed_critical"):
            fixed = policies[comparator]
            ef = fixed["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
            difference = abs(es - ef)
            threshold = (schedule["truth_inclusive_envelope_m"] +
                         fixed["truth_inclusive_envelope_m"])
            comparisons[f"{schedule_name}_vs_{comparator}"] = {
                "rho_tight": ef / es if es > 0 else None,
                "schedule_error_m": es, "fixed_error_m": ef,
                "absolute_error_difference_m": difference,
                "resolution_threshold_m": threshold,
                "resolved": bool(difference > threshold),
                "winner_if_resolved": (
                    (schedule_name if es < ef else comparator)
                    if difference > threshold else None),
            }
    return {
        "sobol_index": index, "name": row["name"],
        "adopted_truth_degree": row["adopted_truth_degree"],
        "design_point": {k: row["design_point"][k] for k in
                         ("hp_km", "ha_km", "incl_deg", "eccentricity")},
        "n_work": row["n_work"], "n_critical": row["n_critical"],
        "truth_self_difference_rms_m": truth_self,
        "truth_status": data[("truth", "tight")][0]["status"],
        "policies": policies, "comparisons": comparisons,
        "trajectory_records": [
            {"policy": p, "level": l,
             "status": data[(p, l)][0]["status"],
             "config_sha256": data[(p, l)][0]["config_sha256"],
             "raw_sha256": data[(p, l)][0]["raw_sha256"],
             "telemetry": data[(p, l)][0]["telemetry"]}
            for p in POLICIES for l in LEVELS],
    }


def summarize(rows: list[dict]) -> dict:
    def stat(values):
        arr = np.asarray([v for v in values if v is not None and np.isfinite(v)])
        if arr.size == 0:
            return None
        return {"n": int(arr.size), "median": float(np.median(arr)),
                "p10": float(np.percentile(arr, 10)),
                "p90": float(np.percentile(arr, 90)),
                "min": float(arr.min()), "max": float(arr.max())}

    out = {"completed_orbits": len(rows)}
    out["truth_self_difference_rms_m"] = stat(
        [r["truth_self_difference_rms_m"] for r in rows])
    out["self_difference_by_policy_m"] = {
        p: stat([r["policies"][p]["self_difference_rms_m"] for r in rows])
        for p in COMPARED}
    decisions = {}
    for key in ("schedule_empirical_vs_fixed_work",
                "schedule_empirical_vs_fixed_critical",
                "schedule_up_vs_fixed_work", "schedule_up_vs_fixed_critical",
                "schedule_down_vs_fixed_work",
                "schedule_down_vs_fixed_critical"):
        cs = [r["comparisons"][key] for r in rows if key in r["comparisons"]]
        schedule_name = key.split("_vs_")[0]
        decisions[key] = {
            "pairs": len(cs),
            "resolved": sum(c["resolved"] for c in cs),
            "unresolved": sum(not c["resolved"] for c in cs),
            "resolved_schedule_wins": sum(
                c["winner_if_resolved"] == schedule_name for c in cs),
            "resolved_fixed_wins": sum(
                c["winner_if_resolved"] not in (None, schedule_name)
                for c in cs),
            "raw_schedule_wins": sum(c["rho_tight"] is not None
                                     and c["rho_tight"] > 1.0 for c in cs),
            "rho": stat([c["rho_tight"] for c in cs]),
        }
    out["decisions"] = decisions
    return out


def write_index(rows, complete, stopped, started, wall0, deadline, smoke,
                failures, workers):
    target = SMOKE_OUTPUT if smoke else OUTPUT
    payload = {
        "schema": "r11_full_convergence_index_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "corrected_baseline_sha256": base.file_hash(CORRECTED),
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "started_utc": started, "updated_utc": base.utc_now(),
        "deadline_utc": deadline.isoformat() if deadline else None,
        "complete": complete, "stopped_for_deadline": stopped,
        "timing_comparable": False,
        "execution": {"mode": "parallel_process_pool", "workers": workers,
                      "note": "accuracy campaign; cost claims use serial R10"},
        "levels": {k: {kk: vv for kk, vv in v.items() if kk != "atol"}
                   for k, v in LEVELS.items()},
        "policies": list(POLICIES),
        "planned_orbits": 64,
        "rows": rows, "failures": failures,
        "summary": summarize(rows),
        "session_wall_ns": time.perf_counter_ns() - wall0,
    }
    if complete or stopped:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(target, payload)


def run(smoke: bool, deadline, workers: int) -> int:
    rows = load_rows()
    if smoke:
        rows = sorted(rows, key=cost_estimate)[:2]
    duration = 2.0 * 3600.0 if smoke else DURATION
    script_sha = base.file_hash(Path(__file__).resolve())
    protocol_sha = base.protocol_payload()["protocol_sha256"]
    corrected_sha = base.file_hash(CORRECTED)

    for row in rows:
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / f"sobolA_{row['sobol_index']:03d}").mkdir(parents=True,
                                                             exist_ok=True)

    # Carry the destination tree explicitly in every task so spawned workers
    # write to this population's tree no matter what module globals they load.
    case_root, raw_root = str(CASE_ROOT), str(RAW_ROOT)
    tasks = []
    for row in sorted(rows, key=cost_estimate, reverse=True):
        for policy in POLICIES:
            for level in LEVELS:
                tasks.append({"row": row, "policy": policy, "level": level,
                              "duration": duration, "smoke": smoke,
                              "script_sha": script_sha,
                              "protocol_sha": protocol_sha,
                              "corrected_sha": corrected_sha,
                              "case_root": case_root, "raw_root": raw_root})

    started = base.utc_now()
    wall0 = time.perf_counter_ns()
    print(f"[r11] orbits={len(rows)} trajectories={len(tasks)} "
          f"workers={workers} deadline={deadline}", flush=True)

    done, failures, stopped = [], [], False
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(worker, t): t for t in tasks}
        for n, future in enumerate(as_completed(futures), start=1):
            try:
                record = future.result()
            except CancelledError:
                # A trajectory cancelled by the deadline guard: skip it and keep
                # draining so the index is still assembled from finished orbits.
                continue
            if record["status"] in ("numerical_failure", "worker_error"):
                failures.append(record)
                print(f"  !! {record['index']:03d} {record['policy']}/"
                      f"{record['level']}: {record.get('message')}", flush=True)
            done.append(record)
            if n % 12 == 0 or n == len(tasks):
                elapsed = time.time() - t_start
                rate = elapsed / n
                print(f"  [{n:4d}/{len(tasks)}] elapsed={elapsed/3600:5.2f}h "
                      f"eta={(len(tasks)-n)*rate/3600:5.2f}h "
                      f"last={record['policy']}/{record['level']} "
                      f"{record['status']} {record['wall_s']:.0f}s", flush=True)
            if remaining(deadline) < 900.0 and not stopped:
                stopped = True
                print("  deadline guard: cancelling queued trajectories",
                      flush=True)
                for pending in futures:
                    pending.cancel()

    summaries = []
    for row in sorted(rows, key=lambda r: r["sobol_index"]):
        try:
            summaries.append(orbit_summary(row, smoke))
        except Exception as exc:
            failures.append({"index": row["sobol_index"], "policy": "summary",
                             "level": "-", "status": "summary_incomplete",
                             "message": f"{type(exc).__name__}: {exc}"})
    complete = len(summaries) == len(rows) and not failures
    write_index(summaries, complete, stopped, started, wall0, deadline, smoke,
                failures, workers)
    print(f"[r11] finished orbits={len(summaries)}/{len(rows)} "
          f"failures={len(failures)} complete={complete}", flush=True)
    return 0 if complete else 3


def status() -> int:
    if not OUTPUT.exists():
        print("no r11 output yet")
        return 0
    data = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({"complete": data["complete"],
                      "stopped_for_deadline": data["stopped_for_deadline"],
                      "orbits": len(data["rows"]),
                      "failures": len(data.get("failures", [])),
                      "summary": data["summary"]}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("run", "smoke", "status"))
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--deadline")
    args = parser.parse_args()
    if args.command == "status":
        return status()
    return run(args.command == "smoke", parse_deadline(args.deadline),
               args.workers)


if __name__ == "__main__":
    raise SystemExit(main())
