"""Selective vector-tolerance convergence for the R10 Sobol A population."""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev7_doe_screening import alt_sched, degree_power, emp_table, kernel_args


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
BASELINE = METRICS / "r10_sobolA_baseline.json"
AUDIT = METRICS / "r10_sobolA_truth_audit.json"
SELECTION = METRICS / "r10_sobolA_convergence_selection.json"
CORRECTED = METRICS / "r10_sobolA_baseline_truth_corrected.json"
OUTPUT = METRICS / "r10_sobolA_convergence.json"
SMOKE_OUTPUT = METRICS / "r10_sobolA_convergence_smoke.json"
ACTIVE = METRICS / "r10_sobolA_convergence_active.json"
CASE_ROOT = METRICS / "r10_cases" / "convergence"
RAW_ROOT = METRICS / "r10_raw" / "convergence"

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
POLICIES = ("truth", "schedule_empirical", "fixed_critical", "fixed_work")


def load_inputs() -> tuple[dict, dict]:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    if not baseline.get("complete") or len(baseline.get("rows", [])) != 64:
        raise RuntimeError("baseline is incomplete")
    if not audit.get("complete") or len(audit.get("rows", [])) != 6:
        raise RuntimeError("truth audit is incomplete")
    return baseline, audit


def corrected_rows(baseline: dict, audit: dict) -> list[dict]:
    audit_rows = {r["sobol_index"]: r for r in audit["rows"]}
    result = []
    for row in baseline["rows"]:
        index = row["sobol_index"]
        adopted = audit_rows.get(index, {}).get("adopted_truth_degree",
                                                row["truth_degree"])
        if adopted == 900:
            policy_errors = audit_rows[index]["policy_errors_against_N900"]
            empirical = policy_errors["schedule_empirical"]["pos_rms_m"]
            rho_work = policy_errors["fixed_work"]["pos_rms_m"] / empirical
            rho_crit = policy_errors["fixed_critical"]["pos_rms_m"] / empirical
            source = "truth_audit_N900"
        else:
            policy_errors = {
                p: row["policies"][p]["error_against_truth"]
                for p in ("schedule_empirical", "fixed_work", "fixed_critical",
                          "schedule_up", "schedule_down")
            }
            rho_work = row["primary_ratios"]["rho_work"]
            rho_crit = row["primary_ratios"]["rho_crit"]
            source = f"baseline_N{row['truth_degree']}"
        result.append({
            "sobol_index": index,
            "name": row["name"],
            "design_point": row["design_point"],
            "original_truth_degree": row["truth_degree"],
            "adopted_truth_degree": adopted,
            "truth_source": source,
            "n_work": row["n_work"],
            "n_critical": row["n_critical"],
            "policy_errors": policy_errors,
            "rho_work": rho_work,
            "rho_crit": rho_crit,
        })
    return result


def selection_parts(rows: list[dict]) -> dict:
    raw = {r["sobol_index"] for r in rows
           if r["rho_work"] > 1.0 or r["rho_crit"] > 1.0}
    close = {r["sobol_index"] for r in rows
             if 0.75 < r["rho_work"] < 1.33
             or 0.75 < r["rho_crit"] < 1.33}
    errors = np.asarray([
        r["policy_errors"]["schedule_empirical"]["pos_rms_m"] for r in rows
    ])
    quantiles = set()
    for q in (0.0, 25.0, 50.0, 75.0, 100.0):
        target = float(np.percentile(errors, q))
        quantiles.add(rows[int(np.argmin(np.abs(errors - target)))]["sobol_index"])
    low_hp = {r["sobol_index"] for r in sorted(
        rows, key=lambda x: x["design_point"]["hp_km"]
    )[:4]}
    high_ecc = {r["sobol_index"] for r in sorted(
        rows, key=lambda x: x["design_point"]["eccentricity"], reverse=True
    )[:4]}
    selected = raw | close | quantiles | low_hp | high_ecc
    regime = lambda i: ("prograde" if i < 60.0 else
                        "high_inclination" if i <= 120.0 else "retrograde")
    additions = set()
    for name in ("prograde", "high_inclination", "retrograde"):
        present = [r for r in rows if r["sobol_index"] in selected
                   and regime(r["design_point"]["incl_deg"]) == name]
        if len(present) < 2:
            candidates = sorted(
                (r for r in rows if r["sobol_index"] not in selected
                 and regime(r["design_point"]["incl_deg"]) == name),
                key=lambda r: r["policy_errors"]["schedule_empirical"]["pos_rms_m"],
            )
            additions.update(r["sobol_index"] for r in candidates[:2-len(present)])
    selected |= additions
    return {
        "raw_win": sorted(raw), "close": sorted(close),
        "error_quantiles": sorted(quantiles), "lowest_hp": sorted(low_hp),
        "highest_eccentricity": sorted(high_ecc),
        "inclination_additions": sorted(additions),
        "selected": sorted(selected),
    }


def freeze_selection() -> dict:
    baseline, audit = load_inputs()
    corrected = corrected_rows(baseline, audit)
    # Pre-audit selection is retained even when a corrected quantile changes.
    pre = corrected_rows(baseline, {"rows": []})
    pre_parts = selection_parts(pre)
    post_parts = selection_parts(corrected)
    final = sorted(set(pre_parts["selected"]) | set(post_parts["selected"]))
    protocol = base.protocol_payload()
    corrected_payload = {
        "schema": "r10_truth_corrected_baseline_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "baseline_sha256": base.file_hash(BASELINE),
        "truth_audit_sha256": base.file_hash(AUDIT),
        "rows": corrected,
    }
    corrected_payload["artifact_sha256"] = base.object_hash(corrected_payload)
    base.atomic_json(CORRECTED, corrected_payload)
    reasons = {}
    labels = {
        "raw_win": "A_raw_win", "close": "B_close",
        "error_quantiles": "C_error_quantile", "lowest_hp": "D_lowest_hp",
        "highest_eccentricity": "E_highest_eccentricity",
        "inclination_additions": "F_inclination_representation",
    }
    for source_name, parts in (("pre_audit", pre_parts),
                               ("post_audit", post_parts)):
        for key, label in labels.items():
            for index in parts[key]:
                reasons.setdefault(str(index), []).append(f"{source_name}:{label}")
    payload = {
        "schema": "r10_sobolA_convergence_selection_v1",
        "created_utc": base.utc_now(),
        "protocol_sha256": protocol["protocol_sha256"],
        "baseline_sha256": base.file_hash(BASELINE),
        "truth_audit_sha256": base.file_hash(AUDIT),
        "corrected_baseline_sha256": corrected_payload["artifact_sha256"],
        "pre_audit": pre_parts,
        "post_audit": post_parts,
        "retention_rule": "union; no pre-audit selected case may be removed",
        "selected_indices": final,
        "selected_count": len(final),
        "selection_reasons": reasons,
    }
    payload["selection_sha256"] = base.object_hash(payload)
    base.atomic_json(SELECTION, payload)
    return payload


def parse_deadline(value: str | None):
    if not value:
        return None
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        raise ValueError("deadline needs UTC offset")
    return result.astimezone(timezone.utc)


def remaining(deadline) -> float:
    return math.inf if deadline is None else (
        deadline - datetime.now(timezone.utc)
    ).total_seconds()


def paths(index: int, policy: str, level: str, smoke: bool):
    suffix = "_smoke" if smoke else ""
    case = CASE_ROOT / f"sobolA_{index:03d}"
    raw = RAW_ROOT / f"sobolA_{index:03d}"
    return (case / f"{policy}_{level}{suffix}.json",
            raw / f"{policy}_{level}{suffix}.npz")


def run_trajectory(row: dict, policy: str, level: str, model, args,
                   schedule, duration: float, smoke: bool):
    tol = LEVELS[level]
    if policy == "truth":
        degree = int(row["adopted_truth_degree"])
        degree_of = lambda t, h, n=degree: n
        spec = {"kind": "fixed_truth", "degree": degree}
    elif policy == "schedule_empirical":
        degree_of = schedule
        spec = {"kind": "frozen_empirical_lookup",
                "lookup_source_degree": row["original_truth_degree"]}
    elif policy == "fixed_critical":
        degree = int(row["n_critical"])
        degree_of = lambda t, h, n=degree: n
        spec = {"kind": "fixed_critical", "degree": degree}
    elif policy == "fixed_work":
        degree = int(row["n_work"])
        degree_of = lambda t, h, n=degree: n
        spec = {"kind": "fixed_work", "degree": degree,
                "source": "scalar-baseline empirical RHS history"}
    else:
        raise ValueError(policy)
    protocol = base.protocol_payload()
    selection = json.loads(SELECTION.read_text(encoding="utf-8"))
    config = {
        "schema": "r10_sobol_convergence_config_v1",
        "protocol_sha256": protocol["protocol_sha256"],
        "selection_sha256": selection["selection_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "sobol_index": row["sobol_index"],
        "initial_state_si": row["design_point"]["initial_state_si"],
        "adopted_truth_degree": row["adopted_truth_degree"],
        "policy": policy, "policy_spec": spec, "level": level,
        "duration_s": duration, "output_step_s": OUTPUT_STEP,
        "integrator": "InstrumentedDOP853",
        "rtol": tol["rtol"],
        "atol_position_m": tol["atol_position_m"],
        "atol_velocity_m_s": tol["atol_velocity_m_s"],
        "max_step_s": MAX_STEP, "timing_comparable": True,
        "source": protocol["provenance"],
    }
    config_sha = base.object_hash(config)
    sidecar, raw = paths(row["sobol_index"], policy, level, smoke)
    if sidecar.exists() and raw.exists() and base.valid_cached(
            sidecar, raw, config_sha, duration):
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        t, y = base.load_raw(raw)
        print(f"    {policy}/{level} cached", flush=True)
        return meta, t, y
    if sidecar.exists() or raw.exists():
        base.preserve_invalid(sidecar); base.preserve_invalid(raw)
    base.atomic_json(ACTIVE, {
        "updated_utc": base.utc_now(), "sobol_index": row["sobol_index"],
        "policy": policy, "level": level, "config_sha256": config_sha,
    })
    grid = np.arange(0.0, duration + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    t, y, status, event, failure, telemetry = base.propagate_event_instrumented(
        model, np.asarray(row["design_point"]["initial_state_si"]), duration,
        grid, degree_of, args, tol["rtol"], tol["atol"], max_step=MAX_STEP,
    )
    if status == "numerical_failure":
        fail = sidecar.with_name(sidecar.stem + f".failure.{int(time.time())}.json")
        base.atomic_json(fail, {"config": config, "config_sha256": config_sha,
                                "status": status, "failure_message": failure,
                                "telemetry": telemetry})
        raise RuntimeError(f"{row['sobol_index']} {policy}/{level}: {failure}")
    arrays = {"t_s": t, "state_si": y}
    if event:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"])
    base.atomic_npz(raw, **arrays)
    meta = {
        "schema": "r10_sobol_convergence_result_v1",
        "created_utc": base.utc_now(), "config": config,
        "config_sha256": config_sha, "status": status, "event": event,
        "failure_message": None, "telemetry": telemetry,
        "raw_path": str(raw.relative_to(ROOT)), "raw_sha256": base.file_hash(raw),
        "n_output_epochs": int(len(t)), "last_output_epoch_s": float(t[-1]),
    }
    base.atomic_json(sidecar, meta)
    print(f"    {policy}/{level} {status} wall={telemetry['total_wall_ns']/1e9:.1f}s", flush=True)
    return meta, t, y


def orbit_summary(row: dict, data: dict) -> dict:
    truth_self = base.common_error(
        data[("truth", "tight")][1], data[("truth", "tight")][2],
        data[("truth", "tighter")][1], data[("truth", "tighter")][2],
    )["pos_rms_m"]
    policies = {}
    for policy in ("schedule_empirical", "fixed_critical", "fixed_work"):
        errors = {}
        for level in LEVELS:
            errors[level] = base.common_error(
                data[(policy, level)][1], data[(policy, level)][2],
                data[("truth", level)][1], data[("truth", level)][2],
            )
        self_diff = base.common_error(
            data[(policy, "tight")][1], data[(policy, "tight")][2],
            data[(policy, "tighter")][1], data[(policy, "tighter")][2],
        )["pos_rms_m"]
        policies[policy] = {
            "errors_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self,
        }
    comparisons = {}
    schedule = policies["schedule_empirical"]
    for comparator in ("fixed_work", "fixed_critical"):
        fixed = policies[comparator]
        es = schedule["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
        ef = fixed["errors_against_same_tolerance_truth"]["tight"]["pos_rms_m"]
        difference = abs(es - ef)
        threshold = (schedule["truth_inclusive_envelope_m"] +
                     fixed["truth_inclusive_envelope_m"])
        comparisons[comparator] = {
            "rho_tight": ef / es,
            "absolute_error_difference_m": difference,
            "resolution_threshold_m": threshold,
            "resolved": bool(difference > threshold),
            "winner_if_resolved": (
                "schedule_empirical" if es < ef else comparator
            ) if difference > threshold else None,
        }
    return {
        "sobol_index": row["sobol_index"], "name": row["name"],
        "adopted_truth_degree": row["adopted_truth_degree"],
        "truth_self_difference_rms_m": truth_self,
        "policies": policies, "comparisons": comparisons,
        "trajectory_records": [
            {"policy": p, "level": l, "status": data[(p, l)][0]["status"],
             "config_sha256": data[(p, l)][0]["config_sha256"],
             "raw_sha256": data[(p, l)][0]["raw_sha256"],
             "telemetry": data[(p, l)][0]["telemetry"]}
            for p in POLICIES for l in LEVELS
        ],
    }


def write_index(rows, selected, complete, stopped, started, wall0, deadline, smoke):
    target = SMOKE_OUTPUT if smoke else OUTPUT
    payload = {
        "schema": "r10_sobolA_convergence_index_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "selection_sha256": selected["selection_sha256"],
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "started_utc": started, "updated_utc": base.utc_now(),
        "deadline_utc": deadline.isoformat() if deadline else None,
        "complete": complete, "stopped_for_deadline": stopped,
        "timing_comparable": True, "selected_count": len(selected["selected_indices"]),
        "rows": rows,
        "summary": {
            "completed_orbits": len(rows),
            "resolved_schedule_wins_vs_work": sum(
                r["comparisons"]["fixed_work"]["winner_if_resolved"] ==
                "schedule_empirical" for r in rows),
            "resolved_schedule_wins_vs_critical": sum(
                r["comparisons"]["fixed_critical"]["winner_if_resolved"] ==
                "schedule_empirical" for r in rows),
        },
        "session_wall_ns": time.perf_counter_ns() - wall0,
    }
    if complete or stopped:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(target, payload)


def estimate_orbit_seconds(row: dict, baseline_by_index: dict,
                           audit_by_index: dict) -> float:
    index = row["sobol_index"]
    old = baseline_by_index[index]
    if row["adopted_truth_degree"] == 900:
        truth_base = audit_by_index[index]["N900_source"]["telemetry"]["total_wall_ns"] / 1e9
    else:
        truth_base = old["policies"]["truth"]["telemetry"]["total_wall_ns"] / 1e9
    policy_base = sum(old["policies"][p]["telemetry"]["total_wall_ns"] / 1e9
                      for p in ("schedule_empirical", "fixed_critical", "fixed_work"))
    return max(600.0, 2.5 * truth_base + 3.2 * policy_base)


def run(smoke: bool, deadline) -> int:
    if base.other_python_processes():
        raise RuntimeError("another Python process is active")
    selected = freeze_selection()
    corrected = json.loads(CORRECTED.read_text(encoding="utf-8"))["rows"]
    row_by_index = {r["sobol_index"]: r for r in corrected}
    baseline, audit = load_inputs()
    baseline_by_index = {r["sobol_index"]: r for r in baseline["rows"]}
    audit_by_index = {r["sobol_index"]: r for r in audit["rows"]}
    indices = selected["selected_indices"]
    # Heavy-first: adopted N900 first, then descending measured baseline cost.
    indices = sorted(indices, key=lambda i: (
        row_by_index[i]["adopted_truth_degree"],
        baseline_by_index[i]["orbit_total_trajectory_wall_s"]), reverse=True)
    if smoke:
        indices = indices[:1]
    needed_truth = sorted({row_by_index[i]["adopted_truth_degree"] for i in indices})
    models = {}
    for degree in needed_truth:
        model = base.load_model(degree); args = base.kernel_args(model); base.warmup(model, args)
        models[degree] = (model, args)
    # Freeze lookup tables from the original N300/N600 policy contract.
    lookup = {}
    for degree in sorted({row_by_index[i]["original_truth_degree"] for i in indices}):
        model = base.load_model(degree)
        lookup[degree] = emp_table(model, degree_power(model))
    target = SMOKE_OUTPUT if smoke else OUTPUT
    previous = []
    if target.exists():
        previous = json.loads(target.read_text(encoding="utf-8")).get("rows", [])
    completed = {r["sobol_index"]: r for r in previous}
    duration = 2.0 * 3600.0 if smoke else DURATION
    started = base.utc_now(); wall0 = time.perf_counter_ns(); stopped = False
    print(f"[convergence] orbits={len(indices)} deadline={deadline}", flush=True)
    for position, index in enumerate(indices, start=1):
        row = row_by_index[index]
        estimate = 300.0 if smoke else estimate_orbit_seconds(
            row, baseline_by_index, audit_by_index
        )
        if remaining(deadline) < estimate + 600.0:
            print(f"  stop before orbit {index}: deadline guard estimate={estimate:.0f}s", flush=True)
            stopped = True
            break
        if base.other_python_processes():
            raise RuntimeError("another Python process appeared during convergence")
        print(f"[{position}/{len(indices)}] sobolA_{index:03d} Ntruth={row['adopted_truth_degree']} estimate={estimate:.0f}s", flush=True)
        model, args = models[row["adopted_truth_degree"]]
        schedule = alt_sched(lookup[row["original_truth_degree"]])
        data = {}
        order = (("truth", "tighter"), ("truth", "tight"),
                 ("schedule_empirical", "tighter"), ("schedule_empirical", "tight"),
                 ("fixed_critical", "tighter"), ("fixed_critical", "tight"),
                 ("fixed_work", "tighter"), ("fixed_work", "tight"))
        for policy, level in order:
            data[(policy, level)] = run_trajectory(
                row, policy, level, model, args, schedule, duration, smoke
            )
        completed[index] = orbit_summary(row, data)
        ordered_rows = [completed[i] for i in indices if i in completed]
        write_index(ordered_rows, selected, False, False, started, wall0,
                    deadline, smoke)
    ordered_rows = [completed[i] for i in indices if i in completed]
    complete = len(ordered_rows) == len(indices)
    write_index(ordered_rows, selected, complete, stopped, started, wall0,
                deadline, smoke)
    if ACTIVE.exists(): ACTIVE.unlink()
    print(f"[convergence] finished complete={complete} orbits={len(ordered_rows)}/{len(indices)}", flush=True)
    return 0 if complete else 3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("freeze", "smoke", "run", "status"))
    parser.add_argument("--deadline")
    args = parser.parse_args()
    if args.command == "freeze":
        print(json.dumps(freeze_selection(), indent=2)); return 0
    if args.command == "status":
        if OUTPUT.exists():
            data = json.loads(OUTPUT.read_text(encoding="utf-8"))
            print(json.dumps({"complete": data["complete"],
                              "completed_orbits": len(data["rows"]),
                              "selected_count": data["selected_count"],
                              "stopped_for_deadline": data["stopped_for_deadline"],
                              "summary": data["summary"]}, indent=2))
        else: print("no convergence output")
        return 0
    return run(args.command == "smoke", parse_deadline(args.deadline))


if __name__ == "__main__":
    raise SystemExit(main())
