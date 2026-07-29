"""Round-17: sixty-day long arcs on a widened geometry set.

The archived month-long stage (R7 stage 3) extended the seven-day comparison to
28 days on four geometries, of which one reached the surface on day 8.2, leaving
three. Two limitations followed and are stated in the manuscript: the sample was
too small to carry a population statement, and it was run at a single scalar
tolerance, so differences smaller than the integration floor could not be ranked.

This campaign addresses both at once.

  * The arc is extended to 60 days.
  * The geometry set is widened from 3 usable cases to 11, by a rule fixed
    before any of them was propagated (see GEOMETRY_RULE).
  * Every configuration is run at two vector-tolerance levels, so each policy
    carries its own numerical envelope and each comparison is resolved under the
    same truth-inclusive rule used everywhere else in the paper.

Scope: only geometries whose perilune is at or above 50 km are included, so the
truth degree is 300 throughout. The LRO-like 30x216 km case of the 28-day stage
needs an N=600 truth, which costs roughly six times as much per arc; it is
deliberately not extended here and its 28-day result stands unchanged. That
exclusion is a cost decision, and because the excluded case is the one geometry
that came closest to favoring a schedule at 28 days, it is recorded as a
restriction on this campaign rather than as a neutral sampling choice.

Schedules, the empirical degree table, quantization and the degree cap are taken
unchanged from the 28-day stage, so the two campaigns differ in arc length,
geometry count and tolerance protocol -- not in policy definition.

Usage:
    python rev17_longarc60.py run --workers 11 --deadline-min 210
    python rev17_longarc60.py smoke --workers 4
    python rev17_longarc60.py summarize
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev3_common import DAY
from rev7_doe_screening import (CANONICAL, CAP, alt_sched, emp_nmin_exact,
                                emp_table, initial_state, kaula_table,
                                perilune_geography)

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r17_longarc60.json"
SMOKE_OUTPUT = METRICS / "r17_longarc60_smoke.json"
CASE_ROOT = METRICS / "r17_cases" / "longarc60"
RAW_ROOT = METRICS / "r17_raw" / "longarc60"

TRUTH_DEGREE = 300
MIN_PERILUNE_KM = 50.0
DURATION = 60.0 * DAY
OUTPUT_STEP = 300.0
MAX_STEP = 120.0
CHECKPOINT_DAYS = (7, 14, 28, 42, 60)

# Same two-level ladder used for the seven-day population reruns.
LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol_position_m": 1.0e-5,
              "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol_position_m": 1.0e-6,
                "atol_velocity_m_s": 1.0e-9},
}

# fixed_work depends on the realized mean squared degree of sched_emp, so it can
# only be built after wave 1 has run.
WAVE1_POLICIES = ("truth", "fixed_crit", "sched_emp", "sched_down", "sched_up")
WAVE2_POLICIES = ("fixed_work",)
COMPARED = ("fixed_crit", "fixed_work", "sched_emp", "sched_down", "sched_up")
SCHEDULES = ("sched_emp", "sched_down", "sched_up")

GEOMETRY_RULE = (
    "The two canonical geometries of the 28-day stage whose perilune is at or "
    "above 50 km (c2_50x300_polar, c3_50x300_i60), in that order, followed by "
    "the orbits of the frozen design-A scrambled-Sobol population in archived "
    "order whose perilune is at or above 50 km, taking the first eight. The rule "
    "was fixed before any 60-day arc was propagated and depends on no measured "
    "outcome; the 50 km floor is the truth-degree boundary, not a selection on "
    "results."
)
N_SOBOL = 8

_MODEL: dict = {}


def atol_vector(level: str) -> np.ndarray:
    tol = LEVELS[level]
    return np.array([tol["atol_position_m"]] * 3
                    + [tol["atol_velocity_m_s"]] * 3, dtype=float)


def model_args():
    """Load the truth-degree model and the three schedule tables once per
    worker process."""
    if "m" not in _MODEL:
        model = base.load_model(TRUTH_DEGREE)
        args = base.kernel_args(model)
        base.warmup(model, args)
        power = base.degree_power(model)
        _MODEL["m"] = (model, args)
        _MODEL["power"] = power
        _MODEL["sched_emp"] = alt_sched(emp_table(model, power))
        _MODEL["sched_down"] = alt_sched(kaula_table(model, "down"))
        _MODEL["sched_up"] = alt_sched(kaula_table(model, "up"))
    return _MODEL["m"]


def geometries() -> list[dict]:
    """Apply GEOMETRY_RULE. Canonical cases carry their stage-3 initial state;
    Sobol cases carry the state archived with the frozen design, so neither is
    reconstructed differently here than where it was first used."""
    model, _ = model_args()
    out = []
    for name, hp, ha, inc, argp, raan in CANONICAL:
        if name in ("c2_50x300_polar", "c3_50x300_i60"):
            orb = {"name": name, "source": "canonical", "hp_km": hp,
                   "ha_km": ha, "incl_deg": inc, "argp_deg": argp,
                   "raan_deg": raan}
            orb["initial_state_si"] = [float(x) for x in
                                       initial_state(model, orb)]
            out.append(orb)
    design = json.loads(
        (METRICS / "r10_sobolA_design.json").read_text(encoding="utf-8"))
    rows = design.get("orbits") or design.get("rows")
    taken = 0
    for row in rows:
        if taken >= N_SOBOL:
            break
        if float(row["hp_km"]) < MIN_PERILUNE_KM:
            continue
        out.append({"name": row["name"], "source": "sobolA",
                    "hp_km": float(row["hp_km"]), "ha_km": float(row["ha_km"]),
                    "incl_deg": float(row["incl_deg"]),
                    "argp_deg": float(row["argp_deg"]),
                    "raan_deg": float(row["raan_deg"]),
                    "initial_state_si": [float(x)
                                         for x in row["initial_state_si"]]})
        taken += 1
    return out


def paths(index: int, policy: str, level: str, smoke: bool) -> tuple[Path, Path]:
    tag = "smoke" if smoke else "run"
    case = CASE_ROOT / f"{tag}_orbit_{index:02d}"
    raw = RAW_ROOT / f"{tag}_orbit_{index:02d}"
    case.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    return case / f"{policy}_{level}.json", raw / f"{policy}_{level}.npz"


def degree_function(policy: str, spec: dict):
    """Return (callable, description). n_crit and n_work arrive through spec so
    a worker never has to reproduce a decision made elsewhere."""
    if policy == "truth":
        return (lambda t, h: TRUTH_DEGREE), {"kind": "fixed",
                                             "degree": TRUTH_DEGREE}
    if policy == "fixed_crit":
        n = int(spec["n_crit"])
        return (lambda t, h: n), {"kind": "fixed", "degree": n,
                                  "basis": "critical-altitude empirical"}
    if policy == "fixed_work":
        n = int(spec["n_work"])
        return (lambda t, h: n), {"kind": "fixed", "degree": n,
                                  "basis": "matched to sched_emp mean N^2"}
    model_args()
    return _MODEL[policy], {"kind": "altitude_schedule", "table": policy}


def worker(task: dict) -> dict:
    index, policy, level = task["index"], task["policy"], task["level"]
    smoke, duration = task["smoke"], task["duration"]
    sidecar, raw = paths(index, policy, level, smoke)
    try:
        model, args = model_args()
        degree_of, spec = degree_function(policy, task["spec"])
        y0 = np.asarray(task["orbit"]["initial_state_si"], dtype=float)
        tol = LEVELS[level]
        config = {
            "schema": "r17_longarc60_config_v1",
            "script_sha256": task["script_sha"],
            "orbit_index": index, "orbit": task["orbit"],
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
            "extends": "r7_doe_stage3_longarc (28 days, 4 geometries, "
                       "scalar tolerance)",
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
            "schema": "r17_longarc60_result_v1",
            "created_utc": base.utc_now(), "config": config,
            "config_sha256": config_sha, "status": status, "event": event,
            "failure_message": None, "telemetry": telemetry,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
        return {"index": index, "policy": policy, "level": level,
                "status": status, "wall_s": telemetry["total_wall_ns"] / 1e9,
                "mean_deg_sq": telemetry.get("mean_degree_sq")}
    except Exception as exc:  # noqa: BLE001
        return {"index": index, "policy": policy, "level": level,
                "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(), "wall_s": 0.0}


def _checkpoints(t, policy_y, truth_y, duration):
    """Cumulative RMS and instantaneous in-track error at each checkpoint that
    the arc actually reached."""
    n = min(policy_y.shape[1], truth_y.shape[1])
    out = {}
    for d in CHECKPOINT_DAYS:
        t_cp = d * DAY
        if t_cp > duration + 1.0 or t_cp > t[n - 1] + 1.0:
            continue
        idx = min(int(np.searchsorted(t[:n], t_cp - 1e-6)), n - 1)
        st = base.err_stats(policy_y[:, :idx + 1], truth_y[:, :idx + 1])
        out[f"d{d}"] = {
            "pos_rms_m": st["pos_rms_m"],
            "in_track_at_checkpoint_m": st["ric_final_m"]["in_track"],
        }
    return out


def orbit_summary(index: int, orbit: dict, spec: dict, smoke: bool,
                  duration: float) -> dict:
    # truth, the critical-altitude comparator and the three schedules are
    # required; the work-matched comparator is optional, because it can only be
    # built after the schedules have run and a campaign stopped at its deadline
    # may not have reached it. An orbit missing only fixed_work still carries
    # every schedule-versus-critical-degree comparison.
    required = ("truth", "fixed_crit") + SCHEDULES
    data = {}
    impact = None
    available = []
    for policy in COMPARED + ("truth",):
        present = True
        for level in LEVELS:
            sidecar, raw = paths(index, policy, level, smoke)
            if not raw.exists():
                if policy in required:
                    return {"orbit_index": index, "orbit": orbit,
                            "status": "incomplete",
                            "missing": f"{policy}_{level}"}
                present = False
                break
            data[(policy, level)] = base.load_raw(raw)
            meta = json.loads(sidecar.read_text(encoding="utf-8"))
            if policy == "truth" and meta.get("event"):
                impact = meta["event"]
        if present and policy != "truth":
            available.append(policy)

    row = {"orbit_index": index, "orbit": orbit, "spec": spec,
           "status": "complete", "policies_available": available,
           "work_matched_comparator_present": "fixed_work" in available}
    row.update(perilune_geography(np.asarray(orbit["initial_state_si"])))

    t_truth, y_truth_tight = data[("truth", "tight")]
    _, y_truth_tighter = data[("truth", "tighter")]
    truth_self = base.common_error(t_truth, y_truth_tight,
                                   t_truth, y_truth_tighter)["pos_rms_m"]
    row["truth_self_difference_rms_m"] = truth_self
    row["arc_end_day"] = float(t_truth[-1] / DAY)
    row["reached_full_arc"] = bool(t_truth[-1] >= duration - 1.0)
    if impact:
        row["truth_impact_day"] = float(impact["epoch_s"] / DAY)

    policies = {}
    for policy in available:
        errors = {}
        for level in LEVELS:
            pt, py = data[(policy, level)]
            tt, ty = data[("truth", level)]
            errors[level] = base.common_error(pt, py, tt, ty)
        pt1, py1 = data[(policy, "tight")]
        pt2, py2 = data[(policy, "tighter")]
        self_diff = base.common_error(pt1, py1, pt2, py2)["pos_rms_m"]
        entry = {
            "errors_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self,
            "checkpoints": _checkpoints(pt2, py2, data[("truth", "tighter")][1],
                                        duration),
        }
        policies[policy] = entry
    row["policies"] = policies

    # Every schedule is compared with both fixed comparators under the same
    # truth-inclusive rule used for the seven-day populations.
    comparisons = {}
    for sched in SCHEDULES:
        e_s = policies[sched]["errors_against_same_tolerance_truth"]["tighter"]["pos_rms_m"]
        env_s = policies[sched]["truth_inclusive_envelope_m"]
        for ref in ("fixed_crit", "fixed_work"):
            if ref not in policies:
                continue
            e_r = policies[ref]["errors_against_same_tolerance_truth"]["tighter"]["pos_rms_m"]
            env_r = policies[ref]["truth_inclusive_envelope_m"]
            gap = abs(e_s - e_r)
            thr = env_s + env_r
            comparisons[f"{sched}_vs_{ref}"] = {
                "schedule_error_m": e_s, "fixed_error_m": e_r,
                "rho": (e_r / e_s) if e_s > 0 else None,
                "absolute_difference_m": gap,
                "resolution_threshold_m": thr,
                "resolved": bool(gap > thr),
                "winner_if_resolved": ((sched if e_s < e_r else ref)
                                       if gap > thr else None)}
    row["comparisons"] = comparisons

    best = min(SCHEDULES, key=lambda s: policies[s][
        "errors_against_same_tolerance_truth"]["tighter"]["pos_rms_m"])
    row["best_schedule_name"] = best
    cps = policies[best]["checkpoints"]
    days = [d for d in CHECKPOINT_DAYS if f"d{d}" in cps]
    it = [abs(cps[f"d{d}"]["in_track_at_checkpoint_m"]) for d in days]
    row["best_schedule_in_track_growth"] = {
        "checkpoint_days": days, "abs_in_track_m": it,
        "d60_over_d7": (it[-1] / it[0]) if len(it) >= 2 and it[0] > 0 else None}
    return row


def summarize(rows: list[dict]) -> dict:
    usable = [r for r in rows if r.get("status") == "complete"
              and r.get("reached_full_arc")]
    impacted = [r for r in rows if r.get("status") == "complete"
                and not r.get("reached_full_arc")]

    def stat(values):
        arr = np.asarray([v for v in values
                          if v is not None and np.isfinite(v)])
        if arr.size == 0:
            return None
        return {"n": int(arr.size), "median": float(np.median(arr)),
                "p10": float(np.percentile(arr, 10)),
                "p90": float(np.percentile(arr, 90)),
                "min": float(arr.min()), "max": float(arr.max())}

    out = {
        "orbits_attempted": len(rows),
        "orbits_reaching_60_days": len(usable),
        "orbits_terminated_early": len(impacted),
        "early_termination_days": [
            {"name": r["orbit"]["name"], "day": r.get("arc_end_day")}
            for r in impacted],
        "truth_self_difference_rms_m": stat(
            [r["truth_self_difference_rms_m"] for r in usable]),
    }
    out["orbits_with_work_matched_comparator"] = sum(
        r.get("work_matched_comparator_present", False) for r in usable)
    out["error_by_policy_m"] = {
        p: stat([r["policies"][p]["errors_against_same_tolerance_truth"]
                 ["tighter"]["pos_rms_m"] for r in usable if p in r["policies"]])
        for p in COMPARED}
    out["envelope_by_policy_m"] = {
        p: stat([r["policies"][p]["truth_inclusive_envelope_m"]
                 for r in usable if p in r["policies"]]) for p in COMPARED}
    out["comparisons"] = {}
    keys = [f"{s}_vs_{ref}" for s in SCHEDULES
            for ref in ("fixed_crit", "fixed_work")]
    for key in keys:
        cs = [r["comparisons"][key] for r in usable if key in r["comparisons"]]
        if not cs:
            continue
        sched = key.split("_vs_")[0]
        out["comparisons"][key] = {
            "n_orbits": len(cs),
            "rho": stat([c["rho"] for c in cs]),
            "resolved": sum(c["resolved"] for c in cs),
            "unresolved": sum(not c["resolved"] for c in cs),
            "resolved_schedule_wins": sum(
                c["winner_if_resolved"] == sched for c in cs),
            "resolved_fixed_wins": sum(
                c["winner_if_resolved"] not in (sched, None) for c in cs)}
    growth = [r["best_schedule_in_track_growth"].get("d60_over_d7")
              for r in usable]
    out["best_schedule_in_track_d60_over_d7"] = stat(growth)
    return out


def run(smoke: bool, workers: int, deadline_min: float) -> int:
    duration = 2.0 * DAY if smoke else DURATION
    out_path = SMOKE_OUTPUT if smoke else OUTPUT
    script_sha = base.file_hash(Path(__file__).resolve())
    orbits = geometries()
    if smoke:
        orbits = orbits[:3]
    model, _ = model_args()
    power = _MODEL["power"]

    specs = []
    for orb in orbits:
        n_crit = int(min(CAP, emp_nmin_exact(power, model.r_ref,
                                             orb["hp_km"] * 1e3)))
        specs.append({"n_crit": n_crit})

    t_start = time.time()
    deadline = t_start + deadline_min * 60.0

    def submit(tasks, label):
        done, failures = [], []
        print(f"[r17] {label}: {len(tasks)} trajectories, workers={workers}",
              flush=True)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(worker, t): t for t in tasks}
            for n, fut in enumerate(as_completed(futures), start=1):
                rec = fut.result()
                if rec["status"] in ("numerical_failure", "worker_error"):
                    failures.append(rec)
                    print(f"  [{n}/{len(tasks)}] FAIL orbit={rec['index']} "
                          f"{rec['policy']}/{rec['level']}: "
                          f"{rec.get('message')}", flush=True)
                else:
                    done.append(rec)
                    el = (time.time() - t_start) / 60.0
                    print(f"  [{n}/{len(tasks)}] orbit={rec['index']:02d} "
                          f"{rec['policy']:11s} {rec['level']:8s} "
                          f"{rec['status']:9s} wall={rec['wall_s']:7.1f}s "
                          f"elapsed={el:5.1f}min", flush=True)
                if time.time() > deadline:
                    print("[r17] deadline reached; cancelling pending work",
                          flush=True)
                    for f in futures:
                        f.cancel()
                    break
        return done, failures

    wave1 = [{"index": i, "orbit": orbits[i], "spec": specs[i],
              "policy": p, "level": l, "duration": duration, "smoke": smoke,
              "script_sha": script_sha}
             for i in range(len(orbits)) for p in WAVE1_POLICIES
             for l in LEVELS]
    done1, fail1 = submit(wave1, "wave 1 (truth, critical, schedules)")

    # fixed_work is matched to the mean squared degree sched_emp actually ran
    # at, so it is only definable once wave 1 is on disk.
    for i, orb in enumerate(orbits):
        sidecar, _ = paths(i, "sched_emp", "tighter", smoke)
        if not sidecar.exists():
            specs[i]["n_work"] = None
            continue
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        mds = meta["telemetry"].get("mean_degree_sq")
        specs[i]["n_work"] = int(round(math.sqrt(mds))) if mds else None

    wave2 = [{"index": i, "orbit": orbits[i], "spec": specs[i],
              "policy": p, "level": l, "duration": duration, "smoke": smoke,
              "script_sha": script_sha}
             for i in range(len(orbits)) if specs[i].get("n_work")
             for p in WAVE2_POLICIES for l in LEVELS]
    done2, fail2 = submit(wave2, "wave 2 (work-matched fixed)")

    rows = [orbit_summary(i, orbits[i], specs[i], smoke, duration)
            for i in range(len(orbits))]
    payload = {
        "schema": "r17_longarc60_v1",
        "created_utc": base.utc_now(),
        "scenario": {
            "purpose": "sixty-day arcs on a widened geometry set with a "
                       "two-level resolution envelope",
            "geometry_rule": GEOMETRY_RULE,
            "excluded": "geometries with perilune below 50 km, which need an "
                        "N=600 truth; the LRO-like 30x216 km case of the "
                        "28-day stage is therefore not extended, and it is "
                        "the one geometry that came closest to favoring a "
                        "schedule at 28 days",
            "duration_days": duration / DAY,
            "output_step_s": OUTPUT_STEP, "max_step_s": MAX_STEP,
            "truth_degree": TRUTH_DEGREE,
            "checkpoint_days": list(CHECKPOINT_DAYS),
            "levels": LEVELS,
            "resolution_rule": "a comparison resolves when the absolute error "
                               "difference exceeds the sum of the two "
                               "truth-inclusive envelopes",
            "policy_definitions": "unchanged from the 28-day stage "
                                  "(rev7_doe_screening tables)",
            "smoke": smoke,
        },
        "provenance": base.provenance(),
        "script_sha256": script_sha,
        "orbits": orbits,
        "rows": rows,
        "summary": summarize(rows),
        "failures": fail1 + fail2,
        "wall_minutes": (time.time() - t_start) / 60.0,
    }
    base.atomic_json(out_path, payload)
    s = payload["summary"]
    print(f"\n[r17] written {out_path.name}  "
          f"wall={payload['wall_minutes']:.1f} min", flush=True)
    print(f"[r17] {s['orbits_reaching_60_days']}/{s['orbits_attempted']} "
          f"orbits reached 60 days; {s['orbits_terminated_early']} terminated "
          f"early", flush=True)
    return 1 if (fail1 or fail2) else 0


def command_summarize(smoke: bool, from_disk: bool) -> int:
    """Rebuild the record. --from-disk reconstructs it from the trajectory
    sidecars alone, so a run stopped by its deadline (or by hand) still yields a
    complete record over the orbits that did finish, rather than nothing."""
    out_path = SMOKE_OUTPUT if smoke else OUTPUT
    if from_disk:
        script_sha = base.file_hash(Path(__file__).resolve())
        orbits = geometries()
        model, _ = model_args()
        power = _MODEL["power"]
        duration = 2.0 * DAY if smoke else DURATION
        rows = []
        for i, orb in enumerate(orbits):
            spec = {"n_crit": int(min(CAP, emp_nmin_exact(
                power, model.r_ref, orb["hp_km"] * 1e3)))}
            sidecar, _ = paths(i, "sched_emp", "tighter", smoke)
            if sidecar.exists():
                meta = json.loads(sidecar.read_text(encoding="utf-8"))
                mds = meta["telemetry"].get("mean_degree_sq")
                spec["n_work"] = int(round(math.sqrt(mds))) if mds else None
            rows.append(orbit_summary(i, orb, spec, smoke, duration))
        payload = {
            "schema": "r17_longarc60_v1",
            "created_utc": base.utc_now(),
            "reconstructed_from_disk": True,
            "scenario": {
                "purpose": "sixty-day arcs on a widened geometry set with a "
                           "two-level resolution envelope",
                "geometry_rule": GEOMETRY_RULE,
                "duration_days": duration / DAY,
                "output_step_s": OUTPUT_STEP, "max_step_s": MAX_STEP,
                "truth_degree": TRUTH_DEGREE,
                "checkpoint_days": list(CHECKPOINT_DAYS),
                "levels": LEVELS,
                "resolution_rule": "a comparison resolves when the absolute "
                                   "error difference exceeds the sum of the "
                                   "two truth-inclusive envelopes",
                "completion_note": "orbits whose trajectory set is incomplete "
                                   "are recorded with status 'incomplete' and "
                                   "excluded from every aggregate",
                "smoke": smoke,
            },
            "provenance": base.provenance(),
            "script_sha256": script_sha,
            "orbits": orbits, "rows": rows,
            "summary": summarize(rows),
        }
    else:
        payload = json.loads(out_path.read_text(encoding="utf-8"))
        payload["summary"] = summarize(payload["rows"])
    base.atomic_json(out_path, payload)
    print(json.dumps(payload["summary"], indent=2))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run", "smoke", "summarize"))
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline-min", type=float, default=210.0)
    ap.add_argument("--from-disk", action="store_true",
                    help="rebuild the record from the trajectory sidecars, for "
                         "a run stopped before it wrote its own output")
    a = ap.parse_args()
    if a.command == "summarize":
        return command_summarize(False, a.from_disk)
    return run(a.command == "smoke", a.workers, a.deadline_min)


if __name__ == "__main__":
    raise SystemExit(main())
