"""O25 Phase B/C: fixed-budget trajectory campaign (R14).

Propagates the budget-calibrated Atallah policy against the constant degree that
spends the same declared per-call gravity-work budget, at the pre-registered
vector tolerances, on the archived Sobol populations.

The Atallah accuracy parameter for each orbit and budget is NOT recalibrated
here: it is read from the frozen Phase-A output (rev14_budget_pareto.py), so the
propagated policy is exactly the one whose deterministic force defect was
measured, and the 10-km binned table is rebuilt from that tolerance
deterministically.

Trajectory reuse. At beta = 1 the fixed comparator N_F(1) = argmin_N |N^2 -
N_crit^2| is the critical-altitude degree itself for every orbit where the
integer rounding is exact. Where that holds, the archived R11 fixed_critical
trajectories are reused unchanged -- identical field, initial state, tolerances,
frame, max step, output grid and code identity -- and the reuse is recorded per
orbit. Where it does not hold, the comparator is propagated fresh.

Truth is always the archived R11 truth at the matching level.

Usage:
    python rev14_budget_trajectory.py run --design A --beta 1.00 --workers 5
    python rev14_budget_trajectory.py status
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PARETO = METRICS / "r14_budget_pareto.json"
PREREG = METRICS / "r14_preregistration.json"

DESIGNS = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "reuse_case": METRICS / "r11_cases" / "convergence",
          "reuse_raw": METRICS / "r11_raw" / "convergence"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "reuse_case": METRICS / "r11_cases" / "designB_convergence",
          "reuse_raw": METRICS / "r11_raw" / "designB_convergence"},
}
CASE_ROOT = METRICS / "r14_cases"
RAW_ROOT = METRICS / "r14_raw"

LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3),
              "atol_position_m": 1.0e-5, "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3),
                "atol_position_m": 1.0e-6, "atol_velocity_m_s": 1.0e-9},
}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
BIN_KM = 10.0
FLOOR = 2

_MODELS: dict[int, tuple] = {}
_GCACHE: dict[int, np.ndarray] = {}


def _model(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def _g(degree: int):
    if degree not in _GCACHE:
        m, _ = _model(degree)
        _GCACHE[degree] = at.precompute_Sn(m, degree)
    return _GCACHE[degree]


def tag_of(beta: float) -> str:
    return f"beta_{beta:.2f}"


def paths(design, beta, index, policy, level):
    sub = f"{design}_{tag_of(beta)}"
    return (CASE_ROOT / sub / f"sobolA_{index:03d}" / f"{policy}_{level}.json",
            RAW_ROOT / sub / f"sobolA_{index:03d}" / f"{policy}_{level}.npz")


def reuse_paths(design, index, policy, level):
    return (DESIGNS[design]["reuse_case"] / f"sobolA_{index:03d}" / f"{policy}_{level}.json",
            DESIGNS[design]["reuse_raw"] / f"sobolA_{index:03d}" / f"{policy}_{level}.npz")


def _propagate(model, args, y0, degree_of, level):
    tol = LEVELS[level]
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    return base.propagate_event_instrumented(
        model, np.asarray(y0), DURATION, grid, degree_of, args,
        tol["rtol"], tol["atol"], max_step=MAX_STEP)


def _save(design, beta, index, policy, level, config, t, y, event, telemetry, status):
    sidecar, raw = paths(design, beta, index, policy, level)
    arrays = {"t_s": t, "state_si": y}
    if event:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"])
    base.atomic_npz(raw, **arrays)
    base.atomic_json(sidecar, {
        "schema": "r14_budget_trajectory_v1", "created_utc": base.utc_now(),
        "config": config, "config_sha256": base.object_hash(config),
        "status": status, "event": event, "telemetry": telemetry,
        "raw_path": str(raw.relative_to(ROOT)), "raw_sha256": base.file_hash(raw),
        "n_output_epochs": int(len(t)), "last_output_epoch_s": float(t[-1])})


def worker(task: dict) -> dict:
    row, design, beta = task["row"], task["design"], task["beta"]
    spec = task["spec"]
    index = int(row["sobol_index"])
    try:
        if spec["censored"]:
            return {"index": index, "status": "censored",
                    "reason": spec["censor_reason"]}
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        g = _g(adopted)

        tol = float(spec["atallah_tol"])
        degree_fn, table = at.atallah_binned_schedule(
            model, g, tol, hp_km, ha_km, floor=FLOOR, cap=adopted, bin_km=BIN_KM)
        n_fixed = int(spec["fixed_degree"])
        reuse_fixed = bool(spec["reuse_fixed_critical"])

        base_cfg = {
            "sobol_index": index, "design": design, "beta_requested": beta,
            "adopted_truth_degree": adopted, "n_critical": n_crit,
            "initial_state_si": [float(v) for v in y0],
            "budget_definition": "beta = W / N_crit^2, per-call quadratic gravity work",
            "atallah_tol_accel_m_s2": tol,
            "atallah_tol_source": ("log-space bisection to <N_A^2> = beta*N_crit^2 on "
                                   "the archived truth epochs; frozen in "
                                   "r14_budget_pareto.json"),
            "atallah_reference": ("Atallah et al. 2022, J.Astronaut.Sci. "
                                  "69(3):745-766, Eq.28/20"),
            "atallah_degree_table": {str(k): int(v) for k, v in table.items()},
            "sampled_beta_achieved_atallah": spec["beta_achieved_atallah"],
            "sampled_beta_achieved_fixed": spec["beta_achieved_fixed"],
            "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
            "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
            "atol_kind": "vector", "timing_comparable": False,
            # resolved once in the parent: one git subprocess per worker
            # exhausts the Windows process table under a full pool
            "source": task["provenance"]}

        telem = {}
        for level in ("tight", "tighter"):
            t, y, st, ev, fail, tel = _propagate(model, args, y0, degree_fn, level)
            if st == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "where": f"atallah_budget/{level}", "message": fail}
            cfg = {**base_cfg, "policy": "atallah_budget", "level": level,
                   "rtol": LEVELS[level]["rtol"],
                   "atol_position_m": LEVELS[level]["atol_position_m"],
                   "atol_velocity_m_s": LEVELS[level]["atol_velocity_m_s"],
                   "policy_spec": {"kind": "budget_calibrated_atallah_radial",
                                   "tol": tol, "bin_km": BIN_KM,
                                   "floor": FLOOR, "cap": adopted}}
            _save(design, beta, index, "atallah_budget", level, cfg, t, y, ev, tel, st)
            telem[f"atallah_{level}"] = tel

        if not reuse_fixed:
            fx = lambda _t, _h, n=n_fixed: n
            for level in ("tight", "tighter"):
                t, y, st, ev, fail, tel = _propagate(model, args, y0, fx, level)
                if st == "numerical_failure":
                    return {"index": index, "status": "numerical_failure",
                            "where": f"fixed_budget/{level}", "message": fail}
                cfg = {**base_cfg, "policy": "fixed_budget", "level": level,
                       "rtol": LEVELS[level]["rtol"],
                       "atol_position_m": LEVELS[level]["atol_position_m"],
                       "atol_velocity_m_s": LEVELS[level]["atol_velocity_m_s"],
                       "policy_spec": {"kind": "fixed_budget_degree",
                                       "degree": n_fixed,
                                       "source": "argmin_N |N^2 - beta*N_crit^2|"}}
                _save(design, beta, index, "fixed_budget", level, cfg, t, y, ev, tel, st)
                telem[f"fixed_{level}"] = tel

        return {"index": index, "status": "complete",
                "atallah_tol": tol, "fixed_degree": n_fixed,
                "reuse_fixed_critical": reuse_fixed,
                "atallah_mean_degree": telem["atallah_tight"].get("mean_degree"),
                "atallah_degree_range": telem["atallah_tight"].get("degree_range"),
                "wall_s": sum(v["total_wall_ns"] for v in telem.values()) / 1e9}
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ------------------------------------------------------------------ parent side
def build_specs(design: str, beta: float) -> dict:
    """Read the frozen Phase-A calibration for this design and budget."""
    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    rows = pareto["designs"][design]["rows"]
    key = tag_of(beta)
    out = {}
    for r in rows:
        b = r["budgets"][key]
        n_fixed = int(b["fixed"]["degree"])
        n_crit = int(r["n_critical"])
        reason = []
        if not b["atallah"]["attainable"]:
            reason.append(f"atallah budget unattainable ({b['atallah']['limit']})")
        if not b["fixed"]["attainable"]:
            reason.append(f"fixed comparator above adopted truth ({b['fixed']['limit']})")
        out[int(r["sobol_index"])] = {
            "atallah_tol": b["atallah"]["tol_accel_m_s2"],
            "fixed_degree": n_fixed,
            "beta_achieved_atallah": b["atallah"]["beta_achieved"],
            "beta_achieved_fixed": b["fixed"]["beta_achieved"],
            "work_mismatch_atallah": b["atallah"]["work_mismatch"],
            # exact-reuse condition: the same integer degree, same everything else
            "reuse_fixed_critical": bool(abs(beta - 1.0) < 1e-12 and n_fixed == n_crit),
            "censored": bool(b["censored"]),
            "censor_reason": "; ".join(reason) or None,
        }
    return out


def _load(path):
    return base.load_raw(path)


def orbit_summary(design, beta, row, spec) -> dict:
    index = int(row["sobol_index"])
    truth = {lv: _load(reuse_paths(design, index, "truth", lv)[1]) for lv in LEVELS}
    atal = {lv: _load(paths(design, beta, index, "atallah_budget", lv)[1])
            for lv in LEVELS}
    if spec["reuse_fixed_critical"]:
        fixed = {lv: _load(reuse_paths(design, index, "fixed_critical", lv)[1])
                 for lv in LEVELS}
    else:
        fixed = {lv: _load(paths(design, beta, index, "fixed_budget", lv)[1])
                 for lv in LEVELS}

    def err(pol, lv):
        return base.common_error(pol[lv][0], pol[lv][1], truth[lv][0], truth[lv][1])

    truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                   truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]

    def envelope(pol):
        sd = base.common_error(pol["tight"][0], pol["tight"][1],
                               pol["tighter"][0], pol["tighter"][1])["pos_rms_m"]
        return sd, sd + truth_self

    entry = {"sobol_index": index, "design": design, "beta": beta,
             "name": row["name"],
             "design_point": {k: row["design_point"][k]
                              for k in ("hp_km", "ha_km", "incl_deg", "eccentricity")},
             "n_critical": int(row["n_critical"]),
             "adopted_truth_degree": int(row["adopted_truth_degree"]),
             "atallah_tol_accel_m_s2": spec["atallah_tol"],
             "fixed_degree": spec["fixed_degree"],
             "reuse_fixed_critical": spec["reuse_fixed_critical"],
             "sampled_beta_achieved": {"atallah": spec["beta_achieved_atallah"],
                                       "fixed": spec["beta_achieved_fixed"]},
             "truth_self_difference_rms_m": truth_self,
             "policies": {}}
    for name, pol in (("atallah_budget", atal), ("fixed_budget", fixed)):
        sd, env = envelope(pol)
        entry["policies"][name] = {
            "error_tight": err(pol, "tight"), "error_tighter": err(pol, "tighter"),
            "self_difference_rms_m": sd, "truth_inclusive_envelope_m": env}

    e_a = entry["policies"]["atallah_budget"]["error_tight"]["pos_rms_m"]
    e_f = entry["policies"]["fixed_budget"]["error_tight"]["pos_rms_m"]
    env_a = entry["policies"]["atallah_budget"]["truth_inclusive_envelope_m"]
    env_f = entry["policies"]["fixed_budget"]["truth_inclusive_envelope_m"]
    diff, thr = abs(e_a - e_f), env_a + env_f
    entry["comparison"] = {
        "atallah_error_m": e_a, "fixed_error_m": e_f,
        "rho_budget": (e_f / e_a) if e_a > 0 else None,
        "absolute_error_difference_m": diff, "resolution_threshold_m": thr,
        "resolution_margin": (diff / thr) if thr > 0 else None,
        "raw_winner": "atallah" if e_a < e_f else "fixed",
        "resolved": bool(diff > thr),
        "resolved_winner": (("atallah" if e_a < e_f else "fixed")
                            if diff > thr else None)}

    # RHS / total-work bookkeeping (section 15 of the protocol)
    def telemetry(policy, reuse):
        p = (reuse_paths(design, index, "fixed_critical", "tight")[0] if reuse
             else paths(design, beta, index, policy, "tight")[0])
        return json.loads(p.read_text(encoding="utf-8"))["telemetry"]

    ta = telemetry("atallah_budget", False)
    tf = telemetry("fixed_budget", spec["reuse_fixed_critical"])
    n_a, n_f = int(ta["n_rhs"]), int(tf["n_rhs"])
    w_a = float(ta["mean_degree_sq"]) * n_a
    w_f = float(spec["fixed_degree"] ** 2) * n_f
    entry["cost"] = {
        "rhs_atallah": n_a, "rhs_fixed": n_f,
        "rhs_ratio_atallah_over_fixed": n_a / n_f,
        "per_call_mean_degree_sq_atallah": float(ta["mean_degree_sq"]),
        "per_call_mean_degree_sq_fixed": float(spec["fixed_degree"] ** 2),
        "per_call_work_ratio": float(ta["mean_degree_sq"]) / (spec["fixed_degree"] ** 2),
        "total_quadratic_work_atallah": w_a,
        "total_quadratic_work_fixed": w_f,
        "total_work_ratio_atallah_over_fixed": w_a / w_f,
        "accepted_steps_atallah": int(ta["n_accepted_steps"]),
        "accepted_steps_fixed": int(tf["n_accepted_steps"]),
        "rejected_trials_atallah": int(ta["n_rejected_trials"]),
        "rejected_trials_fixed": int(tf["n_rejected_trials"]),
        "switch_count_atallah": int(ta["switch_count_at_rhs_samples"]),
        "gravity_kernel_ns_atallah": int(ta["gravity_kernel_ns"]),
        "gravity_kernel_ns_fixed": int(tf["gravity_kernel_ns"]),
        "gravity_kernel_time_comparable": False,
    }
    return entry


def stat(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def summarize(rows: list[dict]) -> dict:
    cs = [r["comparison"] for r in rows]
    return {
        "orbits": len(rows),
        "raw_atallah_wins": sum(c["raw_winner"] == "atallah" for c in cs),
        "raw_fixed_wins": sum(c["raw_winner"] == "fixed" for c in cs),
        "resolved": sum(c["resolved"] for c in cs),
        "unresolved": sum(not c["resolved"] for c in cs),
        "resolved_atallah_wins": sum(c["resolved_winner"] == "atallah" for c in cs),
        "resolved_fixed_wins": sum(c["resolved_winner"] == "fixed" for c in cs),
        "rho_budget": stat([c["rho_budget"] for c in cs]),
        "resolution_margin": stat([c["resolution_margin"] for c in cs]),
        "atallah_error_m": stat([c["atallah_error_m"] for c in cs]),
        "fixed_error_m": stat([c["fixed_error_m"] for c in cs]),
        "rhs_ratio": stat([r["cost"]["rhs_ratio_atallah_over_fixed"] for r in rows]),
        "total_work_ratio": stat([r["cost"]["total_work_ratio_atallah_over_fixed"]
                                  for r in rows]),
        "per_call_work_ratio": stat([r["cost"]["per_call_work_ratio"] for r in rows]),
        "reused_fixed_critical": sum(r["reuse_fixed_critical"] for r in rows),
    }


def parse_deadline(v):
    if not v:
        return None
    d = datetime.fromisoformat(v)
    return d.astimezone(timezone.utc) if d.tzinfo else d.replace(tzinfo=timezone.utc)


def remaining(deadline):
    return math.inf if deadline is None else (
        deadline - datetime.now(timezone.utc)).total_seconds()


def run(design, beta, workers, deadline, limit) -> int:
    rows = json.loads(DESIGNS[design]["rows"].read_text(encoding="utf-8"))["rows"]
    specs = build_specs(design, beta)
    if limit:
        rows = rows[:limit]
    out_path = METRICS / f"r14_trajectory_{design}_{tag_of(beta)}.json"
    prov = base.provenance()
    tasks = [{"row": r, "design": design, "beta": beta, "provenance": prov,
              "spec": specs[int(r["sobol_index"])]} for r in rows]
    live = [t for t in tasks if not t["spec"]["censored"]]
    print(f"[traj] design {design} beta={beta} orbits={len(live)} "
          f"(censored {len(tasks)-len(live)}) workers={workers}", flush=True)
    reused = sum(t["spec"]["reuse_fixed_critical"] for t in live)
    print(f"       fixed comparator reused from R11 fixed_critical on {reused}/{len(live)}",
          flush=True)
    started = base.utc_now()
    done, failures, stopped = [], [], False
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except CancelledError:
                continue
            if rec["status"] not in ("complete", "censored"):
                failures.append(rec)
                print(f"  !! {rec['index']:03d} {rec['status']} "
                      f"{rec.get('where','')}: {rec.get('message','')}", flush=True)
            done.append(rec)
            print(f"  [{n:3d}/{len(tasks)}] idx={rec['index']:03d} {rec['status']} "
                  f"elapsed={(time.time()-t0)/3600:.2f}h", flush=True)
            if remaining(deadline) < 900.0 and not stopped:
                stopped = True
                for p in futs:
                    p.cancel()
                print("  deadline guard: cancelling queued orbits", flush=True)

    summaries, censored = [], []
    for r in sorted(rows, key=lambda x: int(x["sobol_index"])):
        idx = int(r["sobol_index"])
        sp = specs[idx]
        if sp["censored"]:
            censored.append({"sobol_index": idx, "reason": sp["censor_reason"]})
            continue
        if all(paths(design, beta, idx, "atallah_budget", lv)[0].exists()
               for lv in LEVELS):
            try:
                summaries.append(orbit_summary(design, beta, r, sp))
            except Exception as exc:
                failures.append({"index": idx, "status": "summary_error",
                                 "message": f"{type(exc).__name__}: {exc}",
                                 "traceback": traceback.format_exc()})
    payload = {
        "schema": "r14_budget_trajectory_v1",
        "design": design, "beta": beta,
        "protocol": ("budget-calibrated Atallah vs equal-budget constant degree, "
                     "vector tight+tighter, truth-inclusive pairwise resolution"),
        "preregistration_sha256": json.loads(
            PREREG.read_text(encoding="utf-8"))["protocol_sha256"],
        "started_utc": started, "ended_utc": base.utc_now(),
        "complete": len(summaries) + len(censored) == len(rows) and not failures,
        "stopped_for_deadline": stopped, "timing_comparable": False,
        "levels": {k: {kk: vv for kk, vv in v.items() if kk != "atol"}
                   for k, v in LEVELS.items()},
        "censored": censored, "rows": summaries, "failures": failures,
        "summary": summarize(summaries), "source": base.provenance()}
    base.atomic_json(out_path, payload)
    s = payload["summary"]
    print(f"[traj] {out_path.name}: orbits={s['orbits']} censored={len(censored)} "
          f"raw At {s['raw_atallah_wins']}/{s['orbits']} "
          f"resolved At {s['resolved_atallah_wins']} fixed {s['resolved_fixed_wins']} "
          f"unresolved {s['unresolved']} "
          f"rho median={s['rho_budget']['median'] if s['rho_budget'] else float('nan'):.3g}",
          flush=True)
    return 0 if payload["complete"] else 3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run", "status"))
    ap.add_argument("--design", choices=("A", "B"), default="A")
    ap.add_argument("--beta", type=float, default=1.00)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--deadline")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.command == "status":
        for p in sorted(METRICS.glob("r14_trajectory_*.json")):
            d = json.loads(p.read_text(encoding="utf-8"))
            s = d["summary"]
            print(f"{p.name}: complete={d['complete']} orbits={s['orbits']} "
                  f"resolvedAt={s['resolved_atallah_wins']} "
                  f"resolvedFix={s['resolved_fixed_wins']} "
                  f"unresolved={s['unresolved']}")
        return 0
    return run(a.design, a.beta, a.workers, parse_deadline(a.deadline), a.limit)


if __name__ == "__main__":
    raise SystemExit(main())
