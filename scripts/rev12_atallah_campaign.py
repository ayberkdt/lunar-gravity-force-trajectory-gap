"""Direct Atallah (2022) radial-adaptive benchmark on the design-A population (R12).

Phase-0 verified in rev12_atallah.py / r12_atallah_phase0_verification.txt: the
implementation satisfies the paper's own correctness criteria (Eq. 29 bound,
Eq. 31/Fig. 4 error<tol, monotonicity) on the JGGRX field.

Primary definition (per the review): perilune-fidelity-matched. For each orbit
the Atallah accuracy tolerance is frozen so that N_req(r_p) equals the orbit's
critical-altitude empirical degree N_crit. Atallah and the critical fixed degree
then start at the same lowest-altitude fidelity; only Atallah's radial
degree-reduction shape is under test. Atallah degree is capped at N_crit (it is
maximal at perilune and monotone-decreasing outward, so it never exceeds it).

Per orbit, four NEW seven-day trajectories at the R11 vector tolerances:
  atallah / tight, atallah / tighter,
  fixed N_work^A / tight, fixed N_work^A / tighter,
with N_work^A = round(sqrt(<N^2>)) taken from the TIGHT Atallah run's actual RHS
degree history, frozen, and reused for the tighter fixed run. Truth and
fixed_critical are REUSED from the R11 design-A convergence trees (identical
field, tolerances, and grid), so no truth is recomputed.

Two comparisons per orbit, with the R11 resolution rule (a comparison resolves
only when |E_A - E_B| exceeds the sum of the two truth-inclusive envelopes):
  A. Atallah vs critical-altitude fixed
  B. Atallah vs its own work-matched fixed degree

Timing is not comparable under the parallel pool; wall clocks are recorded but
no cost claim is drawn from them (a separate sequential timing panel is run
elsewhere). The measured cost proxy is the gravity-kernel time and <N^2>.

Usage:
    python rev12_atallah_campaign.py pilot --workers 4
    python rev12_atallah_campaign.py run --workers 5 --deadline 2026-07-25T08:30:00+03:00
    python rev12_atallah_campaign.py status
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at
from rev7_doe_screening import CAP

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

# Population selection via environment (inherited by spawned workers), default
# design A. Design B reuses its own r11 truth/fixed_critical trees, giving an
# independent replication of the benchmark.
_POP = os.environ.get("R12_POPULATION", "A").upper()
if _POP == "B":
    ROWS = METRICS / "r11_designB_rows.json"
    REUSE_CASE = METRICS / "r11_cases" / "designB_convergence"
    REUSE_RAW = METRICS / "r11_raw" / "designB_convergence"
    OUTPUT = METRICS / "r12_atallah_campaign_designB.json"
    PILOT_OUTPUT = METRICS / "r12_atallah_pilot_designB.json"
    CASE_ROOT = METRICS / "r12_cases" / "atallah_designB"
    RAW_ROOT = METRICS / "r12_raw" / "atallah_designB"
else:
    ROWS = METRICS / "r10_sobolA_baseline_truth_corrected.json"
    REUSE_CASE = METRICS / "r11_cases" / "convergence"
    REUSE_RAW = METRICS / "r11_raw" / "convergence"
    OUTPUT = METRICS / "r12_atallah_campaign.json"
    PILOT_OUTPUT = METRICS / "r12_atallah_pilot.json"
    CASE_ROOT = METRICS / "r12_cases" / "atallah"
    RAW_ROOT = METRICS / "r12_raw" / "atallah"

LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3),
              "atol_position_m": 1.0e-5, "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3),
                "atol_position_m": 1.0e-6, "atol_velocity_m_s": 1.0e-9},
}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0

# worker-side caches
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


def paths(index: int, policy: str, level: str):
    return (CASE_ROOT / f"sobolA_{index:03d}" / f"{policy}_{level}.json",
            RAW_ROOT / f"sobolA_{index:03d}" / f"{policy}_{level}.npz")


def reuse_paths(index: int, policy: str, level: str):
    return (REUSE_CASE / f"sobolA_{index:03d}" / f"{policy}_{level}.json",
            REUSE_RAW / f"sobolA_{index:03d}" / f"{policy}_{level}.npz")


def _propagate(model, args, y0, degree_of, level):
    tol = LEVELS[level]
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    return base.propagate_event_instrumented(
        model, np.asarray(y0), DURATION, grid, degree_of, args,
        tol["rtol"], tol["atol"], max_step=MAX_STEP)


def _save(index, policy, level, config, t, y, event, telemetry, status):
    sidecar, raw = paths(index, policy, level)
    arrays = {"t_s": t, "state_si": y}
    if event:
        arrays["impact_t_s"] = np.asarray([event["epoch_s"]])
        arrays["impact_state_si"] = np.asarray(event["state_si"])
    base.atomic_npz(raw, **arrays)
    config_sha = base.object_hash(config)
    base.atomic_json(sidecar, {
        "schema": "r12_atallah_trajectory_v1", "created_utc": base.utc_now(),
        "config": config, "config_sha256": config_sha, "status": status,
        "event": event, "telemetry": telemetry,
        "raw_path": str(raw.relative_to(ROOT)), "raw_sha256": base.file_hash(raw),
        "n_output_epochs": int(len(t)), "last_output_epoch_s": float(t[-1])})
    return config_sha


def worker(task: dict) -> dict:
    """Run the four new trajectories for one orbit."""
    row = task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        original = int(row["original_truth_degree"])
        n_crit = int(row["n_critical"])
        model, args = _model(adopted)
        g = _g(adopted)
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        r_p = model.r_ref + hp_km * 1e3
        # Perilune-ACCURACY match: tol = the actual worst-case acceleration
        # truncation error the critical fixed degree already has at perilune
        # (measured against the adopted truth degree over a lat/lon grid), NOT
        # the conservative bound. Atallah then targets that real accuracy.
        tol = at.actual_truncation_error_max(model, args, r_p, n_crit, adopted,
                                             n_lat=25, n_lon=48)
        degree_fn, atal_table = at.atallah_binned_schedule(
            model, g, tol, hp_km, ha_km, floor=2, cap=adopted, bin_km=10.0)

        base_cfg = {
            "sobol_index": index, "adopted_truth_degree": adopted,
            "original_truth_degree": original, "n_critical": n_crit,
            "initial_state_si": [float(v) for v in y0],
            "atallah_tol_accel_m_s2": tol,
            "atallah_reference": "Atallah et al. 2022, J.Astronaut.Sci. 69(3):745-766, Eq.28/20",
            "perilune_match": ("tol = actual worst-case accel truncation error of "
                               "N_crit vs adopted truth at perilune; Atallah N_req "
                               "binned on 10-km altitude bins, capped at adopted"),
            "atallah_degree_table": {str(k): int(v) for k, v in atal_table.items()},
            "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
            "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
            "atol_kind": "vector", "timing_comparable": False,
            "source": base.provenance()}

        # 1. Atallah tight
        t, y, st, ev, fail, tel = _propagate(model, args, y0, degree_fn, "tight")
        if st == "numerical_failure":
            return {"index": index, "status": "numerical_failure",
                    "where": "atallah/tight", "message": fail}
        cfg = {**base_cfg, "policy": "atallah", "level": "tight",
               "rtol": LEVELS["tight"]["rtol"],
               "atol_position_m": 1e-5, "atol_velocity_m_s": 1e-8,
               "policy_spec": {"kind": "atallah_radial_adaptive", "tol": tol}}
        _save(index, "atallah", "tight", cfg, t, y, ev, tel, st)
        n_work = int(round(math.sqrt(tel["mean_degree_sq"])))
        deg_range = tel.get("degree_range")

        # 2. Atallah tighter
        t2, y2, st2, ev2, f2, tel2 = _propagate(model, args, y0, degree_fn, "tighter")
        if st2 == "numerical_failure":
            return {"index": index, "status": "numerical_failure",
                    "where": "atallah/tighter", "message": f2}
        cfg2 = {**base_cfg, "policy": "atallah", "level": "tighter",
                "rtol": LEVELS["tighter"]["rtol"],
                "atol_position_m": 1e-6, "atol_velocity_m_s": 1e-9,
                "policy_spec": {"kind": "atallah_radial_adaptive", "tol": tol}}
        _save(index, "atallah", "tighter", cfg2, t2, y2, ev2, tel2, st2)

        # 3,4. fixed N_work^A tight + tighter
        fw_fn = lambda t, h, n=n_work: n
        for level in ("tight", "tighter"):
            tf, yf, stf, evf, ff, telf = _propagate(model, args, y0, fw_fn, level)
            if stf == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "where": f"fixed_work/{level}", "message": ff}
            cfgf = {**base_cfg, "policy": "fixed_work_atallah", "level": level,
                    "rtol": LEVELS[level]["rtol"],
                    "atol_position_m": LEVELS[level]["atol_position_m"],
                    "atol_velocity_m_s": LEVELS[level]["atol_velocity_m_s"],
                    "policy_spec": {"kind": "fixed_work", "degree": n_work,
                                    "source": "sqrt(<N^2>) of tight Atallah RHS history"}}
            _save(index, "fixed_work_atallah", level, cfgf, tf, yf, evf, telf, stf)

        return {"index": index, "status": "complete", "n_work": n_work,
                "atallah_tol": tol, "atallah_degree_range": deg_range,
                "atallah_mean_degree": tel.get("mean_degree"),
                "wall_s": (tel["total_wall_ns"] + tel2["total_wall_ns"]) / 1e9}
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ------------------------------------------------------------------ parent side
def load_rows():
    payload = json.loads(ROWS.read_text(encoding="utf-8"))
    return payload["rows"]


def _load_reuse(index, policy, level):
    _, raw = reuse_paths(index, policy, level)
    return base.load_raw(raw)


def _load_new(index, policy, level):
    _, raw = paths(index, policy, level)
    return base.load_raw(raw)


def orbit_summary(row: dict) -> dict:
    index = int(row["sobol_index"])
    # same-tolerance truth (reused) and the three compared policies
    truth = {lv: _load_reuse(index, "truth", lv) for lv in LEVELS}
    crit = {lv: _load_reuse(index, "fixed_critical", lv) for lv in LEVELS}
    atal = {lv: _load_new(index, "atallah", lv) for lv in LEVELS}
    work = {lv: _load_new(index, "fixed_work_atallah", lv) for lv in LEVELS}

    def err(pol, lv):
        return base.common_error(pol[lv][0], pol[lv][1], truth[lv][0], truth[lv][1])

    def envelope(pol):
        self_diff = base.common_error(pol["tight"][0], pol["tight"][1],
                                      pol["tighter"][0], pol["tighter"][1])["pos_rms_m"]
        truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                       truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]
        return self_diff, self_diff + truth_self

    entry = {"sobol_index": index, "name": row["name"],
             "adopted_truth_degree": int(row["adopted_truth_degree"]),
             "design_point": {k: row["design_point"][k]
                              for k in ("hp_km", "ha_km", "incl_deg", "eccentricity")},
             "n_critical": int(row["n_critical"]), "policies": {}}
    for name, pol in (("atallah", atal), ("fixed_critical", crit),
                      ("fixed_work_atallah", work)):
        e_tight = err(pol, "tight")
        self_diff, env = envelope(pol)
        entry["policies"][name] = {
            "error_tight": e_tight, "error_tighter": err(pol, "tighter"),
            "self_difference_rms_m": self_diff, "truth_inclusive_envelope_m": env}

    comps = {}
    e_at = entry["policies"]["atallah"]["error_tight"]["pos_rms_m"]
    env_at = entry["policies"]["atallah"]["truth_inclusive_envelope_m"]
    for comp in ("fixed_critical", "fixed_work_atallah"):
        e_f = entry["policies"][comp]["error_tight"]["pos_rms_m"]
        env_f = entry["policies"][comp]["truth_inclusive_envelope_m"]
        diff = abs(e_at - e_f)
        thr = env_at + env_f
        comps[f"atallah_vs_{comp}"] = {
            "atallah_error_m": e_at, "comparator_error_m": e_f,
            "rho": e_f / e_at if e_at > 0 else None,
            "absolute_error_difference_m": diff, "resolution_threshold_m": thr,
            "resolved": bool(diff > thr),
            "winner_if_resolved": (("atallah" if e_at < e_f else comp)
                                   if diff > thr else None)}
    entry["comparisons"] = comps
    return entry


def summarize(rows: list[dict]) -> dict:
    def stat(vals):
        a = np.asarray([v for v in vals if v is not None and np.isfinite(v)])
        if a.size == 0:
            return None
        return {"n": int(a.size), "median": float(np.median(a)),
                "p10": float(np.percentile(a, 10)), "p90": float(np.percentile(a, 90)),
                "min": float(a.min()), "max": float(a.max())}
    out = {"orbits": len(rows), "decisions": {}}
    for key in ("atallah_vs_fixed_critical", "atallah_vs_fixed_work_atallah"):
        cs = [r["comparisons"][key] for r in rows if key in r["comparisons"]]
        out["decisions"][key] = {
            "pairs": len(cs),
            "resolved": sum(c["resolved"] for c in cs),
            "unresolved": sum(not c["resolved"] for c in cs),
            "resolved_atallah_wins": sum(c["winner_if_resolved"] == "atallah" for c in cs),
            "resolved_fixed_wins": sum(
                c["winner_if_resolved"] not in (None, "atallah") for c in cs),
            "raw_atallah_wins": sum(c["rho"] is not None and c["rho"] > 1.0 for c in cs),
            "rho": stat([c["rho"] for c in cs])}
    return out


def parse_deadline(v):
    if not v:
        return None
    d = datetime.fromisoformat(v)
    return d.astimezone(timezone.utc) if d.tzinfo else d.replace(tzinfo=timezone.utc)


def remaining(deadline):
    return math.inf if deadline is None else (deadline - datetime.now(timezone.utc)).total_seconds()


def run(rows, workers, deadline, out_path, tag):
    for r in rows:
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / f"sobolA_{int(r['sobol_index']):03d}").mkdir(parents=True, exist_ok=True)
    tasks = [{"row": r} for r in rows]
    started = base.utc_now(); wall0 = time.perf_counter_ns()
    print(f"[{tag}] orbits={len(rows)} workers={workers} deadline={deadline}", flush=True)
    done, failures, stopped = [], [], False
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except CancelledError:
                continue
            if rec["status"] != "complete":
                failures.append(rec)
                print(f"  !! {rec['index']:03d} {rec['status']} "
                      f"{rec.get('where','')}: {rec.get('message','')}", flush=True)
            done.append(rec)
            el = time.time() - t0
            print(f"  [{n:3d}/{len(tasks)}] idx={rec['index']:03d} {rec['status']} "
                  f"Nwork={rec.get('n_work','-')} elapsed={el/3600:.2f}h", flush=True)
            if remaining(deadline) < 600.0 and not stopped:
                stopped = True
                for p in futs:
                    p.cancel()
                print("  deadline guard: cancelling queued orbits", flush=True)
    # assemble summaries for orbits whose 4 new sidecars all exist
    summaries = []
    for r in sorted(rows, key=lambda x: int(x["sobol_index"])):
        idx = int(r["sobol_index"])
        if all(paths(idx, p, lv)[0].exists()
               for p in ("atallah", "fixed_work_atallah") for lv in LEVELS):
            try:
                summaries.append(orbit_summary(r))
            except Exception as exc:
                failures.append({"index": idx, "status": "summary_error",
                                 "message": f"{type(exc).__name__}: {exc}"})
    payload = {
        "schema": "r12_atallah_campaign_v1",
        "protocol": "perilune-fidelity-matched Atallah, vector tight+tighter",
        "started_utc": started, "ended_utc": base.utc_now(),
        "complete": len(summaries) == len(rows) and not failures,
        "stopped_for_deadline": stopped, "timing_comparable": False,
        "levels": {k: {kk: vv for kk, vv in v.items() if kk != "atol"}
                   for k, v in LEVELS.items()},
        "rows": summaries, "failures": failures, "summary": summarize(summaries),
        "session_wall_ns": time.perf_counter_ns() - wall0}
    base.atomic_json(out_path, payload)
    print(f"[{tag}] done orbits={len(summaries)}/{len(rows)} "
          f"failures={len(failures)} complete={payload['complete']}", flush=True)
    return 0 if payload["complete"] else 3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("pilot", "run", "status"))
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--deadline")
    a = ap.parse_args()
    rows = load_rows()
    if a.command == "status":
        p = OUTPUT if OUTPUT.exists() else PILOT_OUTPUT
        if p.exists():
            d = json.loads(p.read_text())
            print(json.dumps({"file": p.name, "complete": d["complete"],
                              "orbits": len(d["rows"]), "failures": len(d["failures"]),
                              "summary": d["summary"]}, indent=2))
        else:
            print("no r12 output yet")
        return 0
    if a.command == "pilot":
        # lowest perilune, median perilune, high perilune, a retrograde orbit
        by_hp = sorted(rows, key=lambda r: r["design_point"]["hp_km"])
        retro = [r for r in rows if r["design_point"]["incl_deg"] > 120.0]
        picks = [by_hp[0], by_hp[len(by_hp) // 2], by_hp[-1]]
        if retro:
            picks.append(max(retro, key=lambda r: r["design_point"]["incl_deg"]))
        seen = set(); uniq = []
        for r in picks:
            if r["sobol_index"] not in seen:
                seen.add(r["sobol_index"]); uniq.append(r)
        return run(uniq, a.workers, None, PILOT_OUTPUT, "pilot")
    return run(rows, a.workers, parse_deadline(a.deadline), OUTPUT, "atallah")


if __name__ == "__main__":
    raise SystemExit(main())
