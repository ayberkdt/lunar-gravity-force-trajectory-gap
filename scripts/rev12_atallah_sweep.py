"""Accuracy-tolerance sweep of the Atallah benchmark (R12), the paper's own
Figs 5-7 methodology: run Atallah at several physical user tolerances and read
off the cost-accuracy trade-off against a work-matched fixed degree.

Complements the perilune-tolerance-matched primary campaign: instead of one
matched tol per orbit, it sweeps a fixed set of absolute-acceleration tolerances
across the whole population, so the trade-off between Atallah's degree reduction
and its accuracy is characterized across accuracy levels rather than at a single
operating point. Tight level only (the trend, not a resolved verdict; the
primary campaign carries the resolved comparisons). Errors are against the
reused R11 tight-tolerance truth; cost is mean <N^2> and measured kernel time.

Usage:
    python rev12_atallah_sweep.py run --workers 5 --deadline 2026-07-25T09:45:00+03:00
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
ROWS = METRICS / "r10_sobolA_baseline_truth_corrected.json"
REUSE_RAW = METRICS / "r11_raw" / "convergence"
OUTPUT = METRICS / "r12_atallah_sweep.json"
CASE_ROOT = METRICS / "r12_cases" / "atallah_sweep"
RAW_ROOT = METRICS / "r12_raw" / "atallah_sweep"

TOLS = [1.0e-8, 1.0e-10, 1.0e-12]   # absolute acceleration tolerance [m/s^2]
TIGHT = {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3)}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0

_MODELS: dict[int, tuple] = {}
_GCACHE: dict[int, np.ndarray] = {}


def _model(d):
    if d not in _MODELS:
        m = base.load_model(d); a = base.kernel_args(m); base.warmup(m, a)
        _MODELS[d] = (m, a)
    return _MODELS[d]


def _g(d):
    if d not in _GCACHE:
        m, _ = _model(d)
        _GCACHE[d] = at.precompute_Sn(m, d)
    return _GCACHE[d]


def _propagate(model, args, y0, degree_of):
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    return base.propagate_event_instrumented(
        model, np.asarray(y0), DURATION, grid, degree_of, args,
        TIGHT["rtol"], TIGHT["atol"], max_step=MAX_STEP)


def _reuse_truth(index):
    return base.load_raw(REUSE_RAW / f"sobolA_{index:03d}" / "truth_tight.npz")


def _save(index, tag, t, y, tel):
    raw = RAW_ROOT / f"sobolA_{index:03d}" / f"{tag}.npz"
    base.atomic_npz(raw, t_s=t, state_si=y)
    (CASE_ROOT / f"sobolA_{index:03d}" / f"{tag}.json").write_text(
        json.dumps({"telemetry": tel, "raw_sha256": base.file_hash(raw)}),
        encoding="utf-8")


def worker(task):
    row = task["row"]; index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        model, args = _model(adopted); g = _g(adopted)
        hp = float(row["design_point"]["hp_km"]); ha = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], float)
        tt, ty = _reuse_truth(index)
        out = {"index": index, "adopted": adopted, "n_critical": int(row["n_critical"]),
               "design_point": {k: row["design_point"][k]
                                for k in ("hp_km", "ha_km", "incl_deg", "eccentricity")},
               "points": []}
        for tol in TOLS:
            deg_fn, table = at.atallah_binned_schedule(
                model, g, tol, hp, ha, floor=2, cap=adopted, bin_km=10.0)
            t, y, st, ev, fail, tel = _propagate(model, args, y0, deg_fn)
            if st == "numerical_failure":
                out["points"].append({"tol": tol, "status": "numerical_failure"})
                continue
            n = min(y.shape[1], ty.shape[1])
            dp = np.linalg.norm(y[:3, :n] - ty[:3, :n], axis=0)
            e_at = float(np.sqrt(np.mean(dp ** 2)))
            n_work = int(round(math.sqrt(tel["mean_degree_sq"])))
            # work-matched fixed at this tol
            tf, yf, stf, evf, ff, telf = _propagate(model, args, y0,
                                                    lambda t, h, nn=n_work: nn)
            nf = min(yf.shape[1], ty.shape[1])
            dpf = np.linalg.norm(yf[:3, :nf] - ty[:3, :nf], axis=0)
            e_fw = float(np.sqrt(np.mean(dpf ** 2)))
            _save(index, f"atallah_tol{tol:.0e}", t, y, tel)
            out["points"].append({
                "tol": tol, "status": st,
                "atallah_pos_rms_m": e_at, "fixed_work_pos_rms_m": e_fw,
                "rho_work": e_fw / e_at if e_at > 0 else None,
                "n_work": n_work, "atallah_mean_degree": tel.get("mean_degree"),
                "atallah_mean_degree_sq": tel.get("mean_degree_sq"),
                "atallah_degree_range": tel.get("degree_range"),
                "atallah_gravity_kernel_ns": tel.get("gravity_kernel_ns"),
                "atallah_n_rhs": tel.get("n_rhs")})
        return {"index": index, "status": "complete", "record": out}
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def parse_deadline(v):
    if not v:
        return None
    d = datetime.fromisoformat(v)
    return d.astimezone(timezone.utc) if d.tzinfo else d.replace(tzinfo=timezone.utc)


def remaining(dl):
    return math.inf if dl is None else (dl - datetime.now(timezone.utc)).total_seconds()


def run(rows, workers, deadline, smoke_n=None):
    if smoke_n:
        rows = rows[:smoke_n]
    for r in rows:
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / f"sobolA_{int(r['sobol_index']):03d}").mkdir(parents=True, exist_ok=True)
    started = base.utc_now(); wall0 = time.perf_counter_ns()
    print(f"[sweep] orbits={len(rows)} tols={TOLS} workers={workers} "
          f"deadline={deadline}", flush=True)
    records, failures, stopped = [], [], False
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, {"row": r}): r for r in rows}
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                rec = fut.result()
            except CancelledError:
                continue
            if rec["status"] == "complete":
                records.append(rec["record"])
            else:
                failures.append(rec)
                print(f"  !! {rec['index']} {rec['status']}: {rec.get('message')}", flush=True)
            print(f"  [{n}/{len(rows)}] idx={rec['index']} {rec['status']} "
                  f"elapsed={(time.time()-t0)/3600:.2f}h", flush=True)
            if remaining(deadline) < 600.0 and not stopped:
                stopped = True
                for p in futs:
                    p.cancel()
                print("  deadline guard: cancelling", flush=True)
    # aggregate trade-off per tol
    agg = {}
    for tol in TOLS:
        pts = [p for r in records for p in r["points"]
               if p.get("tol") == tol and p.get("status") not in ("numerical_failure",)]
        if not pts:
            continue
        rho = np.array([p["rho_work"] for p in pts if p.get("rho_work")])
        agg[f"{tol:.0e}"] = {
            "n": len(pts),
            "atallah_rms_median_m": float(np.median([p["atallah_pos_rms_m"] for p in pts])),
            "fixed_work_rms_median_m": float(np.median([p["fixed_work_pos_rms_m"] for p in pts])),
            "rho_work_median": float(np.median(rho)) if rho.size else None,
            "rho_work_raw_atallah_wins": int(np.sum(rho > 1.0)) if rho.size else 0,
            "mean_degree_median": float(np.median([p["atallah_mean_degree"] for p in pts])),
            "mean_degree_sq_median": float(np.median([p["atallah_mean_degree_sq"] for p in pts]))}
    payload = {"schema": "r12_atallah_sweep_v1", "tols": TOLS,
               "started_utc": started, "ended_utc": base.utc_now(),
               "complete": len(records) == len(rows) and not failures,
               "stopped_for_deadline": stopped, "timing_comparable": False,
               "records": records, "failures": failures, "aggregate": agg,
               "session_wall_ns": time.perf_counter_ns() - wall0}
    base.atomic_json(OUTPUT, payload)
    print(f"[sweep] done orbits={len(records)}/{len(rows)} complete={payload['complete']}",
          flush=True)
    for k, v in agg.items():
        print(f"  tol={k}: Atallah rms {v['atallah_rms_median_m']:.3f}m, "
              f"fixed-work {v['fixed_work_rms_median_m']:.3f}m, "
              f"rho_med {v['rho_work_median']}, meanN {v['mean_degree_median']:.0f}", flush=True)
    return 0 if payload["complete"] else 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run", "smoke"))
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--deadline")
    a = ap.parse_args()
    rows = json.loads(ROWS.read_text())["rows"]
    return run(rows, a.workers, parse_deadline(a.deadline),
               smoke_n=1 if a.command == "smoke" else None)


if __name__ == "__main__":
    raise SystemExit(main())
