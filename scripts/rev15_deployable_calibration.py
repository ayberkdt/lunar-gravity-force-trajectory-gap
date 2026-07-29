"""Deployable (truth-free) budget calibration (R15-B).

In O25 the tolerance that meets a budget is found by bisecting on the altitude
history of the *archived truth trajectory*. The online degree decision is still a
function of instantaneous radius, but the calibration knows the arc's radial
dwell before the arc is flown, which no operational user does. This replaces that
step with two calibrations that use only information available in advance:

  kepler  -- the two-body altitude history implied by the initial state alone.
             No propagation at all; the cheapest thing a user could do.
  pilot   -- the altitude history of a cheap low-degree, loose-tolerance pilot
             arc (N = 40); measured at a median 8% of the real run's work.

Each calibrated policy is then propagated under the same contract as O25 and
compared with the same fixed comparator. Three things are reported: how far each
deployable calibration lands from the requested budget once measured on the true
arc, how much of that is a calibration error rather than an integrator effect,
and---the question that matters---whether the trajectory verdict changes.

The truth-informed calibration is retained as the oracle upper bound on budget
adherence.

Usage:
    python rev15_deployable_calibration.py run --design A --workers 5
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
import rev12_atallah as at
from rev14_budget_pareto import (DESIGNS, LEVEL, BIN_KM, FLOOR,
                                 calibrate_tolerance, degrees_from_table,
                                 fixed_degree_for, _model, _g)
from rev14_budget_trajectory import LEVELS, MAX_STEP, DURATION, OUTPUT_STEP

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r15_cases"
RAW_ROOT = METRICS / "r15_raw"
BETA = 1.00
PILOT_DEGREE = 40
PILOT_RTOL = 1.0e-9
PILOT_ATOL = np.array([1.0e-3] * 3 + [1.0e-6] * 3)
VARIANTS = ("kepler", "pilot")


def kepler_altitudes(model, y0, t_grid):
    """Two-body altitude history from the initial state alone (no propagation)."""
    mu = model.mu
    r0, v0 = np.asarray(y0[:3], float), np.asarray(y0[3:], float)
    r0n = np.linalg.norm(r0)
    a = 1.0 / (2.0 / r0n - (v0 @ v0) / mu)
    h = np.cross(r0, v0)
    e_vec = np.cross(v0, h) / mu - r0 / r0n
    e = float(np.linalg.norm(e_vec))
    n = math.sqrt(mu / a ** 3)
    # eccentric anomaly at t0 from the r,v pair
    E0 = math.atan2((r0 @ v0) / (a ** 2 * n), 1.0 - r0n / a)
    M0 = E0 - e * math.sin(E0)
    out = np.empty(len(t_grid))
    for i, t in enumerate(t_grid):
        M = M0 + n * float(t)
        E = M if e < 0.8 else math.pi
        for _ in range(60):                     # Newton on Kepler's equation
            f = E - e * math.sin(E) - M
            fp = 1.0 - e * math.cos(E)
            dE = -f / fp
            E += dE
            if abs(dE) < 1e-14:
                break
        out[i] = a * (1.0 - e * math.cos(E))
    return (out - model.r_ref) / 1e3            # altitude [km]


def worker(task: dict) -> dict:
    design, row, spec, prov = (task["design"], task["row"], task["spec"],
                               task["provenance"])
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        g = _g(adopted)
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        target = BETA * n_crit ** 2

        # --- altitude histories available without the truth trajectory
        hist = {"kepler": kepler_altitudes(model, y0, grid)}
        t_p, y_p, st_p, _, fail_p, tel_p = base.propagate_event_instrumented(
            model, y0, DURATION, grid, lambda _t, _h, n=PILOT_DEGREE: n, args,
            PILOT_RTOL, PILOT_ATOL, max_step=MAX_STEP)
        if st_p == "numerical_failure":
            return {"index": index, "status": "pilot_failure", "message": fail_p}
        hist["pilot"] = ((np.linalg.norm(y_p[0:3, :].T, axis=1) - model.r_ref)
                         / 1e3)
        pilot_cost = {"n_rhs": int(tel_p["n_rhs"]),
                      "degree": PILOT_DEGREE,
                      "quadratic_work": float(tel_p["n_rhs"]) * PILOT_DEGREE ** 2}

        # --- truth-informed altitude history, for the realized-budget reference
        t_ref, y_ref = base.load_raw(DESIGNS[design]["r11_raw"]
                                     / f"sobolA_{index:03d}" / f"truth_{LEVEL}.npz")
        h_true = (np.linalg.norm(y_ref[0:3, :].T, axis=1) - model.r_ref) / 1e3

        out = {"index": index, "status": "complete", "n_critical": n_crit,
               "hp_km": hp_km, "adopted_truth_degree": adopted,
               "pilot_cost": pilot_cost, "variants": {}}
        for name in VARIANTS:
            cal = calibrate_tolerance(model, g, hp_km, ha_km, adopted,
                                      hist[name], target)
            deg_fn, table = at.atallah_binned_schedule(
                model, g, cal["tol"], hp_km, ha_km, floor=FLOOR, cap=adopted,
                bin_km=BIN_KM)
            tab = {float(k): int(v) for k, v in table.items()}
            # what that tolerance actually buys once the real arc is flown
            realized = float(np.mean(
                degrees_from_table(tab, h_true).astype(float) ** 2))
            cfg = {"sobol_index": index, "design": design, "beta_requested": BETA,
                   "calibration": name, "adopted_truth_degree": adopted,
                   "n_critical": n_crit,
                   "initial_state_si": [float(v) for v in y0],
                   "atallah_tol_accel_m_s2": cal["tol"],
                   "calibration_source": (
                       "two-body altitude history from the initial state"
                       if name == "kepler" else
                       f"N={PILOT_DEGREE} loose-tolerance pilot arc"),
                   "atallah_degree_table": {str(k): int(v) for k, v in tab.items()},
                   "beta_on_calibration_history": cal["work"] / n_crit ** 2,
                   "beta_realized_on_true_arc": realized / n_crit ** 2,
                   "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
                   "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
                   "atol_kind": "vector", "timing_comparable": False,
                   "source": prov}
            telem = {}
            for level in ("tight", "tighter"):
                tol = LEVELS[level]
                t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                    model, y0, DURATION, grid, deg_fn, args, tol["rtol"],
                    tol["atol"], max_step=MAX_STEP)
                if st == "numerical_failure":
                    return {"index": index, "status": "numerical_failure",
                            "where": f"{name}/{level}", "message": fail}
                sub = f"{design}_{name}"
                raw = RAW_ROOT / sub / f"sobolA_{index:03d}" / f"atallah_{level}.npz"
                base.atomic_npz(raw, t_s=t, state_si=y)
                side = (CASE_ROOT / sub / f"sobolA_{index:03d}"
                        / f"atallah_{level}.json")
                c = {**cfg, "level": level, "rtol": tol["rtol"],
                     "atol_position_m": tol["atol_position_m"],
                     "atol_velocity_m_s": tol["atol_velocity_m_s"]}
                base.atomic_json(side, {
                    "schema": "r15_deployable_trajectory_v1",
                    "created_utc": base.utc_now(), "config": c,
                    "config_sha256": base.object_hash(c), "status": st,
                    "telemetry": tel,
                    "raw_path": str(raw.relative_to(ROOT)),
                    "raw_sha256": base.file_hash(raw),
                    "n_output_epochs": int(len(t))})
                telem[level] = tel
            out["variants"][name] = {
                "tol": cal["tol"],
                "beta_on_calibration_history": cal["work"] / n_crit ** 2,
                "beta_realized_on_true_arc": realized / n_crit ** 2,
                "calibration_error": realized / target - 1.0,
                "rhs_tight": int(telem["tight"]["n_rhs"]),
                "mean_degree_sq_tight": float(telem["tight"]["mean_degree_sq"]),
                "total_work_tight": (float(telem["tight"]["mean_degree_sq"])
                                     * int(telem["tight"]["n_rhs"]))}
        return out
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def summarize(design, rows):
    """Errors against the archived truth and the archived critical comparator."""
    out = []
    for r in rows:
        idx = r["index"]
        truth = {lv: base.load_raw(DESIGNS[design]["r11_raw"]
                                   / f"sobolA_{idx:03d}" / f"truth_{lv}.npz")
                 for lv in LEVELS}
        fixed = {lv: base.load_raw(DESIGNS[design]["reuse_raw_fixed"]
                                   / f"sobolA_{idx:03d}" / f"fixed_critical_{lv}.npz")
                 for lv in LEVELS}
        truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                       truth["tighter"][0],
                                       truth["tighter"][1])["pos_rms_m"]
        sd_f = base.common_error(fixed["tight"][0], fixed["tight"][1],
                                 fixed["tighter"][0], fixed["tighter"][1])["pos_rms_m"]
        e_f = base.common_error(fixed["tight"][0], fixed["tight"][1],
                                truth["tight"][0], truth["tight"][1])["pos_rms_m"]
        entry = {"sobol_index": idx, "hp_km": r["hp_km"],
                 "n_critical": r["n_critical"],
                 "fixed_error_m": e_f, "fixed_envelope_m": sd_f + truth_self,
                 "pilot_cost": r["pilot_cost"], "variants": {}}
        for name in VARIANTS:
            sub = f"{design}_{name}"
            pol = {lv: base.load_raw(RAW_ROOT / sub / f"sobolA_{idx:03d}"
                                     / f"atallah_{lv}.npz") for lv in LEVELS}
            sd_a = base.common_error(pol["tight"][0], pol["tight"][1],
                                     pol["tighter"][0], pol["tighter"][1])["pos_rms_m"]
            e_a = base.common_error(pol["tight"][0], pol["tight"][1],
                                    truth["tight"][0], truth["tight"][1])["pos_rms_m"]
            env_a = sd_a + truth_self
            diff, thr = abs(e_a - e_f), env_a + (sd_f + truth_self)
            v = dict(r["variants"][name])
            v.update({"atallah_error_m": e_a, "rho_budget": (e_f / e_a) if e_a else None,
                      "resolved": bool(diff > thr),
                      "resolution_margin": (diff / thr) if thr else None,
                      "raw_winner": "atallah" if e_a < e_f else "fixed",
                      "resolved_winner": (("atallah" if e_a < e_f else "fixed")
                                          if diff > thr else None)})
            entry["variants"][name] = v
        out.append(entry)
    return out


def stat(v):
    a = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run",))
    ap.add_argument("--design", choices=("A", "B"), default="A")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    DESIGNS["A"]["reuse_raw_fixed"] = METRICS / "r11_raw" / "convergence"
    DESIGNS["B"]["reuse_raw_fixed"] = METRICS / "r11_raw" / "designB_convergence"
    rows = json.loads(DESIGNS[a.design]["rows"].read_text())["rows"]
    if a.limit:
        rows = rows[:a.limit]
    prov = base.provenance()
    tasks = [{"design": a.design, "row": r, "spec": None, "provenance": prov}
             for r in rows]
    print(f"[deployable] design {a.design}: {len(tasks)} orbits x "
          f"{len(VARIANTS)} calibrations at beta={BETA}", flush=True)
    t0 = time.time()
    done, fails = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['index']:03d} {rec.get('message')}", flush=True)
                continue
            done.append(rec)
            if n % 8 == 0 or n == len(tasks):
                print(f"  [{n}/{len(tasks)}] elapsed={(time.time()-t0)/60:.1f}min",
                      flush=True)
    done.sort(key=lambda r: r["index"])
    rows_out = summarize(a.design, done)
    summary = {}
    for name in VARIANTS:
        vs = [r["variants"][name] for r in rows_out]
        summary[name] = {
            "orbits": len(vs),
            "beta_realized": stat([v["beta_realized_on_true_arc"] for v in vs]),
            "abs_calibration_error": stat([abs(v["calibration_error"]) for v in vs]),
            "rho_budget": stat([v["rho_budget"] for v in vs]),
            "raw_atallah_wins": sum(v["raw_winner"] == "atallah" for v in vs),
            "resolved_atallah_wins": sum(v["resolved_winner"] == "atallah" for v in vs),
            "resolved_fixed_wins": sum(v["resolved_winner"] == "fixed" for v in vs),
            "unresolved": sum(not v["resolved"] for v in vs),
        }
    summary["pilot_cost_fraction_of_run"] = stat(
        [r["pilot_cost"]["quadratic_work"]
         / max(r["variants"]["pilot"]["total_work_tight"], 1.0) for r in rows_out])
    payload = {"schema": "r15_deployable_calibration_v1",
               "created_utc": base.utc_now(), "design": a.design, "beta": BETA,
               "variants": list(VARIANTS), "pilot_degree": PILOT_DEGREE,
               "rows": rows_out, "failures": fails, "summary": summary,
               "source": prov}
    base.atomic_json(METRICS / f"r15_deployable_calibration_{a.design}.json", payload)
    for name in VARIANTS:
        s = summary[name]
        print(f"  {name}: realized beta median={s['beta_realized']['median']:.3f} "
              f"|cal err| median={s['abs_calibration_error']['median']:.3f} "
              f"max={s['abs_calibration_error']['max']:.3f} | "
              f"rho median={s['rho_budget']['median']:.3g} "
              f"resolved At {s['resolved_atallah_wins']} fix {s['resolved_fixed_wins']} "
              f"unres {s['unresolved']}", flush=True)
    print("[written] " + f"r15_deployable_calibration_{a.design}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
