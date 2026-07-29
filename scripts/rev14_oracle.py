"""O26: trajectory-informed force-defect allocation bound (R14).

Answers a question the Atallah-versus-fixed comparison cannot: of the achievable
local field accuracy at a given degree budget, how much does the published
radial rule actually capture?

Along each orbit's archived reference trajectory the squared truncation defect

    d(t, N) = || a_N(x_ref(t)) - a_truth(x_ref(t)) ||^2

is tabulated on a degree grid, and the allocation minimizing sum_t d(t, N_t)
subject to sum_t N_t^2 <= B is found through its Lagrangian form

    N_lambda(t) = argmin_N [ d(t, N) + lambda N^2 ],

with lambda bisected to meet the budget. Because d(t, N) is tabulated once per
orbit, every budget on the grid is then free.

This is an ORACLE, not a controller: it reads the reference trajectory and the
local truth field at every epoch, neither of which a propagator has. It bounds
the best instantaneous force defect purchasable with a given degree budget, and
nothing else. It is not flight-realizable, and it carries no implication about
trajectory error -- long-arc displacement is the state-transition-weighted
integral of the defect, so a smaller defect everywhere does not order the
trajectories.

The degree grid is coarse above degree 100. A grid-restricted optimum can only
be worse than the continuous optimum, so the measured Atallah-to-oracle gap is a
conservative (under-)estimate of the true gap.

Usage:
    python rev14_oracle.py --design both --workers 5 --stride 8
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
from rev14_budget_pareto import (DESIGNS, LEVEL, BIN_KM, FLOOR, accel_inertial,
                                 binned_table, degrees_from_table,
                                 calibrate_tolerance, fixed_degree_for, stat,
                                 _model, _g)

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r14_oracle.json"
PREREG = METRICS / "r14_preregistration.json"


def degree_grid(cap: int) -> np.ndarray:
    """Candidate degrees: fine where the defect changes fastest, coarse above."""
    parts = [np.arange(2, min(cap, 40) + 1, 2),
             np.arange(40, min(cap, 100) + 1, 5),
             np.arange(100, cap + 1, 10), np.array([cap])]
    return np.unique(np.concatenate([p for p in parts if p.size]))


def worker(task: dict) -> dict:
    design, row, betas, stride = (task["design"], task["row"], task["betas"],
                                  task["stride"])
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        model, args = _model(adopted)
        g = _g(adopted)
        t_all, y_all = base.load_raw(DESIGNS[design]["r11_raw"]
                                     / f"sobolA_{index:03d}" / f"truth_{LEVEL}.npz")
        sl = slice(None, None, stride)
        t, y = t_all[sl], y_all[:, sl]
        r_all = y[0:3, :].T
        h_km = (np.linalg.norm(r_all, axis=1) - model.r_ref) / 1e3
        n_ep = len(t)
        w_crit = float(n_crit ** 2)

        # ---- calibrate every policy first, so the degree grid can contain the
        # exact degrees they use: the policies are then evaluated at their own
        # degrees, never snapped onto a coarse grid, and only the oracle is
        # grid-restricted (which can only understate its advantage).
        cals, fixes = {}, {}
        for beta in betas:
            cals[beta] = calibrate_tolerance(model, g, hp_km, ha_km, adopted,
                                             h_km, beta * w_crit)
            fixes[beta] = fixed_degree_for(beta, n_crit, adopted)
        policy_degrees = np.unique(np.concatenate(
            [c["degrees"] for c in cals.values()]
            + [np.array([f[0] for f in fixes.values()])]))
        grid = np.unique(np.concatenate([degree_grid(adopted), policy_degrees]))

        # ---- tabulate d(t, N) once
        d2 = np.empty((n_ep, len(grid)))
        for j in range(n_ep):
            rj, tj = r_all[j], float(t[j])
            a_truth = accel_inertial(rj, tj, adopted, args)
            for i, n in enumerate(grid):
                dv = accel_inertial(rj, tj, int(n), args) - a_truth
                d2[j, i] = float(dv @ dv)
        cost = grid.astype(float) ** 2

        def rms_of(deg_idx):
            return float(np.sqrt(np.mean(d2[np.arange(n_ep), deg_idx])))

        def exact_idx(deg: np.ndarray):
            """Index of each degree in the grid; every policy degree is in it."""
            i = np.searchsorted(grid, deg)
            if not np.array_equal(grid[i], np.asarray(deg)):
                raise RuntimeError("policy degree missing from the tabulated grid")
            return i

        def oracle_for(target_work):
            """Bisect lambda so the Lagrangian allocation meets the budget."""
            def alloc(lam):
                idx = np.argmin(d2 + lam * cost[None, :], axis=1)
                return idx, float(np.mean(grid[idx].astype(float) ** 2))
            lo, hi = -30.0, 10.0            # log10 lambda
            i_hi, w_hi = alloc(10.0 ** hi)   # heavy cost penalty -> lowest degree
            i_lo, w_lo = alloc(10.0 ** lo)
            if target_work >= w_lo:
                return i_lo, w_lo, False
            if target_work <= w_hi:
                return i_hi, w_hi, False
            best = None
            for _ in range(120):
                mid = 0.5 * (lo + hi)
                idx, w = alloc(10.0 ** mid)
                err = abs(w / target_work - 1.0)
                if best is None or err < best[2]:
                    best = (idx, w, err)
                if w > target_work:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-12:
                    break
            return best[0], best[1], bool(best[2] < 0.05)

        out = {}
        for beta in betas:
            target = beta * w_crit
            cal = cals[beta]
            n_f, f_cens = fixes[beta]
            i_or, w_or, ok = oracle_for(target)
            i_at = exact_idx(cal["degrees"])
            i_fx = exact_idx(np.full(n_ep, n_f))
            e_or, e_at, e_fx = rms_of(i_or), rms_of(i_at), rms_of(i_fx)
            out[f"beta_{beta:.2f}"] = {
                "beta_requested": beta,
                "censored": bool(f_cens or not cal["attainable"] or not ok),
                "oracle": {"defect_rms_m_s2": e_or,
                           "achieved_work": w_or, "beta_achieved": w_or / w_crit,
                           "mean_degree": float(grid[i_or].mean()),
                           "degree_range": [int(grid[i_or].min()), int(grid[i_or].max())],
                           "budget_met": ok},
                "atallah": {"defect_rms_m_s2": e_at,
                            "achieved_work": cal["work"],
                            "beta_achieved": cal["work"] / w_crit,
                            "mean_degree": float(cal["degrees"].mean())},
                "fixed": {"defect_rms_m_s2": e_fx, "degree": n_f,
                          "achieved_work": float(n_f ** 2),
                          "beta_achieved": (n_f ** 2) / w_crit},
                # eta < 1: the oracle achieves a smaller defect at the same budget
                "eta_atallah": (e_or / e_at) if e_at > 0 else None,
                "eta_fixed": (e_or / e_fx) if e_fx > 0 else None,
                "atallah_penalty_over_oracle": (e_at / e_or) if e_or > 0 else None,
                "fixed_penalty_over_oracle": (e_fx / e_or) if e_or > 0 else None,
            }
        return {"design": design, "sobol_index": index, "status": "complete",
                "adopted_truth_degree": adopted, "n_critical": n_crit,
                "hp_km": hp_km, "ha_km": ha_km, "n_epochs": n_ep,
                "stride": stride, "degree_grid_size": int(len(grid)),
                "budgets": out}
    except Exception as exc:
        return {"design": design, "sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def summarize(rows, key):
    live = [r["budgets"][key] for r in rows if not r["budgets"][key]["censored"]]
    cens = sum(r["budgets"][key]["censored"] for r in rows)
    if not live:
        return {"orbits": 0, "censored": cens}
    return {
        "orbits": len(live), "censored": cens,
        "beta_requested": live[0]["beta_requested"],
        "oracle_defect_rms_m_s2": stat([e["oracle"]["defect_rms_m_s2"] for e in live]),
        "atallah_defect_rms_m_s2": stat([e["atallah"]["defect_rms_m_s2"] for e in live]),
        "fixed_defect_rms_m_s2": stat([e["fixed"]["defect_rms_m_s2"] for e in live]),
        "atallah_penalty_over_oracle": stat(
            [e["atallah_penalty_over_oracle"] for e in live]),
        "fixed_penalty_over_oracle": stat([e["fixed_penalty_over_oracle"] for e in live]),
        "oracle_beta_achieved": stat([e["oracle"]["beta_achieved"] for e in live]),
        "atallah_closer_than_fixed": int(sum(
            e["atallah_penalty_over_oracle"] < e["fixed_penalty_over_oracle"]
            for e in live)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B", "both"), default="both")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--stride", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    betas = list(prereg["budget_grid"])
    designs = ["A", "B"] if a.design == "both" else [a.design]
    payload = {"schema": "r14_oracle_v1", "created_utc": base.utc_now(),
               "reference_level": LEVEL, "budget_grid": betas,
               "epoch_stride": a.stride,
               "preregistration_sha256": prereg["protocol_sha256"],
               "nature": ("trajectory-informed force-defect allocation bound; an "
                          "oracle, not a deployable controller, and not a "
                          "trajectory-optimality statement"),
               "source": base.provenance(), "designs": {}}
    if OUTPUT.exists():
        try:
            payload["designs"] = json.loads(OUTPUT.read_text())["designs"]
        except Exception:
            payload["designs"] = {}
    for d in designs:
        rows = json.loads(DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]
        if a.limit:
            rows = rows[:a.limit]
        tasks = [{"design": d, "row": r, "betas": betas, "stride": a.stride}
                 for r in rows]
        print(f"[oracle] design {d}: {len(tasks)} orbits", flush=True)
        t0 = time.time()
        done, fails = [], []
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(worker, t): t for t in tasks}
            for n, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                if rec["status"] != "complete":
                    fails.append(rec)
                    print(f"  !! {rec['sobol_index']:03d} {rec['message']}", flush=True)
                    continue
                done.append(rec)
                if n % 8 == 0 or n == len(tasks):
                    print(f"  [{n:3d}/{len(tasks)}] elapsed={(time.time()-t0)/60:.1f}min",
                          flush=True)
        done.sort(key=lambda r: r["sobol_index"])
        keys = [f"beta_{b:.2f}" for b in betas]
        payload["designs"][d] = {"rows": done, "failures": fails,
                                 "summary": {k: summarize(done, k) for k in keys}}
        for k in keys:
            e = payload["designs"][d]["summary"][k]
            if not e.get("orbits"):
                continue
            print(f"  {k}: n={e['orbits']} At/oracle median="
                  f"{e['atallah_penalty_over_oracle']['median']:.3g} "
                  f"fixed/oracle median={e['fixed_penalty_over_oracle']['median']:.3g} "
                  f"At closer {e['atallah_closer_than_fixed']}/{e['orbits']}", flush=True)
    base.atomic_json(OUTPUT, payload)
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
