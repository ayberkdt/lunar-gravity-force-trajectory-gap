"""O25 Phase A: integration-noise-free fixed-budget Pareto sweep (R14).

For every orbit of both 64-point scrambled-Sobol populations, and for every
normalized gravity-work budget beta = W / N_crit^2 on the pre-registered grid,
this compares two ways of spending the same per-call budget:

  * budget-calibrated Atallah -- the published radial mapping N_req(r; eps_A)
    with its accuracy parameter eps_A calibrated by log-space bisection so that
    the 10-km-binned degree history consumes <N_A^2> = beta * N_crit^2;
  * the constant degree N_F(beta) = argmin_N |N^2 - beta * N_crit^2|.

Both are evaluated as the deterministic truncation acceleration defect along the
orbit's archived truth trajectory,

    Delta_a_P(t) = a_{N_P(t)}(x_ref(t)) - a_{N_truth}(x_ref(t)),

so the comparison carries no integration noise and reruns bit-identically.

Each orbit's archived accuracy-target Atallah operating point is carried as an
additional non-grid budget at its own measured beta, never approximated.

Budget points needing degrees above the adopted truth degree are censored: they
are flagged, counted, and excluded from that budget's population statistics,
because a defect measured against a truth the policy has reached is not a
truncation error.

Usage:
    python rev14_budget_pareto.py --design both --workers 5
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
from rev3_common import OMEGA_MOON
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r14_budget_pareto.json"
PREREG = METRICS / "r14_preregistration.json"

DESIGNS = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "r12_case": METRICS / "r12_cases" / "atallah",
          "r11_raw": METRICS / "r11_raw" / "convergence"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "r12_case": METRICS / "r12_cases" / "atallah_designB",
          "r11_raw": METRICS / "r11_raw" / "designB_convergence"},
}
LEVEL = "tighter"
BIN_KM = 10.0
FLOOR = 2
WORK_TOLERANCE = 0.01          # pre-registered |W_A/(beta W_crit) - 1| target
BISECTION_STEPS = 200

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


def accel_inertial(r_vec, t, n, args):
    """Inertial acceleration of the uniformly rotating body-fixed N-field,
    identical to the propagator's right-hand-side convention."""
    x, y, z = r_vec
    th = OMEGA_MOON * t
    c, s = math.cos(th), math.sin(th)
    xb = c * x + s * y
    yb = -s * x + c * y
    axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *args)
    return np.array([c * axb - s * ayb, s * axb + c * ayb, azb])


# ---------------------------------------------------------------- degree tables
def binned_table(model, g, tol, hp_km, ha_km, cap):
    _, table = at.atallah_binned_schedule(model, g, tol, hp_km, ha_km,
                                          floor=FLOOR, cap=cap, bin_km=BIN_KM)
    return {float(k): int(v) for k, v in table.items()}


def degrees_from_table(table: dict, h_km: np.ndarray) -> np.ndarray:
    hmin, hmax = min(table), max(table)
    keys = np.array(sorted(table))
    vals = np.array([table[k] for k in keys], dtype=int)
    hb = np.clip(BIN_KM * np.floor(h_km / BIN_KM), hmin, hmax)
    idx = np.searchsorted(keys, hb - 1e-9)
    idx = np.clip(idx, 0, len(keys) - 1)
    return vals[idx]


def calibrate_tolerance(model, g, hp_km, ha_km, cap, h_km, target_work):
    """Log-space bisection on eps_A so that <N_A^2> over the sampled epochs
    matches ``target_work``. Degree work is monotone non-increasing in eps."""
    def work_of(tol):
        table = binned_table(model, g, tol, hp_km, ha_km, cap)
        deg = degrees_from_table(table, h_km)
        return float(np.mean(deg.astype(float) ** 2)), table, deg

    lo_log, hi_log = -18.0, 2.0          # eps_A [m/s^2]
    w_hi, _, _ = work_of(10.0 ** hi_log)  # loosest tolerance -> lowest degree
    w_lo, _, _ = work_of(10.0 ** lo_log)  # tightest tolerance -> highest degree
    if target_work > w_lo:
        w, table, deg = work_of(10.0 ** lo_log)
        return {"tol": 10.0 ** lo_log, "table": table, "degrees": deg, "work": w,
                "attainable": False, "limit": "cap"}
    if target_work < w_hi:
        w, table, deg = work_of(10.0 ** hi_log)
        return {"tol": 10.0 ** hi_log, "table": table, "degrees": deg, "work": w,
                "attainable": False, "limit": "floor"}
    best = None
    for _ in range(BISECTION_STEPS):
        mid_log = 0.5 * (lo_log + hi_log)
        w, table, deg = work_of(10.0 ** mid_log)
        err = abs(w / target_work - 1.0)
        if best is None or err < best["mismatch"]:
            best = {"tol": 10.0 ** mid_log, "table": table, "degrees": deg,
                    "work": w, "mismatch": err}
        if w > target_work:            # too much work -> loosen tolerance
            lo_log = mid_log
        else:
            hi_log = mid_log
        if hi_log - lo_log < 1e-12:
            break
    best["attainable"] = bool(best["mismatch"] < WORK_TOLERANCE)
    best["limit"] = None if best["attainable"] else "integer_degree_discretization"
    return best


def fixed_degree_for(beta: float, n_crit: int, cap: int) -> tuple[int, bool]:
    """N_F(beta) = argmin_N |N^2 - beta N_crit^2|, no interpolation."""
    target = beta * n_crit ** 2
    n_real = math.sqrt(target)
    cands = [max(1, int(math.floor(n_real))), max(1, int(math.ceil(n_real)))]
    n = min(cands, key=lambda k: abs(k ** 2 - target))
    return (min(n, cap), n > cap)


# ------------------------------------------------------------------- diagnostics
def ric_axes(r, v):
    radial = r / np.linalg.norm(r)
    cross = np.cross(r, v)
    cross /= np.linalg.norm(cross)
    in_track = np.cross(cross, radial)
    return radial, in_track, cross


def allocation_stats(deg: np.ndarray, h_km: np.ndarray, t: np.ndarray,
                     cap: int) -> dict:
    d = deg.astype(float)
    vals, counts = np.unique(deg, return_counts=True)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    hb = BIN_KM * np.floor(h_km / BIN_KM)
    bvals, bcounts = np.unique(hb, return_counts=True)
    return {
        "min_degree": int(deg.min()), "max_degree": int(deg.max()),
        "mean_degree": float(d.mean()), "rms_degree": float(np.sqrt(np.mean(d ** 2))),
        "mean_degree_sq": float(np.mean(d ** 2)),
        "fraction_at_cap": float(np.mean(deg >= cap)),
        "fraction_at_floor": float(np.mean(deg <= FLOOR)),
        "degree_histogram": {str(int(k)): int(c) for k, c in zip(vals, counts)},
        "time_in_altitude_bin_s": {f"{b:.0f}": float(c * dt)
                                   for b, c in zip(bvals, bcounts)},
        "switch_count_binned": int(np.count_nonzero(np.diff(deg))),
    }


def defect_stats(d: np.ndarray, proj: dict, t: np.ndarray) -> dict:
    mag = np.linalg.norm(d, axis=1)
    T = float(t[-1])
    weight = T - t
    return {
        "defect_rms_m_s2": float(np.sqrt(np.mean(mag ** 2))),
        "defect_max_m_s2": float(mag.max()),
        "defect_p95_m_s2": float(np.percentile(mag, 95)),
        "radial_rms_m_s2": float(np.sqrt(np.mean(proj["R"] ** 2))),
        "in_track_rms_m_s2": float(np.sqrt(np.mean(proj["I"] ** 2))),
        "cross_track_rms_m_s2": float(np.sqrt(np.mean(proj["C"] ** 2))),
        "in_track_mean_m_s2": float(np.mean(proj["I"])),
        "in_track_coherence": float(abs(np.mean(proj["I"]))
                                    / np.sqrt(np.mean(proj["I"] ** 2))
                                    if np.any(proj["I"]) else 0.0),
        "in_track_impulse_m_s": float(np.trapezoid(proj["I"], t)),
        "in_track_displacement_proxy_m": float(np.trapezoid(weight * proj["I"], t)),
    }


def worker(task: dict) -> dict:
    design, row, betas = task["design"], task["row"], task["betas"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        model, args = _model(adopted)
        g = _g(adopted)
        t, y = base.load_raw(DESIGNS[design]["r11_raw"] / f"sobolA_{index:03d}"
                             / f"truth_{LEVEL}.npz")
        r_all = y[0:3, :].T
        v_all = y[3:6, :].T
        h_km = (np.linalg.norm(r_all, axis=1) - model.r_ref) / 1e3
        n_ep = len(t)
        w_crit = float(n_crit ** 2)

        # ---- policies: budget grid + the archived accuracy-target point
        policies = {}
        for beta in betas:
            cal = calibrate_tolerance(model, g, hp_km, ha_km, adopted, h_km,
                                      beta * w_crit)
            n_f, f_censored = fixed_degree_for(beta, n_crit, adopted)
            policies[f"beta_{beta:.2f}"] = {
                "beta_requested": beta,
                "atallah": {"kind": "budget_calibrated_atallah",
                            "tol_accel_m_s2": cal["tol"],
                            "degrees": cal["degrees"],
                            "table": cal["table"],
                            "achieved_work": cal["work"],
                            "work_mismatch": cal["work"] / (beta * w_crit) - 1.0,
                            "beta_achieved": cal["work"] / w_crit,
                            "attainable": cal["attainable"],
                            "limit": cal["limit"]},
                "fixed": {"kind": "fixed_degree", "degree": n_f,
                          "degrees": np.full(n_ep, n_f, dtype=int),
                          "achieved_work": float(n_f ** 2),
                          "work_mismatch": (n_f ** 2) / (beta * w_crit) - 1.0,
                          "beta_achieved": (n_f ** 2) / w_crit,
                          "attainable": not f_censored,
                          "limit": "cap" if f_censored else None},
                "censored": bool(f_censored or not cal["attainable"]),
            }

        # archived accuracy-target operating point (measured beta, never assumed)
        cfg_dir = DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
        at_cfg = json.loads((cfg_dir / "atallah_tight.json").read_text())["config"]
        fw_cfg = json.loads((cfg_dir / "fixed_work_atallah_tight.json"
                             ).read_text())["config"]
        orig_table = {float(k): int(v) for k, v
                      in at_cfg["atallah_degree_table"].items()}
        orig_deg = degrees_from_table(orig_table, h_km)
        orig_work = float(np.mean(orig_deg.astype(float) ** 2))
        n_work_orig = int(fw_cfg["policy_spec"]["degree"])
        policies["original"] = {
            "beta_requested": None,
            "atallah": {"kind": "accuracy_target_atallah",
                        "tol_accel_m_s2": float(at_cfg["atallah_tol_accel_m_s2"]),
                        "degrees": orig_deg, "table": orig_table,
                        "achieved_work": orig_work, "work_mismatch": None,
                        "beta_achieved": orig_work / w_crit,
                        "attainable": True, "limit": None},
            "fixed": {"kind": "fixed_work_archived", "degree": n_work_orig,
                      "degrees": np.full(n_ep, n_work_orig, dtype=int),
                      "achieved_work": float(n_work_orig ** 2),
                      "work_mismatch": None,
                      "beta_achieved": (n_work_orig ** 2) / w_crit,
                      "attainable": True, "limit": None},
            "censored": False,
        }

        # ---- evaluate defects, caching one kernel call per (epoch, degree)
        deg_matrix = {name: {k: p[k]["degrees"] for k in ("atallah", "fixed")}
                      for name, p in policies.items()}
        acc = {name: {k: np.empty((n_ep, 3)) for k in ("atallah", "fixed")}
               for name in policies}
        proj = {name: {k: {c: np.empty(n_ep) for c in "RIC"}
                       for k in ("atallah", "fixed")} for name in policies}

        for j in range(n_ep):
            rj, vj, tj = r_all[j], v_all[j], float(t[j])
            a_truth = accel_inertial(rj, tj, adopted, args)
            cache: dict[int, np.ndarray] = {}
            axes = ric_axes(rj, vj)
            for name in policies:
                for k in ("atallah", "fixed"):
                    n = int(deg_matrix[name][k][j])
                    if n not in cache:
                        cache[n] = accel_inertial(rj, tj, n, args) - a_truth
                    d = cache[n]
                    acc[name][k][j] = d
                    for c, ax in zip("RIC", axes):
                        proj[name][k][c][j] = ax @ d

        out_policies = {}
        for name, p in policies.items():
            entry = {"beta_requested": p["beta_requested"],
                     "censored": p["censored"]}
            for k in ("atallah", "fixed"):
                deg = p[k]["degrees"]
                entry[k] = {
                    "kind": p[k]["kind"],
                    "tol_accel_m_s2": p[k].get("tol_accel_m_s2"),
                    "degree": p[k].get("degree"),
                    "achieved_work": p[k]["achieved_work"],
                    "beta_achieved": p[k]["beta_achieved"],
                    "work_mismatch": p[k]["work_mismatch"],
                    "attainable": p[k]["attainable"],
                    "limit": p[k]["limit"],
                    "allocation": allocation_stats(deg, h_km, t, adopted),
                    "defect": defect_stats(acc[name][k], proj[name][k], t),
                }
            a_d, f_d = entry["atallah"]["defect"], entry["fixed"]["defect"]
            entry["ratios"] = {
                "R_a_defect_rms": (f_d["defect_rms_m_s2"] / a_d["defect_rms_m_s2"]
                                   if a_d["defect_rms_m_s2"] > 0 else None),
                "R_a_defect_max": (f_d["defect_max_m_s2"] / a_d["defect_max_m_s2"]
                                   if a_d["defect_max_m_s2"] > 0 else None),
                "R_in_track_rms": (f_d["in_track_rms_m_s2"] / a_d["in_track_rms_m_s2"]
                                   if a_d["in_track_rms_m_s2"] > 0 else None),
                "displacement_proxy_ratio": (
                    abs(f_d["in_track_displacement_proxy_m"])
                    / abs(a_d["in_track_displacement_proxy_m"])
                    if a_d["in_track_displacement_proxy_m"] != 0.0 else None),
                "work_ratio_atallah_over_fixed": (entry["atallah"]["achieved_work"]
                                                  / entry["fixed"]["achieved_work"]),
                "atallah_smaller_defect": bool(
                    a_d["defect_rms_m_s2"] < f_d["defect_rms_m_s2"]),
                "atallah_smaller_displacement_proxy": bool(
                    abs(a_d["in_track_displacement_proxy_m"])
                    < abs(f_d["in_track_displacement_proxy_m"])),
            }
            out_policies[name] = entry

        return {"design": design, "sobol_index": index, "status": "complete",
                "adopted_truth_degree": adopted, "n_critical": n_crit,
                "w_crit": w_crit, "hp_km": hp_km, "ha_km": ha_km,
                "incl_deg": float(row["design_point"]["incl_deg"]),
                "eccentricity": float(row["design_point"]["eccentricity"]),
                "n_epochs": n_ep, "arc_s": float(t[-1]),
                "budgets": out_policies}
    except Exception as exc:
        return {"design": design, "sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ---------------------------------------------------------------------- summary
def stat(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p50": float(np.median(a)),
            "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def summarize_budget(rows: list[dict], key: str) -> dict:
    live = [r for r in rows if key in r["budgets"]
            and not r["budgets"][key]["censored"]]
    censored = [r for r in rows if key in r["budgets"]
                and r["budgets"][key]["censored"]]
    if not live:
        return {"orbits": 0, "censored": len(censored)}
    ent = [r["budgets"][key] for r in live]
    ra = [e["ratios"]["R_a_defect_rms"] for e in ent]
    return {
        "orbits": len(live), "censored": len(censored),
        "beta_requested": ent[0]["beta_requested"],
        "beta_achieved_atallah": stat([e["atallah"]["beta_achieved"] for e in ent]),
        "beta_achieved_fixed": stat([e["fixed"]["beta_achieved"] for e in ent]),
        "work_match_mismatch_abs": stat(
            [abs(e["atallah"]["work_mismatch"]) for e in ent
             if e["atallah"]["work_mismatch"] is not None]),
        "R_a_defect_rms": stat(ra),
        "R_in_track_rms": stat([e["ratios"]["R_in_track_rms"] for e in ent]),
        "displacement_proxy_ratio": stat(
            [e["ratios"]["displacement_proxy_ratio"] for e in ent]),
        "atallah_defect_rms_m_s2": stat([e["atallah"]["defect"]["defect_rms_m_s2"]
                                         for e in ent]),
        "fixed_defect_rms_m_s2": stat([e["fixed"]["defect"]["defect_rms_m_s2"]
                                       for e in ent]),
        "atallah_wins_defect": int(sum(e["ratios"]["atallah_smaller_defect"]
                                       for e in ent)),
        "fixed_wins_defect": int(sum(not e["ratios"]["atallah_smaller_defect"]
                                     for e in ent)),
        "atallah_wins_displacement_proxy": int(sum(
            e["ratios"]["atallah_smaller_displacement_proxy"] for e in ent)),
        "atallah_mean_degree": stat([e["atallah"]["allocation"]["mean_degree"]
                                     for e in ent]),
        "atallah_fraction_at_cap": stat([e["atallah"]["allocation"]["fraction_at_cap"]
                                         for e in ent]),
        "fixed_degree": stat([float(e["fixed"]["degree"]) for e in ent]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B", "both"), default="both")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    betas = list(prereg["budget_grid"])
    designs = ["A", "B"] if a.design == "both" else [a.design]

    payload = {"schema": "r14_budget_pareto_v1", "created_utc": base.utc_now(),
               "reference_level": LEVEL, "budget_grid": betas,
               "preregistration_sha256": prereg["protocol_sha256"],
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
        tasks = [{"design": d, "row": r, "betas": betas} for r in rows]
        print(f"[pareto] design {d}: {len(tasks)} orbits x "
              f"{len(betas) + 1} budget points", flush=True)
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
        keys = [f"beta_{b:.2f}" for b in betas] + ["original"]
        payload["designs"][d] = {
            "rows": done, "failures": fails,
            "summary": {k: summarize_budget(done, k) for k in keys},
            "beta_original": stat([r["budgets"]["original"]["atallah"]["beta_achieved"]
                                   for r in done]),
        }
        s = payload["designs"][d]["summary"]
        for k in keys:
            e = s[k]
            if not e.get("orbits"):
                print(f"  {k}: all {e['censored']} orbits censored", flush=True)
                continue
            print(f"  {k}: n={e['orbits']} censored={e['censored']} "
                  f"R_a median={e['R_a_defect_rms']['median']:.3g} "
                  f"At wins {e['atallah_wins_defect']}/{e['orbits']} "
                  f"displ wins {e['atallah_wins_displacement_proxy']}/{e['orbits']}",
                  flush=True)
        bo = payload["designs"][d]["beta_original"]
        print(f"  original operating point: beta median={bo['median']:.3f} "
              f"[{bo['min']:.3f}, {bo['max']:.3f}]", flush=True)

    base.atomic_json(OUTPUT, payload)
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
