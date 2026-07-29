"""Forced-variational calibration of the force-defect ranking (R13).

rev13_force_defect.py measures the truncation acceleration defect of each policy
along the archived truth trajectory without any integration noise, but converts
it into meters only through a first-order free-propagation proxy. This script
does the conversion properly, with the manuscript's own forced variational
system (Eq. 11 of the main text):

    d/dt dr = dv,
    d/dt dv = G(t) dr + Delta_a_P(t),      Delta_a_P = a_{N_P} - a_{N_truth},

integrated along the reference trajectory for three policies at once (the
Atallah degree history, its work-matched fixed degree, and the critical fixed
degree). All three share one reference and one gravity gradient, so a single
augmented integration gives all three linear predictions.

The critical fixed degree is the calibration channel: its measured seven-day
error is a large, tolerance-stable signal (a median relative change of 1.8%
under tolerance refinement), so agreement between its linear prediction and its
measured error certifies the predictions for the two near-truth policies, whose
measured errors are noise-limited.

The gravity gradient is evaluated at degree 120 by central differences; the
gradient is dominated by the central term and the lowest degrees, and the linear
response is insensitive to that choice at the reported precision.

Usage:
    python rev13_variational_check.py --orbits 8 --workers 4
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
from scipy.integrate import solve_ivp

import rev10_sobol_confirmatory as base
from rev3_common import OMEGA_MOON
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r13_variational_check.json"
TABLE = METRICS / "r13_variational_check_table.tex"

DESIGNS = {
    "A": {"r12_case": METRICS / "r12_cases" / "atallah",
          "campaign": METRICS / "r12_atallah_campaign.json"},
    "B": {"r12_case": METRICS / "r12_cases" / "atallah_designB",
          "campaign": METRICS / "r12_atallah_campaign_designB.json"},
}
GRADIENT_DEGREE = 120
FD_STEP_M = 1.0
RTOL = 1.0e-11
ATOL = 1.0e-9
MAX_STEP_S = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0

_MODELS: dict[int, tuple] = {}


def _model(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def accel_inertial(r_vec, t, n, args):
    x, y, z = r_vec
    th = OMEGA_MOON * t
    c, s = math.cos(th), math.sin(th)
    xb = c * x + s * y
    yb = -s * x + c * y
    axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *args)
    return np.array([c * axb - s * ayb, s * axb + c * ayb, azb])


def gradient(r_vec, t, n, args, h=FD_STEP_M):
    g = np.empty((3, 3))
    for j in range(3):
        rp = r_vec.copy(); rp[j] += h
        rm = r_vec.copy(); rm[j] -= h
        g[:, j] = (accel_inertial(rp, t, n, args)
                   - accel_inertial(rm, t, n, args)) / (2.0 * h)
    return g


def binned_lookup(table: dict, bin_km: float = 10.0):
    tab = {float(k): int(v) for k, v in table.items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(h_m):
        hb = min(hmax, max(hmin, bin_km * math.floor(h_m / 1e3 / bin_km)))
        return tab[hb]

    return degree_of


def ric_axes(r, v):
    radial = r / np.linalg.norm(r)
    cross = np.cross(r, v)
    cross /= np.linalg.norm(cross)
    return radial, np.cross(cross, radial), cross


def worker(task: dict) -> dict:
    design, index = task["design"], task["index"]
    try:
        cfg_dir = DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
        at_cfg = json.loads((cfg_dir / "atallah_tight.json").read_text())["config"]
        fw_cfg = json.loads((cfg_dir / "fixed_work_atallah_tight.json").read_text())["config"]
        adopted = int(at_cfg["adopted_truth_degree"])
        n_work = int(fw_cfg["policy_spec"]["degree"])
        n_crit = int(at_cfg["n_critical"])
        y0 = np.asarray(at_cfg["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        gmodel, gargs = _model(min(GRADIENT_DEGREE, adopted))
        deg_at = binned_lookup(at_cfg["atallah_degree_table"])
        r_ref = model.r_ref
        policies = ("atallah", "fixed_work", "fixed_critical")

        def rhs(t, Y):
            r = Y[0:3]
            v = Y[3:6]
            a_ref = accel_inertial(r, t, adopted, args)
            G = gradient(r, t, min(GRADIENT_DEGREE, adopted), gargs)
            alt = float(np.linalg.norm(r)) - r_ref
            degs = (deg_at(alt), n_work, n_crit)
            dY = np.empty_like(Y)
            dY[0:3] = v
            dY[3:6] = a_ref
            for k, n in enumerate(degs):
                off = 6 + 6 * k
                da = accel_inertial(r, t, n, args) - a_ref
                dY[off:off + 3] = Y[off + 3:off + 6]
                dY[off + 3:off + 6] = G @ Y[off:off + 3] + da
            return dY

        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        Y0 = np.zeros(6 + 6 * len(policies))
        Y0[0:6] = y0
        t0 = time.time()
        sol = solve_ivp(rhs, (0.0, DURATION), Y0, method="DOP853", t_eval=grid,
                        rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S)
        if not sol.success:
            return {"design": design, "index": index, "status": "integration_failed",
                    "message": sol.message}
        ref = sol.y[0:6]
        out = {"design": design, "sobol_index": index, "status": "complete",
               "adopted_truth_degree": adopted, "n_work": n_work,
               "n_critical": n_crit, "wall_s": time.time() - t0,
               "n_rhs": int(sol.nfev), "policies": {}}
        for k, name in enumerate(policies):
            off = 6 + 6 * k
            d = sol.y[off:off + 3]
            ric = np.empty((d.shape[1], 3))
            for j in range(d.shape[1]):
                axes = ric_axes(ref[0:3, j], ref[3:6, j])
                ric[j] = [ax @ d[:, j] for ax in axes]
            mag = np.linalg.norm(d, axis=0)
            out["policies"][name] = {
                "predicted_pos_rms_m": float(np.sqrt(np.mean(mag ** 2))),
                "predicted_pos_final_m": float(mag[-1]),
                "predicted_in_track_rms_m": float(np.sqrt(np.mean(ric[:, 1] ** 2))),
                "predicted_in_track_final_m": float(ric[-1, 1]),
                "predicted_radial_rms_m": float(np.sqrt(np.mean(ric[:, 0] ** 2))),
            }
        return out
    except Exception as exc:
        return {"design": design, "index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def measured(design: str, index: int) -> dict:
    camp = json.loads(DESIGNS[design]["campaign"].read_text())
    for r in camp["rows"]:
        if int(r["sobol_index"]) == index:
            return {"atallah": r["policies"]["atallah"]["error_tight"]["pos_rms_m"],
                    "fixed_work": r["policies"]["fixed_work_atallah"]["error_tight"]["pos_rms_m"],
                    "fixed_critical": r["policies"]["fixed_critical"]["error_tight"]["pos_rms_m"],
                    "envelope_atallah": r["policies"]["atallah"]["truth_inclusive_envelope_m"],
                    "envelope_fixed_work": r["policies"]["fixed_work_atallah"]["truth_inclusive_envelope_m"]}
    return {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orbits", type=int, default=8)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--design", choices=("A", "B"), default="A")
    a = ap.parse_args()
    camp = json.loads(DESIGNS[a.design]["campaign"].read_text())
    rows = sorted(camp["rows"], key=lambda r: r["design_point"]["hp_km"])
    picks = [int(rows[int(i)]["sobol_index"])
             for i in np.linspace(0, len(rows) - 1, a.orbits).round()]
    tasks = [{"design": a.design, "index": i} for i in sorted(set(picks))]
    print(f"[variational] {len(tasks)} orbits, gradient degree {GRADIENT_DEGREE}",
          flush=True)
    results, fails = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['index']:03d} {rec.get('message')}", flush=True)
                continue
            rec["measured_tight_m"] = measured(rec["design"], rec["sobol_index"])
            m = rec["measured_tight_m"]
            p = rec["policies"]
            rec["calibration"] = {
                "critical_predicted_over_measured": (
                    p["fixed_critical"]["predicted_pos_rms_m"] / m["fixed_critical"]),
                "predicted_ratio_fixed_over_atallah": (
                    p["fixed_work"]["predicted_pos_rms_m"]
                    / p["atallah"]["predicted_pos_rms_m"]),
                "predicted_gap_m": abs(p["fixed_work"]["predicted_pos_rms_m"]
                                       - p["atallah"]["predicted_pos_rms_m"]),
                "measured_threshold_m": (m["envelope_atallah"]
                                         + m["envelope_fixed_work"]),
            }
            rec["calibration"]["predicted_gap_over_threshold"] = (
                rec["calibration"]["predicted_gap_m"]
                / rec["calibration"]["measured_threshold_m"])
            results.append(rec)
            c = rec["calibration"]
            print(f"  [{n:2d}/{len(tasks)}] idx={rec['sobol_index']:03d} "
                  f"crit pred/meas {c['critical_predicted_over_measured']:.2f}  "
                  f"pred ratio fix/At {c['predicted_ratio_fixed_over_atallah']:.1f}  "
                  f"pred gap {c['predicted_gap_m']:.3f} m vs thr "
                  f"{c['measured_threshold_m']:.3f} m  "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    results.sort(key=lambda r: r["sobol_index"])

    def stat(v):
        x = np.asarray([q for q in v if q is not None and np.isfinite(q)], float)
        return None if x.size == 0 else {
            "n": int(x.size), "median": float(np.median(x)),
            "min": float(x.min()), "max": float(x.max())}

    payload = {"schema": "r13_variational_check_v1", "created_utc": base.utc_now(),
               "gradient_degree": GRADIENT_DEGREE, "rtol": RTOL, "atol": ATOL,
               "rows": results, "failures": fails,
               "summary": {
                   "orbits": len(results),
                   "critical_predicted_over_measured": stat(
                       [r["calibration"]["critical_predicted_over_measured"]
                        for r in results]),
                   "predicted_ratio_fixed_over_atallah": stat(
                       [r["calibration"]["predicted_ratio_fixed_over_atallah"]
                        for r in results]),
                   "predicted_gap_over_threshold": stat(
                       [r["calibration"]["predicted_gap_over_threshold"]
                        for r in results]),
                   "atallah_predicted_smaller": int(sum(
                       r["policies"]["atallah"]["predicted_pos_rms_m"]
                       < r["policies"]["fixed_work"]["predicted_pos_rms_m"]
                       for r in results))}}
    base.atomic_json(OUTPUT, payload)
    body = "\n".join(
        f"    {r['sobol_index']} & {r['n_critical']} & {r['n_work']} & "
        f"{r['policies']['fixed_critical']['predicted_pos_rms_m']:.1f} & "
        f"{r['measured_tight_m']['fixed_critical']:.1f} & "
        f"{r['calibration']['critical_predicted_over_measured']:.2f} & "
        f"{r['policies']['atallah']['predicted_pos_rms_m']:.4f} & "
        f"{r['policies']['fixed_work']['predicted_pos_rms_m']:.4f} & "
        f"{r['calibration']['predicted_ratio_fixed_over_atallah']:.1f} & "
        f"{r['calibration']['measured_threshold_m']:.2f}\\\\" for r in results)
    TABLE.write_text(f"""% auto-generated by rev13_variational_check.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Forced-variational calibration of the matched-work comparison. For
  each orbit one augmented integration along the reference trajectory gives the
  linear position response to the measured truncation force defect of three
  policies. The critical fixed degree is the calibration channel: its predicted
  and measured seven-day errors are compared directly, since that error is a
  large, tolerance-stable signal. The Atallah and work-matched predictions are
  then read on the same footing, and their predicted gap is compared with the
  measured resolution threshold of the trajectory experiment.}}
  \\label{{tab:variational-check}}
  \\begin{{tabular}}{{r r r r r r r r r r}}
    \\toprule
    & & & \\multicolumn{{3}}{{c}}{{critical fixed degree [m]}} &
      \\multicolumn{{3}}{{c}}{{matched work, predicted [m]}} & thr.\\\\
    \\cmidrule(lr){{4-6}}\\cmidrule(lr){{7-9}}
    idx & $N_{{\\mathrm{{crit}}}}$ & $N_{{\\mathrm{{work}}}}$ & pred & meas & ratio &
      Atallah & fixed & ratio & [m]\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
""", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))
    print(f"[written] {OUTPUT.name}, {TABLE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
