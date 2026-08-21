"""R35: the RMS-form control on the second rung of the instrument ladder.

The instrument ladder of Section 7.2 scores its second rung with a terminal,
in-track, free-propagation proxy

    D_I = |int_0^T (T - tau) Delta a_I(tau) d tau|,

whereas the trajectory metric the ladder is being compared against is the RMS
over the arc of the three-dimensional position error.  The two therefore differ
in form as well as in content, which leaves open the reading that the second
rung fails only because it is a terminal scalar and not because signed
remaining-horizon weighting is insufficient.

This script closes that gap without any new propagation.  Along the same
archived reference trajectory, and for the same two policies at beta = 1, it
integrates the free-propagation (gradient-free) displacement

    d(t) = int_0^t (t - tau) Delta a(tau) d tau

and reports it in the metric's own form:

  * D_I_rms   = sqrt( (1/T) int_0^T d_I(t)^2 dt ), the in-track component
                projected exactly as the terminal proxy projects it;
  * D_3d_rms  = sqrt( (1/T) int_0^T |d(t)|^2 dt ), the inertial three-vector,
                which is the free-propagation analog of the propagated
                metric itself.

Both are still free-propagation statistics: they carry the sign and the
remaining-horizon weight but no gravity-gradient coupling.  The terminal proxy
D_I is recomputed alongside them and checked against the archived R14 value, so
a drift in the calibration or the reference records shows up as a failure here
rather than as a silently different number.

Every quantity is a deterministic function of the field, the calibrated degree
histories and the archived reference states; reruns reproduce it bit for bit.

Usage:
    python rev35_proxy_rms.py --design both --workers 5
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_pareto as bp

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r35_proxy_rms.json"
TABLE = METRICS / "r35_proxy_rms_table.tex"
LADDER = METRICS / "r34_instrument_ladder.json"

BETA = 1.0
BETA_KEY = "beta_1.00"
LEVEL = bp.LEVEL
PROXY_RTOL = 1e-6          # agreement required against the archived R14 proxy


def cum_trapz(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral of y over t, starting at zero."""
    out = np.zeros_like(y, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))
    return out


def free_displacement(a: np.ndarray, t: np.ndarray) -> np.ndarray:
    """d(t) = int_0^t (t - tau) a(tau) d tau, by double cumulative integration.

    Equivalent to integrating the impulse J(t) = int_0^t a d tau once more,
    which is what a body under forcing a and no restoring gradient does.
    """
    if a.ndim == 1:
        return cum_trapz(cum_trapz(a, t), t)
    return np.column_stack([cum_trapz(cum_trapz(a[:, k], t), t)
                            for k in range(a.shape[1])])


def time_rms(x: np.ndarray, t: np.ndarray) -> float:
    """sqrt( (1/T) int_0^T x(t)^2 dt ), the arc-RMS in the metric's own form."""
    T = float(t[-1] - t[0])
    return float(np.sqrt(np.trapezoid(x ** 2, t) / T))


def worker(task: dict) -> dict:
    design, row = task["design"], task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        model, args = bp._model(adopted)
        g = bp._g(adopted)
        t, y = base.load_raw(bp.DESIGNS[design]["r11_raw"] / f"sobolA_{index:03d}"
                             / f"truth_{LEVEL}.npz")
        r_all = y[0:3, :].T
        v_all = y[3:6, :].T
        h_km = (np.linalg.norm(r_all, axis=1) - model.r_ref) / 1e3
        n_ep = len(t)
        w_crit = float(n_crit ** 2)

        # Same two policies as the ladder's beta = 1 column, recalibrated by the
        # same deterministic bisection rather than copied.
        cal = bp.calibrate_tolerance(model, g, hp_km, ha_km, adopted, h_km,
                                     BETA * w_crit)
        n_f, f_censored = bp.fixed_degree_for(BETA, n_crit, adopted)
        censored = bool(f_censored or not cal["attainable"])
        degrees = {"atallah": cal["degrees"],
                   "fixed": np.full(n_ep, n_f, dtype=int)}

        defect = {k: np.empty((n_ep, 3)) for k in degrees}
        in_track = {k: np.empty(n_ep) for k in degrees}
        for j in range(n_ep):
            rj, vj, tj = r_all[j], v_all[j], float(t[j])
            a_truth = bp.accel_inertial(rj, tj, adopted, args)
            axes = bp.ric_axes(rj, vj)
            cache: dict[int, np.ndarray] = {}
            for k, deg in degrees.items():
                n = int(deg[j])
                if n not in cache:
                    cache[n] = bp.accel_inertial(rj, tj, n, args) - a_truth
                defect[k][j] = cache[n]
                in_track[k][j] = axes[1] @ cache[n]

        T = float(t[-1])
        weight = T - t
        out = {}
        for k in degrees:
            d_it = free_displacement(in_track[k], t)          # scalar, in-track
            d_3d = free_displacement(defect[k], t)            # inertial vector
            out[k] = {
                "degree": (int(n_f) if k == "fixed" else None),
                "mean_degree_sq": float(np.mean(degrees[k].astype(float) ** 2)),
                # terminal proxy, recomputed for the archived-value check
                "in_track_displacement_proxy_m":
                    float(np.trapezoid(weight * in_track[k], t)),
                # RMS-form controls
                "in_track_displacement_rms_m": time_rms(d_it, t),
                "displacement_rms_3d_m": time_rms(
                    np.linalg.norm(d_3d, axis=1), t),
                "in_track_displacement_final_m": float(d_it[-1]),
            }

        a_o, f_o = out["atallah"], out["fixed"]
        ratios = {}
        for key, name in (("in_track_displacement_proxy_m", "proxy_terminal"),
                          ("in_track_displacement_rms_m", "proxy_in_track_rms"),
                          ("displacement_rms_3d_m", "proxy_3d_rms")):
            va, vf = abs(a_o[key]), abs(f_o[key])
            ratios[name] = {
                "ratio_fixed_over_radial": (vf / va if va > 0.0 else None),
                "radial_smaller": bool(va < vf),
            }

        return {"design": design, "sobol_index": index, "status": "complete",
                "adopted_truth_degree": adopted, "n_critical": n_crit,
                "hp_km": hp_km, "ha_km": ha_km, "n_epochs": n_ep, "arc_s": T,
                "beta_achieved_radial": cal["work"] / w_crit,
                "beta_achieved_fixed": (n_f ** 2) / w_crit,
                "censored": censored,
                "policies": out, "ratios": ratios}
    except Exception as exc:
        return {"design": design, "sobol_index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def check_against_archive(rows: list[dict], design: str) -> dict:
    """The terminal proxy recomputed here must reproduce the archived R14 value."""
    pareto = json.loads((METRICS / "r14_budget_pareto.json").read_text())
    archived = {int(r["sobol_index"]): r for r in pareto["designs"][design]["rows"]}
    worst = 0.0
    checked = 0
    failures = []
    for row in rows:
        ref = archived.get(row["sobol_index"])
        if ref is None:
            continue
        entry = ref.get("budgets", {}).get(BETA_KEY)
        if entry is None:
            continue
        for k in ("atallah", "fixed"):
            old = entry.get(k, {}).get("defect", {}).get(
                "in_track_displacement_proxy_m")
            if old is None:
                continue
            new = row["policies"][k]["in_track_displacement_proxy_m"]
            rel = abs(new - old) / abs(old) if old != 0.0 else abs(new)
            checked += 1
            worst = max(worst, rel)
            if rel > PROXY_RTOL:
                failures.append({"sobol_index": row["sobol_index"],
                                 "policy": k, "archived": old, "recomputed": new,
                                 "rel": rel})
    return {"checked": checked, "worst_rel": worst,
            "tolerance": PROXY_RTOL, "passed": not failures,
            "failures": failures}


def stat(vals):
    x = np.asarray([v for v in vals if v is not None and np.isfinite(v)],
                   dtype=float)
    if x.size == 0:
        return None
    return {"n": int(x.size), "median": float(np.median(x)),
            "p10": float(np.percentile(x, 10)), "p90": float(np.percentile(x, 90)),
            "min": float(x.min()), "max": float(x.max())}


def summarize(rows: list[dict]) -> dict:
    used = [r for r in rows if not r["censored"]]
    out = {"orbits": len(used), "censored": len(rows) - len(used)}
    for name in ("proxy_terminal", "proxy_in_track_rms", "proxy_3d_rms"):
        out[name] = {
            "radial_smaller": int(sum(r["ratios"][name]["radial_smaller"]
                                      for r in used)),
            "ratio_fixed_over_radial": stat(
                [r["ratios"][name]["ratio_fixed_over_radial"] for r in used]),
        }
    return out


def build_table(payload: dict) -> str:
    """Three forms of the same free-propagation statistic, against the
    propagated arcs read from the ladder record so the two cannot drift."""
    a = payload["designs"]["A"]["summary"]
    b = payload["designs"]["B"]["summary"]
    rows = [
        (r"Terminal, in-track",
         r"$\lvert\int_0^T(T-\tau)\Delta a_I\,\mathrm d\tau\rvert$",
         "proxy_terminal"),
        (r"Arc-RMS, in-track",
         r"$\langle d_I(t)^2\rangle_t^{1/2}$", "proxy_in_track_rms"),
        (r"Arc-RMS, three-vector",
         r"$\langle\lVert\mathbf d(t)\rVert^2\rangle_t^{1/2}$",
         "proxy_3d_rms"),
    ]
    lines = [r"\begin{tabular}{@{}l l cc cc@{}}", r"\toprule",
             r" & & \multicolumn{2}{c}{Orbits favoring radial}"
             r" & \multicolumn{2}{c}{Median ratio} \\",
             r"\cmidrule(lr){3-4}\cmidrule(l){5-6}",
             r"Form & Statistic & A & B & A & B \\", r"\midrule"]
    for label, expr, key in rows:
        lines.append(
            f"{label} & {expr} & "
            f"{a[key]['radial_smaller']}/{a['orbits']} & "
            f"{b[key]['radial_smaller']}/{b['orbits']} & "
            f"{a[key]['ratio_fixed_over_radial']['median']:.2f} & "
            f"{b[key]['ratio_fixed_over_radial']['median']:.2f}" + r" \\")
    if LADDER.exists():
        lad = json.loads(LADDER.read_text())
        pa = lad["designs"]["A"]["propagated"]
        pb = lad["designs"]["B"]["propagated"]
        lines += [r"\midrule",
                  r"Propagated arcs & (reference) & "
                  f"{pa['raw_radial_better']}/{pa['n_orbits']} & "
                  f"{pb['raw_radial_better']}/{pb['n_orbits']} & "
                  f"{pa['median_ratio']:.2f} & {pb['median_ratio']:.2f}"
                  + r" \\"]
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B", "both"), default="both")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0, help="first N orbits only")
    a = ap.parse_args()
    designs = ["A", "B"] if a.design == "both" else [a.design]

    payload = {"schema": "r35_proxy_rms_v1", "created_utc": base.utc_now(),
               "beta": BETA, "reference_level": LEVEL,
               "note": ("Free-propagation displacement proxies in the arc-RMS "
                        "form of the trajectory metric; no propagation is run, "
                        "and the terminal proxy is checked against R14."),
               "designs": {}}
    if OUTPUT.exists():
        payload["designs"] = json.loads(OUTPUT.read_text()).get("designs", {})

    for d in designs:
        rows_in = json.loads(
            bp.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]
        if a.limit:
            rows_in = rows_in[:a.limit]
        tasks = [{"design": d, "row": r} for r in rows_in]
        print(f"[proxy-rms] design {d}: {len(tasks)} orbits", flush=True)
        t0 = time.time()
        rows, fails = [], []
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(worker, t): t for t in tasks}
            for n, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                if rec["status"] != "complete":
                    fails.append(rec)
                    print(f"  !! {rec['sobol_index']:03d} {rec['message']}",
                          flush=True)
                    continue
                rows.append(rec)
                if n % 8 == 0:
                    print(f"  [{n:3d}/{len(tasks)}] "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
        rows.sort(key=lambda r: r["sobol_index"])
        payload["designs"][d] = {"rows": rows, "failures": fails,
                                 "summary": summarize(rows),
                                 "archive_check": check_against_archive(rows, d)}
        s = payload["designs"][d]["summary"]
        chk = payload["designs"][d]["archive_check"]
        print(f"[design {d}] over {s['orbits']} uncensored orbits, radial "
              f"favored by: terminal proxy "
              f"{s['proxy_terminal']['radial_smaller']}, in-track RMS "
              f"{s['proxy_in_track_rms']['radial_smaller']}, 3-D RMS "
              f"{s['proxy_3d_rms']['radial_smaller']}", flush=True)
        print(f"           archive check {'PASS' if chk['passed'] else 'FAIL'} "
              f"({chk['checked']} values, worst rel {chk['worst_rel']:.2e})",
              flush=True)

    base.atomic_json(OUTPUT, payload)
    if set(payload["designs"]) >= {"A", "B"}:
        TABLE.write_text(build_table(payload), encoding="utf-8")
        print(f"[written] {OUTPUT.name}, {TABLE.name}")
    else:
        print(f"[written] {OUTPUT.name} (table needs both designs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
