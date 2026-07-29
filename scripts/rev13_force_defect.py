"""Noise-free force-defect ranking of Atallah against its work-matched fixed
degree, over both 64-orbit populations (R13).

The seven-day trajectory comparison at matched work is limited by the numerical
resolution of the experiment, not by the difference between the two force models
(rev13_resolution_diagnosis.py). The quantity that separates them without any
integration noise is the truncation acceleration defect itself, evaluated along
the common archived truth trajectory:

    Delta_a_P(t) = a_{N_P(t)}(x_truth(t)) - a_{N_truth}(x_truth(t)),

for P in {Atallah (binned degree history), work-matched fixed degree}. This is a
deterministic function of the field, the degree history, and the reference
states: rerunning it reproduces bit-identical numbers.

Reported per orbit:
  * RMS and maximum of |Delta_a| over the arc,
  * the signed in-track projection, whose coherent part drives the t^2 in-track
    growth identified by the manuscript's variational analysis,
  * the in-track velocity impulse J = int Delta_a_it dt and the first-order
    in-track displacement proxy D = int (T - tau) Delta_a_it(tau) d tau, which
    is the free-propagation estimate of the accumulated displacement,
  * the sampled quadratic work of each degree history, as a check that the two
    policies really are matched on the sampled epochs.

The ranking is then read from the ratios; a full forced-variational check on a
subset is run separately by rev13_variational_check.py.

Usage:
    python rev13_force_defect.py --design both --workers 5
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
from rev3_common import OMEGA_MOON
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r13_force_defect.json"
TABLE = METRICS / "r13_force_defect_table.tex"

DESIGNS = {
    "A": {"r12_case": METRICS / "r12_cases" / "atallah",
          "r11_raw": METRICS / "r11_raw" / "convergence"},
    "B": {"r12_case": METRICS / "r12_cases" / "atallah_designB",
          "r11_raw": METRICS / "r11_raw" / "designB_convergence"},
}
LEVEL = "tighter"   # reference trajectory level (the more converged of the two)

_MODELS: dict[int, tuple] = {}


def _model(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def accel_inertial(r_vec, t, n, args):
    """Inertial acceleration of the uniformly rotating body-fixed N-field,
    identical to the convention of the propagator's right-hand side."""
    x, y, z = r_vec
    th = OMEGA_MOON * t
    c, s = math.cos(th), math.sin(th)
    xb = c * x + s * y
    yb = -s * x + c * y
    axb, ayb, azb = sh_accel_fixed_numba(xb, yb, z, n, *args)
    return np.array([c * axb - s * ayb, s * axb + c * ayb, azb])


def binned_lookup(table: dict, bin_km: float = 10.0):
    tab = {float(k): int(v) for k, v in table.items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(h_m):
        hb = min(hmax, max(hmin, bin_km * math.floor(h_m / 1e3 / bin_km)))
        return tab[hb]

    return degree_of


def in_track_axis(r, v):
    radial = r / np.linalg.norm(r)
    cross = np.cross(r, v)
    cross /= np.linalg.norm(cross)
    return np.cross(cross, radial)


def worker(task: dict) -> dict:
    design, index = task["design"], task["index"]
    try:
        cfg_dir = DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
        at_cfg = json.loads((cfg_dir / "atallah_tight.json").read_text())["config"]
        fw_cfg = json.loads((cfg_dir / "fixed_work_atallah_tight.json").read_text())["config"]
        adopted = int(at_cfg["adopted_truth_degree"])
        n_work = int(fw_cfg["policy_spec"]["degree"])
        model, args = _model(adopted)
        t, y = base.load_raw(DESIGNS[design]["r11_raw"] / f"sobolA_{index:03d}"
                             / f"truth_{LEVEL}.npz")
        deg_of = binned_lookup(at_cfg["atallah_degree_table"])
        n_ep = len(t)
        d_at = np.empty((n_ep, 3))
        d_fw = np.empty((n_ep, 3))
        it_at = np.empty(n_ep)
        it_fw = np.empty(n_ep)
        deg_hist = np.empty(n_ep, dtype=int)
        for k in range(n_ep):
            r = y[0:3, k]
            v = y[3:6, k]
            tk = float(t[k])
            n_at = deg_of(float(np.linalg.norm(r)) - model.r_ref)
            deg_hist[k] = n_at
            a_truth = accel_inertial(r, tk, adopted, args)
            d_at[k] = accel_inertial(r, tk, n_at, args) - a_truth
            d_fw[k] = accel_inertial(r, tk, n_work, args) - a_truth
            axis = in_track_axis(r, v)
            it_at[k] = axis @ d_at[k]
            it_fw[k] = axis @ d_fw[k]

        T = float(t[-1])
        weight = T - t                      # first-order propagation kernel

        def summarize(d, it):
            mag = np.linalg.norm(d, axis=1)
            return {
                "defect_rms_m_s2": float(np.sqrt(np.mean(mag ** 2))),
                "defect_max_m_s2": float(mag.max()),
                "in_track_rms_m_s2": float(np.sqrt(np.mean(it ** 2))),
                "in_track_mean_m_s2": float(np.mean(it)),
                "in_track_impulse_m_s": float(np.trapezoid(it, t)),
                "in_track_displacement_proxy_m": float(np.trapezoid(weight * it, t)),
            }

        return {
            "design": design, "sobol_index": index, "status": "complete",
            "adopted_truth_degree": adopted, "n_work": n_work,
            "n_epochs": n_ep, "arc_s": T,
            "atallah": {**summarize(d_at, it_at),
                        "sampled_mean_N2": float(np.mean(deg_hist.astype(float) ** 2)),
                        "mean_degree": float(deg_hist.mean()),
                        "degree_range": [int(deg_hist.min()), int(deg_hist.max())]},
            "fixed_work": {**summarize(d_fw, it_fw),
                           "sampled_mean_N2": float(n_work ** 2),
                           "mean_degree": float(n_work),
                           "degree_range": [n_work, n_work]},
        }
    except Exception as exc:
        return {"design": design, "index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def ratios(row: dict) -> dict:
    a, f = row["atallah"], row["fixed_work"]
    def r(key, invert=True):
        va, vf = abs(a[key]), abs(f[key])
        if va == 0.0:
            return None
        return vf / va          # >1 favors Atallah (smaller defect)
    return {
        "defect_rms_ratio": r("defect_rms_m_s2"),
        "in_track_rms_ratio": r("in_track_rms_m_s2"),
        "in_track_impulse_ratio": r("in_track_impulse_m_s"),
        "displacement_proxy_ratio": r("in_track_displacement_proxy_m"),
        "sampled_work_ratio_atallah_over_fixed": (a["sampled_mean_N2"]
                                                  / f["sampled_mean_N2"]),
        "atallah_favored_by_defect_rms": bool(
            a["defect_rms_m_s2"] < f["defect_rms_m_s2"]),
        "atallah_favored_by_displacement_proxy": bool(
            abs(a["in_track_displacement_proxy_m"])
            < abs(f["in_track_displacement_proxy_m"])),
    }


def stat(v):
    x = np.asarray([q for q in v if q is not None and np.isfinite(q)], dtype=float)
    if x.size == 0:
        return None
    return {"n": int(x.size), "median": float(np.median(x)),
            "p10": float(np.percentile(x, 10)), "p90": float(np.percentile(x, 90)),
            "min": float(x.min()), "max": float(x.max())}


def summarize(rows: list[dict]) -> dict:
    rs = [r["ratios"] for r in rows]
    return {
        "orbits": len(rows),
        "defect_rms_ratio": stat([x["defect_rms_ratio"] for x in rs]),
        "in_track_rms_ratio": stat([x["in_track_rms_ratio"] for x in rs]),
        "displacement_proxy_ratio": stat([x["displacement_proxy_ratio"] for x in rs]),
        "sampled_work_ratio": stat([x["sampled_work_ratio_atallah_over_fixed"]
                                    for x in rs]),
        "atallah_smaller_defect_rms": int(sum(x["atallah_favored_by_defect_rms"]
                                              for x in rs)),
        "atallah_smaller_displacement_proxy": int(sum(
            x["atallah_favored_by_displacement_proxy"] for x in rs)),
        "atallah_defect_rms_m_s2": stat([r["atallah"]["defect_rms_m_s2"] for r in rows]),
        "fixed_defect_rms_m_s2": stat([r["fixed_work"]["defect_rms_m_s2"] for r in rows]),
        "atallah_displacement_proxy_m": stat(
            [abs(r["atallah"]["in_track_displacement_proxy_m"]) for r in rows]),
        "fixed_displacement_proxy_m": stat(
            [abs(r["fixed_work"]["in_track_displacement_proxy_m"]) for r in rows]),
    }


def build_table(payload: dict) -> str:
    def sci(v, digits=2):
        e = int(math.floor(math.log10(abs(v))))
        return f"${v / 10.0 ** e:.{digits}f}\\times 10^{{{e}}}$"

    def row(d):
        s = payload["designs"][d]["summary"]
        return (f"    {d} & {s['orbits']} & "
                f"{s['sampled_work_ratio']['median']:.3f} & "
                f"{sci(s['atallah_defect_rms_m_s2']['median'])} & "
                f"{sci(s['fixed_defect_rms_m_s2']['median'])} & "
                f"{s['defect_rms_ratio']['median']:.1f} "
                f"[{s['defect_rms_ratio']['p10']:.1f}, {s['defect_rms_ratio']['p90']:.1f}] & "
                f"{s['atallah_smaller_defect_rms']}/{s['orbits']} & "
                f"{s['displacement_proxy_ratio']['median']:.1f} & "
                f"{s['atallah_smaller_displacement_proxy']}/{s['orbits']}\\\\")
    body = "\n".join(row(d) for d in payload["designs"])
    return f"""% auto-generated by rev13_force_defect.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Integration-noise-free comparison at matched work. Along each
  orbit's archived truth trajectory, the truncation acceleration defect
  $\\Delta a = a_{{N_P}}-a_{{N_{{\\mathrm{{truth}}}}}}$ is evaluated for the Atallah
  degree history and for its work-matched fixed degree at the same 5041 epochs.
  $w$ is the ratio of the sampled $\\langle N^2\\rangle$ of the two histories (a
  check that the work match holds on the sampled epochs). The defect ratio is
  fixed over Atallah, so values above unity favor Atallah. The displacement
  proxy is the first-order in-track accumulation
  $\\int_0^T (T-\\tau)\\,\\Delta a_{{\\mathrm{{it}}}}(\\tau)\\,d\\tau$. Every quantity in
  this table is a deterministic function of the field and the degree histories,
  so none of it carries integration noise.}}
  \\label{{tab:force-defect}}
  \\begin{{tabular}}{{l r r r r r r r r}}
    \\toprule
    Design & orbits & $w$ & $|\\Delta a|_{{\\mathrm{{At}}}}$ & $|\\Delta a|_{{\\mathrm{{fix}}}}$
      & defect ratio & At.\\ smaller & displ.\\ ratio & At.\\ smaller\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B", "both"), default="both")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0, help="first N orbits only")
    a = ap.parse_args()
    designs = ["A", "B"] if a.design == "both" else [a.design]
    payload = {"schema": "r13_force_defect_v1", "created_utc": base.utc_now(),
               "reference_level": LEVEL, "designs": {}}
    if OUTPUT.exists():
        payload["designs"] = json.loads(OUTPUT.read_text()).get("designs", {})
    for d in designs:
        campaign = json.loads((METRICS / ("r12_atallah_campaign.json" if d == "A"
                               else "r12_atallah_campaign_designB.json")).read_text())
        idx = [int(r["sobol_index"]) for r in campaign["rows"]]
        if a.limit:
            idx = idx[:a.limit]
        tasks = [{"design": d, "index": i} for i in idx]
        t0 = time.time()
        rows, fails = [], []
        print(f"[force-defect] design {d}: {len(tasks)} orbits", flush=True)
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(worker, t): t for t in tasks}
            for n, fut in enumerate(as_completed(futs), 1):
                rec = fut.result()
                if rec["status"] != "complete":
                    fails.append(rec)
                    print(f"  !! {rec['index']:03d} {rec['message']}", flush=True)
                    continue
                rec["ratios"] = ratios(rec)
                rows.append(rec)
                if n % 8 == 0:
                    print(f"  [{n:3d}/{len(tasks)}] "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
        rows.sort(key=lambda r: r["sobol_index"])
        payload["designs"][d] = {"rows": rows, "failures": fails,
                                 "summary": summarize(rows)}
        s = payload["designs"][d]["summary"]
        print(f"[design {d}] defect ratio median {s['defect_rms_ratio']['median']:.2f}, "
              f"Atallah smaller in {s['atallah_smaller_defect_rms']}/{s['orbits']}; "
              f"displacement proxy median {s['displacement_proxy_ratio']['median']:.2f}, "
              f"Atallah smaller in {s['atallah_smaller_displacement_proxy']}/{s['orbits']}; "
              f"sampled work ratio median {s['sampled_work_ratio']['median']:.3f}",
              flush=True)
    base.atomic_json(OUTPUT, payload)
    TABLE.write_text(build_table(payload), encoding="utf-8")
    print(f"[written] {OUTPUT.name}, {TABLE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
