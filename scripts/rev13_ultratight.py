"""Targeted ultra-tighter reruns of the contested Atallah comparisons (R13).

The diagnosis (rev13_resolution_diagnosis.py) selects the orbits worth another
decade of tolerance: every matched-work comparison that already resolves (as a
validation set) and every one with 0.5 < M_res <= 1. For those orbits the truth,
the Atallah rule, and its work-matched fixed degree are propagated at a third
tolerance level,

    ultra:  rtol 3e-14, atol [1e-7 m] * 3 + [1e-10 m/s] * 3,

one decade below the ``tighter`` level of the campaign. The maximum step is kept
at the campaign's 60 s so that the tighter-to-ultra difference isolates the
tolerance change; step sensitivity is controlled separately elsewhere.

The comparison is then repeated one level up: errors are taken at the tighter
level against the tighter truth (already archived), and each policy's numerical
envelope becomes its tighter-to-ultra self-difference plus the truth's
tighter-to-ultra self-difference. If the gap is a physical truncation signal it
survives; if it is integration noise it shrinks with the envelope and the
comparison stays unresolved.

Usage:
    python rev13_ultratight.py pilot --workers 3
    python rev13_ultratight.py run --workers 5
    python rev13_ultratight.py aggregate
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

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SELECTION = METRICS / "r13_ultratight_selection.json"
OUTPUT = METRICS / "r13_ultratight.json"
TABLE = METRICS / "r13_ultratight_table.tex"
CASE_ROOT = METRICS / "r13_cases" / "ultratight"
RAW_ROOT = METRICS / "r13_raw" / "ultratight"

DESIGNS = {
    "A": {"r12_case": METRICS / "r12_cases" / "atallah",
          "r11_case": METRICS / "r11_cases" / "convergence",
          "r11_raw": METRICS / "r11_raw" / "convergence",
          "r12_raw": METRICS / "r12_raw" / "atallah"},
    "B": {"r12_case": METRICS / "r12_cases" / "atallah_designB",
          "r11_case": METRICS / "r11_cases" / "designB_convergence",
          "r11_raw": METRICS / "r11_raw" / "designB_convergence",
          "r12_raw": METRICS / "r12_raw" / "atallah_designB"},
}

ULTRA = {"rtol": 3.0e-14,
         "atol": np.array([1.0e-7] * 3 + [1.0e-10] * 3),
         "atol_position_m": 1.0e-7, "atol_velocity_m_s": 1.0e-10}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
POLICIES = ("truth", "atallah", "fixed_work_atallah")

_MODELS: dict[int, tuple] = {}


def _model(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def paths(design: str, index: int, policy: str):
    stem = f"{policy}_ultra"
    return (CASE_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.json",
            RAW_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.npz")


def binned_lookup(table: dict[str, int], bin_km: float = 10.0):
    """Rebuild the archived binned Atallah selector from its stored table, with
    the same clamping as rev12_atallah.atallah_binned_schedule."""
    tab = {float(k): int(v) for k, v in table.items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, bin_km * math.floor(h_m / 1e3 / bin_km)))
        return tab[hb]

    return degree_of


def worker(task: dict) -> dict:
    design, index = task["design"], task["index"]
    try:
        cfg_dir = DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
        at_cfg = json.loads((cfg_dir / "atallah_tight.json").read_text())["config"]
        fw_cfg = json.loads((cfg_dir / "fixed_work_atallah_tight.json").read_text())["config"]
        adopted = int(at_cfg["adopted_truth_degree"])
        n_work = int(fw_cfg["policy_spec"]["degree"])
        y0 = np.asarray(at_cfg["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        selectors = {
            "truth": (lambda t, h, n=adopted: n),
            "atallah": binned_lookup(at_cfg["atallah_degree_table"]),
            "fixed_work_atallah": (lambda t, h, n=n_work: n),
        }
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        out = {"design": design, "index": index, "status": "complete",
               "adopted_truth_degree": adopted, "n_work": n_work, "telemetry": {}}
        for policy in POLICIES:
            sidecar, raw = paths(design, index, policy)
            if sidecar.exists() and raw.exists():
                out["telemetry"][policy] = json.loads(sidecar.read_text())["telemetry"]
                continue
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, selectors[policy], args,
                ULTRA["rtol"], ULTRA["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"design": design, "index": index,
                        "status": "numerical_failure",
                        "where": f"{policy}/ultra", "message": fail}
            cfg = {"design": design, "sobol_index": index,
                   "adopted_truth_degree": adopted,
                   "n_critical": int(at_cfg["n_critical"]),
                   "initial_state_si": [float(v) for v in y0],
                   "atallah_tol_accel_m_s2": at_cfg["atallah_tol_accel_m_s2"],
                   "atallah_degree_table": at_cfg["atallah_degree_table"],
                   "policy": policy, "level": "ultra",
                   "policy_spec": ({"kind": "fixed_truth", "degree": adopted}
                                   if policy == "truth" else
                                   {"kind": "fixed_work", "degree": n_work}
                                   if policy == "fixed_work_atallah" else
                                   {"kind": "atallah_radial_adaptive_binned",
                                    "tol": at_cfg["atallah_tol_accel_m_s2"]}),
                   "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
                   "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
                   "atol_kind": "vector", "rtol": ULTRA["rtol"],
                   "atol_position_m": ULTRA["atol_position_m"],
                   "atol_velocity_m_s": ULTRA["atol_velocity_m_s"],
                   "timing_comparable": False,
                   "purpose": ("third tolerance level for the targeted "
                               "resolution retest of the matched-work "
                               "Atallah comparison"),
                   "source": base.provenance()}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r13_ultratight_trajectory_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": base.object_hash(cfg), "status": st,
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
            out["telemetry"][policy] = tel
        return out
    except Exception as exc:
        return {"design": design, "index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ------------------------------------------------------------------ aggregation
def _r11_raw(design, index, policy, level):
    return DESIGNS[design]["r11_raw"] / f"sobolA_{index:03d}" / f"{policy}_{level}.npz"


def _r12_raw(design, index, policy, level):
    return DESIGNS[design]["r12_raw"] / f"sobolA_{index:03d}" / f"{policy}_{level}.npz"


def orbit_summary(design: str, index: int) -> dict:
    truth_tighter = base.load_raw(_r11_raw(design, index, "truth", "tighter"))
    truth_ultra = base.load_raw(paths(design, index, "truth")[1])
    res = {"design": design, "sobol_index": index, "policies": {}}

    def diff(a, b):
        return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]

    truth_self = diff(truth_tighter, truth_ultra)
    res["truth_self_difference_tighter_to_ultra_m"] = truth_self
    for policy in ("atallah", "fixed_work_atallah"):
        p_tighter = base.load_raw(_r12_raw(design, index, policy, "tighter"))
        p_ultra = base.load_raw(paths(design, index, policy)[1])
        self_diff = diff(p_tighter, p_ultra)
        res["policies"][policy] = {
            "error_tighter_m": diff(p_tighter, truth_tighter),
            "error_ultra_m": diff(p_ultra, truth_ultra),
            "self_difference_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self}
    a = res["policies"]["atallah"]
    f = res["policies"]["fixed_work_atallah"]
    gap = abs(a["error_tighter_m"] - f["error_tighter_m"])
    thr = a["truth_inclusive_envelope_m"] + f["truth_inclusive_envelope_m"]
    gap_u = abs(a["error_ultra_m"] - f["error_ultra_m"])
    res["comparison"] = {
        "level": "tighter errors, tighter-to-ultra envelopes",
        "atallah_error_m": a["error_tighter_m"],
        "comparator_error_m": f["error_tighter_m"],
        "rho": (f["error_tighter_m"] / a["error_tighter_m"]
                if a["error_tighter_m"] > 0 else None),
        "absolute_error_difference_m": gap,
        "resolution_threshold_m": thr,
        "m_res": gap / thr,
        "resolved": bool(gap > thr),
        "winner_if_resolved": (("atallah" if a["error_tighter_m"] < f["error_tighter_m"]
                                else "fixed_work_atallah") if gap > thr else None),
        "gap_at_ultra_level_m": gap_u,
        "gap_sign_stable_tighter_to_ultra": bool(
            np.sign(a["error_tighter_m"] - f["error_tighter_m"])
            == np.sign(a["error_ultra_m"] - f["error_ultra_m"]))}
    return res


def previous_margin(design: str, index: int) -> dict:
    d = json.loads((METRICS / "r13_resolution_diagnosis.json").read_text())
    for r in d["designs"][design]["rows"]:
        if r["sobol_index"] == index:
            return r
    return {}


def aggregate() -> int:
    sel = json.loads(SELECTION.read_text())
    rows, missing = [], []
    for design in ("A", "B"):
        for index in sel[design]["all"]:
            if all(paths(design, index, p)[0].exists() for p in POLICIES):
                r = orbit_summary(design, index)
                prev = previous_margin(design, index)
                r["previous"] = {
                    "m_res": prev.get("m_res"),
                    "resolved": prev.get("resolved"),
                    "threshold_m": prev.get("threshold_m"),
                    "group": ("resolved" if index in sel[design]["resolved"]
                              else "borderline")}
                rows.append(r)
            else:
                missing.append({"design": design, "index": index})
    summ = {}
    for group in ("resolved", "borderline", "all"):
        sub = [r for r in rows
               if group == "all" or r["previous"]["group"] == group]
        if not sub:
            continue
        m = np.array([r["comparison"]["m_res"] for r in sub])
        m0 = np.array([r["previous"]["m_res"] for r in sub])
        thr = np.array([r["comparison"]["resolution_threshold_m"] for r in sub])
        thr0 = np.array([r["previous"]["threshold_m"] for r in sub])
        summ[group] = {
            "orbits": len(sub),
            "resolved_now": int(sum(r["comparison"]["resolved"] for r in sub)),
            "atallah_wins": int(sum(r["comparison"]["winner_if_resolved"] == "atallah"
                                    for r in sub)),
            "fixed_wins": int(sum(r["comparison"]["winner_if_resolved"]
                                  == "fixed_work_atallah" for r in sub)),
            "m_res_median_before": float(np.median(m0)),
            "m_res_median_after": float(np.median(m)),
            "threshold_median_before_m": float(np.median(thr0)),
            "threshold_median_after_m": float(np.median(thr)),
            "envelope_shrink_factor_median": float(np.median(thr0 / thr)),
            "gap_sign_stable": int(sum(
                r["comparison"]["gap_sign_stable_tighter_to_ultra"] for r in sub))}
    payload = {"schema": "r13_ultratight_v1",
               "level": {k: v for k, v in ULTRA.items() if k != "atol"},
               "max_step_s": MAX_STEP,
               "selection_rule": sel["rule"],
               "rows": rows, "missing": missing, "summary": summ}
    base.atomic_json(OUTPUT, payload)
    print(json.dumps(summ, indent=2))
    if missing:
        print(f"  missing orbits: {missing}")
    return 0


def run(tasks, workers) -> int:
    for t in tasks:
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / t["design"] / f"sobolA_{t['index']:03d}").mkdir(parents=True,
                                                                  exist_ok=True)
    t0 = time.time()
    print(f"[ultratight] orbits={len(tasks)} workers={workers}", flush=True)
    fails = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['design']}{rec['index']:03d} {rec['status']} "
                      f"{rec.get('where','')}: {rec.get('message','')}", flush=True)
            print(f"  [{n:2d}/{len(tasks)}] {rec['design']}{rec['index']:03d} "
                  f"{rec['status']} elapsed={(time.time()-t0)/60:.1f}min", flush=True)
    print(f"[ultratight] done failures={len(fails)} "
          f"wall={(time.time()-t0)/60:.1f}min", flush=True)
    return 0 if not fails else 3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("pilot", "run", "aggregate"))
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    if a.command == "aggregate":
        return aggregate()
    sel = json.loads(SELECTION.read_text())
    tasks = [{"design": d, "index": i} for d in ("A", "B") for i in sel[d]["all"]]
    if a.command == "pilot":
        # one already-resolved and two borderline orbits, design A
        pick = [t for t in tasks if t["design"] == "A"][:3]
        return run(pick, a.workers)
    return run(tasks, a.workers)


if __name__ == "__main__":
    raise SystemExit(main())
