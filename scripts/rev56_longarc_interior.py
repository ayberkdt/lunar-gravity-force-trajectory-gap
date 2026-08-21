"""R56 (O56): the interior member at sixty days, recalibrated on the sixty-day
arc and matched at the scoring tolerance.

(O30) propagated the seven-day allocation family for sixty days and found the
interior member no longer ahead. That result carries two confounds that it
cannot separate on its own:

  calibration  the k = 0.5 degree table is the frozen seven-day one, reused
               verbatim. Over sixty days the member spends a median 1.22 times
               the constant endpoint's per-call budget, so the two policies are
               no longer on the same budget at the horizon they are scored at.

  accounting   the comparison is the nominal per-call one. (O42), (O53), (O54)
               and (O55) all showed that where a crossing sits depends on
               whether realized work is matched at the tolerance the errors are
               scored at.

So "the interior advantage does not survive the horizon" and "a seven-day
allocation is the wrong allocation for a sixty-day problem" are not currently
distinguishable. This campaign separates them: the member is recalibrated so
that <N_k^2> equals the budget over the *sixty-day* reference epochs, and its
comparator is a constant degree matched on realized total quadratic work read
at the tighter level.

Scope is deliberately minimal. The panel is the eight Design-A orbits that
already carry a sixty-day reference; no population is selected here. One
budget, beta = 1. Two policies, k = 0 and k = 0.5. Those three choices are
fixed in the registration and are not revisited after the result is seen.

Usage:
    python rev56_longarc_interior.py plan
    python rev56_longarc_interior.py member     --workers 8
    python rev56_longarc_interior.py comparator --workers 8
    python rev56_longarc_interior.py summarize
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
import rev12_atallah as at
import rev14_budget_trajectory as r14
import rev17_longarc60 as r17
import rev18_span_sweep as r18
import rev20_span_longarc as r20

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r56_cases"
RAW_ROOT = METRICS / "r56_raw"
OUT = METRICS / "r56_longarc_interior.json"

BETA = 1.0
K_MEMBER = 0.50
DESIGN = "A"
LEVELS = r17.LEVELS
DURATION = r17.DURATION
OUTPUT_STEP = r17.OUTPUT_STEP
MAX_STEP = r17.MAX_STEP


def paths(index: int, which: str, level: str):
    return (CASE_ROOT / which / f"sobolA_{index:03d}" / f"arc_{level}.json",
            RAW_ROOT / which / f"sobolA_{index:03d}" / f"arc_{level}.npz")


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def sixtyday_altitudes(r17_index: int) -> np.ndarray | None:
    """Altitude history of the sixty-day reference, tight level.

    (O25) calibrates on the tight-level archived reference; using the same
    level here keeps the calibration convention identical, so the only thing
    that changes between the two campaigns is the arc.
    """
    _, raw = r20.truth_paths(r17_index, "tight")
    if not raw.exists():
        return None
    _, y = _load(raw)
    model, _ = r14._model(300)
    r = np.linalg.norm(y[:3], axis=0)
    return (r - model.r_ref) / 1e3


def member_table(entry: dict):
    """The k = 0.5 table with its scale bisected on the sixty-day epochs."""
    row = entry["r14_row"]
    adopted = int(row["adopted_truth_degree"])
    n_crit = int(row["n_critical"])
    n0 = int(row["fixed_degree"])
    model, margs = r14._model(adopted)
    g = r14._g(adopted)
    _, table_a = at.atallah_binned_schedule(
        model, g, float(row["atallah_tol_accel_m_s2"]),
        float(row["design_point"]["hp_km"]), float(row["design_point"]["ha_km"]),
        floor=r18.FLOOR, cap=adopted, bin_km=r18.BIN_KM)
    table_a = {float(a): int(b) for a, b in table_a.items()}
    h_km = sixtyday_altitudes(entry["r17_index"])
    if h_km is None:
        return None
    cal = r18.calibrate_scale(table_a, n0, K_MEMBER, adopted, h_km,
                              BETA * n_crit ** 2)
    return {"table": cal["table"], "scale": cal["scale"],
            "work_achieved": cal["work"], "work_mismatch": cal["mismatch"],
            "adopted": adopted, "n_crit": n_crit, "n0": n0}


def _work_tighter(sidecar: Path) -> tuple[float | None, int | None]:
    """Realized total quadratic work at the tighter level, from telemetry."""
    if not sidecar.exists():
        return None, None
    d = json.loads(sidecar.read_text(encoding="utf-8"))
    tel = d.get("telemetry") or {}
    n_rhs = tel.get("n_rhs")
    msq = tel.get("mean_degree_sq")
    if not n_rhs or not msq:
        return None, None
    return float(msq) * float(n_rhs), int(n_rhs)


def constant_endpoint_work(index: int) -> float | None:
    """W_0 at the tighter level from the archived (O30) constant endpoint."""
    sidecar, _ = r20.paths(index, 0.0, "tighter")
    w, _ = _work_tighter(sidecar)
    return w


def propagate(y0, adopted, degree_fn, cfg, index, which):
    model, margs = r14._model(adopted)
    for level in ("tight", "tighter"):
        sidecar, raw = paths(index, which, level)
        if sidecar.exists() and raw.exists():
            prev = json.loads(sidecar.read_text(encoding="utf-8"))
            if prev.get("config_sha256") == base.object_hash(cfg) \
                    and prev.get("status") == "ok":
                continue
        tol = LEVELS[level]
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        t, y, st, ev, fail, tel = base.propagate_event_instrumented(
            model, y0, DURATION, grid, degree_fn, margs,
            tol["rtol"], r17.atol_vector(level), max_step=MAX_STEP)
        if st == "numerical_failure":
            return f"numerical_failure: {fail}"
        base.atomic_npz(raw, t_s=t, state_si=y)
        base.atomic_json(sidecar, {
            "schema": "r56_longarc_interior_v1", "created_utc": base.utc_now(),
            "config": cfg, "config_sha256": base.object_hash(cfg),
            "status": "ok", "event": ev, "telemetry": tel,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
    return None


def member_worker(task: dict) -> dict:
    idx = int(task["index"])
    try:
        cal = task["cal"]
        table = {float(a): int(b) for a, b in cal["table"].items()}
        cfg = {"sobol_index": idx, "design": DESIGN, "policy": "interior",
               "k": K_MEMBER, "beta": BETA, "arc_days": 60,
               "calibrated_on": "sixty-day reference epochs, tight level",
               "degree_table": {f"{a}": b for a, b in table.items()},
               "scale": cal["scale"], "work_achieved": cal["work_achieved"],
               "work_mismatch": cal["work_mismatch"],
               "adopted_truth_degree": cal["adopted"],
               "n_critical": cal["n_crit"],
               "constant_degree_endpoint": cal["n0"],
               "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
               "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
               "atol_kind": "vector", "source": task["provenance"]}
        err = propagate(np.asarray(task["y0"], dtype=float), cal["adopted"],
                        r18.degree_fn_of(table), cfg, idx, "member")
        if err:
            return {"index": idx, "status": err}
        return {"index": idx, "status": "ok", "scale": cal["scale"],
                "mismatch": cal["work_mismatch"]}
    except Exception as exc:                                   # noqa: BLE001
        return {"index": idx, "status": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def comparator_worker(task: dict) -> dict:
    idx = int(task["index"])
    try:
        degree = int(task["degree"])
        cfg = {"sobol_index": idx, "design": DESIGN, "policy": "constant",
               "beta": BETA, "arc_days": 60, "degree": degree,
               "matched_on": "realized total quadratic work at the tighter "
                             "level against the recalibrated k = 0.5 member",
               "target_work_tighter": task["target_work"],
               "constant_endpoint_degree": task["n0"],
               "constant_endpoint_work_tighter": task["w0"],
               "adopted_truth_degree": task["adopted"],
               "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
               "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
               "atol_kind": "vector", "source": task["provenance"]}

        def degree_of(t, h_m):
            return degree

        err = propagate(np.asarray(task["y0"], dtype=float), task["adopted"],
                        degree_of, cfg, idx, "comparator")
        if err:
            return {"index": idx, "status": err}
        return {"index": idx, "status": "ok", "degree": degree}
    except Exception as exc:                                   # noqa: BLE001
        return {"index": idx, "status": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


PROV = {"driver": "rev56_longarc_interior.py",
        "reuses": "r17 sixty-day references (both levels), r20 constant-"
                  "endpoint telemetry, r14/r18 frozen configuration",
        "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"}


def build_member_tasks():
    tasks, missing = [], []
    for e in r20.eligible_orbits():
        cal = member_table(e)
        if cal is None:
            missing.append(e["sobol_index"])
            continue
        tasks.append({"index": e["sobol_index"], "cal": cal,
                      "y0": [float(v) for v in e["r14_row"]["initial_state_si"]],
                      "provenance": PROV})
    return tasks, missing


def build_comparator_tasks():
    tasks, missing, censored = [], [], []
    for e in r20.eligible_orbits():
        idx = e["sobol_index"]
        row = e["r14_row"]
        adopted = int(row["adopted_truth_degree"])
        n0 = int(row["fixed_degree"])
        wk, _ = _work_tighter(paths(idx, "member", "tighter")[0])
        w0 = constant_endpoint_work(idx)
        if wk is None or w0 is None:
            missing.append(idx)
            continue
        degree = int(round(n0 * (wk / w0) ** 0.5))
        if degree >= adopted:
            censored.append({"sobol_index": idx, "requested_degree": degree,
                             "adopted_truth_degree": adopted})
            continue
        tasks.append({"index": idx, "degree": degree, "n0": n0, "w0": w0,
                      "target_work": wk, "adopted": adopted,
                      "y0": [float(v) for v in row["initial_state_si"]],
                      "provenance": PROV})
    return tasks, missing, censored


def _drive(tasks, worker, label, workers, deadline_min):
    print(f"[r56-{label}] {len(tasks)} orbits, workers={workers}", flush=True)
    t0 = time.time()
    deadline = t0 + deadline_min * 60.0
    done = fail = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for fut in as_completed(futs):
            res = fut.result()
            done += 1
            if res["status"] != "ok":
                fail += 1
                print(f"  [FAIL] orbit={res['index']:02d} {res['status']}",
                      flush=True)
            else:
                print(f"  [{done}/{len(tasks)}] orbit={res['index']:02d} "
                      f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
            if time.time() > deadline:
                print(f"[r56-{label}] deadline; stopping", flush=True)
                break
    print(f"[r56-{label}] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min", flush=True)
    return fail


def plan(args) -> int:
    tasks, missing = build_member_tasks()
    print(f"[r56-plan] {len(tasks)} orbits eligible, {len(missing)} missing")
    for t in tasks:
        c = t["cal"]
        degs = list(c["table"].values())
        print(f"  orbit {t['index']:02d}: N0={c['n0']:3d} "
              f"scale={c['scale']:.4f} span={max(degs)/max(1,min(degs)):.2f} "
              f"work mismatch={c['work_mismatch']:+.4f}")
    return 0


def member(args) -> int:
    tasks, missing = build_member_tasks()
    if missing:
        print(f"[r56] no sixty-day reference for {missing}")
    return _drive(tasks, member_worker, "member", args.workers,
                  args.deadline_min)


def comparator(args) -> int:
    tasks, missing, censored = build_comparator_tasks()
    if missing:
        print(f"[r56] member telemetry not yet available for {missing}")
    if censored:
        (METRICS / "r56_censored.json").write_text(
            json.dumps(censored, indent=2), encoding="utf-8")
        print(f"[r56] censored (comparator at or above reference): "
              f"{[c['sobol_index'] for c in censored]}")
    for t in tasks:
        print(f"  orbit {t['index']:02d}: N0={t['n0']:3d} -> N*={t['degree']:3d}"
              f"  (Wk/W0 tighter = {t['target_work']/t['w0']:.4f})")
    return _drive(tasks, comparator_worker, "comparator", args.workers,
                  args.deadline_min)


def summarize(args) -> int:
    rows = []
    for e in r20.eligible_orbits():
        idx = e["sobol_index"]
        truth = {}
        ok = True
        for lv in ("tight", "tighter"):
            _, raw = r20.truth_paths(e["r17_index"], lv)
            if not raw.exists():
                ok = False
                break
            truth[lv] = _load(raw)
        if not ok:
            continue
        truth_self = base.common_error(
            truth["tight"][0], truth["tight"][1],
            truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]

        got = {}
        for which in ("member", "comparator"):
            side = {}
            for lv in ("tight", "tighter"):
                sc, raw = paths(idx, which, lv)
                if not (sc.exists() and raw.exists()):
                    side = None
                    break
                side[lv] = _load(raw)
            if side is None:
                got = None
                break
            got[which] = side
        if got is None:
            continue

        rec = {"sobol_index": idx, "name": e["name"],
               "hp_km": float(e["r14_row"]["design_point"]["hp_km"]),
               "truth_self_difference_rms_m": truth_self}
        for which in ("member", "comparator"):
            err = base.common_error(
                got[which]["tighter"][0], got[which]["tighter"][1],
                truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]
            self_diff = base.common_error(
                got[which]["tight"][0], got[which]["tight"][1],
                got[which]["tighter"][0], got[which]["tighter"][1])["pos_rms_m"]
            rec[f"{which}_error_m"] = err
            rec[f"{which}_self_difference_m"] = self_diff
            sc, _ = paths(idx, which, "tighter")
            cfg = json.loads(sc.read_text(encoding="utf-8"))["config"]
            w, _ = _work_tighter(sc)
            rec[f"{which}_work_tighter"] = w
            if which == "comparator":
                rec["comparator_degree"] = cfg["degree"]
            else:
                rec["member_scale"] = cfg["scale"]

        env = (rec["member_self_difference_m"]
               + rec["comparator_self_difference_m"] + 2 * truth_self)
        diff = rec["comparator_error_m"] - rec["member_error_m"]
        rec["resolution_threshold_m"] = env
        rec["resolved"] = bool(abs(diff) > env)
        rec["winner"] = ("interior" if diff > env
                         else "fixed" if -diff > env else None)
        rec["rho"] = (rec["comparator_error_m"] / rec["member_error_m"]
                      if rec["member_error_m"] else None)
        rec["achieved_work_ratio_tighter"] = (
            rec["comparator_work_tighter"] / rec["member_work_tighter"]
            if rec["member_work_tighter"] else None)
        rows.append(rec)

    if not rows:
        print("[r56] no complete orbits yet")
        return 1
    res = [r for r in rows if r["resolved"]]
    interior = sum(1 for r in res if r["winner"] == "interior")
    fixed = sum(1 for r in res if r["winner"] == "fixed")
    rhos = [r["rho"] for r in rows if r["rho"]]
    ratios = [r["achieved_work_ratio_tighter"] for r in rows
              if r["achieved_work_ratio_tighter"]]
    summary = {"orbits": len(rows), "resolved": len(res),
               "resolved_interior_wins": interior,
               "resolved_fixed_wins": fixed,
               "unresolved": len(rows) - len(res),
               "median_rho": float(np.median(rhos)) if rhos else None,
               "achieved_work_ratio_tighter": {
                   "median": float(np.median(ratios)) if ratios else None,
                   "min": float(min(ratios)) if ratios else None,
                   "max": float(max(ratios)) if ratios else None}}
    OUT.write_text(json.dumps(
        {"schema": "r56_longarc_interior_v1", "created_utc": base.utc_now(),
         "design": DESIGN, "beta": BETA, "k_member": K_MEMBER, "arc_days": 60,
         "what_is_held_equal": "realized total quadratic work at the tighter "
                               "level; the member is recalibrated so that "
                               "<N_k^2> meets the budget over the sixty-day "
                               "reference epochs",
         "summary": summary, "rows": rows}, indent=2), encoding="utf-8")
    print(f"[r56] written {OUT.name}: {summary['orbits']} orbits; "
          f"resolved {summary['resolved']}: interior {interior}, "
          f"fixed {fixed}, unresolved {summary['unresolved']}; "
          f"median rho {summary['median_rho']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("plan", plan), ("member", member),
                     ("comparator", comparator), ("summarize", summarize)):
        s = sub.add_parser(name)
        if name in ("member", "comparator"):
            s.add_argument("--workers", type=int, default=8)
            s.add_argument("--deadline-min", type=float, default=240.0)
        s.set_defaults(func=fn)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
