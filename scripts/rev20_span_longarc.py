"""R20: does the interior optimum survive a sixty-day arc?

R18 finds an interior member of the constant-to-radial interpolation family
better than either endpoint at equal nominal per-call budget, over seven days.
This paper's own long-arc evidence says truncation error accumulates
superlinearly, so the seven-day ranking is not automatically the sixty-day one:
a policy that wins on a week could lose on two months if its defect is the more
coherent of the two.

This campaign propagates the same family, at the same budget, for 60 days on
the orbits where a 60-day truth already exists. Nothing about the family is
re-derived: the per-orbit degree tables are rebuilt from the frozen R18
configuration, so the only thing that changes is the arc.

Truths are reused from R17 at both tolerance levels, which fixes the comparison
contract; the endpoints are propagated here rather than read from R14, because
R14 is a seven-day campaign and has no sixty-day endpoint to reuse.

Usage:
    python rev20_span_longarc.py run --workers 11 --deadline-min 200
    python rev20_span_longarc.py summarize
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev17_longarc60 as r17
import rev18_span_sweep as r18

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r20_cases"
RAW_ROOT = METRICS / "r20_raw"
OUT = METRICS / "r20_span_longarc.json"

K_ALL = (0.0, 0.25, 0.50, 0.75, 1.0)
BETA = 1.0
DESIGN = "A"

LEVELS = r17.LEVELS
DURATION = r17.DURATION
OUTPUT_STEP = r17.OUTPUT_STEP
MAX_STEP = r17.MAX_STEP
CHECKPOINTS = (7, 14, 28, 42, 60)


def paths(index: int, k: float, level: str):
    return (CASE_ROOT / f"k_{k:.2f}" / f"sobolA_{index:03d}" / f"span_{level}.json",
            RAW_ROOT / f"k_{k:.2f}" / f"sobolA_{index:03d}" / f"span_{level}.npz")


def truth_paths(r17_orbit_index: int, level: str):
    return (METRICS / "r17_cases" / "longarc60"
            / f"run_orbit_{r17_orbit_index:02d}" / f"truth_{level}.json",
            METRICS / "r17_raw" / "longarc60"
            / f"run_orbit_{r17_orbit_index:02d}" / f"truth_{level}.npz")


def eligible_orbits() -> list:
    """Orbits carrying both a frozen R18 span record and an R17 60-day truth."""
    span = json.loads(
        (METRICS / f"r18_span_sweep_{DESIGN}_beta_1.00.json"
         ).read_text(encoding="utf-8"))
    long60 = json.loads((METRICS / "r17_longarc60.json").read_text(
        encoding="utf-8"))
    by_name = {r["orbit"]["name"]: r for r in long60["rows"]
               if r.get("status") == "complete" and r.get("reached_full_arc")}
    rows14 = {int(r["sobol_index"]): r for r in r18.load_rows(DESIGN, BETA)}
    out = []
    for r in span["rows"]:
        name = r["name"]
        if name not in by_name:
            continue
        idx = int(r["sobol_index"])
        if idx not in rows14:
            continue
        out.append({"span_row": r, "r14_row": rows14[idx],
                    "r17_index": int(by_name[name]["orbit_index"]),
                    "sobol_index": idx, "name": name})
    return out


def degree_table_for(entry: dict, k: float) -> dict | None:
    """Rebuild the k-member's binned table from the frozen R18 sidecar."""
    if k == 0.0:
        n0 = entry["span_row"].get("constant_degree")
        return {0.0: int(n0)} if n0 else None
    sidecar, _ = r18.paths(DESIGN, BETA, entry["sobol_index"], k, "tighter")
    if k in (0.25, 0.50, 0.75):
        if not sidecar.exists():
            return None
        cfg = json.loads(sidecar.read_text(encoding="utf-8"))["config"]
        return {float(a): int(b) for a, b in cfg["degree_table"].items()}
    # k = 1: the budget-calibrated radial endpoint, rebuilt from its tolerance
    row = entry["r14_row"]
    adopted = int(row["adopted_truth_degree"])
    model, _ = r14._model(adopted)
    g = r14._g(adopted)
    import rev12_atallah as at
    _, table = at.atallah_binned_schedule(
        model, g, float(row["atallah_tol_accel_m_s2"]),
        float(row["design_point"]["hp_km"]), float(row["design_point"]["ha_km"]),
        floor=r18.FLOOR, cap=adopted, bin_km=r18.BIN_KM)
    return {float(a): int(b) for a, b in table.items()}


def degree_fn_of(table: dict):
    if len(table) == 1:
        only = next(iter(table.values()))

        def constant(t, h_m):
            return only
        return constant
    hmin, hmax = min(table), max(table)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, r18.BIN_KM * math.floor(h_m / 1e3 / r18.BIN_KM)))
        return table[hb]
    return degree_of


def worker(task: dict) -> dict:
    idx, k = int(task["sobol_index"]), float(task["k"])
    try:
        row = task["r14_row"]
        adopted = int(row["adopted_truth_degree"])
        y0 = np.asarray(row["initial_state_si"], dtype=float)
        table = task["table"]
        model, args = r14._model(adopted)

        degs = list(table.values())
        cfg = {
            "sobol_index": idx, "design": DESIGN, "k": k, "beta": BETA,
            "arc": "60 days", "duration_s": DURATION,
            "degree_table": {str(a): int(b) for a, b in table.items()},
            "degree_span": max(degs) / max(1, min(degs)),
            "table_source": ("frozen R18 seven-day configuration, reused "
                             "verbatim; only the arc length differs"),
            "adopted_truth_degree": adopted,
            "initial_state_si": [float(v) for v in y0],
            "output_step_s": OUTPUT_STEP, "max_step_s": MAX_STEP,
            "integrator": "InstrumentedDOP853", "atol_kind": "vector",
            "timing_comparable": False, "source": task["provenance"]}

        degree_fn = degree_fn_of(table)
        for level in ("tight", "tighter"):
            sidecar, raw = paths(idx, k, level)
            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == base.object_hash(cfg) \
                        and prev.get("status") == "ok":
                    continue
            tol = LEVELS[level]
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            # R17 stores the component tolerances separately; build the same
            # vector its own driver builds rather than assuming an "atol" key
            atol = r17.atol_vector(level)
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, degree_fn, args,
                tol["rtol"], atol, max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"sobol_index": idx, "k": k,
                        "status": "numerical_failure", "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r20_span_longarc_v1", "created_utc": base.utc_now(),
                "config": cfg, "config_sha256": base.object_hash(cfg),
                "status": "ok", "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
        return {"sobol_index": idx, "k": k, "status": "ok",
                "span": cfg["degree_span"]}
    except Exception as exc:                                   # noqa: BLE001
        return {"sobol_index": idx, "k": k, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def run(args) -> int:
    orbits = eligible_orbits()
    prov = {"driver": "rev20_span_longarc.py",
            "reuses": "r18 span tables (frozen), r17 60-day truths, "
                      "r14 initial states and tolerances",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"}
    tasks = []
    for e in orbits:
        for k in K_ALL:
            table = degree_table_for(e, k)
            if table is None:
                continue
            tasks.append({**e, "k": k, "table": table, "provenance": prov})
    print(f"[r20] {len(orbits)} orbits with a 60-day truth, {len(tasks)} "
          f"trajectories ({len(K_ALL)} k values x 2 levels), "
          f"workers={args.workers}")
    if args.limit:
        tasks = tasks[:args.limit]

    t0 = time.time()
    deadline = t0 + args.deadline_min * 60.0
    done = fail = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        try:
            for fut in as_completed(futs):
                res = fut.result()
                done += 1
                if res["status"] != "ok":
                    fail += 1
                    print(f"  [FAIL] orbit={res['sobol_index']:02d} "
                          f"k={res['k']:.2f} {res['status']}: "
                          f"{res.get('detail')}")
                else:
                    print(f"  [{done}/{len(tasks)}] "
                          f"orbit={res['sobol_index']:02d} k={res['k']:.2f} "
                          f"span={res['span']:.2f} "
                          f"elapsed={(time.time()-t0)/60:.1f}min")
                if time.time() > deadline:
                    print("[r20] deadline; cancelling pending")
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r20] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace())


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def summarize(args) -> int:
    rows = []
    for e in eligible_orbits():
        idx = e["sobol_index"]
        truth = {}
        ok = True
        for lv in ("tight", "tighter"):
            _, raw = truth_paths(e["r17_index"], lv)
            if not raw.exists():
                ok = False
                break
            truth[lv] = _load(raw)
        if not ok:
            continue
        truth_self = base.common_error(
            truth["tight"][0], truth["tight"][1],
            truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]

        entries = {}
        complete = True
        for k in K_ALL:
            got = {}
            for lv in ("tight", "tighter"):
                sc, raw = paths(idx, k, lv)
                if not (sc.exists() and raw.exists()):
                    complete = False
                    break
                got[lv] = _load(raw)
            if not complete:
                break
            err = base.common_error(got["tighter"][0], got["tighter"][1],
                                    truth["tighter"][0], truth["tighter"][1]
                                    )["pos_rms_m"]
            self_diff = base.common_error(got["tight"][0], got["tight"][1],
                                          got["tighter"][0], got["tighter"][1]
                                          )["pos_rms_m"]
            entries[f"{k:.2f}"] = {"error_m": err,
                                   "self_difference_rms_m": self_diff,
                                   "envelope_m": self_diff + truth_self}
        if not complete:
            continue

        best = min(entries, key=lambda kk: entries[kk]["error_m"])

        def beats(a, b):
            ea, eb = entries[a], entries[b]
            return (eb["error_m"] - ea["error_m"]) > (ea["envelope_m"]
                                                      + eb["envelope_m"])

        rows.append({
            "sobol_index": idx, "name": e["name"],
            "hp_km": e["span_row"]["hp_km"],
            "truth_self_difference_rms_m": truth_self,
            "entries": entries, "best_k": best,
            "interior_best": best not in ("0.00", "1.00"),
            "beats_constant_resolved": beats(best, "0.00"),
            "beats_radial_resolved": beats(best, "1.00"),
            "beats_both_resolved": beats(best, "0.00") and beats(best, "1.00"),
        })

    if not rows:
        print("[r20] no complete orbits yet")
        return 1
    med = {}
    for k in K_ALL:
        key = f"{k:.2f}"
        v = [r["entries"][key]["error_m"] for r in rows if key in r["entries"]]
        med[key] = float(np.median(v)) if v else None
    payload = {
        "schema": "r20_span_longarc_v1", "created_utc": base.utc_now(),
        "design": DESIGN, "beta": BETA, "arc_days": 60,
        "summary": {
            "orbits": len(rows),
            "interior_best": sum(1 for r in rows if r["interior_best"]),
            "beats_constant_resolved": sum(1 for r in rows
                                           if r["beats_constant_resolved"]),
            "beats_both_resolved": sum(1 for r in rows
                                       if r["beats_both_resolved"]),
            "median_error_by_k": med},
        "rows": rows}
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    s = payload["summary"]
    print(f"[r20] written {OUT.name}: {s['orbits']} orbits")
    print(f"  interior best {s['interior_best']}, resolved vs constant "
          f"{s['beats_constant_resolved']}, vs both {s['beats_both_resolved']}")
    for k in K_ALL:
        print(f"    k={k:.2f}  median 60-day error = {med[f'{k:.2f}']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=200.0)
    r.add_argument("--limit", type=int, default=0)
    r.set_defaults(func=run)
    s = sub.add_parser("summarize")
    s.set_defaults(func=summarize)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
