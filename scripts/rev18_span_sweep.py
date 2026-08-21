"""R18: is there an interior optimum between constant and radial allocation?

The budget campaign (R14) compares two endpoints of a spectrum: a constant
degree, which spends its budget uniformly, and the budget-recalibrated radial
rule, which concentrates it at perilune. The radial endpoint loses, and R14 also
shows the loss correlates with the width of the degree span the rule runs
(r = -0.37 and -0.27 against log rho). A constant degree is span 1; the radial
rule at beta = 1 runs a median span of 5.7. Nothing in the paper tests between
them, so "radial allocation spends the budget badly" and "this profile is too
aggressive at fixed budget" are not separated. This campaign separates them.

Family. For each orbit, let N_0 be the equal-budget constant degree and
N_A(h) the budget-calibrated Atallah table, both taken unchanged from the frozen
R14 records. Interpolate geometrically in degree,

    N_k(h) = round( s_k * N_0^(1-k) * N_A(h)^k ),   k in {0, 0.25, 0.5, 0.75, 1}

clipped to [FLOOR, adopted truth degree]. Geometric interpolation is the natural
choice because it makes the span multiplicative: span(k) = span(1)^k, so k moves
the aggressiveness on a scale where the endpoints are 1 and the rule's own span.
The scale s_k is bisected so that <N_k^2> over the archived truth epochs equals
beta * N_crit^2, i.e. every k spends the same declared per-call budget. At k = 0
and k = 1 the construction returns the two archived endpoints, which is used as
a correctness check rather than re-propagated.

What this can and cannot answer. It samples one interpolation path between two
policies at one budget; an interior minimum on this path is evidence that the
radial endpoint is over-aggressive, and a monotone path is evidence that
concentration itself is what costs. Neither statement is a claim about the
optimal allocation, which the O26 benchmark bounds separately.

Usage:
    python rev18_span_sweep.py run --design A --workers 11 --deadline-min 200
    python rev18_span_sweep.py summarize [--from-disk]
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
import rev12_atallah as at
import rev14_budget_trajectory as r14

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r18_cases"
RAW_ROOT = METRICS / "r18_raw"
def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def out_path(design: str, beta: float) -> Path:
    return METRICS / f"r18_span_sweep_{design}_{beta_tag(beta)}.json"

# Endpoints are archived, so only the interior of the path is propagated.
K_INTERIOR = (0.25, 0.50, 0.75)
K_ALL = (0.0, 0.25, 0.50, 0.75, 1.0)

LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP
BIN_KM = r14.BIN_KM
FLOOR = r14.FLOOR
WORK_TOLERANCE = 0.01     # same 1% work-match target as R14


def tag(k: float) -> str:
    return f"k_{k:.2f}"


def paths(design: str, beta: float, index: int, k: float, level: str):
    sub = f"{design}_{beta_tag(beta)}_{tag(k)}"
    return (CASE_ROOT / sub / f"sobolA_{index:03d}" / f"span_{level}.json",
            RAW_ROOT / sub / f"sobolA_{index:03d}" / f"span_{level}.npz")


def interpolated_table(table_a: dict, n0: int, k: float, scale: float,
                       cap: int) -> dict:
    """Geometric blend of the constant degree and the Atallah table."""
    out = {}
    for h, na in table_a.items():
        v = scale * (float(n0) ** (1.0 - k)) * (float(na) ** k)
        out[float(h)] = int(min(cap, max(FLOOR, round(v))))
    return out


def work_of(table: dict, h_km: np.ndarray) -> float:
    deg = r14_degrees_from_table(table, h_km)
    return float(np.mean(deg.astype(float) ** 2))


def r14_degrees_from_table(table: dict, h_km: np.ndarray) -> np.ndarray:
    hmin, hmax = min(table), max(table)
    keys = np.array(sorted(table))
    vals = np.array([table[key] for key in keys], dtype=int)
    hb = np.clip(BIN_KM * np.floor(h_km / BIN_KM), hmin, hmax)
    idx = np.clip(np.searchsorted(keys, hb - 1e-9), 0, len(keys) - 1)
    return vals[idx]


def calibrate_scale(table_a: dict, n0: int, k: float, cap: int,
                    h_km: np.ndarray, target_work: float) -> dict:
    """Bisect the multiplicative scale so <N_k^2> matches the budget.

    Work is monotone non-decreasing in the scale, but integer rounding makes it
    a staircase, so the bisection keeps the best-so-far rather than assuming it
    can drive the mismatch to zero.
    """
    lo, hi = 0.2, 5.0
    best = None
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        tbl = interpolated_table(table_a, n0, k, mid, cap)
        w = work_of(tbl, h_km)
        err = abs(w / target_work - 1.0)
        if best is None or err < best["mismatch"]:
            best = {"scale": mid, "table": tbl, "work": w, "mismatch": err}
        if w > target_work:
            hi = mid
        else:
            lo = mid
        if hi / lo - 1.0 < 1e-12:
            break
    best["attainable"] = bool(best["mismatch"] < WORK_TOLERANCE)
    return best


def degree_fn_of(table: dict):
    hmin, hmax = min(table), max(table)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, BIN_KM * math.floor(h_m / 1e3 / BIN_KM)))
        return table[hb]

    return degree_of


def worker(task: dict) -> dict:
    index, k = int(task["index"]), float(task["k"])
    design = task["design"]
    beta = float(task["beta"])
    try:
        row = task["row"]
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        n0 = int(row["fixed_degree"])
        tol_a = float(row["atallah_tol_accel_m_s2"])
        hp_km = float(row["design_point"]["hp_km"])
        ha_km = float(row["design_point"]["ha_km"])
        y0 = np.asarray(row["initial_state_si"], dtype=float)

        model, args = r14._model(adopted)
        g = r14._g(adopted)
        _, table_a = at.atallah_binned_schedule(
            model, g, tol_a, hp_km, ha_km, floor=FLOOR, cap=adopted,
            bin_km=BIN_KM)
        table_a = {float(a): int(b) for a, b in table_a.items()}

        # Altitude history of the archived truth: the same epochs R14 calibrated
        # its own budget on, so the budgets are defined identically.
        h_km = np.asarray(task["truth_alt_km"], dtype=float)
        target = beta * n_crit ** 2
        cal = calibrate_scale(table_a, n0, k, adopted, h_km, target)
        table_k = cal["table"]
        degs = list(table_k.values())
        span = max(degs) / max(1, min(degs))

        cfg = {
            "sobol_index": index, "design": design, "k": k, "beta": beta,
            "adopted_truth_degree": adopted, "n_critical": n_crit,
            "constant_degree_endpoint": n0,
            "atallah_tol_accel_m_s2": tol_a,
            "family": "N_k(h) = round(scale * N_0^(1-k) * N_A(h)^k), "
                      "clipped to [floor, truth]",
            "scale": cal["scale"], "work_achieved": cal["work"],
            "work_target": target,
            "work_mismatch": cal["work"] / target - 1.0,
            "work_attainable": cal["attainable"],
            "degree_table": {str(a): int(b) for a, b in table_k.items()},
            "degree_span": span,
            # span of the k = 1 endpoint, recorded here so the summary need not
            # rebuild the Atallah table (which would need the loaded model)
            "atallah_span": (max(table_a.values())
                             / max(1, min(table_a.values()))),
            "initial_state_si": [float(v) for v in y0],
            "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
            "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
            "atol_kind": "vector", "timing_comparable": False,
            "source": task["provenance"]}

        degree_fn = degree_fn_of(table_k)
        for level in ("tight", "tighter"):
            sidecar, raw = paths(design, beta, index, k, level)
            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == base.object_hash(cfg) \
                        and prev.get("status") == "ok":
                    continue
            tol = LEVELS[level]
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, degree_fn, args,
                tol["rtol"], tol["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"index": index, "k": k, "status": "numerical_failure",
                        "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r18_span_sweep_v1", "created_utc": base.utc_now(),
                "config": cfg, "config_sha256": base.object_hash(cfg),
                "status": "ok", "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
        return {"index": index, "k": k, "status": "ok",
                "span": span, "work_mismatch": cfg["work_mismatch"]}
    except Exception as exc:                                  # noqa: BLE001
        return {"index": index, "k": k, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def load_rows(design: str, beta: float):
    """Join the frozen R14 trajectory record with the design source, which is
    where the initial states live."""
    src = json.loads((METRICS / f"r14_trajectory_{design}_{beta_tag(beta)}.json"
                      ).read_text(encoding="utf-8"))
    design_src = json.loads(
        r14.DESIGNS[design]["rows"].read_text(encoding="utf-8"))
    drows = (design_src["rows"] if isinstance(design_src, dict)
             else design_src)
    by_index = {int(d.get("sobol_index", d.get("index", -1))): d for d in drows}
    out = []
    for r in src["rows"]:
        if r.get("censored"):
            continue
        d = by_index.get(int(r["sobol_index"]))
        if d is None:
            continue
        state = (d.get("initial_state_si")
                 or d.get("design_point", {}).get("initial_state_si")
                 or d.get("state_si"))
        if state is None:
            continue
        merged = dict(r)
        merged["initial_state_si"] = [float(v) for v in state]
        out.append(merged)
    return out


def truth_altitudes(design: str, index: int) -> np.ndarray | None:
    """Altitude history of the archived truth, tight level."""
    _, raw = r14.reuse_paths(design, index, "truth", "tight")
    if not raw.exists():
        return None
    d = np.load(raw)
    y = d["state_si"]
    model, _ = r14._model(300)
    r = np.linalg.norm(y[:3], axis=0)
    return (r - model.r_ref) / 1e3


def run(args) -> int:
    design, beta = args.design, float(args.beta)
    rows = load_rows(design, beta)
    if args.limit:
        rows = rows[:args.limit]
    prov = {"driver": "rev18_span_sweep.py",
            "reuses": f"r14_trajectory_{design}_{beta_tag(beta)}.json (tolerances, "
                      "comparator degree, truth), r11 convergence truths",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"}
    tasks = []
    skipped = []
    for row in rows:
        idx = int(row["sobol_index"])
        alt = truth_altitudes(design, idx)
        if alt is None:
            skipped.append(idx)
            continue
        for k in K_INTERIOR:
            tasks.append({"index": idx, "k": k, "row": row,
                          "design": design, "beta": beta,
                          "truth_alt_km": alt.tolist(),
                          "provenance": prov})
    print(f"[r18-{design}-{beta_tag(beta)}] {len(rows)} orbits, "
          f"{len(tasks)} trajectories "
          f"({len(K_INTERIOR)} interior k values x 2 levels each), "
          f"workers={args.workers}")
    if skipped:
        print(f"[r18] skipped (no archived truth): {skipped}")

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
                    print(f"  [FAIL] orbit={res['index']:02d} k={res['k']:.2f} "
                          f"{res['status']}: {res.get('detail')}")
                else:
                    print(f"  [{done}/{len(tasks)}] orbit={res['index']:02d} "
                          f"k={res['k']:.2f} span={res['span']:.2f} "
                          f"dW={res['work_mismatch']:+.4f} "
                          f"elapsed={(time.time()-t0)/60:.1f}min")
                if time.time() > deadline:
                    print(f"[r18] deadline at {args.deadline_min} min; "
                          f"cancelling pending")
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r18] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace(from_disk=True, design=design,
                                        beta=beta))


# ------------------------------------------------------------------ summary
def _load(path):
    d = np.load(path)
    return d["t_s"], d["state_si"]


def endpoint_record(row: dict, which: str, design: str, beta: float) -> dict:
    """k = 0 and k = 1 come from the frozen R14 record, not re-propagated."""
    pol = row["policies"]["fixed_budget" if which == "k0" else "atallah_budget"]
    return {"error_m": pol["error_tighter"]["pos_rms_m"],
            "envelope_m": pol["truth_inclusive_envelope_m"],
            "source": f"r14_trajectory_{design}_{beta_tag(beta)} (frozen)"}


def orbit_summary(design: str, beta: float, row: dict) -> dict | None:
    index = int(row["sobol_index"])
    truth = {}
    for lv in ("tight", "tighter"):
        _, raw = r14.reuse_paths(design, index, "truth", lv)
        if not raw.exists():
            return None
        truth[lv] = _load(raw)
    truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                   truth["tighter"][0], truth["tighter"][1]
                                   )["pos_rms_m"]

    entries = {}
    # archived endpoints; the k = 1 span is read from an interior sidecar,
    # which recorded the Atallah table it interpolated from
    entries["0.00"] = {**endpoint_record(row, "k0", design, beta), "span": 1.0}
    entries["1.00"] = {**endpoint_record(row, "k1", design, beta), "span": None}

    for k in K_INTERIOR:
        got = {}
        for lv in ("tight", "tighter"):
            sidecar, raw = paths(design, beta, index, k, lv)
            if not (sidecar.exists() and raw.exists()):
                return None
            t, y = _load(raw)
            got[lv] = (t, y, json.loads(sidecar.read_text(encoding="utf-8")))
        err = {lv: base.common_error(got[lv][0], got[lv][1],
                                     truth[lv][0], truth[lv][1])["pos_rms_m"]
               for lv in ("tight", "tighter")}
        self_diff = base.common_error(got["tight"][0], got["tight"][1],
                                      got["tighter"][0], got["tighter"][1]
                                      )["pos_rms_m"]
        cfg = got["tighter"][2]["config"]
        # Cost bookkeeping is read from the tight level only, because that is
        # the level R14 recorded its endpoint costs at; mixing levels would
        # compare step counts taken under different tolerances.
        tel_t = got["tight"][2].get("telemetry") or {}
        n_rhs = tel_t.get("n_rhs")
        mds = tel_t.get("mean_degree_sq")
        entries[f"{k:.2f}"] = {
            "error_m": err["tighter"], "error_tight_m": err["tight"],
            "self_difference_rms_m": self_diff,
            "envelope_m": self_diff + truth_self,
            "span": cfg["degree_span"], "scale": cfg["scale"],
            "work_mismatch": cfg["work_mismatch"],
            "work_attainable": cfg["work_attainable"],
            "n_rhs_tight": n_rhs,
            "total_quadratic_work_tight": (mds * n_rhs
                                           if (mds and n_rhs) else None),
            "source": "propagated"}
        entries["1.00"]["span"] = cfg.get("atallah_span")

    # endpoint costs come from the frozen R14 record, same level
    cost = row.get("cost", {})
    entries["0.00"]["n_rhs_tight"] = cost.get("rhs_fixed")
    entries["0.00"]["total_quadratic_work_tight"] = cost.get(
        "total_quadratic_work_fixed")
    entries["1.00"]["n_rhs_tight"] = cost.get("rhs_atallah")
    entries["1.00"]["total_quadratic_work_tight"] = cost.get(
        "total_quadratic_work_atallah")
    w0 = entries["0.00"]["total_quadratic_work_tight"]
    r0 = entries["0.00"]["n_rhs_tight"]
    for e in entries.values():
        e["total_work_ratio_vs_constant"] = (
            e["total_quadratic_work_tight"] / w0
            if (w0 and e.get("total_quadratic_work_tight")) else None)
        e["rhs_ratio_vs_constant"] = (
            e["n_rhs_tight"] / r0 if (r0 and e.get("n_rhs_tight")) else None)

    def resolved_better(ka: str, kb: str) -> bool:
        """Is ka better than kb by more than the summed envelope?"""
        a, b = entries[ka], entries[kb]
        if a.get("error_m") is None or b.get("error_m") is None:
            return False
        thr = (a.get("envelope_m") or 0.0) + (b.get("envelope_m") or 0.0)
        return (b["error_m"] - a["error_m"]) > thr

    best_k = min(entries, key=lambda kk: entries[kk]["error_m"])
    e0, eb = entries["0.00"]["error_m"], entries[best_k]["error_m"]
    thr = entries["0.00"]["envelope_m"] + entries[best_k]["envelope_m"]
    interior = best_k not in ("0.00", "1.00")

    # The location of the optimum inside the interior is a separate and much
    # weaker claim than its existence: at the best k the error often falls to
    # the numerical floor, so neighboring interior points need not be
    # distinguishable even when both clearly beat the endpoints.
    interior_ks = [f"{k:.2f}" for k in K_INTERIOR]
    distinguishable = {
        f"{best_k}_vs_{other}": resolved_better(best_k, other)
        for other in interior_ks if other != best_k}
    return {
        "sobol_index": index, "name": row["name"],
        "hp_km": row["design_point"]["hp_km"],
        "ha_km": row["design_point"]["ha_km"],
        "incl_deg": row["design_point"].get("incl_deg"),
        "n_critical": row["n_critical"],
        "constant_degree": row["fixed_degree"],
        "truth_self_difference_rms_m": truth_self,
        "entries": entries,
        "best_k": best_k,
        "interior_optimum": interior,
        "best_beats_constant_resolved": bool(
            interior and (e0 - eb) > thr and eb < e0),
        "best_beats_radial_resolved": bool(
            interior and resolved_better(best_k, "1.00")),
        "best_beats_both_endpoints_resolved": bool(
            interior and resolved_better(best_k, "0.00")
            and resolved_better(best_k, "1.00")),
        "best_distinguishable_from_other_interior": distinguishable,
    }


def summarize(args) -> int:
    design = getattr(args, "design", "A")
    beta = float(getattr(args, "beta", 1.0))
    src = json.loads(
        (METRICS / f"r14_trajectory_{design}_{beta_tag(beta)}.json"
         ).read_text(encoding="utf-8"))
    rows = [orbit_summary(design, beta, r) for r in src["rows"]
            if not r.get("censored")]
    rows = [r for r in rows if r]
    if not rows:
        print("[r18] no complete orbits yet")
        return 1
    counts = {}
    for r in rows:
        counts[r["best_k"]] = counts.get(r["best_k"], 0) + 1
    interior = [r for r in rows if r["interior_optimum"]]
    resolved = [r for r in interior if r["best_beats_constant_resolved"]]
    both = [r for r in interior if r["best_beats_both_endpoints_resolved"]]
    # how often the optimum's location inside the interior is itself resolved
    located = [r for r in interior
               if all(r["best_distinguishable_from_other_interior"].values())]
    med = {}
    for k in K_ALL:
        key = f"{k:.2f}"
        vals = [r["entries"][key]["error_m"] for r in rows
                if key in r["entries"]]
        spans = [r["entries"][key]["span"] for r in rows
                 if key in r["entries"] and r["entries"][key]["span"]]
        wr = [r["entries"][key]["total_work_ratio_vs_constant"] for r in rows
              if key in r["entries"]
              and r["entries"][key].get("total_work_ratio_vs_constant")]
        rr = [r["entries"][key]["rhs_ratio_vs_constant"] for r in rows
              if key in r["entries"]
              and r["entries"][key].get("rhs_ratio_vs_constant")]
        med[key] = {"n": len(vals),
                    "median_error_m": float(np.median(vals)) if vals else None,
                    "median_span": float(np.median(spans)) if spans else None,
                    "median_total_work_ratio_vs_constant":
                        float(np.median(wr)) if wr else None,
                    "median_rhs_ratio_vs_constant":
                        float(np.median(rr)) if rr else None}
    payload = {
        "schema": "r18_span_sweep_v1", "created_utc": base.utc_now(),
        "design": design, "beta": beta,
        "family": "N_k(h) = round(scale * N_0^(1-k) * N_A(h)^k); "
                  "scale bisected to <N_k^2> = beta * N_crit^2",
        "k_values": list(K_ALL),
        "endpoints_reused": f"k=0 and k=1 read from r14_trajectory_{design}_{beta_tag(beta)}",
        "summary": {
            "orbits": len(rows),
            "best_k_counts": counts,
            "orbits_with_interior_best": len(interior),
            "interior_best_resolved_against_constant": len(resolved),
            "interior_best_resolved_against_radial": sum(
                1 for r in interior if r["best_beats_radial_resolved"]),
            "interior_best_resolved_against_both": len(both),
            "interior_best_location_resolved": len(located),
            "by_k": med},
        "rows": rows}
    out = out_path(design, beta)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r18-{design}-{beta_tag(beta)}] written {out.name}: "
          f"{len(rows)} orbits")
    print(f"  best-k counts: {counts}")
    print(f"  interior best: {len(interior)}; resolved vs constant "
          f"{len(resolved)}, vs radial "
          f"{sum(1 for r in interior if r['best_beats_radial_resolved'])}, "
          f"vs both {len(both)}; location resolved {len(located)}")
    for k in K_ALL:
        m = med[f"{k:.2f}"]
        w = m["median_total_work_ratio_vs_constant"]
        print(f"    k={k:.2f}  n={m['n']:2d}  "
              f"err={m['median_error_m']:.3f}  "
              f"span={m['median_span']:.2f}  "
              f"work x{w:.3f}" if w else
              f"    k={k:.2f}  n={m['n']:2d}  err={m['median_error_m']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--design", default="A")
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=200.0)
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--beta", type=float, default=1.0)
    r.set_defaults(func=run)
    s = sub.add_parser("summarize")
    s.add_argument("--from-disk", action="store_true")
    s.add_argument("--design", default="A")
    s.add_argument("--beta", type=float, default=1.0)
    s.set_defaults(func=summarize)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
