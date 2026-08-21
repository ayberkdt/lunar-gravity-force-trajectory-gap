"""R58: the budget-calibrated radial *endpoint* against a constant degree
matched on realized total quadratic work.

Why this exists
---------------
R19 already runs an equal-realized-work comparison, but only for the interior
member k=0.50 (its ``K_TARGET`` is a module constant).  The main text's
wide-elliptic result is about the *endpoint*, k=1.00, and it is reported at
equal *nominal per-call* budget while the endpoint spends a median 1.31 times
the comparator's realized quadratic work.  A reader is entitled to ask whether
the endpoint wins there because of the geometry or because the adaptive
integrator hands it more work.  This driver answers that by giving the constant
comparator the endpoint's realized work and re-scoring.

Nothing here writes into an R19 path.  The record, the case tree and the raw
tree all carry an r58 prefix, so the sealed R19 campaign is untouched and the
two experiments can be compared side by side.

K_TARGET is a module constant rather than a flag on purpose: the pool workers
are spawned on Windows and re-import this module, so a value patched at runtime
in the parent would not reach them and the children would silently write the
wrong comparator into the wrong tree.

Usage
-----
    python rev58_endpoint_equal_work.py run --design OE --beta 1.00 --workers 3
    python rev58_endpoint_equal_work.py summarize --design OE --beta 1.00
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev18_span_sweep as r18
import rev19_equal_total_work as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r58_cases"
RAW_ROOT = METRICS / "r58_raw"

K_TARGET = "1.00"          # the budget-calibrated radial endpoint

# The wide-elliptic populations are not in rev14's static registry; the strata
# driver installs them from argv at import time. Doing it that way here would
# be a spawn trap, because a pool worker re-imports this module without the
# parent's argv and would resolve the wrong rows. Both keys are therefore
# registered unconditionally, so parent and child agree by construction. The
# file names are the ones revJ2_fullforce.py already validated against the
# archived trajectories.
_WIDE_ELLIPTIC = {
    "OE":  ("r31_operational_elliptical_rows.json",
            "stratum_operational_elliptical_convergence"),
    "OEU": ("r38_operational_elliptical_uncapped_rows.json",
            "stratum_operational_elliptical_uncapped_convergence"),
}
for _key, (_rows, _conv) in _WIDE_ELLIPTIC.items():
    r14.DESIGNS.setdefault(_key, {
        "rows": METRICS / _rows,
        "reuse_case": METRICS / "r11_cases" / _conv,
        "reuse_raw": METRICS / "r11_raw" / _conv,
    })

LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def out_path(design: str, beta: float) -> Path:
    return METRICS / f"r58_endpoint_equal_work_{design}_{r18.beta_tag(beta)}.json"


def paths(design: str, beta: float, index: int, level: str):
    sub = f"{design}_endpoint_workmatched_{r18.beta_tag(beta)}"
    return (CASE_ROOT / sub / f"sobolA_{index:03d}" / f"fixed_{level}.json",
            RAW_ROOT / sub / f"sobolA_{index:03d}" / f"fixed_{level}.npz")


def worker(task: dict) -> dict:
    index = int(task["index"])
    design = task["design"]
    beta = float(task["beta"])
    try:
        row = task["row"]
        adopted = int(row["adopted_truth_degree"])
        degree = int(task["degree"])
        y0 = np.asarray(row["initial_state_si"], dtype=float)
        model, args = r14._model(adopted)

        cfg = {
            "sobol_index": index, "design": design, "beta": beta,
            "policy": "constant degree matched on realized total quadratic work",
            "matched_to": f"span-sweep member k={K_TARGET} "
                          f"(budget-calibrated radial endpoint)",
            "degree": degree,
            "target_total_quadratic_work": task["target_work"],
            "constant_endpoint_degree": int(row["fixed_degree"]),
            "constant_endpoint_total_work": task["work_constant"],
            "first_estimate_rule": "N* = round(N_0 * sqrt(W_k / W_0)); the "
                                   "achieved ratio is measured, not assumed",
            "adopted_truth_degree": adopted,
            "initial_state_si": [float(v) for v in y0],
            "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
            "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
            "atol_kind": "vector", "timing_comparable": False,
            "source": task["provenance"]}

        def degree_of(t, h_m):
            return degree

        for level in ("tight", "tighter"):
            sidecar, raw = paths(design, beta, index, level)
            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == base.object_hash(cfg) \
                        and prev.get("status") == "ok":
                    continue
            sidecar.parent.mkdir(parents=True, exist_ok=True)
            raw.parent.mkdir(parents=True, exist_ok=True)
            tol = LEVELS[level]
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, degree_of, args,
                tol["rtol"], tol["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"status": "numerical_failure", "index": index,
                        "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r58_endpoint_equal_work_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": base.object_hash(cfg), "status": "ok",
                "level": level, "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
        return {"status": "ok", "index": index, "degree": degree}
    except Exception as exc:                       # noqa: BLE001
        import traceback
        return {"status": "error", "index": index,
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def build_tasks(design: str, beta: float) -> list:
    tag = r18.beta_tag(beta)
    span = json.loads(
        (METRICS / f"r18_span_sweep_{design}_{tag}.json"
         ).read_text(encoding="utf-8"))
    rows14 = {int(r["sobol_index"]): r for r in r18.load_rows(design, beta)}
    prov = {"driver": "rev58_endpoint_equal_work.py",
            "reuses": f"r18_span_sweep_{design}_{tag}.json (realized work), "
                      f"r14_trajectory_{design}_{tag}.json, r11 truths",
            "matched_member": K_TARGET}
    tasks, censored, identical = [], [], []
    for r in span["rows"]:
        idx = int(r["sobol_index"])
        e_k = r["entries"].get(K_TARGET, {})
        e_0 = r["entries"].get("0.00", {})
        wk = e_k.get("total_quadratic_work_tight")
        w0 = e_0.get("total_quadratic_work_tight")
        n0 = r.get("constant_degree")
        if not (wk and w0 and n0) or idx not in rows14:
            continue
        degree = int(round(n0 * (wk / w0) ** 0.5))
        adopted = int(rows14[idx]["adopted_truth_degree"])
        if degree >= adopted:
            # A comparator at or above the reference degree cannot be scored:
            # its error would vanish by construction. Clamping would hand the
            # comparison to the constant side, so the orbit is censored and the
            # censoring is written out rather than absorbed.
            censored.append({"sobol_index": idx, "requested_degree": degree,
                             "adopted_truth_degree": adopted,
                             "reason": "work-matched degree at or above "
                                       "reference"})
            continue
        if degree == n0:
            identical.append(idx)      # work gap under half an integer degree
            continue
        tasks.append({"index": idx, "design": design, "beta": beta,
                      "row": rows14[idx], "degree": degree,
                      "target_work": wk, "work_constant": w0,
                      "provenance": prov})
    if censored:
        (METRICS / f"r58_censored_{design}_{tag}.json").write_text(
            json.dumps(censored, indent=2), encoding="utf-8")
        print(f"[r58-{design}-{tag}] censored {len(censored)} orbits "
              f"(work-matched degree at or above reference)", flush=True)
    if identical:
        print(f"[r58-{design}-{tag}] {len(identical)} orbits keep the constant "
              f"endpoint as comparator (work gap under one degree)", flush=True)
    return tasks


def summarize(args) -> int:
    design = args.design
    beta = float(args.beta)
    tag = r18.beta_tag(beta)
    span = json.loads(
        (METRICS / f"r18_span_sweep_{design}_{tag}.json"
         ).read_text(encoding="utf-8"))
    rows = []
    for r in span["rows"]:
        idx = int(r["sobol_index"])
        e_k = r["entries"].get(K_TARGET, {})
        e_0 = r["entries"].get("0.00", {})
        if e_k.get("error_m") is None:
            continue

        truth, ok = {}, True
        for lv in ("tight", "tighter"):
            _, raw = r14.reuse_paths(design, idx, "truth", lv)
            if not raw.exists():
                ok = False
                break
            truth[lv] = r19._load(raw)
        if not ok:
            continue
        truth_self = base.common_error(
            truth["tight"][0], truth["tight"][1],
            truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]

        sidecar_t, raw_t = paths(design, beta, idx, "tight")
        _, raw_g = paths(design, beta, idx, "tighter")
        if sidecar_t.exists() and raw_t.exists() and raw_g.exists():
            st = json.loads(sidecar_t.read_text(encoding="utf-8"))
            tel = st.get("telemetry") or {}
            deg = int(st["config"]["degree"])
            w_fixed = (tel.get("mean_degree_sq") or deg ** 2) * tel["n_rhs"]
            got_t, got_g = r19._load(raw_t), r19._load(raw_g)
            err = base.common_error(got_g[0], got_g[1],
                                    truth["tighter"][0], truth["tighter"][1]
                                    )["pos_rms_m"]
            self_diff = base.common_error(got_t[0], got_t[1],
                                          got_g[0], got_g[1])["pos_rms_m"]
            source = "propagated"
        else:
            n0 = r.get("constant_degree")
            wk = e_k.get("total_quadratic_work_tight")
            w0 = e_0.get("total_quadratic_work_tight")
            if not (n0 and wk and w0):
                continue
            if int(round(n0 * (wk / w0) ** 0.5)) != n0:
                continue                      # pending, not equal by accident
            deg, w_fixed = n0, w0
            err = e_0.get("error_m")
            self_diff = max((e_0.get("envelope_m") or 0.0) - truth_self, 0.0)
            source = "constant endpoint (work gap below one degree)"

        if err is None or w_fixed is None:
            continue
        env_k = e_k.get("envelope_m") or 0.0
        env_f = self_diff + truth_self
        diff = err - e_k["error_m"]
        thr = env_k + env_f
        rows.append({
            "sobol_index": idx, "name": r["name"], "hp_km": r["hp_km"],
            "ha_km": r.get("ha_km"), "n_critical": r["n_critical"],
            "constant_endpoint_degree": r.get("constant_degree"),
            "work_matched_degree": deg,
            "radial_endpoint_error_m": e_k["error_m"],
            "work_matched_error_m": err,
            "achieved_total_work_ratio": (
                e_k.get("total_quadratic_work_tight") / w_fixed
                if w_fixed else None),
            "rho_workmatched": (err / e_k["error_m"]
                                if e_k["error_m"] else None),
            "resolution_threshold_m": thr,
            "resolved": bool(abs(diff) > thr),
            "winner": ("radial" if diff > thr else
                       ("constant" if -diff > thr else None)),
            "comparator_source": source,
        })

    if not rows:
        print("[r58] nothing to summarize yet")
        return 1
    res = [r for r in rows if r["resolved"]]
    wins = sum(1 for r in res if r["winner"] == "radial")
    ratios = [r["achieved_total_work_ratio"] for r in rows
              if r["achieved_total_work_ratio"]]
    rho = [r["rho_workmatched"] for r in rows if r["rho_workmatched"]]
    payload = {
        "schema": "r58_endpoint_equal_work_v1", "created_utc": base.utc_now(),
        "design": design, "beta": beta, "matched_member_k": K_TARGET,
        "what_is_held_equal": "realized total quadratic work at the tight "
                              "level, not nominal per-call work",
        "compare_with": f"r14_trajectory_{design}_{tag}.json holds the same "
                        f"pair at equal nominal per-call budget",
        "summary": {
            "orbits": len(rows),
            "resolved": len(res),
            "resolved_radial_wins": wins,
            "resolved_constant_wins": len(res) - wins,
            "unresolved": len(rows) - len(res),
            "median_rho": float(np.median(rho)) if rho else None,
            "achieved_work_ratio": {
                "median": float(np.median(ratios)),
                "min": float(np.min(ratios)),
                "max": float(np.max(ratios))} if ratios else None},
        "rows": rows}
    out_path(design, beta).write_text(json.dumps(payload, indent=2),
                                      encoding="utf-8")
    s = payload["summary"]
    print(f"[r58-{design}-{tag}] written {out_path(design, beta).name}: "
          f"{s['orbits']} orbits", flush=True)
    print(f"  resolved {s['resolved']}: radial {s['resolved_radial_wins']}, "
          f"constant {s['resolved_constant_wins']}, "
          f"unresolved {s['unresolved']}", flush=True)
    print(f"  median rho (constant/radial) = {s['median_rho']}", flush=True)
    if s["achieved_work_ratio"]:
        a = s["achieved_work_ratio"]
        print(f"  achieved total-work ratio: median {a['median']:.4f} "
              f"[{a['min']:.4f}, {a['max']:.4f}]", flush=True)
    return 0


def run(args) -> int:
    design = args.design
    beta = float(args.beta)
    tag = r18.beta_tag(beta)
    tasks = build_tasks(design, beta)
    print(f"[r58-{design}-{tag}] {len(tasks)} orbits need a distinct "
          f"work-matched comparator, workers={args.workers}", flush=True)
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
                    print(f"  [FAIL] orbit={res['index']:02d} "
                          f"{res['status']}: {res.get('detail')}", flush=True)
                else:
                    print(f"  [{done}/{len(tasks)}] orbit={res['index']:02d} "
                          f"N*={res['degree']:3d} "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
                if time.time() > deadline:
                    print(f"[r58-{design}] deadline; cancelling pending",
                          flush=True)
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r58-{design}-{tag}] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min", flush=True)
    return summarize(argparse.Namespace(design=design, beta=beta))


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--design", default="OE")
    r.add_argument("--beta", type=float, default=1.00)
    r.add_argument("--workers", type=int, default=3)
    r.add_argument("--deadline-min", type=float, default=120.0)
    r.set_defaults(func=run)
    s = sub.add_parser("summarize")
    s.add_argument("--design", default="OE")
    s.add_argument("--beta", type=float, default=1.00)
    s.set_defaults(func=summarize)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
