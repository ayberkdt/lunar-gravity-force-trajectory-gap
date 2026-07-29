"""R19: the interior policy against a constant degree matched on REALIZED work.

R18 holds the *nominal per-call* budget equal, which is the same quantity R14
shows the integrator does not preserve: at beta = 1 the interior member at
k = 0.5 ends up consuming about 13% more total quadratic work than the constant
degree it is compared with. Using that budget definition to condemn the radial
endpoint in one section and to license an interior optimum in the next is not
consistent, and the interior result is the weaker of the two for it.

This campaign removes the objection by changing what is held equal. For each
orbit it reads the realized total quadratic work of the k = 0.5 policy,

    W_k = <N_k^2>_calls * n_RHS(k)          (tight level, as archived by R18)

and propagates the constant degree N* whose realized total work matches it.
Because a constant degree's RHS count varies only weakly with the degree, the
first estimate

    N* = round( N_0 * sqrt(W_k / W_0) )

lands close; the achieved ratio is then *measured* from the run and reported
rather than assumed, so a reader sees how equal the comparison actually was.

If k = 0.5 still wins here, the constructive claim no longer depends on the
per-call abstraction at all. If it does not, the honest statement is that the
interior member buys accuracy with extra realized work, and the paper says so.

The campaign is run at more than one budget scale. beta = 1 is the anchor the
main text reports; beta = 0.5 is the scale at which Section 8.3 makes its
strongest constructive claim, and that claim has never been tested under this
accounting. Passing --beta selects the scale; the beta = 1 record paths are
left byte-identical to the original run so the archived campaign is untouched.

Usage:
    python rev19_equal_total_work.py run --design A --workers 11
    python rev19_equal_total_work.py run --design A --beta 0.5 --workers 6
    python rev19_equal_total_work.py summarize --design A --beta 0.5
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev18_span_sweep as r18

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r19_cases"
RAW_ROOT = METRICS / "r19_raw"

K_TARGET = "0.50"          # the interior member the main text reports
ANCHOR_BETA = 1.0          # the scale the original campaign ran at
LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def _suffix(beta: float) -> str:
    """Empty at the anchor scale, so the beta = 1 archive keeps its paths."""
    return "" if abs(beta - ANCHOR_BETA) < 1e-12 else f"_{r18.beta_tag(beta)}"


def out_path(design: str, beta: float = ANCHOR_BETA) -> Path:
    return METRICS / f"r19_equal_total_work_{design}{_suffix(beta)}.json"


def paths(design: str, beta: float, index: int, level: str):
    sub = f"{design}_workmatched{_suffix(beta)}"
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
            "matched_to": f"span-sweep member k={K_TARGET}",
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
            tol = LEVELS[level]
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, degree_of, args,
                tol["rtol"], tol["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r19_equal_total_work_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": base.object_hash(cfg), "status": "ok",
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
        return {"index": index, "status": "ok", "degree": degree}
    except Exception as exc:                                   # noqa: BLE001
        return {"index": index, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def build_tasks(design: str, beta: float = ANCHOR_BETA) -> list:
    """Read the archived R18 record for the realized-work targets."""
    tag = r18.beta_tag(beta)
    span = json.loads(
        (METRICS / f"r18_span_sweep_{design}_{tag}.json"
         ).read_text(encoding="utf-8"))
    rows14 = {int(r["sobol_index"]): r
              for r in r18.load_rows(design, beta)}
    prov = {"driver": "rev19_equal_total_work.py",
            "reuses": f"r18_span_sweep_{design}_{tag}.json (realized work), "
                      f"r14_trajectory_{design}_{tag}.json, r11 truths",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"}
    tasks = []
    censored = []
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
            # A comparator at or above the truth degree cannot be scored: its
            # error would be zero by construction. Clamping would silently
            # hand the comparison to the fixed side, so the orbit is censored
            # and the censoring is reported.
            censored.append({"sobol_index": idx, "requested_degree": degree,
                             "adopted_truth_degree": adopted,
                             "reason": "work-matched degree at or above truth"})
            continue
        if degree == n0:
            # the work gap is under half an integer degree; the constant
            # endpoint already is the work-matched comparator
            continue
        tasks.append({"index": idx, "design": design, "beta": beta,
                      "row": rows14[idx],
                      "degree": degree, "target_work": wk,
                      "work_constant": w0, "provenance": prov})
    if censored:
        (METRICS / f"r19_censored_{design}{_suffix(beta)}.json").write_text(
            json.dumps(censored, indent=2), encoding="utf-8")
        print(f"[r19-{design}-{tag}] censored {len(censored)} orbits "
              f"(work-matched degree at or above truth)")
    return tasks


def run(args) -> int:
    design = args.design
    beta = float(getattr(args, "beta", ANCHOR_BETA))
    tasks = build_tasks(design, beta)
    print(f"[r19-{design}-{r18.beta_tag(beta)}] {len(tasks)} orbits need a "
          f"distinct work-matched comparator, workers={args.workers}")
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
                    print(f"  [FAIL] orbit={res['index']:02d} {res['status']}: "
                          f"{res.get('detail')}")
                else:
                    print(f"  [{done}/{len(tasks)}] orbit={res['index']:02d} "
                          f"N*={res['degree']:3d} "
                          f"elapsed={(time.time()-t0)/60:.1f}min")
                if time.time() > deadline:
                    print(f"[r19-{design}] deadline; cancelling pending")
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r19-{design}-{r18.beta_tag(beta)}] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace(design=design, beta=beta))


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def summarize(args) -> int:
    design = args.design
    beta = float(getattr(args, "beta", ANCHOR_BETA))
    span = json.loads(
        (METRICS / f"r18_span_sweep_{design}_{r18.beta_tag(beta)}.json"
         ).read_text(encoding="utf-8"))
    rows = []
    for r in span["rows"]:
        idx = int(r["sobol_index"])
        e_k = r["entries"].get(K_TARGET, {})
        e_0 = r["entries"].get("0.00", {})
        if e_k.get("error_m") is None:
            continue

        truth = {}
        ok = True
        for lv in ("tight", "tighter"):
            _, raw = r14.reuse_paths(design, idx, "truth", lv)
            if not raw.exists():
                ok = False
                break
            truth[lv] = _load(raw)
        if not ok:
            continue
        truth_self = base.common_error(
            truth["tight"][0], truth["tight"][1],
            truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]

        sidecar_t, raw_t = paths(design, beta, idx, "tight")
        sidecar_g, raw_g = paths(design, beta, idx, "tighter")
        if sidecar_t.exists() and raw_t.exists() and raw_g.exists():
            st = json.loads(sidecar_t.read_text(encoding="utf-8"))
            tel = st.get("telemetry") or {}
            deg = int(st["config"]["degree"])
            w_fixed = (tel.get("mean_degree_sq") or deg ** 2) * tel["n_rhs"]
            got_t, got_g = _load(raw_t), _load(raw_g)
            err = base.common_error(got_g[0], got_g[1],
                                    truth["tighter"][0], truth["tighter"][1]
                                    )["pos_rms_m"]
            self_diff = base.common_error(got_t[0], got_t[1],
                                          got_g[0], got_g[1])["pos_rms_m"]
            source = "propagated"
        else:
            # The constant endpoint may legitimately *be* the work-matched
            # comparator, but only when the required degree rounds back to it.
            # For any other orbit the absence of a sidecar means the run has
            # not happened yet, and silently substituting the endpoint would
            # score the wrong comparator.
            n0 = r.get("constant_degree")
            wk = e_k.get("total_quadratic_work_tight")
            w0 = e_0.get("total_quadratic_work_tight")
            if not (n0 and wk and w0):
                continue
            if int(round(n0 * (wk / w0) ** 0.5)) != n0:
                continue                      # pending, not equal by accident
            deg = n0
            w_fixed = w0
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
            "n_critical": r["n_critical"],
            "constant_endpoint_degree": r.get("constant_degree"),
            "work_matched_degree": deg,
            "interior_error_m": e_k["error_m"],
            "work_matched_error_m": err,
            "achieved_total_work_ratio": (
                e_k.get("total_quadratic_work_tight") / w_fixed
                if w_fixed else None),
            "rho_workmatched": (err / e_k["error_m"]
                                if e_k["error_m"] else None),
            "resolution_threshold_m": thr,
            "resolved": bool(abs(diff) > thr),
            "winner": ("interior" if diff > thr else
                       ("fixed" if -diff > thr else None)),
            "comparator_source": source,
        })

    if not rows:
        print("[r19] nothing to summarize yet")
        return 1
    res = [r for r in rows if r["resolved"]]
    wins = sum(1 for r in res if r["winner"] == "interior")
    ratios = [r["achieved_total_work_ratio"] for r in rows
              if r["achieved_total_work_ratio"]]
    rho = [r["rho_workmatched"] for r in rows if r["rho_workmatched"]]
    payload = {
        "schema": "r19_equal_total_work_v1", "created_utc": base.utc_now(),
        "design": design, "beta": beta, "interior_member_k": K_TARGET,
        "what_is_held_equal": "realized total quadratic work at the tight "
                              "level, not nominal per-call work",
        "summary": {
            "orbits": len(rows),
            "resolved": len(res),
            "resolved_interior_wins": wins,
            "resolved_fixed_wins": len(res) - wins,
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
    print(f"[r19-{design}-{r18.beta_tag(beta)}] written "
          f"{out_path(design, beta).name}: {s['orbits']} orbits")
    print(f"  resolved {s['resolved']}: interior {s['resolved_interior_wins']}, "
          f"fixed {s['resolved_fixed_wins']}, unresolved {s['unresolved']}")
    print(f"  median rho (fixed/interior) = {s['median_rho']}")
    if s["achieved_work_ratio"]:
        a = s["achieved_work_ratio"]
        print(f"  achieved total-work ratio: median {a['median']:.4f} "
              f"[{a['min']:.4f}, {a['max']:.4f}]")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--design", default="A")
    r.add_argument("--beta", type=float, default=ANCHOR_BETA)
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=120.0)
    r.set_defaults(func=run)
    s = sub.add_parser("summarize")
    s.add_argument("--design", default="A")
    s.add_argument("--beta", type=float, default=ANCHOR_BETA)
    s.set_defaults(func=summarize)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
