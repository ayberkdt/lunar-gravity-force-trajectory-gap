"""R44: the equal-realized-work comparison re-matched at the TIGHTER level.

Council round 13/14 (Sert Elestirmen, MAJOR): R19 matches the comparator on the
k = 0.5 member's realized total quadratic work read at the TIGHT level, while
every error in the comparison is read at the TIGHTER level. On the calls the
tighter propagation actually makes, the interior member outspends its
"work-matched" comparator by 9.5%/10.1% (A/B at beta = 1), so the printed
match is level-inconsistent with the metric it licenses.

This campaign closes the item by matching where the errors are read:

    W_k^tighter = mean(N_k^2)_calls * n_RHS(k)     (tighter level, from the
                                                    archived R18 telemetry)
    N* = round( N_0 * sqrt(W_k^tighter / W_0^tighter) )

with W_0^tighter = N_0^2 * n_RHS(fixed, tighter) taken from the archived
R11/R14 comparator telemetry. The comparator is propagated fresh at BOTH
tolerance levels (the envelope needs the pair), and the achieved work ratio is
measured at both levels and reported, so the published table can print the
convention next to its number instead of implying it.

Nothing in the R18/R19 archives is touched: this campaign has its own case
and raw trees, and its leaf directories carry the design in their name
(sobolA_/sobolB_), closing the naming defect the council flagged in R18.

Usage:
    python rev44_equal_work_tighter.py plan --design A --beta 1.00
    python rev44_equal_work_tighter.py run  --design A --beta 1.00 --workers 10
    python rev44_equal_work_tighter.py summarize --design A --beta 1.00
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
CASE_ROOT = METRICS / "r44_cases"
RAW_ROOT = METRICS / "r44_raw"
R18_CASES = METRICS / "r18_cases"

K_TARGET = "0.50"
LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def out_path(design: str, beta: float) -> Path:
    return METRICS / (f"r44_equal_work_tighter_{design}_"
                      f"{r18.beta_tag(beta)}.json")


def paths(design: str, beta: float, index: int, level: str):
    sub = f"{design}_workmatched_tighter_{r18.beta_tag(beta)}"
    leaf = f"sobol{design}_{index:03d}"
    return (CASE_ROOT / sub / leaf / f"fixed_{level}.json",
            RAW_ROOT / sub / leaf / f"fixed_{level}.npz")


def _r18_member_sidecar(design: str, beta: float, index: int, level: str):
    """Archived R18 member sidecar; leaf name is sobolA_ for every design,
    a recorded archive convention (config.design is authoritative)."""
    sub = f"{design}_{r18.beta_tag(beta)}_k_{K_TARGET}"
    return R18_CASES / sub / f"sobolA_{index:03d}" / f"span_{level}.json"


def _work_from_sidecar(p: Path):
    if not p.exists():
        return None, None
    d = json.loads(p.read_text(encoding="utf-8"))
    tel = d.get("telemetry") or {}
    n_rhs = tel.get("n_rhs")
    msq = tel.get("mean_degree_sq")
    if not n_rhs or not msq:
        return None, None
    return float(msq) * float(n_rhs), int(n_rhs)


def _fixed_tighter_work(design: str, beta: float, index: int, row14: dict,
                        n0: int, level: str = "tighter"):
    """Realized total quadratic work of the archived constant comparator at
    `level`: N_0^2 x n_RHS, with n_RHS from the archived sidecar telemetry.
    At beta = 1 the comparator is the reused R11 fixed_critical run; at other
    budgets it is the R14 fixed_budget run."""
    policies = (["fixed_critical", "fixed_budget"]
                if row14.get("reuse_fixed_critical") or abs(beta - 1.0) < 1e-12
                else ["fixed_budget", "fixed_critical"])
    for pol in policies:
        try:
            sidecar, _ = r14.reuse_paths(design, index, pol, level)
        except Exception:                                      # noqa: BLE001
            continue
        if not sidecar.exists():
            continue
        d = json.loads(sidecar.read_text(encoding="utf-8"))
        tel = d.get("telemetry") or {}
        if tel.get("n_rhs"):
            return float(n0) ** 2 * float(tel["n_rhs"]), pol
    return None, None


def build_tasks(design: str, beta: float, verbose: bool = False):
    tag = r18.beta_tag(beta)
    span = json.loads((METRICS / f"r18_span_sweep_{design}_{tag}.json"
                       ).read_text(encoding="utf-8"))
    rows14 = {int(r["sobol_index"]): r for r in r18.load_rows(design, beta)}
    prov = {"driver": "rev44_equal_work_tighter.py",
            "reuses": (f"r18_span_sweep_{design}_{tag}.json + r18 tighter "
                       "telemetry (member work), archived comparator "
                       "telemetry (fixed work), r11 truths"),
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"}
    tasks, censored, degenerate, missing = [], [], [], []
    for r in span["rows"]:
        idx = int(r["sobol_index"])
        if idx not in rows14:
            continue
        n0 = r.get("constant_degree")
        e_k = r["entries"].get(K_TARGET, {})
        if not n0 or e_k.get("error_m") is None:
            continue
        wk, _ = _work_from_sidecar(
            _r18_member_sidecar(design, beta, idx, "tighter"))
        w0, pol = _fixed_tighter_work(design, beta, idx, rows14[idx], int(n0))
        if wk is None or w0 is None:
            missing.append(idx)
            continue
        degree = int(round(n0 * (wk / w0) ** 0.5))
        adopted = int(rows14[idx]["adopted_truth_degree"])
        if degree >= adopted:
            censored.append({"sobol_index": idx, "requested_degree": degree,
                             "adopted_truth_degree": adopted,
                             "reason": "work-matched degree at or above truth"})
            continue
        if degree == int(n0):
            degenerate.append(idx)
            continue
        tasks.append({"index": idx, "design": design, "beta": beta,
                      "row": rows14[idx], "degree": degree,
                      "n0": int(n0), "target_work_tighter": wk,
                      "work_constant_tighter": w0,
                      "fixed_source_policy": pol, "provenance": prov})
        if verbose:
            print(f"  orbit {idx:02d}: N0={n0:3d} -> N*={degree:3d} "
                  f"(Wk/W0 tighter = {wk / w0:.4f}, fixed from {pol})")
    if censored:
        (METRICS / f"r44_censored_{design}_{tag}.json").write_text(
            json.dumps(censored, indent=2), encoding="utf-8")
    return tasks, censored, degenerate, missing


def plan(args) -> int:
    """Dry run: build every task, print the plan, propagate nothing."""
    design, beta = args.design, float(args.beta)
    tasks, censored, degenerate, missing = build_tasks(
        design, beta, verbose=True)
    per_run_s = 100.0                       # 7-day tight/tighter arc, measured
    total_h = len(tasks) * 2 * per_run_s / 3600.0
    print(f"[r44-plan {design}@{r18.beta_tag(beta)}] "
          f"{len(tasks)} orbits to propagate x 2 levels, "
          f"{len(degenerate)} degenerate (N*=N0, reuse endpoint), "
          f"{len(censored)} censored, {len(missing)} missing telemetry")
    print(f"  est. {total_h:.1f} h single-thread, "
          f"~{total_h / 10 * 60:.0f} min at 10 workers")
    return 0 if not missing else 1


def worker(task: dict) -> dict:
    index = int(task["index"])
    design = task["design"]
    beta = float(task["beta"])
    try:
        row = task["row"]
        adopted = int(row["adopted_truth_degree"])
        degree = int(task["degree"])
        y0 = np.asarray(row["initial_state_si"], dtype=float)
        model, margs = r14._model(adopted)
        cfg = {
            "sobol_index": index, "design": design, "beta": beta,
            "policy": ("constant degree matched on realized total quadratic "
                       "work at the tighter level"),
            "matched_to": f"span-sweep member k={K_TARGET}",
            "match_level": "tighter",
            "degree": degree,
            "target_total_quadratic_work_tighter": task["target_work_tighter"],
            "constant_endpoint_degree": task["n0"],
            "constant_endpoint_total_work_tighter":
                task["work_constant_tighter"],
            "fixed_work_source_policy": task["fixed_source_policy"],
            "first_estimate_rule": ("N* = round(N_0 * sqrt(W_k/W_0)) with "
                                    "both works at the tighter level; the "
                                    "achieved ratio is measured, not assumed"),
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
                model, y0, DURATION, grid, degree_of, margs,
                tol["rtol"], tol["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r44_equal_work_tighter_v1",
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


def run(args) -> int:
    design, beta = args.design, float(args.beta)
    tasks, _, _, missing = build_tasks(design, beta)
    if missing:
        print(f"[r44-{design}] ABORT: {len(missing)} orbits lack telemetry: "
              f"{missing}")
        return 1
    print(f"[r44-{design}-{r18.beta_tag(beta)}] {len(tasks)} orbits, "
          f"workers={args.workers}")
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
                          f"{res['status']}: {res.get('detail')}")
                else:
                    print(f"  [{done}/{len(tasks)}] orbit={res['index']:02d} "
                          f"N*={res['degree']:3d} "
                          f"elapsed={(time.time() - t0) / 60:.1f}min")
                if time.time() > deadline:
                    print(f"[r44-{design}] deadline; cancelling pending")
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r44-{design}-{r18.beta_tag(beta)}] {done} finished, "
          f"{fail} failed, wall={(time.time() - t0) / 60:.1f} min")
    return summarize(argparse.Namespace(design=design, beta=beta))


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def summarize(args) -> int:
    design, beta = args.design, float(args.beta)
    tag = r18.beta_tag(beta)
    span = json.loads((METRICS / f"r18_span_sweep_{design}_{tag}.json"
                       ).read_text(encoding="utf-8"))
    rows = []
    for r in span["rows"]:
        idx = int(r["sobol_index"])
        e_k = r["entries"].get(K_TARGET, {})
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

        wk_tighter, _ = _work_from_sidecar(
            _r18_member_sidecar(design, beta, idx, "tighter"))
        wk_tight, _ = _work_from_sidecar(
            _r18_member_sidecar(design, beta, idx, "tight"))

        sidecar_t, raw_t = paths(design, beta, idx, "tight")
        sidecar_g, raw_g = paths(design, beta, idx, "tighter")
        if not (sidecar_g.exists() and raw_g.exists() and raw_t.exists()):
            continue
        sg = json.loads(sidecar_g.read_text(encoding="utf-8"))
        st = json.loads(sidecar_t.read_text(encoding="utf-8"))
        deg = int(sg["config"]["degree"])
        telg, telt = sg.get("telemetry") or {}, st.get("telemetry") or {}
        wf_tighter = (telg.get("mean_degree_sq") or deg ** 2) * telg["n_rhs"]
        wf_tight = (telt.get("mean_degree_sq") or deg ** 2) * telt["n_rhs"]
        got_t, got_g = _load(raw_t), _load(raw_g)
        err = base.common_error(got_g[0], got_g[1],
                                truth["tighter"][0], truth["tighter"][1]
                                )["pos_rms_m"]
        self_diff = base.common_error(got_t[0], got_t[1],
                                      got_g[0], got_g[1])["pos_rms_m"]
        env_k = e_k.get("envelope_m") or 0.0
        env_f = self_diff + truth_self
        diff = err - e_k["error_m"]
        thr = env_k + env_f
        rows.append({
            "sobol_index": idx, "name": r.get("name"), "hp_km": r.get("hp_km"),
            "constant_endpoint_degree": r.get("constant_degree"),
            "work_matched_degree": deg,
            "interior_error_m": e_k["error_m"],
            "work_matched_error_m": err,
            "achieved_work_ratio_tighter": (wk_tighter / wf_tighter
                                            if wk_tighter and wf_tighter
                                            else None),
            "achieved_work_ratio_tight": (wk_tight / wf_tight
                                          if wk_tight and wf_tight else None),
            "rho_workmatched": (err / e_k["error_m"]
                                if e_k["error_m"] else None),
            "resolution_threshold_m": thr,
            "resolved": bool(abs(diff) > thr),
            "winner": ("interior" if diff > thr else
                       ("fixed" if -diff > thr else None))})

    if not rows:
        print("[r44] nothing to summarize yet")
        return 1
    res = [r for r in rows if r["resolved"]]
    wins = sum(1 for r in res if r["winner"] == "interior")
    rt = [r["achieved_work_ratio_tighter"] for r in rows
          if r["achieved_work_ratio_tighter"]]
    rho = [r["rho_workmatched"] for r in rows if r["rho_workmatched"]]
    payload = {
        "schema": "r44_equal_work_tighter_v1", "created_utc": base.utc_now(),
        "design": design, "beta": beta, "interior_member_k": K_TARGET,
        "what_is_held_equal": ("realized total quadratic work at the TIGHTER "
                               "level, the level every error in the "
                               "comparison is read at"),
        "summary": {
            "orbits": len(rows), "resolved": len(res),
            "resolved_interior_wins": wins,
            "resolved_fixed_wins": len(res) - wins,
            "unresolved": len(rows) - len(res),
            "median_rho": float(np.median(rho)) if rho else None,
            "achieved_work_ratio_tighter": {
                "median": float(np.median(rt)),
                "min": float(np.min(rt)),
                "max": float(np.max(rt))} if rt else None},
        "rows": rows}
    out_path(design, beta).write_text(json.dumps(payload, indent=2),
                                      encoding="utf-8")
    s = payload["summary"]
    print(f"[r44-{design}-{tag}] written {out_path(design, beta).name}: "
          f"{s['orbits']} orbits; resolved {s['resolved']}: "
          f"interior {s['resolved_interior_wins']}, "
          f"fixed {s['resolved_fixed_wins']}, "
          f"unresolved {s['unresolved']}; median rho {s['median_rho']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("plan", plan), ("run", run), ("summarize", summarize)):
        s = sub.add_parser(name)
        s.add_argument("--design", default="A")
        s.add_argument("--beta", type=float, default=1.0)
        if name == "run":
            s.add_argument("--workers", type=int, default=10)
            s.add_argument("--deadline-min", type=float, default=180.0)
        s.set_defaults(func=fn)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
