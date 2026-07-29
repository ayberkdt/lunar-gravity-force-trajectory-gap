"""R23-B: the fixed oracle applied to the constructive claim.

R15-A represented the constant-degree family by its best member under the
budget, F_oracle, and used that stronger comparator to harden the paper's
negative result: the radial rule loses to the family, not merely to one
nominated member. The same control was never applied to the constructive
result, which is still scored against the single budget-saturating degree.

That asymmetry is load-bearing. On this panel the oracle is a median 2.75 times
more accurate than the saturating degree, while the interior member's entire
reported margin is a median near 2.56. Two numbers of the same size cannot both
be read as decisive, so the comparison has to be run.

Convention
----------
The R15 ladder was scored at the *tight* level against the tight truth, while
R18 reports the interior member at the *tighter* level against the tighter
truth. Those are different statistics and must not be divided into each other.
This campaign therefore re-propagates the two constant comparators at the
tighter level, which is what makes them commensurable with the interior member
and, separately, what supplies the self-difference each needs before the
paper's resolution rule can be applied at all.

Two comparators are propagated, not one:

  N_sat     the budget-saturating degree the paper already nominates. Its
            tighter-level error is archived in R18 as the k = 0 endpoint, so
            re-propagating it reproduces a number the archive already contains
            and checks this pipeline against it before any new claim rests on it.
  N_oracle  the best member of the ladder, from the frozen R15 offsets.

Usage:
    python rev23_oracle_vs_interior.py run --workers 11
    python rev23_oracle_vs_interior.py summarize
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
ORACLE_RECORD = METRICS / "r15_fixed_oracle.json"
R15_RAW = METRICS / "r15_raw" / "fixed_oracle"
CASE_ROOT = METRICS / "r23_cases" / "oracle_vs_interior"
RAW_ROOT = METRICS / "r23_raw" / "oracle_vs_interior"
OUTPUT = METRICS / "r23_oracle_vs_interior.json"

BETA = 1.00
K_TARGET = "0.50"
LEVEL = "tighter"
# how far the reproduction of the archived k = 0 endpoint may drift before the
# comparison is abandoned rather than patched
REPRODUCTION_TOLERANCE = 0.05


def paths(design: str, index: int, kind: str, degree: int):
    stem = f"{kind}_{degree}_{LEVEL}"
    return (CASE_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.json",
            RAW_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.npz")


def archived_tight_raw(design: str, index: int, degree: int) -> Path:
    return R15_RAW / design / f"sobolA_{index:03d}" / f"fixed_{degree}_tight.npz"


def load_panel() -> list[dict]:
    rec = json.loads(ORACLE_RECORD.read_text(encoding="utf-8"))
    return rec["rows"]


def _design_rows(design: str) -> dict:
    """R18's loader is used because it is what attaches the initial states."""
    return {int(r["sobol_index"]): r for r in r18.load_rows(design, BETA)}


def worker(task: dict) -> dict:
    design, index = task["design"], int(task["index"])
    kind, degree = task["kind"], int(task["degree"])
    try:
        sidecar, raw = paths(design, index, kind, degree)
        adopted = int(task["adopted_truth_degree"])
        y0 = np.asarray(task["initial_state_si"], dtype=float)
        cfg = {
            "sobol_index": index, "design": design, "beta": BETA,
            "policy": f"constant degree ({kind})", "degree": degree,
            "comparator_role": task["role"],
            "ladder_source": "r15_fixed_oracle.json, offsets frozen in R15-A",
            "adopted_truth_degree": adopted,
            "initial_state_si": [float(v) for v in y0],
            "duration_s": r14.DURATION, "output_step_s": r14.OUTPUT_STEP,
            "integrator": "InstrumentedDOP853", "max_step_s": r14.MAX_STEP,
            "level": LEVEL, "atol_kind": "vector", "timing_comparable": False,
            "source": task["provenance"]}

        if sidecar.exists() and raw.exists():
            prev = json.loads(sidecar.read_text(encoding="utf-8"))
            if prev.get("config_sha256") == base.object_hash(cfg) \
                    and prev.get("status") == "ok":
                return {"index": index, "design": design, "kind": kind,
                        "status": "cached", "degree": degree}

        model, args = r14._model(adopted)
        tol = r14.LEVELS[LEVEL]
        grid = np.arange(0.0, r14.DURATION + 0.5 * r14.OUTPUT_STEP,
                         r14.OUTPUT_STEP)
        t, y, st, ev, fail, tel = base.propagate_event_instrumented(
            model, y0, r14.DURATION, grid, lambda _t, _h, k=degree: k, args,
            tol["rtol"], tol["atol"], max_step=r14.MAX_STEP)
        if st == "numerical_failure":
            return {"index": index, "design": design, "kind": kind,
                    "status": "numerical_failure", "detail": fail}
        base.atomic_npz(raw, t_s=t, state_si=y)
        base.atomic_json(sidecar, {
            "schema": "r23_oracle_vs_interior_v1",
            "created_utc": base.utc_now(), "config": cfg,
            "config_sha256": base.object_hash(cfg), "status": "ok",
            "event": ev, "telemetry": tel,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
        return {"index": index, "design": design, "kind": kind,
                "status": "ok", "degree": degree}
    except Exception as exc:                                    # noqa: BLE001
        return {"index": index, "design": design, "kind": kind,
                "status": "error", "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def build_tasks() -> list:
    prov = {"driver": "rev23_oracle_vs_interior.py",
            "reuses": "r15_fixed_oracle.json (frozen ladder and its argmin), "
                      "r14_trajectory_{design}_beta_1.00.json (states, truths), "
                      "r11 truths",
            "preregistration": "r23_preregistration.json"}
    rows14 = {d: _design_rows(d) for d in ("A", "B")}
    tasks = []
    for row in load_panel():
        design, index = row["design"], int(row["sobol_index"])
        src = rows14[design].get(index)
        if src is None:
            continue
        common = {
            "design": design, "index": index,
            "adopted_truth_degree": int(src["adopted_truth_degree"]),
            "initial_state_si": src["initial_state_si"],
            "provenance": prov}
        tasks.append({**common, "kind": "sat", "degree": int(row["n_sat"]),
                      "role": "the budget-saturating degree the paper nominates"})
        if int(row["n_oracle"]) != int(row["n_sat"]):
            tasks.append({
                **common, "kind": "oracle", "degree": int(row["n_oracle"]),
                "role": "post-hoc best member of the frozen ladder under the "
                        "budget; a lower envelope, not a selectable policy"})
    return tasks


def run(args) -> int:
    tasks = build_tasks()
    print(f"[r23b] {len(tasks)} comparator propagations at the {LEVEL} level, "
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
                if res["status"] not in ("ok", "cached"):
                    fail += 1
                    print(f"  [FAIL] {res['design']}{res['index']:03d} "
                          f"{res['kind']}: {res.get('detail')}")
                else:
                    print(f"  [{done}/{len(tasks)}] {res['design']}"
                          f"{res['index']:03d} {res['kind']} N={res['degree']} "
                          f"{res['status']} "
                          f"elapsed={(time.time()-t0)/60:.1f}min")
                if time.time() > deadline:
                    print("[r23b] deadline; cancelling pending")
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r23b] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace())


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def _err(a, b) -> float:
    return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]


def summarize(args) -> int:
    panel = load_panel()
    spans = {d: json.loads(
        (METRICS / f"r18_span_sweep_{d}_beta_1.00.json"
         ).read_text(encoding="utf-8")) for d in ("A", "B")}
    span_rows = {d: {int(r["sobol_index"]): r for r in spans[d]["rows"]}
                 for d in ("A", "B")}

    rows, drift = [], []
    for prow in panel:
        design, index = prow["design"], int(prow["sobol_index"])
        sr = span_rows[design].get(index)
        if sr is None:
            continue
        interior = sr["entries"].get(K_TARGET, {})
        endpoint = sr["entries"].get("0.00", {})
        if interior.get("error_m") is None or endpoint.get("error_m") is None:
            continue

        _, truth_tight = r14.reuse_paths(design, index, "truth", "tight")
        _, truth_tighter = r14.reuse_paths(design, index, "truth", "tighter")
        if not (truth_tight.exists() and truth_tighter.exists()):
            continue
        tt, tg = _load(truth_tight), _load(truth_tighter)
        truth_self = _err(tt, tg)

        comp = {}
        ok = True
        for kind, degree in (("sat", int(prow["n_sat"])),
                             ("oracle", int(prow["n_oracle"]))):
            if kind == "oracle" and degree == int(prow["n_sat"]):
                comp["oracle"] = dict(comp["sat"])
                comp["oracle"]["is_sat"] = True
                continue
            _, raw = paths(design, index, kind, degree)
            arch = archived_tight_raw(design, index, degree)
            if not (raw.exists() and arch.exists()):
                ok = False
                break
            got = _load(raw)
            self_diff = _err(_load(arch), got)
            comp[kind] = {
                "degree": degree,
                "error_m": _err(got, tg),
                "envelope_m": self_diff + truth_self,
                "self_difference_rms_m": self_diff,
                "is_sat": kind == "sat"}
        if not ok:
            continue

        # The k = 0 endpoint of R18 is this same saturating degree at this same
        # level. Reproducing it is the check that this pipeline is on the
        # archive's footing before any new comparison is read off it.
        rel = abs(comp["sat"]["error_m"] - endpoint["error_m"]) \
            / endpoint["error_m"]
        drift.append(rel)

        def verdict(a: dict, b: dict) -> tuple[str | None, bool]:
            """a versus b under the paper's truth-inclusive resolution rule."""
            diff = b["error_m"] - a["error_m"]
            thr = (a.get("envelope_m") or 0.0) + (b.get("envelope_m") or 0.0)
            if abs(diff) <= thr:
                return None, False
            return ("a" if diff > 0 else "b"), True

        int_side = {"error_m": interior["error_m"],
                    "envelope_m": interior.get("envelope_m") or 0.0}
        w_or, r_or = verdict(int_side, comp["oracle"])
        w_sat, r_sat = verdict(int_side, comp["sat"])
        rows.append({
            "design": design, "sobol_index": index, "hp_km": prow["hp_km"],
            "n_critical": prow["n_critical"],
            "n_sat": comp["sat"]["degree"],
            "n_oracle": comp["oracle"]["degree"],
            "interior_error_m": interior["error_m"],
            "interior_envelope_m": int_side["envelope_m"],
            "sat_error_m": comp["sat"]["error_m"],
            "sat_envelope_m": comp["sat"]["envelope_m"],
            "oracle_error_m": comp["oracle"]["error_m"],
            "oracle_envelope_m": comp["oracle"]["envelope_m"],
            "archived_endpoint_error_m": endpoint["error_m"],
            "endpoint_reproduction_rel_diff": rel,
            "rho_interior_vs_sat": (comp["sat"]["error_m"]
                                    / interior["error_m"]),
            "rho_interior_vs_oracle": (comp["oracle"]["error_m"]
                                       / interior["error_m"]),
            "oracle_gain_over_sat": (comp["sat"]["error_m"]
                                     / comp["oracle"]["error_m"]),
            "resolved_vs_sat": r_sat,
            "winner_vs_sat": ({"a": "interior", "b": "sat"}.get(w_sat)),
            "resolved_vs_oracle": r_or,
            "winner_vs_oracle": ({"a": "interior", "b": "oracle"}.get(w_or)),
        })

    if not rows:
        print("[r23b] nothing to summarize yet")
        return 1

    def med(vals):
        return float(np.median(vals)) if vals else None

    res_or = [r for r in rows if r["resolved_vs_oracle"]]
    res_sat = [r for r in rows if r["resolved_vs_sat"]]
    worst_drift = max(drift) if drift else 0.0
    payload = {
        "schema": "r23_oracle_vs_interior_v1", "created_utc": base.utc_now(),
        "beta": BETA, "interior_member_k": K_TARGET, "level": LEVEL,
        "preregistration": "r23_preregistration.json",
        "ladder_source": ("frozen R15-A offsets [0,1,2,3,4,6,8,12,16,24] "
                          "below the budget-saturating degree"),
        "what_is_compared": (
            "the interior span member against the best constant degree under "
            "the same nominal per-call budget, on the R15-A panel, with both "
            "comparators propagated fresh at the tighter level so the interior "
            "member and the constant family are scored by the same statistic"),
        "endpoint_reproduction_check": {
            "statement": ("the saturating degree re-propagated here is the same "
                          "object as the R18 k = 0 endpoint; relative difference "
                          "in seven-day position RMS"),
            "worst_relative_difference": worst_drift,
            "tolerance": REPRODUCTION_TOLERANCE,
            "passes": bool(worst_drift <= REPRODUCTION_TOLERANCE)},
        "summary": {
            "orbits": len(rows),
            "orbits_where_oracle_is_sat": sum(1 for r in rows
                                              if r["n_oracle"] == r["n_sat"]),
            "vs_oracle": {
                "resolved": len(res_or),
                "interior_wins": sum(1 for r in res_or
                                     if r["winner_vs_oracle"] == "interior"),
                "oracle_wins": sum(1 for r in res_or
                                   if r["winner_vs_oracle"] == "oracle"),
                "unresolved": len(rows) - len(res_or),
                "median_rho": med([r["rho_interior_vs_oracle"] for r in rows])},
            "vs_saturating": {
                "resolved": len(res_sat),
                "interior_wins": sum(1 for r in res_sat
                                     if r["winner_vs_sat"] == "interior"),
                "sat_wins": sum(1 for r in res_sat
                                if r["winner_vs_sat"] == "sat"),
                "unresolved": len(rows) - len(res_sat),
                "median_rho": med([r["rho_interior_vs_sat"] for r in rows])},
            "median_oracle_gain_over_sat": med(
                [r["oracle_gain_over_sat"] for r in rows])},
        "rows": rows}
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    s = payload["summary"]
    print(f"[r23b] written {OUTPUT.name}: {s['orbits']} orbits")
    chk = payload["endpoint_reproduction_check"]
    print(f"  endpoint reproduction: worst {chk['worst_relative_difference']:.2%}"
          f" -> {'PASS' if chk['passes'] else 'FAIL'}")
    for label, key in (("vs oracle", "vs_oracle"),
                       ("vs saturating", "vs_saturating")):
        d = s[key]
        won = d.get("interior_wins")
        lost = d.get("oracle_wins", d.get("sat_wins"))
        print(f"  interior {label}: resolved {d['resolved']} "
              f"(interior {won}, comparator {lost}), "
              f"unresolved {d['unresolved']}, median rho {d['median_rho']:.3f}")
    if not chk["passes"]:
        print("[r23b] reproduction check FAILED; per the registration this "
              "comparison is abandoned rather than patched")
        return 2
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=90.0)
    r.set_defaults(func=run)
    s = sub.add_parser("summarize")
    s.set_defaults(func=summarize)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
