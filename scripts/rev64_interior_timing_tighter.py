"""R64 (O57): the measured-kernel-time comparator, matched at the tolerance
the errors are scored at.

(O48) matched the interior member's comparator on kernel time measured at the
*tight* level while every error in the comparison is read at the *tighter*
level. That is the same level inconsistency (O42) removed for realized work,
and it is why (O48) bounds cost without superseding the operation-count
comparison. This campaign repeats (O48) with the timing moved to the scoring
level.

What changes, and only this: the member's contention-free re-run, the
first-pass comparator and the refined comparator are all timed at the tighter
level, and the first-pass degree is inverted on the member's tighter call
histogram rather than its tight one. The refined comparator is still
propagated at both levels, because the envelope needs the pair, but the tight
run is no longer a timed one.

What does not change: the panel is (O48)'s fourteen orbits read from its
selection record, the member is k = 0.5, the budget is beta = 1, and no orbit
is selected here.

Timing comparability requires an idle machine: a pipeline invocation refuses
to start while other python processes are alive. The check runs once, at the
head of the invocation, and does not re-arm; the after-the-fact quiet audit is
what covers the rest of the run.

Usage:
    python rev64_interior_timing_tighter.py pipeline [--cutoff ...]
    python rev64_interior_timing_tighter.py aggregate
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev13_timing_match as r13
import rev14_budget_trajectory as r14
import rev48_interior_timing as r48

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r64_cases" / "interior_timing_tighter"
RAW_ROOT = METRICS / "r64_raw" / "interior_timing_tighter"
OUT = METRICS / "r64_interior_timing_tighter.json"
STATE = METRICS / "r64_interior_timing_state.json"

K_TARGET = r48.K_TARGET
COST_EXPONENT = r48.COST_EXPONENT
LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def case_dir(design: str, index: int) -> Path:
    return CASE_ROOT / design / f"sobol{design}_{index:03d}"


def raw_dir(design: str, index: int) -> Path:
    return RAW_ROOT / design / f"sobol{design}_{index:03d}"


def _propagate(cfg_extra: dict, degree_of, y0, adopted: int, level: str,
               sidecar: Path, raw: Path, timed: bool) -> dict | None:
    if sidecar.exists() and raw.exists():
        return json.loads(sidecar.read_text(encoding="utf-8"))
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    raw.parent.mkdir(parents=True, exist_ok=True)
    model = base.load_model(adopted)
    args = base.kernel_args(model)
    base.warmup(model, args)
    tol = LEVELS[level]
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    t, y, st, ev, fail, tel = base.propagate_event_instrumented(
        model, y0, DURATION, grid, degree_of, args,
        tol["rtol"], tol["atol"], max_step=MAX_STEP)
    if st == "numerical_failure":
        print(f"  !! numerical failure: {fail}", flush=True)
        return None
    base.atomic_npz(raw, t_s=t, state_si=y)
    payload = {"schema": "r64_interior_timing_tighter_case_v1",
               "created_utc": base.utc_now(),
               "config": dict(cfg_extra, level=level,
                              timing_comparable=timed,
                              adopted_truth_degree=adopted,
                              source=base.provenance()),
               "status": st, "telemetry": tel,
               "raw_path": str(raw.relative_to(ROOT)),
               "raw_sha256": base.file_hash(raw)}
    base.atomic_json(sidecar, payload)
    return payload


def panel() -> dict:
    """(O48)'s selection, reused verbatim. No orbit is chosen here."""
    sel = json.loads(r48.SELECTION.read_text(encoding="utf-8"))
    deg_tab, ns_tab = r13.cost_curve()
    out = {}
    for design, entries in sel["designs"].items():
        rows = []
        for e in entries:
            index = int(e["sobol_index"])
            # First-pass degree from the member's TIGHTER call histogram, the
            # level this campaign matches and scores at.
            side = r48.member_sidecar(design, index, "tighter")
            counts = {int(k): int(v) for k, v in
                      side["telemetry"]["degree_counts"].items()}
            n_rhs = sum(counts.values())
            mean_cost = sum(r13.cost_of([k], deg_tab, ns_tab)[0] * v
                            for k, v in counts.items()) / n_rhs
            rows.append({
                "sobol_index": index, "hp_km": e["hp_km"],
                "adopted_truth_degree": int(e["adopted_truth_degree"]),
                "n_time_tighter_histogram":
                    int(r13.inverse_cost(mean_cost, deg_tab, ns_tab)),
                "n_time_tight_histogram_r48": int(e["n_time_measured_cost"]),
                "member_mean_per_call_ns_tighter": float(mean_cost)})
        out[design] = rows
    return out


def pipeline(args) -> int:
    cutoff = datetime.fromisoformat(args.cutoff) if args.cutoff else None
    state = (json.loads(STATE.read_text(encoding="utf-8"))
             if STATE.exists() else {"designs": panel()})
    if not r48.idle_or_die():
        return 2
    t0 = time.time()
    for design, entries in state["designs"].items():
        for e in entries:
            if cutoff and datetime.now() > cutoff:
                print("[r64] cutoff reached; stopping cleanly", flush=True)
                STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
                return aggregate(args)
            index = int(e["sobol_index"])
            adopted = int(e["adopted_truth_degree"])
            degree_of, mcfg = r48._member_degree_of(design, index)
            y0 = np.asarray(mcfg["initial_state_si"], dtype=float)
            cd, rd = case_dir(design, index), raw_dir(design, index)

            p = _propagate(
                {"design": design, "sobol_index": index,
                 "policy": f"span member k={K_TARGET} (serial re-run)",
                 "purpose": "contention-free member kernel time at the "
                            "scoring tolerance"},
                degree_of, y0, adopted, "tighter",
                cd / "member_serial_tighter.json",
                rd / "member_serial_tighter.npz", timed=True)
            if p is None:
                continue
            e["member_serial_kernel_ns"] = p["telemetry"]["gravity_kernel_ns"]

            n1 = int(e["n_time_tighter_histogram"])

            def const1(tt, h_m, _n=n1):
                return _n
            p1 = _propagate(
                {"design": design, "sobol_index": index,
                 "policy": "fixed matched on mean per-call measured cost",
                 "degree": n1}, const1, y0, adopted, "tighter",
                cd / "fixed_time_tighter.json",
                rd / "fixed_time_tighter.npz", timed=True)
            if p1 is None:
                continue
            ratio = (p1["telemetry"]["gravity_kernel_ns"]
                     / e["member_serial_kernel_ns"])
            e["first_pass_time_ratio"] = ratio

            n2 = int(round(n1 * (1.0 / ratio) ** (1.0 / COST_EXPONENT)))
            n2 = max(2, min(n2, adopted))
            e["n_time_refined"] = n2

            def const2(tt, h_m, _n=n2):
                return _n
            cfg2 = {"design": design, "sobol_index": index,
                    "policy": "fixed matched on measured total kernel time "
                              "at the scoring tolerance",
                    "degree": n2, "first_pass_degree": n1}
            p2 = _propagate(cfg2, const2, y0, adopted, "tighter",
                            cd / "fixed_time2_tighter.json",
                            rd / "fixed_time2_tighter.npz", timed=True)
            if p2 is None:
                continue
            e["refined_time_ratio"] = (p2["telemetry"]["gravity_kernel_ns"]
                                       / e["member_serial_kernel_ns"])
            # The envelope needs the pair; this one is not a timed run.
            _propagate(cfg2, const2, y0, adopted, "tight",
                       cd / "fixed_time2_tight.json",
                       rd / "fixed_time2_tight.npz", timed=False)
            print(f"  [{(time.time() - t0) / 60:5.1f} min] {design}{index:03d} "
                  f"N1={n1} ratio={ratio:.2f} -> N2={n2} "
                  f"achieved={e['refined_time_ratio']:.2f}", flush=True)
    STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    return aggregate(args)


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def aggregate(args) -> int:
    if not STATE.exists():
        print("[r64] no state yet")
        return 1
    state = json.loads(STATE.read_text(encoding="utf-8"))
    rows_out = []
    for design, entries in state["designs"].items():
        span = {int(r["sobol_index"]): r for r in r48.span_rows(design)}
        for e in entries:
            index = int(e["sobol_index"])
            e_k = span[index]["entries"][K_TARGET]
            cd, rd = case_dir(design, index), raw_dir(design, index)
            rg, rt = (rd / "fixed_time2_tighter.npz",
                      rd / "fixed_time2_tight.npz")
            if not (rg.exists() and rt.exists()):
                continue
            truth = {}
            ok = True
            for lv in ("tight", "tighter"):
                _, raw = r14.reuse_paths(design, index, "truth", lv)
                if not raw.exists():
                    ok = False
                    break
                truth[lv] = _load(raw)
            if not ok:
                continue
            truth_self = base.common_error(
                truth["tight"][0], truth["tight"][1],
                truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]
            got_g, got_t = _load(rg), _load(rt)
            err = base.common_error(
                got_g[0], got_g[1], truth["tighter"][0],
                truth["tighter"][1])["pos_rms_m"]
            self_diff = base.common_error(
                got_t[0], got_t[1], got_g[0], got_g[1])["pos_rms_m"]
            thr = (e_k.get("envelope_m") or 0.0) + self_diff + truth_self
            diff = err - e_k["error_m"]
            rows_out.append({
                "design": design, "sobol_index": index, "hp_km": e["hp_km"],
                "member_error_m": e_k["error_m"],
                "comparator_degree": e.get("n_time_refined"),
                "comparator_degree_r48_tight_match": e.get(
                    "n_time_tight_histogram_r48"),
                "comparator_error_m": err,
                "achieved_time_ratio": e.get("refined_time_ratio"),
                "rho_fixed_over_member": (err / e_k["error_m"]
                                          if e_k["error_m"] else None),
                "resolution_threshold_m": thr,
                "resolved": bool(abs(diff) > thr),
                "winner": ("interior" if diff > thr else
                           ("fixed" if -diff > thr else None))})
    if not rows_out:
        print("[r64] nothing to aggregate yet")
        return 1
    res = [r for r in rows_out if r["resolved"]]
    interior = sum(1 for r in res if r["winner"] == "interior")
    fixed = sum(1 for r in res if r["winner"] == "fixed")
    ratios = [r["achieved_time_ratio"] for r in rows_out
              if r.get("achieved_time_ratio")]
    rhos = [r["rho_fixed_over_member"] for r in rows_out
            if r.get("rho_fixed_over_member")]
    summary = {"orbits": len(rows_out), "resolved": len(res),
               "interior_wins": interior, "fixed_wins": fixed,
               "unresolved": len(rows_out) - len(res),
               "median_rho": float(np.median(rhos)) if rhos else None,
               "achieved_time_ratio": {
                   "median": float(np.median(ratios)) if ratios else None,
                   "min": float(min(ratios)) if ratios else None,
                   "max": float(max(ratios)) if ratios else None}}
    OUT.write_text(json.dumps(
        {"schema": "r64_interior_timing_tighter_v1",
         "created_utc": base.utc_now(),
         "protocol": "(O48) protocol with every timed stage and the "
                     "first-pass histogram moved to the tighter level, the "
                     "level the errors are scored at",
         "summary": summary, "rows": rows_out}, indent=2), encoding="utf-8")
    print(f"[r64] written {OUT.name}: {summary['orbits']} orbits; "
          f"resolved {len(res)}: interior {interior}, fixed {fixed}, "
          f"unresolved {summary['unresolved']}; "
          f"median rho {summary['median_rho']}; "
          f"time ratio median {summary['achieved_time_ratio']['median']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("pipeline", pipeline), ("aggregate", aggregate)):
        s = sub.add_parser(name)
        if name == "pipeline":
            s.add_argument("--cutoff", default=None)
        s.set_defaults(func=fn)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
