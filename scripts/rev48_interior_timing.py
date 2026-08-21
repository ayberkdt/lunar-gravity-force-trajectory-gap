"""O48/R48: measured-kernel-time-matched comparator for the INTERIOR member.

The R13 panel established the measured-time comparison for the radial
endpoint; the interior member (k = 0.5) has only ever been compared under
operation-count proxies, and the paper itself records that the quadratic
proxy's low-degree flattening biases those proxies in the interior member's
favor. This campaign closes that gap with the same protocol as R13, sources
swapped to the R18 span-sweep member at beta = 1:

  1. select   7 orbits per design spread over perilune (extremes retained),
              N_time from the call-weighted mean per-call measured cost of the
              member's archived degree histogram, inverted on the measured
              cost curve (r12_kernel_cost_curve.json);
  2. serial   member schedule re-run serially at tight for a contention-free
              kernel time (the archived campaign ran with concurrent workers);
  3. first    comparator at N_time, serial, tight;
  4. refine   second-pass degree from the measured total-kernel-time ratio,
              c(N) ~ N^2 locally, capped at the orbit's adopted truth degree;
  5. second   refined comparator at tight (serial, timed) and tighter (for
              the envelope);
  6. aggregate errors at the tighter level against the archived truths under
              the campaign envelope rule; achieved time ratios reported.

Timing comparability requires an idle machine: a pipeline invocation refuses
to start while other python processes are alive. The check runs once, at the
head of the invocation, and does not re-arm, so a machine joined by other work
mid-run is caught by the after-the-fact quiet audit rather than by this guard.

Usage:
    python rev48_interior_timing.py pipeline [--cutoff 2026-08-10T06:45]
    python rev48_interior_timing.py aggregate
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

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SELECTION = METRICS / "r48_interior_timing_selection.json"
OUTPUT = METRICS / "r48_interior_timing.json"
CASE_ROOT = METRICS / "r48_cases" / "interior_timing"
RAW_ROOT = METRICS / "r48_raw" / "interior_timing"

K_TARGET = "0.50"
BETA = 1.0
N_PER_DESIGN = 7
COST_EXPONENT = 2.0
LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def member_sidecar(design: str, index: int, level: str) -> dict:
    p = (METRICS / "r18_cases" / f"{design}_beta_1.00_k_{K_TARGET}"
         / f"sobolA_{index:03d}" / f"span_{level}.json")
    return json.loads(p.read_text(encoding="utf-8"))


def span_rows(design: str) -> list:
    d = json.loads((METRICS / f"r18_span_sweep_{design}_beta_1.00.json"
                    ).read_text(encoding="utf-8"))
    return d["rows"]


def case_dir(design: str, index: int) -> Path:
    return CASE_ROOT / design / f"sobol{design}_{index:03d}"


def raw_dir(design: str, index: int) -> Path:
    return RAW_ROOT / design / f"sobol{design}_{index:03d}"


def idle_or_die() -> bool:
    others = base.other_python_processes()
    if others:
        print(f"!! {len(others)} other python processes running; timing "
              f"stages need an idle machine — aborting this stage")
        return False
    return True


def select() -> int:
    if SELECTION.exists():
        print(f"[r48] selection exists, keeping {SELECTION.name}")
        return 0
    deg_tab, ns_tab = r13.cost_curve()
    out = {"schema": "r48_interior_timing_selection_v1",
           "created_utc": base.utc_now(),
           "rule": (f"{N_PER_DESIGN} orbits per design spread over perilune "
                    "with the extremes retained, the R13 panel rule applied "
                    "to the R18 beta=1 population"),
           "cost_curve": "r12_kernel_cost_curve.json", "designs": {}}
    for design in ("A", "B"):
        rows = [r for r in span_rows(design)
                if r["entries"].get(K_TARGET, {}).get("error_m") is not None]
        rows.sort(key=lambda r: r["hp_km"])
        pick = [rows[int(i)] for i in
                np.linspace(0, len(rows) - 1, N_PER_DESIGN).round()]
        entries = []
        for r in pick:
            index = int(r["sobol_index"])
            side = member_sidecar(design, index, "tight")
            counts = {int(k): int(v) for k, v in
                      side["telemetry"]["degree_counts"].items()}
            n_rhs = sum(counts.values())
            mean_cost = sum(r13.cost_of([k], deg_tab, ns_tab)[0] * v
                            for k, v in counts.items()) / n_rhs
            n_time = r13.inverse_cost(mean_cost, deg_tab, ns_tab)
            entries.append({
                "sobol_index": index, "hp_km": r["hp_km"],
                "adopted_truth_degree":
                    int(side["config"]["adopted_truth_degree"]),
                "n_work_proxy": int(round(math.sqrt(
                    side["telemetry"]["mean_degree_sq"]))),
                "n_time_measured_cost": int(n_time),
                "member_mean_per_call_ns": float(mean_cost)})
            print(f"  {design}{index:03d} hp={r['hp_km']:6.1f} km "
                  f"N_work={entries[-1]['n_work_proxy']:3d} -> "
                  f"N_time={n_time:3d}")
        out["designs"][design] = entries
    SELECTION.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[written] {SELECTION.name}")
    return 0


def _propagate(cfg_extra: dict, degree_of, y0, adopted: int, level: str,
               sidecar: Path, raw: Path, timed: bool) -> dict | None:
    if sidecar.exists() and (not timed or raw.exists()):
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
    payload = {"schema": "r48_interior_timing_case_v1",
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


def _member_degree_of(design: str, index: int):
    cfg = member_sidecar(design, index, "tight")["config"]
    tab = {float(k): int(v) for k, v in cfg["degree_table"].items()}
    lo, hi = min(tab), max(tab)

    def degree_of(tt, h_m, _tab=tab, _lo=lo, _hi=hi):
        hb = min(_hi, max(_lo, 10.0 * math.floor(h_m / 1e4)))
        return _tab[hb]
    return degree_of, cfg


def pipeline(args) -> int:
    cutoff = (datetime.fromisoformat(args.cutoff)
              if args.cutoff else None)
    select()
    sel = json.loads(SELECTION.read_text(encoding="utf-8"))
    deg_tab, ns_tab = r13.cost_curve()
    if not idle_or_die():
        return 2
    t0 = time.time()
    for design, entries in sel["designs"].items():
        for e in entries:
            if cutoff and datetime.now() > cutoff:
                print(f"[r48] cutoff reached; stopping cleanly")
                SELECTION.write_text(json.dumps(sel, indent=2),
                                     encoding="utf-8")
                return aggregate(args)
            index = int(e["sobol_index"])
            adopted = int(e["adopted_truth_degree"])
            degree_of, mcfg = _member_degree_of(design, index)
            y0 = np.asarray(mcfg["initial_state_si"], dtype=float)
            cd, rd = case_dir(design, index), raw_dir(design, index)

            # 2. member serial baseline (tight, timed)
            p = _propagate(
                {"design": design, "sobol_index": index,
                 "policy": f"span member k={K_TARGET} (serial re-run)",
                 "purpose": "contention-free member kernel time"},
                degree_of, y0, adopted, "tight",
                cd / "member_serial_tight.json",
                rd / "member_serial_tight.npz", timed=True)
            if p is None:
                continue
            e["member_serial_kernel_ns"] = p["telemetry"]["gravity_kernel_ns"]

            # 3. first-pass comparator (tight, timed)
            n1 = int(e["n_time_measured_cost"])

            def const1(tt, h_m, _n=n1):
                return _n
            p1 = _propagate(
                {"design": design, "sobol_index": index,
                 "policy": "fixed matched on mean per-call measured cost",
                 "degree": n1}, const1, y0, adopted, "tight",
                cd / "fixed_time_tight.json",
                rd / "fixed_time_tight.npz", timed=True)
            if p1 is None:
                continue
            ratio = (p1["telemetry"]["gravity_kernel_ns"]
                     / e["member_serial_kernel_ns"])
            e["first_pass_time_ratio"] = ratio

            # 4. refine on the measured total-time ratio
            n2 = int(round(n1 * (1.0 / ratio) ** (1.0 / COST_EXPONENT)))
            n2 = max(2, min(n2, adopted))
            e["n_time_refined"] = n2

            # 5. refined comparator: tight (timed) + tighter (envelope)
            def const2(tt, h_m, _n=n2):
                return _n
            p2 = _propagate(
                {"design": design, "sobol_index": index,
                 "policy": "fixed matched on measured total kernel time",
                 "degree": n2, "first_pass_degree": n1},
                const2, y0, adopted, "tight",
                cd / "fixed_time2_tight.json",
                rd / "fixed_time2_tight.npz", timed=True)
            if p2 is None:
                continue
            e["refined_time_ratio"] = (p2["telemetry"]["gravity_kernel_ns"]
                                       / e["member_serial_kernel_ns"])
            _propagate(
                {"design": design, "sobol_index": index,
                 "policy": "fixed matched on measured total kernel time",
                 "degree": n2, "first_pass_degree": n1},
                const2, y0, adopted, "tighter",
                cd / "fixed_time2_tighter.json",
                rd / "fixed_time2_tighter.npz", timed=False)
            print(f"  [{(time.time() - t0) / 60:5.1f} min] {design}{index:03d} "
                  f"N1={n1} ratio={ratio:.2f} -> N2={n2} "
                  f"achieved={e['refined_time_ratio']:.2f}", flush=True)
    SELECTION.write_text(json.dumps(sel, indent=2), encoding="utf-8")
    return aggregate(args)


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def aggregate(args) -> int:
    sel = json.loads(SELECTION.read_text(encoding="utf-8"))
    rows_out = []
    for design, entries in sel["designs"].items():
        span = {int(r["sobol_index"]): r for r in span_rows(design)}
        for e in entries:
            index = int(e["sobol_index"])
            e_k = span[index]["entries"][K_TARGET]
            cd, rd = case_dir(design, index), raw_dir(design, index)
            s2t = cd / "fixed_time2_tight.json"
            r2t, r2g = (rd / "fixed_time2_tight.npz",
                        rd / "fixed_time2_tighter.npz")
            if not (s2t.exists() and r2t.exists() and r2g.exists()):
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
            got_t, got_g = _load(r2t), _load(r2g)
            err = base.common_error(
                got_g[0], got_g[1], truth["tighter"][0],
                truth["tighter"][1])["pos_rms_m"]
            self_diff = base.common_error(
                got_t[0], got_t[1], got_g[0], got_g[1])["pos_rms_m"]
            env_k = e_k.get("envelope_m") or 0.0
            thr = env_k + self_diff + truth_self
            diff = err - e_k["error_m"]
            rows_out.append({
                "design": design, "sobol_index": index,
                "hp_km": e["hp_km"],
                "member_error_m": e_k["error_m"],
                "comparator_degree": e.get("n_time_refined"),
                "comparator_error_m": err,
                "achieved_time_ratio": e.get("refined_time_ratio"),
                "rho_fixed_over_member": (err / e_k["error_m"]
                                          if e_k["error_m"] else None),
                "resolution_threshold_m": thr,
                "resolved": bool(abs(diff) > thr),
                "winner": ("interior" if diff > thr else
                           ("fixed" if -diff > thr else None))})
    if not rows_out:
        print("[r48] nothing to aggregate yet")
        return 1
    res = [r for r in rows_out if r["resolved"]]
    wins = sum(1 for r in res if r["winner"] == "interior")
    ratios = [r["achieved_time_ratio"] for r in rows_out
              if r["achieved_time_ratio"]]
    payload = {
        "schema": "r48_interior_timing_v1", "created_utc": base.utc_now(),
        "protocol": "R13 measured-time panel protocol, interior member",
        "summary": {
            "orbits": len(rows_out), "resolved": len(res),
            "interior_wins": wins, "fixed_wins": len(res) - wins,
            "unresolved": len(rows_out) - len(res),
            "achieved_time_ratio": {
                "median": float(np.median(ratios)),
                "min": float(np.min(ratios)),
                "max": float(np.max(ratios))} if ratios else None},
        "rows": rows_out}
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    s = payload["summary"]
    print(f"[r48] written {OUTPUT.name}: {s['orbits']} orbits, "
          f"resolved {s['resolved']} (interior {s['interior_wins']}, "
          f"fixed {s['fixed_wins']}), unresolved {s['unresolved']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("pipeline")
    pl.add_argument("--cutoff", default=None,
                    help="local ISO datetime; stop starting new orbits after")
    pl.set_defaults(func=pipeline)
    ag = sub.add_parser("aggregate")
    ag.set_defaults(func=aggregate)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
