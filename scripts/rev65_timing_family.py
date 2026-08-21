"""R65 (O58): the sampled interior family under level-consistent measured time.

(O48) and (O57) both compare one member, k = 0.5, against a constant degree.
That member was chosen because it was the most frequent sampled minimum of the
seven-day nominal sweep, not because the family has an optimum there. So
(O57)'s negative result is a statement about k = 0.5 and not about interior
allocation: nothing in it rules out a different concentration holding a
measured-time advantage.

This campaign asks the family-level question on the same fourteen orbits, the
same budget and the same level-consistent timing: each of k = 0.25, 0.50 and
0.75 against its \\emph{own} constant degree matched on measured kernel time at
the tighter tolerance, the level the errors are scored at. The k = 0.5 column
is reused from (O57); only 0.25 and 0.75 are propagated here.

Two things (O57) left open are fixed rather than repeated. Its comparator
degree came from a single refinement step and two of fourteen cells missed the
0.90-1.10 timing band; here the degree is refined until the band is met or the
integer step cannot improve it, and a cell that still misses is reported as a
timing-match miss rather than quietly kept.

All three members are reported. Selecting the favourable one after the fact
would turn a family sweep into an oracle, and the registration forbids it: no
per-orbit argmin over k is taken, and no single deployable k is claimed here.

Timing comparability requires an idle machine. The pipeline refuses to start
while other python processes are alive, but that check runs once, before the
first cell, and does not re-arm: a machine joined by other work mid-run is not
caught by it, which is how the first attempt at this campaign was lost.
probe_kernel_quiet.py guards the start and rev65_quiet_audit.py is the
after-the-fact test that does catch a mid-run change.

Usage:
    python rev65_timing_family.py plan
    python rev65_timing_family.py pipeline [--cutoff 2026-08-20T10:00]
    python rev65_timing_family.py aggregate
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
CASE_ROOT = METRICS / "r65_cases" / "timing_family"
RAW_ROOT = METRICS / "r65_raw" / "timing_family"
OUT = METRICS / "r65_timing_family.json"
STATE = METRICS / "r65_timing_family_state.json"

K_NEW = ("0.25", "0.75")
K_ALL = ("0.25", "0.50", "0.75")
BAND_LO, BAND_HI = 0.90, 1.10
MAX_REFINE = 4
COST_EXPONENT = r48.COST_EXPONENT
LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def member_sidecar(design: str, index: int, k: str, level: str) -> dict:
    p = (METRICS / "r18_cases" / f"{design}_beta_1.00_k_{k}"
         / f"sobolA_{index:03d}" / f"span_{level}.json")
    return json.loads(p.read_text(encoding="utf-8"))


def member_degree_fn(design: str, index: int, k: str):
    cfg = member_sidecar(design, index, k, "tight")["config"]
    tab = {float(a): int(b) for a, b in cfg["degree_table"].items()}
    lo, hi = min(tab), max(tab)

    def degree_of(tt, h_m, _t=tab, _lo=lo, _hi=hi):
        hb = min(_hi, max(_lo, 10.0 * math.floor(h_m / 1e4)))
        return _t[hb]
    return degree_of, cfg


def case_dir(design: str, index: int, k: str) -> Path:
    return CASE_ROOT / f"k_{k}" / design / f"sobol{design}_{index:03d}"


def raw_dir(design: str, index: int, k: str) -> Path:
    return RAW_ROOT / f"k_{k}" / design / f"sobol{design}_{index:03d}"


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
    base.atomic_json(sidecar, {
        "schema": "r65_timing_family_case_v1", "created_utc": base.utc_now(),
        "config": dict(cfg_extra, level=level, timing_comparable=timed,
                       adopted_truth_degree=adopted,
                       source=base.provenance()),
        "status": st, "telemetry": tel,
        "raw_path": str(raw.relative_to(ROOT)),
        "raw_sha256": base.file_hash(raw)})
    return json.loads(sidecar.read_text(encoding="utf-8"))


def panel() -> dict:
    """(O57)'s panel, itself (O48)'s. No orbit is selected here."""
    sel = json.loads(r48.SELECTION.read_text(encoding="utf-8"))
    deg_tab, ns_tab = r13.cost_curve()
    out = {}
    for design, entries in sel["designs"].items():
        rows = []
        for e in entries:
            index = int(e["sobol_index"])
            per_k = {}
            for k in K_NEW:
                side = member_sidecar(design, index, k, "tighter")
                counts = {int(a): int(b) for a, b in
                          side["telemetry"]["degree_counts"].items()}
                n_rhs = sum(counts.values())
                mean_cost = sum(r13.cost_of([a], deg_tab, ns_tab)[0] * b
                                for a, b in counts.items()) / n_rhs
                per_k[k] = {"n_time_first_pass":
                            int(r13.inverse_cost(mean_cost, deg_tab, ns_tab))}
            rows.append({"sobol_index": index, "hp_km": e["hp_km"],
                         "adopted_truth_degree": int(e["adopted_truth_degree"]),
                         "k": per_k})
        out[design] = rows
    return out


def run_cell(design: str, index: int, k: str, adopted: int, first: int,
             rec: dict) -> dict:
    """One (orbit, k) cell: member timing, then comparator refined to band."""
    degree_of, mcfg = member_degree_fn(design, index, k)
    y0 = np.asarray(mcfg["initial_state_si"], dtype=float)
    cd, rd = case_dir(design, index, k), raw_dir(design, index, k)

    p = _propagate(
        {"design": design, "sobol_index": index, "k": k,
         "policy": f"span member k={k} (serial re-run)",
         "purpose": "contention-free member kernel time at the scoring "
                    "tolerance"},
        degree_of, y0, adopted, "tighter",
        cd / "member_tighter.json", rd / "member_tighter.npz", timed=True)
    if p is None:
        return {"status": "member_failed"}
    t_member = p["telemetry"]["gravity_kernel_ns"]
    rec["member_kernel_ns"] = t_member

    tried, degree, best = [], int(first), None
    for step in range(MAX_REFINE):
        def const(tt, h_m, _n=degree):
            return _n
        pc = _propagate(
            {"design": design, "sobol_index": index, "k": k,
             "policy": "constant matched on measured kernel time at the "
                       "scoring tolerance",
             "degree": degree, "refine_step": step},
            const, y0, adopted, "tighter",
            cd / f"fixed_{degree}_tighter.json",
            rd / f"fixed_{degree}_tighter.npz", timed=True)
        if pc is None:
            return {"status": "comparator_failed"}
        ratio = pc["telemetry"]["gravity_kernel_ns"] / t_member
        tried.append({"degree": degree, "ratio": ratio})
        if best is None or abs(ratio - 1.0) < abs(best["ratio"] - 1.0):
            best = {"degree": degree, "ratio": ratio}
        if BAND_LO <= ratio <= BAND_HI:
            break
        nxt = int(round(degree * (1.0 / ratio) ** (1.0 / COST_EXPONENT)))
        nxt = max(2, min(nxt, adopted))
        if nxt == degree or any(t["degree"] == nxt for t in tried):
            break
        degree = nxt

    rec["refinements"] = tried
    rec["comparator_degree"] = best["degree"]
    rec["achieved_time_ratio"] = best["ratio"]
    rec["timing_match_miss"] = not (BAND_LO <= best["ratio"] <= BAND_HI)

    # The envelope needs the pair; this run is not timed.
    def constb(tt, h_m, _n=best["degree"]):
        return _n
    _propagate(
        {"design": design, "sobol_index": index, "k": k,
         "policy": "constant matched on measured kernel time at the scoring "
                   "tolerance", "degree": best["degree"], "envelope_run": True},
        constb, y0, adopted, "tight",
        cd / f"fixed_{best['degree']}_tight.json",
        rd / f"fixed_{best['degree']}_tight.npz", timed=False)
    rec["status"] = "ok"
    return rec


def plan(args) -> int:
    st = panel()
    n = sum(len(v) for v in st.values())
    print(f"[r65-plan] {n} orbits x {len(K_NEW)} new members "
          f"(k = {', '.join(K_NEW)}); k=0.50 reused from (O57)")
    for design, rows in st.items():
        for e in rows:
            fp = ", ".join(f"k={k}:N1={e['k'][k]['n_time_first_pass']}"
                           for k in K_NEW)
            print(f"  {design}{e['sobol_index']:03d} hp={e['hp_km']:6.1f}  {fp}")
    return 0


def pipeline(args) -> int:
    cutoff = datetime.fromisoformat(args.cutoff) if args.cutoff else None
    state = (json.loads(STATE.read_text(encoding="utf-8"))
             if STATE.exists() else {"designs": panel()})
    if not r48.idle_or_die():
        return 2
    t0 = time.time()
    for design, rows in state["designs"].items():
        for e in rows:
            for k in K_NEW:
                if cutoff and datetime.now() > cutoff:
                    print("[r65] cutoff reached; stopping cleanly", flush=True)
                    STATE.write_text(json.dumps(state, indent=2),
                                     encoding="utf-8")
                    return aggregate(args)
                rec = e["k"][k]
                if rec.get("status") == "ok":
                    continue
                out = run_cell(design, int(e["sobol_index"]), k,
                               int(e["adopted_truth_degree"]),
                               int(rec["n_time_first_pass"]), rec)
                if out.get("status") != "ok":
                    print(f"  !! {design}{e['sobol_index']:03d} k={k}: "
                          f"{out.get('status')}", flush=True)
                    continue
                miss = " MISS" if rec["timing_match_miss"] else ""
                print(f"  [{(time.time()-t0)/60:5.1f} min] "
                      f"{design}{e['sobol_index']:03d} k={k} "
                      f"N={rec['comparator_degree']} "
                      f"ratio={rec['achieved_time_ratio']:.2f} "
                      f"({len(rec['refinements'])} passes){miss}", flush=True)
                STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
    return aggregate(args)


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def _score(design, index, k, comparator_degree, cd, rd):
    span = {int(r["sobol_index"]): r for r in r48.span_rows(design)}
    e_k = span[index]["entries"][k]
    truth = {}
    for lv in ("tight", "tighter"):
        _, raw = r14.reuse_paths(design, index, "truth", lv)
        if not raw.exists():
            return None
        truth[lv] = _load(raw)
    s_ref = base.common_error(truth["tight"][0], truth["tight"][1],
                              truth["tighter"][0],
                              truth["tighter"][1])["pos_rms_m"]
    rg = rd / f"fixed_{comparator_degree}_tighter.npz"
    rt = rd / f"fixed_{comparator_degree}_tight.npz"
    if not (rg.exists() and rt.exists()):
        return None
    got_g, got_t = _load(rg), _load(rt)
    err = base.common_error(got_g[0], got_g[1], truth["tighter"][0],
                            truth["tighter"][1])["pos_rms_m"]
    s_fix = base.common_error(got_t[0], got_t[1], got_g[0],
                              got_g[1])["pos_rms_m"]
    thr = (e_k.get("envelope_m") or 0.0) + s_fix + s_ref
    diff = err - e_k["error_m"]
    return {"member_error_m": e_k["error_m"], "comparator_error_m": err,
            "resolution_threshold_m": thr,
            "M_res": abs(diff) / thr if thr else None,
            "resolved": bool(abs(diff) > thr),
            "winner": ("interior" if diff > thr else
                       ("fixed" if -diff > thr else None)),
            "rho_fixed_over_member": (err / e_k["error_m"]
                                      if e_k["error_m"] else None)}


def aggregate(args) -> int:
    rows = []
    if STATE.exists():
        state = json.loads(STATE.read_text(encoding="utf-8"))
        for design, entries in state["designs"].items():
            for e in entries:
                idx = int(e["sobol_index"])
                for k in K_NEW:
                    rec = e["k"][k]
                    if rec.get("status") != "ok":
                        continue
                    sc = _score(design, idx, k, rec["comparator_degree"],
                                case_dir(design, idx, k),
                                raw_dir(design, idx, k))
                    if sc is None:
                        continue
                    rows.append(dict(
                        {"design": design, "sobol_index": idx, "k": k,
                         "hp_km": e["hp_km"],
                         "comparator_degree": rec["comparator_degree"],
                         "achieved_time_ratio": rec["achieved_time_ratio"],
                         "timing_match_miss": rec["timing_match_miss"],
                         "refinement_passes": len(rec["refinements"])}, **sc))
    # k = 0.50 reused from (O57), reported in the same table
    p57 = METRICS / "r64_interior_timing_tighter.json"
    if p57.exists():
        for r in json.loads(p57.read_text(encoding="utf-8"))["rows"]:
            thr = r["resolution_threshold_m"]
            gap = abs(r["comparator_error_m"] - r["member_error_m"])
            rows.append({
                "design": r["design"], "sobol_index": r["sobol_index"],
                "k": "0.50", "hp_km": r["hp_km"],
                "comparator_degree": r["comparator_degree"],
                "achieved_time_ratio": r["achieved_time_ratio"],
                "timing_match_miss": not (BAND_LO <= r["achieved_time_ratio"]
                                          <= BAND_HI),
                "refinement_passes": 1,
                "member_error_m": r["member_error_m"],
                "comparator_error_m": r["comparator_error_m"],
                "resolution_threshold_m": thr,
                "M_res": gap / thr if thr else None,
                "resolved": r["resolved"], "winner": r["winner"],
                "rho_fixed_over_member": r["rho_fixed_over_member"],
                "source": "reused from (O57)"})
    if not rows:
        print("[r65] nothing to aggregate yet")
        return 1
    by_k = {}
    for k in K_ALL:
        sel = [r for r in rows if r["k"] == k]
        if not sel:
            continue
        res = [r for r in sel if r["resolved"]]
        rhos = [r["rho_fixed_over_member"] for r in sel
                if r["rho_fixed_over_member"]]
        by_k[k] = {"orbits": len(sel), "resolved": len(res),
                   "interior_wins": sum(1 for r in res
                                        if r["winner"] == "interior"),
                   "fixed_wins": sum(1 for r in res
                                     if r["winner"] == "fixed"),
                   "unresolved": len(sel) - len(res),
                   "median_rho": float(np.median(rhos)) if rhos else None,
                   "timing_match_misses": sum(1 for r in sel
                                              if r["timing_match_miss"])}
    OUT.write_text(json.dumps(
        {"schema": "r65_timing_family_v1", "created_utc": base.utc_now(),
         "protocol": "each sampled interior member against its own constant "
                     "degree matched on measured kernel time at the tighter "
                     "tolerance, the level the errors are scored at; all "
                     "three members reported, no per-orbit argmin taken",
         "timing_band": [BAND_LO, BAND_HI], "by_k": by_k, "rows": rows},
        indent=2), encoding="utf-8")
    print(f"[r65] written {OUT.name}")
    for k, s in by_k.items():
        print(f"  k={k}: resolved {s['resolved']}/{s['orbits']} "
              f"(interior {s['interior_wins']}, fixed {s['fixed_wins']}), "
              f"median rho {s['median_rho']:.3f}, "
              f"timing misses {s['timing_match_misses']}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("plan", plan), ("pipeline", pipeline),
                     ("aggregate", aggregate)):
        s = sub.add_parser(name)
        if name == "pipeline":
            s.add_argument("--cutoff", default=None)
        s.set_defaults(func=fn)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
