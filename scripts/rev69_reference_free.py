"""R69 (O61): the interior candidate calibrated without a reference trajectory.

Why this exists
---------------
(O60) shows the exploratory k = 0.5 candidate keeping its advantage at equal
measured kernel time on both full designs. Its budget, though, is still set the
way every policy in this paper sets it: by bisecting on the altitude history of
an already-propagated high-fidelity reference arc. The manuscript says so in
its limitations, and a practitioner's next question follows from it. If the
reference propagation has to be run first, when is the schedule usable?

This campaign removes that foreknowledge. Nothing about the policy family
changes; only where the two calibrations read their altitude history:

    archived (O60):   h(t) from the propagated N = 300/600 reference arc
    here    (O61):    h(E) from the two-body orbit of the initial osculating
                      elements, r(E) = a0 (1 - e0 cos E), sampled on the same
                      uniform time grid through Kepler's equation

Both calibrations are re-run on that predicted history: the radial rule's
tolerance eps_A, which fixes the k = 1 endpoint's degree table, and the family
scale s_k, which fixes the k = 0.5 member's. The resulting schedule is frozen
before the arc is propagated and reads nothing from any reference trajectory.

What stays reference-free, and what does not
--------------------------------------------
N_crit is the empirical minimum degree at the orbit's own perilune altitude, a
design parameter, so it needs no propagation. The constant comparator degree
N_0 follows from N_crit and beta. The reference trajectory is still used for
one thing, scoring: errors are read against the same archived reference as
every other campaign, because the question here is whether the schedule can be
*built* without foreknowledge, not whether it can be *scored* without a truth.

Two comparisons are reported. The member's realized budget against the target
it was calibrated to, which is what the Kepler prediction can get wrong; and
the member against a constant degree matched on measured kernel time under the
(O60) protocol, which is what tells us whether the ordering survives.

Usage:
    python rev69_reference_free.py plan
    python rev69_reference_free.py pipeline [--cutoff ISO8601]
    python rev69_reference_free.py aggregate
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rev10_sobol_confirmatory as base            # noqa: E402
import rev13_timing_match as r13                   # noqa: E402
import rev14_budget_pareto as pareto               # noqa: E402
import rev14_budget_trajectory as r14              # noqa: E402
import rev15_deployable_calibration as r15         # noqa: E402
import rev18_span_sweep as r18                     # noqa: E402
import rev48_interior_timing as r48                # noqa: E402
import rev68_timing_full as r68                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r69_cases"
RAW_ROOT = METRICS / "r69_raw"
OUT = METRICS / "r69_reference_free.json"
STATE = METRICS / "r69_reference_free_state.json"

K_TARGET = 0.5
BETA = 1.0
BETA_TAG = "beta_1.00"
DESIGNS = ("A", "B")
BAND_LO, BAND_HI = r68.BAND_LO, r68.BAND_HI
MAX_REFINE = r68.MAX_REFINE
BUDGET_MISS_GATE = 0.05
DISK_FLOOR_GB = 6.0

LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


# ------------------------------------------------------- Kepler prediction

def osculating(y0: np.ndarray, mu: float) -> dict:
    r = np.asarray(y0[:3], dtype=float)
    v = np.asarray(y0[3:], dtype=float)
    rn, vn = float(np.linalg.norm(r)), float(np.linalg.norm(v))
    a = 1.0 / (2.0 / rn - vn * vn / mu)
    evec = ((vn * vn - mu / rn) * r - float(np.dot(r, v)) * v) / mu
    e = float(np.linalg.norm(evec))
    return {"a_m": float(a), "e": e,
            "period_s": 2.0 * math.pi * math.sqrt(a ** 3 / mu),
            "rp_m": float(a * (1.0 - e)), "ra_m": float(a * (1.0 + e))}


# The two-body altitude history comes from R15, which already solves it
# for the deployable-calibration route of Section S6. Two Kepler solvers
# in one package is one more than a reader should have to reconcile; the
# one dropped here agreed with R15's to 1e-11 km on all 128 orbits, so
# adopting R15's leaves every frozen schedule unchanged.


# ------------------------------------------------------- schedule building

def archived_rows(design: str) -> list:
    p = METRICS / f"r14_trajectory_{design}_{BETA_TAG}.json"
    return json.loads(p.read_text(encoding="utf-8"))["rows"]


def span_rows(design: str) -> dict:
    p = METRICS / f"r18_span_sweep_{design}_{BETA_TAG}.json"
    return {int(r["sobol_index"]): r
            for r in json.loads(p.read_text(encoding="utf-8"))["rows"]}


def build_schedule(row: dict) -> dict:
    """Calibrate both stages on the Kepler-predicted altitude history."""
    adopted = int(row["adopted_truth_degree"])
    n_crit = int(row["n_critical"])
    n0 = int(row["fixed_degree"])
    hp_km = float(row["design_point"]["hp_km"])
    ha_km = float(row["design_point"]["ha_km"])
    # The trajectory record carries the design point but not the state; the
    # state lives in the archived case sidecar, which is where R68 reads it.
    side = r68.member_sidecar("endpoint", row["design"],
                              int(row["sobol_index"]), "tighter")
    y0 = np.asarray(side["config"]["initial_state_si"], dtype=float)

    model, _ = r14._model(adopted)
    g = r14._g(adopted)
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    h_km = np.asarray(r15.kepler_altitudes(model, y0, grid), dtype=float)
    target = BETA * n_crit ** 2

    tol = pareto.calibrate_tolerance(model, g, hp_km, ha_km, adopted,
                                     h_km, target)
    table_a = {float(a): int(b) for a, b in tol["table"].items()}
    cal = r18.calibrate_scale(table_a, n0, K_TARGET, adopted, h_km, target)
    table_k = cal["table"]

    el = osculating(y0, model.mu)
    return {
        "adopted_truth_degree": adopted, "n_critical": n_crit,
        "constant_degree_endpoint": n0,
        "work_target": target,
        "kepler_elements": el,
        "kepler_alt_km": {"min": float(h_km.min()), "max": float(h_km.max()),
                          "samples": int(h_km.size)},
        "atallah_tol_accel_m_s2": float(tol["tol"]),
        "atallah_tol_attainable": bool(tol["attainable"]),
        "scale": float(cal["scale"]),
        "work_predicted": float(cal["work"]),
        "work_predicted_mismatch": float(cal["work"] / target - 1.0),
        "work_attainable": bool(cal["attainable"]),
        "degree_table": {str(a): int(b) for a, b in table_k.items()},
        "degree_span": max(table_k.values()) / max(1, min(table_k.values())),
        "initial_state_si": [float(v) for v in y0],
    }


def degree_fn(table: dict):
    tab = {float(a): int(b) for a, b in table.items()}
    lo, hi = min(tab), max(tab)

    def degree_of(tt, h_m, _t=tab, _lo=lo, _hi=hi):
        hb = min(_hi, max(_lo, 10.0 * math.floor(h_m / 1e4)))
        return _t[hb]
    return degree_of


# ------------------------------------------------------------- propagation

def case_dir(design: str, index: int) -> Path:
    return CASE_ROOT / design / f"sobol{design}_{index:03d}"


def raw_dir(design: str, index: int) -> Path:
    return RAW_ROOT / design / f"sobol{design}_{index:03d}"


def _propagate(cfg_extra, degree_of, y0, adopted, level, sidecar, raw, timed):
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
        "schema": "r69_reference_free_case_v1", "created_utc": base.utc_now(),
        "config": dict(cfg_extra, level=level, timing_comparable=timed,
                       adopted_truth_degree=adopted,
                       source=base.provenance()),
        "status": st, "event": ev, "telemetry": tel,
        "raw_path": str(raw.relative_to(ROOT)),
        "raw_sha256": base.file_hash(raw)})
    return json.loads(sidecar.read_text(encoding="utf-8"))


def realized_budget(design: str, index: int, sched: dict) -> dict | None:
    """What the frozen schedule actually spent, on the arc it actually flew."""
    raw = raw_dir(design, index) / "member_tighter.npz"
    if not raw.exists():
        return None
    d = np.load(raw)
    # state_si is stored components-by-epoch, (6, n), not (n, 6).
    y = d["state_si"]
    model, _ = r14._model(int(sched["adopted_truth_degree"]))
    h_km = (np.linalg.norm(y[:3, :], axis=0) - model.r_ref) / 1000.0
    tab = {float(a): int(b) for a, b in sched["degree_table"].items()}
    deg = r18.r14_degrees_from_table(tab, h_km)
    sampled = float(np.mean(deg.astype(float) ** 2))
    target = sched["work_target"]
    return {"work_sampled_on_flown_arc": sampled,
            "budget_miss_sampled": sampled / target - 1.0,
            "flown_alt_km": {"min": float(h_km.min()),
                             "max": float(h_km.max())}}


def run_cell(design: str, e: dict) -> dict:
    index = int(e["sobol_index"])
    sched = e["schedule"]
    adopted = int(sched["adopted_truth_degree"])
    y0 = np.asarray(sched["initial_state_si"], dtype=float)
    cd, rd = case_dir(design, index), raw_dir(design, index)
    rec = {"refinements": []}

    p = _propagate(
        {"design": design, "sobol_index": index, "k": K_TARGET,
         "policy": "k=0.5 member calibrated on the Kepler-predicted arc",
         "scale": sched["scale"],
         "atallah_tol_accel_m_s2": sched["atallah_tol_accel_m_s2"],
         "degree_table": sched["degree_table"],
         "purpose": "reference-free budget calibration, serial for timing"},
        degree_fn(sched["degree_table"]), y0, adopted, "tighter",
        cd / "member_tighter.json", rd / "member_tighter.npz", timed=True)
    if p is None:
        return {"status": "member_failed"}
    t_member = p["telemetry"]["gravity_kernel_ns"]
    rec["member_kernel_ns"] = t_member
    rec["member_n_rhs"] = p["telemetry"]["n_rhs"]
    rec["member_mean_degree_sq"] = p["telemetry"]["mean_degree_sq"]
    rec["budget_miss_calls"] = (p["telemetry"]["mean_degree_sq"]
                                / sched["work_target"] - 1.0)
    rec.update(realized_budget(design, index, sched) or {})

    deg_tab, ns_tab = r13.cost_curve()
    degree, rule = r68.first_guess(t_member, e.get("archived_constant_n_rhs"),
                                   p["telemetry"]["n_rhs"], adopted,
                                   deg_tab, ns_tab)
    rec["first_guess_rule"] = rule
    tried, best = [], None
    for step in range(MAX_REFINE):
        def const(tt, h_m, _n=degree):
            return _n
        pc = _propagate(
            {"design": design, "sobol_index": index,
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
        nxt = int(round(degree * (1.0 / ratio) ** (1.0 / r68._exponent(tried))))
        nxt = max(2, min(nxt, adopted))
        if nxt == degree or any(t["degree"] == nxt for t in tried):
            break
        degree = nxt

    rec["refinements"] = tried
    rec["comparator_degree"] = best["degree"]
    rec["achieved_time_ratio"] = best["ratio"]
    rec["timing_match_miss"] = not (BAND_LO <= best["ratio"] <= BAND_HI)
    rec["comparator_at_ceiling"] = bool(best["degree"] >= adopted)

    def constb(tt, h_m, _n=best["degree"]):
        return _n
    _propagate(
        {"design": design, "sobol_index": index,
         "policy": "constant matched on measured kernel time at the scoring "
                   "tolerance", "degree": best["degree"], "envelope_run": True},
        constb, y0, adopted, "tight",
        cd / f"fixed_{best['degree']}_tight.json",
        rd / f"fixed_{best['degree']}_tight.npz", timed=False)
    rec["status"] = "ok"
    return rec


# ------------------------------------------------------------- scoring

def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def score_cell(design: str, e: dict) -> dict | None:
    rec = e["cell"]
    if rec.get("status") != "ok":
        return None
    index = int(e["sobol_index"])
    truth = {}
    for lv in ("tight", "tighter"):
        _, raw = r14.reuse_paths(design, index, "truth", lv)
        if not raw.exists():
            return None
        truth[lv] = _load(raw)
    s_ref = base.common_error(truth["tight"][0], truth["tight"][1],
                              truth["tighter"][0],
                              truth["tighter"][1])["pos_rms_m"]
    rd = raw_dir(design, index)
    mem = rd / "member_tighter.npz"
    n = rec["comparator_degree"]
    rg, rt = rd / f"fixed_{n}_tighter.npz", rd / f"fixed_{n}_tight.npz"
    if not (mem.exists() and rg.exists() and rt.exists()):
        return None
    got_m, got_g, got_t = _load(mem), _load(rg), _load(rt)
    err_m = base.common_error(got_m[0], got_m[1], truth["tighter"][0],
                              truth["tighter"][1])["pos_rms_m"]
    err_f = base.common_error(got_g[0], got_g[1], truth["tighter"][0],
                              truth["tighter"][1])["pos_rms_m"]
    s_fix = base.common_error(got_t[0], got_t[1], got_g[0],
                              got_g[1])["pos_rms_m"]
    # The member has no second level here, so its envelope is taken from the
    # archived member's, the same policy family at the same budget. It is
    # recorded as reused rather than measured.
    env_member = e.get("archived_member_envelope_m") or 0.0
    thr = env_member + s_fix + s_ref
    diff = err_f - err_m
    return {
        "sobol_index": index, "design": design,
        "hp_km": e["hp_km"], "ha_km": e.get("ha_km"),
        "n_critical": e["schedule"]["n_critical"],
        "comparator_degree": n,
        "achieved_time_ratio": rec["achieved_time_ratio"],
        "timing_match_miss": rec["timing_match_miss"],
        "comparator_at_ceiling": rec["comparator_at_ceiling"],
        "scale": e["schedule"]["scale"],
        "archived_scale": e.get("archived_scale"),
        "budget_miss_predicted": e["schedule"]["work_predicted_mismatch"],
        "budget_miss_sampled": rec.get("budget_miss_sampled"),
        "budget_miss_calls": rec.get("budget_miss_calls"),
        "member_error_m": err_m,
        "archived_member_error_m": e.get("archived_member_error_m"),
        "comparator_error_m": err_f,
        "resolution_threshold_m": thr,
        "envelope_source": "archived k=0.5 member at the same budget",
        "resolved": bool(abs(diff) > thr),
        "winner": ("interior" if diff > thr else
                   ("constant" if -diff > thr else None)),
        "rho_constant_over_member": err_f / err_m if err_m else None}


def panel() -> dict:
    out = {}
    for design in DESIGNS:
        span = span_rows(design)
        rows = []
        for row in archived_rows(design):
            index = int(row["sobol_index"])
            sr = span.get(index, {})
            entry = sr.get("entries", {}).get("0.50", {})
            rows.append({
                "sobol_index": index,
                "hp_km": row["design_point"]["hp_km"],
                "ha_km": row["design_point"]["ha_km"],
                "archived_constant_n_rhs":
                    r68.archived_constant_rhs(design, index),
                "archived_member_error_m": entry.get("error_m"),
                "archived_member_envelope_m": entry.get("envelope_m"),
                "archived_scale": entry.get("scale"),
                "schedule": build_schedule(row),
                "cell": {"status": "pending"}})
        out[design] = rows
    return out


def aggregate(args) -> int:
    st = json.loads(STATE.read_text(encoding="utf-8"))
    by_design, rows = {}, []
    for design, entries in st["designs"].items():
        drows = [r for r in (score_cell(design, e) for e in entries)
                 if r is not None]
        rows.extend(drows)
        scored = [r for r in drows if not r["comparator_at_ceiling"]]
        res = [r for r in scored if r["resolved"]]
        miss = sorted(abs(r["budget_miss_sampled"]) for r in scored
                      if r["budget_miss_sampled"] is not None)
        rho = sorted(r["rho_constant_over_member"] for r in scored
                     if r["rho_constant_over_member"])
        by_design[design] = {
            "orbits_in_design": len(entries),
            "cells_scored": len(drows),
            "ceiling_censored": sum(1 for r in drows
                                    if r["comparator_at_ceiling"]),
            "timing_match_misses": sum(1 for r in scored
                                       if r["timing_match_miss"]),
            "resolved": len(res),
            "unresolved": len(scored) - len(res),
            "resolved_interior_wins":
                sum(1 for r in res if r["winner"] == "interior"),
            "resolved_constant_wins":
                sum(1 for r in res if r["winner"] == "constant"),
            "median_abs_budget_miss_sampled":
                (float(np.median(miss)) if miss else None),
            "max_abs_budget_miss_sampled": (miss[-1] if miss else None),
            "median_rho_constant_over_member":
                (float(np.median(rho)) if rho else None),
            "median_achieved_time_ratio":
                float(np.median([r["achieved_time_ratio"] for r in scored]))
                if scored else None}
    majority = all(d["resolved_interior_wins"] > d["resolved_constant_wins"]
                   for d in by_design.values())
    budget_ok = all((d["median_abs_budget_miss_sampled"] or 1.0)
                    <= BUDGET_MISS_GATE for d in by_design.values())
    lost = [k for k, d in by_design.items()
            if d["resolved_interior_wins"] <= d["resolved_constant_wins"]]
    outcome = ("A" if majority and budget_ok else
               "B" if majority else
               "C" if len(lost) == 1 else "D")
    payload = {
        "schema": "r69_reference_free_v1", "created_utc": base.utc_now(),
        "k": K_TARGET, "beta": BETA,
        "protocol": (
            "the radial tolerance and the family scale are both bisected on "
            "the two-body altitude history of the initial osculating "
            "elements, not on a propagated reference arc; the frozen schedule "
            "is then flown and compared with a constant degree matched on "
            "measured kernel time at the tighter level, the (O60) protocol "
            "unchanged"),
        "timing_band": [BAND_LO, BAND_HI],
        "budget_gate": BUDGET_MISS_GATE,
        "registered_outcome": outcome,
        "by_design": by_design, "rows": rows, "source": base.provenance()}
    base.atomic_json(OUT, payload)
    print(f"[written] {OUT.name}  registered outcome {outcome}")
    for design, s in by_design.items():
        print(f"  {design}: resolved {s['resolved']} "
              f"({s['resolved_interior_wins']} interior / "
              f"{s['resolved_constant_wins']} constant), "
              f"median |budget miss| "
              f"{s['median_abs_budget_miss_sampled']}, "
              f"median rho {s['median_rho_constant_over_member']}, "
              f"misses {s['timing_match_misses']}")
    return 0


def plan(args) -> int:
    st = panel()
    for design, rows in st.items():
        pred = [abs(r["schedule"]["work_predicted_mismatch"]) for r in rows]
        scales = [r["schedule"]["scale"] for r in rows]
        arch = [r["archived_scale"] for r in rows if r["archived_scale"]]
        print(f"[r69-plan] {design}: {len(rows)} orbits, "
              f"predicted work mismatch median "
              f"{float(np.median(pred)):.2e}, scale median "
              f"{float(np.median(scales)):.3f}"
              + (f" against archived {float(np.median(arch)):.3f}"
                 if arch else ""))
    STATE.parent.mkdir(parents=True, exist_ok=True)
    base.atomic_json(STATE, {"schema": "r69_reference_free_state_v1",
                             "created_utc": base.utc_now(), "designs": st})
    print(f"[written] {STATE.name}")
    return 0


def pipeline(args) -> int:
    cutoff = datetime.fromisoformat(args.cutoff) if args.cutoff else None
    if not STATE.exists():
        plan(args)
    state = json.loads(STATE.read_text(encoding="utf-8"))
    if not r48.idle_or_die():
        return 2
    t0 = time.time()
    for design, entries in state["designs"].items():
        for e in entries:
            if e["cell"].get("status") == "ok":
                continue
            if cutoff and datetime.now() > cutoff:
                print("[r69] cutoff reached; stopping cleanly", flush=True)
                base.atomic_json(STATE, state)
                return aggregate(args)
            free_gb = shutil.disk_usage(METRICS).free / 1e9
            if free_gb < DISK_FLOOR_GB:
                print(f"!! {free_gb:.1f} GB free below the floor; stopping",
                      flush=True)
                base.atomic_json(STATE, state)
                return 3
            rec = run_cell(design, e)
            e["cell"] = rec
            base.atomic_json(STATE, state)
            if rec.get("status") != "ok":
                print(f"  !! {design}{e['sobol_index']:03d}: "
                      f"{rec.get('status')}", flush=True)
                continue
            flags = (" MISS" if rec["timing_match_miss"] else "")
            flags += (" CEIL" if rec["comparator_at_ceiling"] else "")
            print(f"  [{(time.time()-t0)/60:6.1f} min] "
                  f"{design}{e['sobol_index']:03d} "
                  f"hp={e['hp_km']:6.1f} N={rec['comparator_degree']:3d} "
                  f"ratio={rec['achieved_time_ratio']:.3f} "
                  f"budget {rec.get('budget_miss_sampled', float('nan')):+.3f}"
                  f"{flags}", flush=True)
    base.atomic_json(STATE, state)
    print(f"[r69] complete in {(time.time()-t0)/3600:.2f} h", flush=True)
    return aggregate(args)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "pipeline", "aggregate"):
        p = sub.add_parser(name)
        if name == "pipeline":
            p.add_argument("--cutoff", default=None)
    args = ap.parse_args()
    return {"plan": plan, "pipeline": pipeline,
            "aggregate": aggregate}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
