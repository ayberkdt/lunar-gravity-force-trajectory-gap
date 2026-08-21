"""R68 (O59/O60): the full-design comparison at equal *measured* kernel time.

Why this exists
---------------
The headline population comparison is matched on a per-call quadratic proxy
calibrated on the reference arc. Measured machine time has only ever been
matched on declared panels of fourteen orbits -- (O13) for the radial
endpoint, (O48)/(O57)/(O58) for the interior member -- and those panels are
explicitly not population tallies. A reader is therefore entitled to say that
the paper is framed around compute allocation while its principal comparison
is not matched on compute as the machine spends it.

This campaign runs the same measured-time construction on the whole of both
coverage designs, 64 orbits each, at beta = 1:

  endpoint arm (O59): the budget-calibrated radial endpoint, k = 1.00, against
                      a constant degree matched on its measured gravity-kernel
                      time;
  interior arm (O60): the sampled interior member k = 0.50 against its own
                      measured-time-matched constant degree.

Both arms are scored where the errors are read, at the tighter level, which is
the convention R44 established and R48/R65 follow. The member is re-propagated
serially here even though an archived tighter run exists: the archived one was
produced with concurrent workers and its kernel time is not comparable.

What is new against (O58)
-------------------------
Only the population and the first-guess rule. The band is tightened from
0.90--1.10 to 0.95--1.05 because a population tally is more sensitive to a
systematic band edge than a fourteen-cell panel is, and the first comparator
degree is taken from

    N_1 = c^{-1}( T_member_measured / n_RHS(constant, archived) )

rather than from the member's mean per-call cost. The quantity being matched
is *total* kernel time, and the constant-degree run makes fewer calls than the
member does, so inverting the mean per-call cost alone starts systematically
low: in (O58) the first pass landed near 0.6 and every cell paid for a second
propagation. n_RHS for a constant degree varies weakly with the degree, so the
archived count is a good enough predictor to start from.

Nothing is selected after the fact. Every orbit of both designs enters, every
cell is reported, and a cell that cannot be matched inside the band after
MAX_REFINE passes is reported as a timing-match miss rather than dropped. A
cell whose matched degree reaches the orbit's reference degree cannot be
scored against that reference and is reported as ceiling-censored, not as a
win for either policy.

Timing comparability requires an idle machine. The guard runs once per
invocation and does not re-arm; rev68_quiet_audit.py is the after-the-fact
test that catches a machine joined mid-run.

Usage:
    python rev68_timing_full.py plan --arm endpoint
    python rev68_timing_full.py pipeline --arm endpoint [--cutoff ISO8601]
    python rev68_timing_full.py aggregate --arm endpoint
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
import rev14_budget_trajectory as r14              # noqa: E402
import rev48_interior_timing as r48                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r68_cases"
RAW_ROOT = METRICS / "r68_raw"

ARMS = {"endpoint": "1.00", "interior": "0.50"}
DESIGNS = ("A", "B")
BETA_TAG = "beta_1.00"
BAND_LO, BAND_HI = 0.95, 1.05
MAX_REFINE = 5
COST_EXPONENT = r48.COST_EXPONENT
DISK_FLOOR_GB = 6.0

LEVELS = r14.LEVELS
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


def out_path(arm: str) -> Path:
    return METRICS / f"r68_timing_full_{arm}.json"


def state_path(arm: str) -> Path:
    return METRICS / f"r68_timing_full_{arm}_state.json"


def case_dir(arm: str, design: str, index: int) -> Path:
    return CASE_ROOT / arm / design / f"sobol{design}_{index:03d}"


def raw_dir(arm: str, design: str, index: int) -> Path:
    return RAW_ROOT / arm / design / f"sobol{design}_{index:03d}"


def span_rows(design: str) -> list:
    p = METRICS / f"r18_span_sweep_{design}_{BETA_TAG}.json"
    return json.loads(p.read_text(encoding="utf-8"))["rows"]


def member_sidecar(arm: str, design: str, index: int, level: str) -> dict:
    """Archived member run. Leaf names are sobolA_ for every design, an
    archive convention recorded in R18; config.design is authoritative."""
    if arm == "endpoint":
        p = (METRICS / "r14_cases" / f"{design}_{BETA_TAG}"
             / f"sobolA_{index:03d}" / f"atallah_budget_{level}.json")
    else:
        p = (METRICS / "r18_cases" / f"{design}_{BETA_TAG}_k_{ARMS[arm]}"
             / f"sobolA_{index:03d}" / f"span_{level}.json")
    return json.loads(p.read_text(encoding="utf-8"))


def member_degree_fn(arm: str, design: str, index: int):
    """The archived degree table, read back as the same step function the
    campaign propagated: 10 km bins, clamped at both ends."""
    cfg = member_sidecar(arm, design, index, "tighter")["config"]
    key = "atallah_degree_table" if arm == "endpoint" else "degree_table"
    tab = {float(a): int(b) for a, b in cfg[key].items()}
    lo, hi = min(tab), max(tab)

    def degree_of(tt, h_m, _t=tab, _lo=lo, _hi=hi):
        hb = min(_hi, max(_lo, 10.0 * math.floor(h_m / 1e4)))
        return _t[hb]
    return degree_of, cfg


def archived_constant_rhs(design: str, index: int) -> int | None:
    """n_RHS of the archived constant-degree run at the tighter level."""
    sidecar, _ = r14.reuse_paths(design, index, "fixed_critical", "tighter")
    if not sidecar.exists():
        return None
    tel = json.loads(sidecar.read_text(encoding="utf-8")).get("telemetry") or {}
    n = tel.get("n_rhs")
    return int(n) if n else None


def panel(arm: str) -> dict:
    """Every orbit of both designs that carries a scored member. No selection."""
    k = ARMS[arm]
    out = {}
    for design in DESIGNS:
        rows = []
        for r in span_rows(design):
            e = r["entries"].get(k, {})
            if e.get("error_m") is None:
                continue
            index = int(r["sobol_index"])
            side = member_sidecar(arm, design, index, "tighter")
            rows.append({
                "sobol_index": index,
                "hp_km": r["hp_km"],
                "ha_km": r.get("ha_km"),
                "n_critical": int(r["n_critical"]),
                "constant_degree": r.get("constant_degree"),
                "adopted_truth_degree":
                    int(side["config"]["adopted_truth_degree"]),
                "archived_constant_n_rhs": archived_constant_rhs(design, index),
                "member_error_m": e["error_m"],
                "member_envelope_m": e.get("envelope_m") or 0.0,
                "cell": {"status": "pending"}})
        out[design] = rows
    return out


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
        "schema": "r68_timing_full_case_v1", "created_utc": base.utc_now(),
        "config": dict(cfg_extra, level=level, timing_comparable=timed,
                       adopted_truth_degree=adopted,
                       source=base.provenance()),
        "status": st, "event": ev, "telemetry": tel,
        "raw_path": str(raw.relative_to(ROOT)),
        "raw_sha256": base.file_hash(raw)})
    return json.loads(sidecar.read_text(encoding="utf-8"))


def first_guess(t_member_ns: float, n_rhs_constant: int | None,
                n_rhs_member: int, adopted: int,
                deg_tab, ns_tab) -> tuple[int, str]:
    """Total-time inversion; falls back to the member's own call count."""
    if n_rhs_constant:
        rule = "total time over the archived constant-degree call count"
        calls = n_rhs_constant
    else:
        rule = "total time over the member's own call count"
        calls = n_rhs_member
    per_call = t_member_ns / max(calls, 1)
    n = int(r13.inverse_cost(per_call, deg_tab, ns_tab))
    return max(2, min(n, adopted)), rule


def _exponent(tried: list) -> float:
    """Local slope of measured time against degree, from the run's own passes.

    The first pass has nothing to fit and uses the quadratic assumption the
    panel campaigns used. From the second on, two measured points give the
    exponent directly, which matters because the kernel is super-quadratic at
    high degree and the comparator's call count is not the member's: assuming
    2.0 throughout costs a third propagation on most cells.
    """
    if len(tried) < 2:
        return COST_EXPONENT
    (d1, r1), (d2, r2) = ((t["degree"], t["ratio"]) for t in tried[-2:])
    if d1 == d2 or r1 <= 0 or r2 <= 0:
        return COST_EXPONENT
    p = math.log(r2 / r1) / math.log(d2 / d1)
    return min(3.5, max(1.2, p))


def run_cell(arm: str, design: str, e: dict, deg_tab, ns_tab) -> dict:
    """One orbit: member kernel time, comparator refined to the band, envelope."""
    index = int(e["sobol_index"])
    adopted = int(e["adopted_truth_degree"])
    rec = {"refinements": []}
    degree_of, mcfg = member_degree_fn(arm, design, index)
    y0 = np.asarray(mcfg["initial_state_si"], dtype=float)
    cd, rd = case_dir(arm, design, index), raw_dir(arm, design, index)

    p = _propagate(
        {"design": design, "sobol_index": index, "arm": arm, "k": ARMS[arm],
         "policy": f"member k={ARMS[arm]} (serial re-run)",
         "purpose": "contention-free member kernel time at the scoring "
                    "tolerance"},
        degree_of, y0, adopted, "tighter",
        cd / "member_tighter.json", rd / "member_tighter.npz", timed=True)
    if p is None:
        return {"status": "member_failed"}
    t_member = p["telemetry"]["gravity_kernel_ns"]
    rec["member_kernel_ns"] = t_member
    rec["member_n_rhs"] = p["telemetry"]["n_rhs"]
    rec["member_mean_degree_sq"] = p["telemetry"]["mean_degree_sq"]

    degree, rule = first_guess(t_member, e.get("archived_constant_n_rhs"),
                               p["telemetry"]["n_rhs"], adopted,
                               deg_tab, ns_tab)
    rec["first_guess_rule"] = rule
    rec["first_guess_degree"] = degree

    tried, best = [], None
    for step in range(MAX_REFINE):
        def const(tt, h_m, _n=degree):
            return _n
        pc = _propagate(
            {"design": design, "sobol_index": index, "arm": arm,
             "policy": "constant matched on measured kernel time at the "
                       "scoring tolerance",
             "degree": degree, "refine_step": step},
            const, y0, adopted, "tighter",
            cd / f"fixed_{degree}_tighter.json",
            rd / f"fixed_{degree}_tighter.npz", timed=True)
        if pc is None:
            return {"status": "comparator_failed"}
        ratio = pc["telemetry"]["gravity_kernel_ns"] / t_member
        tried.append({"degree": degree, "ratio": ratio,
                      "kernel_ns": pc["telemetry"]["gravity_kernel_ns"]})
        if best is None or abs(ratio - 1.0) < abs(best["ratio"] - 1.0):
            best = {"degree": degree, "ratio": ratio}
        if BAND_LO <= ratio <= BAND_HI:
            break
        nxt = int(round(degree * (1.0 / ratio) ** (1.0 / _exponent(tried))))
        nxt = max(2, min(nxt, adopted))
        if nxt == degree or any(t["degree"] == nxt for t in tried):
            break
        degree = nxt

    rec["refinements"] = tried
    rec["comparator_degree"] = best["degree"]
    rec["achieved_time_ratio"] = best["ratio"]
    rec["timing_match_miss"] = not (BAND_LO <= best["ratio"] <= BAND_HI)
    rec["comparator_at_ceiling"] = bool(best["degree"] >= adopted)

    # The two-level envelope needs the pair; this run is not timed.
    def constb(tt, h_m, _n=best["degree"]):
        return _n
    _propagate(
        {"design": design, "sobol_index": index, "arm": arm,
         "policy": "constant matched on measured kernel time at the scoring "
                   "tolerance", "degree": best["degree"],
         "envelope_run": True},
        constb, y0, adopted, "tight",
        cd / f"fixed_{best['degree']}_tight.json",
        rd / f"fixed_{best['degree']}_tight.npz", timed=False)
    rec["status"] = "ok"
    return rec


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def score_cell(arm: str, design: str, e: dict) -> dict | None:
    index = int(e["sobol_index"])
    rec = e["cell"]
    if rec.get("status") != "ok":
        return None
    truth = {}
    for lv in ("tight", "tighter"):
        _, raw = r14.reuse_paths(design, index, "truth", lv)
        if not raw.exists():
            return None
        truth[lv] = _load(raw)
    s_ref = base.common_error(truth["tight"][0], truth["tight"][1],
                              truth["tighter"][0],
                              truth["tighter"][1])["pos_rms_m"]
    rd = raw_dir(arm, design, index)
    n = rec["comparator_degree"]
    rg, rt = rd / f"fixed_{n}_tighter.npz", rd / f"fixed_{n}_tight.npz"
    if not (rg.exists() and rt.exists()):
        return None
    got_g, got_t = _load(rg), _load(rt)
    err = base.common_error(got_g[0], got_g[1], truth["tighter"][0],
                            truth["tighter"][1])["pos_rms_m"]
    s_fix = base.common_error(got_t[0], got_t[1], got_g[0],
                              got_g[1])["pos_rms_m"]
    thr = e["member_envelope_m"] + s_fix + s_ref
    diff = err - e["member_error_m"]
    member_name = "radial" if arm == "endpoint" else "interior"
    return {
        "sobol_index": index, "design": design,
        "hp_km": e["hp_km"], "ha_km": e.get("ha_km"),
        "n_critical": e["n_critical"],
        "constant_endpoint_degree": e.get("constant_degree"),
        "adopted_truth_degree": e["adopted_truth_degree"],
        "comparator_degree": n,
        "achieved_time_ratio": rec["achieved_time_ratio"],
        "timing_match_miss": rec["timing_match_miss"],
        "comparator_at_ceiling": rec["comparator_at_ceiling"],
        "refinement_passes": len(rec["refinements"]),
        "member_kernel_ns": rec["member_kernel_ns"],
        "member_error_m": e["member_error_m"],
        "comparator_error_m": err,
        "resolution_threshold_m": thr,
        "resolution_margin": abs(diff) / thr if thr else None,
        "resolved": bool(abs(diff) > thr),
        "winner": (member_name if diff > thr else
                   ("constant" if -diff > thr else None)),
        "rho_constant_over_member": (err / e["member_error_m"]
                                     if e["member_error_m"] else None)}


def aggregate(args) -> int:
    arm = args.arm
    st = json.loads(state_path(arm).read_text(encoding="utf-8"))
    member_name = "radial" if arm == "endpoint" else "interior"
    by_design, rows = {}, []
    for design, entries in st["designs"].items():
        drows = [r for r in (score_cell(arm, design, e) for e in entries)
                 if r is not None]
        rows.extend(drows)
        scored = [r for r in drows if not r["comparator_at_ceiling"]]
        res = [r for r in scored if r["resolved"]]
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
            f"resolved_{member_name}_wins":
                sum(1 for r in res if r["winner"] == member_name),
            "resolved_constant_wins":
                sum(1 for r in res if r["winner"] == "constant"),
            f"raw_{member_name}_wins":
                sum(1 for r in scored
                    if r["comparator_error_m"] > r["member_error_m"]),
            "median_rho_constant_over_member":
                (float(np.median(rho)) if rho else None),
            "median_achieved_time_ratio":
                float(np.median([r["achieved_time_ratio"] for r in scored]))
                if scored else None,
            "mean_refinement_passes":
                float(np.mean([r["refinement_passes"] for r in scored]))
                if scored else None}
    payload = {
        "schema": "r68_timing_full_v1", "created_utc": base.utc_now(),
        "arm": arm, "member_k": ARMS[arm], "beta": 1.0,
        "protocol": (
            "every orbit of both coverage designs; the member is re-propagated "
            "serially for a contention-free gravity-kernel time at the tighter "
            "level, and a constant degree is refined until its measured kernel "
            "time falls in the band; errors are read at the tighter level, the "
            "level the match is made at"),
        "timing_band": [BAND_LO, BAND_HI],
        "max_refine": MAX_REFINE,
        "compare_with": (f"r14_trajectory_*_{BETA_TAG}.json (endpoint) and "
                         f"r18_span_sweep_*_{BETA_TAG}.json (interior) hold "
                         "the same pairs at equal nominal per-call budget"),
        "by_design": by_design, "rows": rows, "source": base.provenance()}
    base.atomic_json(out_path(arm), payload)
    print(f"[written] {out_path(arm).name}")
    for design, s in by_design.items():
        print(f"  {design}: scored {s['cells_scored']}, "
              f"resolved {s['resolved']} "
              f"({s[f'resolved_{member_name}_wins']} {member_name} / "
              f"{s['resolved_constant_wins']} constant), "
              f"median rho {s['median_rho_constant_over_member']}, "
              f"misses {s['timing_match_misses']}, "
              f"ceiling {s['ceiling_censored']}")
    return 0


def plan(args) -> int:
    st = panel(args.arm)
    deg_tab, ns_tab = r13.cost_curve()
    total = sum(len(v) for v in st.values())
    print(f"[r68-plan] arm={args.arm} k={ARMS[args.arm]}  {total} orbits "
          f"(designs {', '.join(DESIGNS)}), band "
          f"{BAND_LO:.2f}-{BAND_HI:.2f}")
    for design, rows in st.items():
        miss = sum(1 for r in rows if not r["archived_constant_n_rhs"])
        print(f"  {design}: {len(rows)} orbits, "
              f"{miss} without an archived constant call count")
    return 0


def pipeline(args) -> int:
    arm = args.arm
    cutoff = datetime.fromisoformat(args.cutoff) if args.cutoff else None
    sp = state_path(arm)
    state = (json.loads(sp.read_text(encoding="utf-8")) if sp.exists()
             else {"schema": "r68_timing_full_state_v1",
                   "arm": arm, "created_utc": base.utc_now(),
                   "designs": panel(arm)})
    if not r48.idle_or_die():
        return 2
    deg_tab, ns_tab = r13.cost_curve()
    t0 = time.time()
    done = 0
    for design, entries in state["designs"].items():
        for e in entries:
            if e["cell"].get("status") == "ok":
                continue
            if cutoff and datetime.now() > cutoff:
                print("[r68] cutoff reached; stopping cleanly", flush=True)
                base.atomic_json(sp, state)
                return aggregate(args)
            free_gb = shutil.disk_usage(METRICS).free / 1e9
            if free_gb < DISK_FLOOR_GB:
                print(f"!! {free_gb:.1f} GB free below the {DISK_FLOOR_GB} GB "
                      f"floor; stopping before this cell", flush=True)
                base.atomic_json(sp, state)
                return 3
            rec = run_cell(arm, design, e, deg_tab, ns_tab)
            e["cell"] = rec
            base.atomic_json(sp, state)
            done += 1
            if rec.get("status") != "ok":
                print(f"  !! {design}{e['sobol_index']:03d}: "
                      f"{rec.get('status')}", flush=True)
                continue
            flags = ("" if not rec["timing_match_miss"] else " MISS")
            flags += ("" if not rec["comparator_at_ceiling"] else " CEIL")
            print(f"  [{(time.time()-t0)/60:6.1f} min] "
                  f"{design}{e['sobol_index']:03d} "
                  f"hp={e['hp_km']:6.1f} N={rec['comparator_degree']:3d} "
                  f"ratio={rec['achieved_time_ratio']:.3f} "
                  f"({len(rec['refinements'])} passes){flags}", flush=True)
    base.atomic_json(sp, state)
    print(f"[r68] arm {arm} complete: {done} cells this session, "
          f"{(time.time()-t0)/3600:.2f} h", flush=True)
    return aggregate(args)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "pipeline", "aggregate"):
        p = sub.add_parser(name)
        p.add_argument("--arm", choices=sorted(ARMS), required=True)
        if name == "pipeline":
            p.add_argument("--cutoff", default=None,
                           help="ISO 8601; stop cleanly before the next cell")
    args = ap.parse_args()
    return {"plan": plan, "pipeline": pipeline,
            "aggregate": aggregate}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
