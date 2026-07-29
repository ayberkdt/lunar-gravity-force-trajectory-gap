"""R23-D: one cost curve, measured under both archived protocols at once.

The manuscript quotes two measured cost curves for the same kernel on the same
machine and calls them one curve. Section 5 anchors on ``r1_kernel_timing``
(N=60 -> 25.9 us, N=120 -> 79.0, N=300 -> 423.6); the timing-matched comparator
degrees are produced by inverting ``r12_kernel_cost_curve`` (18.5 / 61.5 /
366.5), which is 1.40, 1.28 and 1.16 times faster at the same degrees. Neither
protocol is described in the text, and a reader cannot tell which result rests
on which curve.

Reading the two scripts shows they differ in two ways at once, and either could
produce a gap of this size:

  timing protocol   r1 times a *block* of consecutive calls and divides by the
                    block length, over 9 blocks, with the block length falling
                    from 400 calls at low degree to 11 at degree 600. r12 times
                    every call separately and takes the median over all of them,
                    so each sample carries one timer-call overhead and the full
                    per-call jitter.
  model degree      r12 evaluates with the coefficient arrays of a degree-300
                    model, r13's high-degree extension with a degree-900 model.
                    The arrays differ in size by an order of magnitude, so the
                    same evaluation degree does not see the same cache.

This script measures both protocols, at both model degrees, on one degree
ladder, in one session on an idle machine. That turns an unexplained
discrepancy into a measured decomposition: whatever the gap is, it is
attributable here rather than left as two numbers that disagree.

It must run alone. Nothing else may be on the machine.

Usage:
    python rev23_cost_curve_unified.py --target-s 0.25
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np

from rev3_common import load_model, kernel_args, warmup, OUT
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

# the degrees the manuscript actually quotes or inverts, plus enough of the
# ladder to see the shape; 400-900 only exist under the degree-900 model
LADDER_LOW = [10, 20, 40, 60, 80, 100, 120, 160, 200, 250, 300]
LADDER_HIGH = [400, 600, 900]
BLOCKS = 9                      # r1's block count
OUTPUT = Path(OUT) / "r23_cost_curve_unified.json"


def _percall(x, y, z, degree, args, repeats: int) -> np.ndarray:
    """r12's protocol: every call timed on its own."""
    s = np.empty(repeats)
    for i in range(repeats):
        t0 = time.perf_counter_ns()
        sh_accel_fixed_numba(x, y, z, degree, *args)
        s[i] = time.perf_counter_ns() - t0
    return s


def _block(x, y, z, degree, args, reps: int, blocks: int) -> np.ndarray:
    """r1's protocol: a block of consecutive calls, divided by the block."""
    s = np.empty(blocks)
    for b in range(blocks):
        t0 = time.perf_counter_ns()
        for _ in range(reps):
            sh_accel_fixed_numba(x, y, z, degree, *args)
        s[b] = (time.perf_counter_ns() - t0) / reps
    return s


def _calibrate_reps(x, y, z, degree, args, target_s: float) -> int:
    """How many calls fit in the per-cell time budget, measured not guessed."""
    t0 = time.perf_counter_ns()
    for _ in range(5):
        sh_accel_fixed_numba(x, y, z, degree, *args)
    per_call_s = (time.perf_counter_ns() - t0) / 5 / 1e9
    if per_call_s <= 0:
        return 400
    return int(min(4000, max(11, target_s / per_call_s)))


def stats(ns: np.ndarray) -> dict:
    return {"per_call_ns_median": float(np.median(ns)),
            "per_call_ns_p10": float(np.percentile(ns, 10)),
            "per_call_ns_p90": float(np.percentile(ns, 90)),
            "per_call_ns_min": float(ns.min()),
            "samples": int(ns.size),
            "spread_p90_over_p10": float(
                np.percentile(ns, 90) / np.percentile(ns, 10))}


def measure_model(model_degree: int, ladder: list[int],
                  target_s: float) -> list[dict]:
    model = load_model(model_degree)
    args = kernel_args(model)
    warmup(model, args)
    r_m = model.r_ref + 100e3
    x, y, z = r_m, 0.0, 0.0
    rows = []
    for n in ladder:
        for _ in range(50):                       # warm this specific degree
            sh_accel_fixed_numba(x, y, z, n, *args)
        reps = _calibrate_reps(x, y, z, n, args, target_s)
        pc = _percall(x, y, z, n, args, reps)
        bl = _block(x, y, z, n, args, reps, BLOCKS)
        row = {"model_degree": model_degree, "degree": n,
               "reps_per_block": reps, "blocks": BLOCKS,
               "eval_radius_m": r_m,
               "per_call_protocol": stats(pc),
               "block_protocol": stats(bl)}
        row["block_over_percall"] = (row["block_protocol"]["per_call_ns_median"]
                                     / row["per_call_protocol"]
                                     ["per_call_ns_median"])
        rows.append(row)
        print(f"  model {model_degree:3d}  N={n:4d}  "
              f"per-call {pc.mean()*0+np.median(pc)/1e3:8.2f} us   "
              f"block {np.median(bl)/1e3:8.2f} us   "
              f"block/per-call {row['block_over_percall']:.3f}", flush=True)
    return rows


def archived_comparison(rows: list[dict]) -> dict:
    """What the two archived curves say at the degrees the text quotes."""
    out = {}
    for name, fname, key in (
            ("r1_kernel_timing", "r1_kernel_timing.json", "degree_sweep"),
            ("r12_kernel_cost_curve", "r12_kernel_cost_curve.json", "rows")):
        p = Path(OUT) / fname
        if not p.exists():
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        got = {}
        for r in d.get(key, []):
            if "median_us" in r:
                got[int(r["degree"])] = float(r["median_us"]) * 1e3
            elif "per_call_ns_median" in r:
                got[int(r["degree"])] = float(r["per_call_ns_median"])
        out[name] = got
    table = []
    for n in (60, 120, 300):
        cell = {"degree": n}
        for name, got in out.items():
            if n in got:
                cell[name + "_ns"] = got[n]
        for r in rows:
            if r["degree"] != n:
                continue
            cell[f"remeasured_model{r['model_degree']}_block_ns"] = \
                r["block_protocol"]["per_call_ns_median"]
            cell[f"remeasured_model{r['model_degree']}_percall_ns"] = \
                r["per_call_protocol"]["per_call_ns_median"]
        table.append(cell)
    return {"archived": out, "at_quoted_degrees": table}


def session_spread(sessions: list[list[dict]]) -> list[dict]:
    """Session-to-session reproducibility, which is the quantity actually in
    dispute. A tighter median inside one session says nothing about whether a
    second session would land on the same number; the archived repeatability
    check found factors of 1.78, 0.88 and 1.33 doing exactly this comparison."""
    if len(sessions) < 2:
        return []
    keyed: dict[tuple, list[dict]] = {}
    for s in sessions:
        for r in s:
            keyed.setdefault((r["model_degree"], r["degree"]), []).append(r)
    out = []
    for (md, n), rs in sorted(keyed.items()):
        row = {"model_degree": md, "degree": n, "sessions": len(rs)}
        for proto in ("per_call_protocol", "block_protocol"):
            v = np.array([r[proto]["per_call_ns_median"] for r in rs])
            row[proto] = {
                "median_ns": float(np.median(v)),
                "min_ns": float(v.min()), "max_ns": float(v.max()),
                "max_over_min": float(v.max() / v.min()),
                "rel_spread": float((v.max() - v.min()) / np.median(v))}
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-s", type=float, default=0.25,
                    help="measurement budget per protocol per degree")
    ap.add_argument("--sessions", type=int, default=1,
                    help="repeat the whole panel this many times, so "
                         "session-to-session reproducibility is measured "
                         "rather than assumed")
    ap.add_argument("--out", default=None,
                    help="output filename under metrics/")
    a = ap.parse_args()

    out_path = Path(OUT) / a.out if a.out else OUTPUT
    t0 = time.time()
    print(f"[r23d] unified cost curve, target {a.target_s}s per cell, "
          f"{a.sessions} session(s)", flush=True)
    sessions = []
    for s in range(a.sessions):
        print(f"[r23d] === session {s + 1}/{a.sessions} ===", flush=True)
        srows = []
        print("[r23d] model degree 300", flush=True)
        srows += measure_model(300, LADDER_LOW, a.target_s)
        print("[r23d] model degree 900", flush=True)
        srows += measure_model(900, LADDER_LOW + LADDER_HIGH, a.target_s)
        for r in srows:
            r["session"] = s + 1
        sessions.append(srows)
    rows = sessions[-1]

    by = {(r["model_degree"], r["degree"]): r for r in rows}
    model_effect = []
    for n in LADDER_LOW:
        lo, hi = by.get((300, n)), by.get((900, n))
        if lo and hi:
            model_effect.append({
                "degree": n,
                "model900_over_model300_percall": (
                    hi["per_call_protocol"]["per_call_ns_median"]
                    / lo["per_call_protocol"]["per_call_ns_median"]),
                "model900_over_model300_block": (
                    hi["block_protocol"]["per_call_ns_median"]
                    / lo["block_protocol"]["per_call_ns_median"])})

    r900 = by.get((900, 900))
    r300 = by.get((900, 300))
    ratio = {}
    if r900 and r300:
        for proto in ("block_protocol", "per_call_protocol"):
            ratio[proto] = (r900[proto]["per_call_ns_median"]
                            / r300[proto]["per_call_ns_median"])

    payload = {
        "schema": "r23_cost_curve_unified_v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "machine": {"platform": platform.platform(),
                    "processor": platform.processor(),
                    "single_threaded": True, "exclusive_machine": True},
        "why": ("the manuscript quotes two measured curves for one kernel; "
                "this measures both archived protocols at both archived model "
                "degrees in one session so the discrepancy is decomposed "
                "rather than restated"),
        "protocols": {
            "per_call": "each call timed separately, median over all samples "
                        "(the r12 protocol)",
            "block": f"{BLOCKS} blocks of consecutive calls, block time divided "
                     f"by block length, median over blocks (the r1 protocol)"},
        "target_s_per_cell": a.target_s,
        "n_sessions": a.sessions,
        "rows": rows,
        "all_sessions": [r for s in sessions for r in s]
        if a.sessions > 1 else None,
        "session_to_session_spread": session_spread(sessions),
        "ratio_900_over_300_per_session": [
            {"session": i + 1,
             "block": (by_s.get((900, 900), {}).get("block_protocol", {})
                       .get("per_call_ns_median", float("nan"))
                       / by_s.get((900, 300), {}).get("block_protocol", {})
                       .get("per_call_ns_median", float("nan"))),
             "per_call": (by_s.get((900, 900), {}).get("per_call_protocol", {})
                          .get("per_call_ns_median", float("nan"))
                          / by_s.get((900, 300), {})
                          .get("per_call_protocol", {})
                          .get("per_call_ns_median", float("nan")))}
            for i, by_s in enumerate(
                {(r["model_degree"], r["degree"]): r for r in s}
                for s in sessions)],
        "model_degree_effect": model_effect,
        "ratio_900_over_300_remeasured": ratio,
        "quadratic_expectation_900_over_300": 9.0,
        "archived_curve_comparison": archived_comparison(rows),
        "wall_s": time.time() - t0,
    }
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r23d] written {out_path.name} in {(time.time()-t0)/60:.1f} min")
    for r in payload["session_to_session_spread"]:
        if r["degree"] in (60, 120, 300, 900):
            print(f"  session spread N={r['degree']} model {r['model_degree']}:"
                  f" per-call max/min {r['per_call_protocol']['max_over_min']:.2f}"
                  f", block {r['block_protocol']['max_over_min']:.2f}")
    for r in payload["ratio_900_over_300_per_session"]:
        print(f"  session {r['session']}: c(900)/c(300) block {r['block']:.2f}, "
              f"per-call {r['per_call']:.2f}")
    for cell in payload["archived_curve_comparison"]["at_quoted_degrees"]:
        print(f"  N={cell['degree']}: " + ", ".join(
            f"{k.replace('_ns','')} {v/1e3:.1f}us"
            for k, v in cell.items() if k != "degree"))
    if ratio:
        print(f"  c(900)/c(300) remeasured: " + ", ".join(
            f"{k} {v:.2f}" for k, v in ratio.items())
            + " (pure quadratic 9.00)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
