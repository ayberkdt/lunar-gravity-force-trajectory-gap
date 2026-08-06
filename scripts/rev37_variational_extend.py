"""R37: enlarge the forced-variational panel toward both full designs.

Section 7.2 states an asymmetry it does not repair. The force-trajectory
reversal is a 128-orbit result and the two force-level rungs of the instrument
ladder are evaluated on all 128, but the rung that identifies the mechanism ---
the only instrument carrying the gravity gradient G --- is an eight-orbit
panel, because each augmented solve costs a propagation of its own. This script
enlarges that panel along the same deterministic rule it was drawn with.

Nothing here is a new method. The per-orbit solve is imported verbatim from
rev14_variational_budget, so an orbit computed here and an orbit computed there
are the same computation; the only additions are the level chain, a wall clock,
and a checkpoint after every orbit.

  selection   for each design the uncensored beta = 1 rows are sorted by
              perilune and n_per points taken at linspace(0, len-1, n_per).
              The levels 4, 7, 13, 25, 40, 56, 64 per design are nested, each a
              superset of the last, and level 1 is exactly the archived panel.
              Every level therefore spans the perilune range, so a run stopped
              part-way leaves a spanning panel and not a low-perilune prefix.

  clock       levels are attempted in order and the run stops submitting at a
              deadline fixed before it started. What the clock leaves
              incomplete is reported as incomplete, with the completion
              fraction of the truncated level.

  reuse       the eight archived orbits are carried over rather than recomputed,
              which is only legitimate if the archived record and the current
              source still describe the same computation. Before any new orbit
              is accepted, one archived orbit is recomputed and its predicted
              ratio compared with the archived value; a disagreement above the
              registered threshold aborts without writing.

The registration is metrics/r37_preregistration.json. This script writes only
metrics/r37_variational_extension.json; metrics/r14_variational_budget.json is
sealed under the R14 manifest and is read, never touched.

Usage:
    python rev37_variational_extend.py --workers 8 --deadline 2026-08-04T08:45
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_variational_budget as vb
from rev13_variational_check import GRADIENT_DEGREE

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PREREG = METRICS / "r37_preregistration.json"
ARCHIVED = METRICS / "r14_variational_budget.json"
OUTPUT = METRICS / "r37_variational_extension.json"

BETA_KEY = f"beta_{vb.BETA:.2f}"


def selection(pareto, design, n_per):
    """The archived rule, unchanged: perilune-sorted, linspace-sampled."""
    pr = [r for r in pareto["designs"][design]["rows"]
          if not r["budgets"][BETA_KEY]["censored"]]
    pr.sort(key=lambda r: r["hp_km"])
    idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
    return [pr[i] for i in idx]


def build_levels(pareto, levels):
    """Ordered (level, design, pareto_row) tasks, each orbit appearing once."""
    seen = set()
    out = []
    for n_per in levels:
        for design in ("A", "B"):
            for p in selection(pareto, design, n_per):
                key = (design, int(p["sobol_index"]))
                if key in seen:
                    continue
                seen.add(key)
                out.append({"level": n_per, "design": design,
                            "sobol_index": int(p["sobol_index"]),
                            "spec": p["budgets"][BETA_KEY],
                            "hp_km": float(p["hp_km"])})
    return out


def make_task(entry, rows_by_design):
    row = rows_by_design[entry["design"]][entry["sobol_index"]]
    return {"design": entry["design"], "row": row, "spec": entry["spec"]}


def finish(rec):
    """The archived driver's post-processing, applied identically."""
    rec["measured"] = vb.measured(rec["design"], rec["sobol_index"])
    if rec["measured"]:
        mf = rec["measured"]["fixed_budget"]
        pf = rec["policies"]["fixed_budget"]["predicted_pos_rms_m"]
        rec["calibration_ratio_fixed"] = (pf / mf) if mf > 0 else None
    return rec


def sign_agrees(rec):
    pred = rec.get("predicted_ratio_fixed_over_atallah")
    meas = rec.get("measured")
    if pred is None or not meas:
        return None
    mr = meas["fixed_budget"] / meas["atallah_budget"]
    return bool((pred < 1.0) == (mr < 1.0))


def summarize(rows, levels, level_counts):
    done = {(r["design"], r["sobol_index"]) for r in rows}
    complete = [n for n in levels
                if level_counts[n] and level_counts[n] <= len(done)
                and all(k in done for k in level_counts["members"][n])]
    agrees = [sign_agrees(r) for r in rows]
    scored = [a for a in agrees if a is not None]
    cal = [r.get("calibration_ratio_fixed") for r in rows
           if r.get("calibration_ratio_fixed") is not None]
    pred_radial = sum(1 for r in rows
                      if (r.get("predicted_ratio_fixed_over_atallah") or 0) > 1.0)
    return {
        "orbits": len(rows),
        "highest_complete_level_per_design": max(complete) if complete else None,
        "highest_complete_level_orbits": (2 * max(complete)) if complete else 0,
        "sign_agreement": sum(1 for a in scored if a),
        "sign_scored": len(scored),
        "sign_disagreements": [
            {"design": r["design"], "sobol_index": r["sobol_index"],
             "hp_km": r["hp_km"],
             "predicted": r["predicted_ratio_fixed_over_atallah"],
             "measured": (r["measured"]["fixed_budget"]
                          / r["measured"]["atallah_budget"])}
            for r in rows if sign_agrees(r) is False],
        "predicted_favors_radial": pred_radial,
        "calibration_ratio_fixed": (
            {"n": len(cal), "median": float(np.median(cal)),
             "min": float(np.min(cal)), "max": float(np.max(cal)),
             "outside_archived_band": [
                 {"design": r["design"], "sobol_index": r["sobol_index"],
                  "hp_km": r["hp_km"],
                  "calibration_ratio_fixed": r["calibration_ratio_fixed"]}
                 for r in rows
                 if r.get("calibration_ratio_fixed") is not None
                 and not (0.94 <= r["calibration_ratio_fixed"] <= 1.04)]}
            if cal else None),
    }


def checkpoint(rows, levels, level_counts, meta):
    payload = {"schema": "r37_variational_extension_v1",
               "created_utc": base.utc_now(),
               "beta": vb.BETA, "gradient_degree": GRADIENT_DEGREE,
               **meta,
               "rows": sorted(rows, key=lambda r: (r["design"],
                                                   r["sobol_index"])),
               "summary": summarize(rows, levels, level_counts)}
    base.atomic_json(OUTPUT, payload)
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--deadline", required=True,
                    help="local ISO time after which no new orbit is submitted")
    a = ap.parse_args()

    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    levels = prereg["selection_rule"]["levels_per_design"]
    deadline = datetime.fromisoformat(a.deadline).timestamp()

    pareto = json.loads(vb.PARETO.read_text(encoding="utf-8"))
    rows_by_design = {
        d: {int(r["sobol_index"]): r
            for r in json.loads(vb.ROWS[d].read_text())["rows"]}
        for d in ("A", "B")}

    entries = build_levels(pareto, levels)
    level_counts = {"members": {}}
    for n_per in levels:
        members = set()
        for d in ("A", "B"):
            for p in selection(pareto, d, n_per):
                members.add((d, int(p["sobol_index"])))
        level_counts["members"][n_per] = members
        level_counts[n_per] = len(members)

    archived = {(r["design"], r["sobol_index"]): r
                for r in json.loads(ARCHIVED.read_text(encoding="utf-8"))["rows"]}
    print(f"[r37] {len(entries)} orbits declared over levels {levels}; "
          f"{len(archived)} archived, deadline {a.deadline}", flush=True)

    # ---- admissibility self-check: recompute one archived orbit
    check_name = prereg["admissibility_self_check"]["recomputed_orbit"]
    cd, ci = check_name[0], int(check_name[1:])
    ref = archived[(cd, ci)]
    spec = next(e["spec"] for e in entries
                if e["design"] == cd and e["sobol_index"] == ci)
    print(f"[r37] self-check: recomputing {check_name}", flush=True)
    t0 = time.time()
    got = vb.worker(make_task({"design": cd, "sobol_index": ci, "spec": spec},
                              rows_by_design))
    if got["status"] != "complete":
        print(f"[r37] ABORT: self-check orbit failed: {got.get('message')}")
        return 1
    old = ref["predicted_ratio_fixed_over_atallah"]
    new = got["predicted_ratio_fixed_over_atallah"]
    rel = abs(new - old) / abs(old)
    thr = float(prereg["admissibility_self_check"]["abort_threshold_rel"])
    print(f"[r37] self-check {check_name}: archived {old:.10g}, recomputed "
          f"{new:.10g}, rel {rel:.2e} (threshold {thr:g}), "
          f"{time.time()-t0:.0f}s", flush=True)
    if rel > thr:
        print("[r37] ABORT: archived record and current source disagree; "
              "nothing written.")
        return 1

    meta = {"registration_sha256": prereg.get("preregistration_sha256"),
            "deadline_local": a.deadline,
            "levels_per_design": levels,
            "self_check": {"orbit": check_name, "archived": old,
                           "recomputed": new, "rel": rel, "threshold": thr,
                           "passed": True},
            "reuse_note": ("the eight archived orbits are carried over from "
                           "r14_variational_budget.json unchanged; every other "
                           "orbit here was solved by this run with the same "
                           "imported worker."),
            "completed_levels_note": (
                "highest_complete_level_orbits is the panel; any orbits beyond "
                "it belong to a level the wall clock truncated and are reported "
                "as such.")}

    rows = [dict(archived[k], reused_from_r14=True) for k in archived]
    pending = [e for e in entries
               if (e["design"], e["sobol_index"]) not in archived]
    checkpoint(rows, levels, level_counts, meta)
    print(f"[r37] {len(pending)} orbits to solve", flush=True)

    t0 = time.time()
    submitted, completed, skipped = 0, 0, 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {}
        it = iter(pending)
        # keep the pool full, but stop submitting once the clock is spent
        for _ in range(min(a.workers, len(pending))):
            e = next(it, None)
            if e is None:
                break
            futs[pool.submit(vb.worker, make_task(e, rows_by_design))] = e
            submitted += 1
        while futs:
            for fut in as_completed(list(futs)):
                e = futs.pop(fut)
                rec = fut.result()
                if rec["status"] == "complete":
                    rec = finish(rec)
                    rec["level"] = e["level"]
                    rec["reused_from_r14"] = False
                    rows.append(rec)
                    completed += 1
                    agr = sign_agrees(rec)
                    print(f"  [{completed}/{len(pending)}] L{e['level']:02d} "
                          f"{e['design']}{e['sobol_index']:03d} "
                          f"hp={e['hp_km']:6.1f} "
                          f"pred={rec['predicted_ratio_fixed_over_atallah']:.4g} "
                          f"sign={'ok' if agr else 'DISAGREE' if agr is False else 'n/a'} "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
                else:
                    print(f"  !! {e['design']}{e['sobol_index']:03d} "
                          f"{rec.get('message')}", flush=True)
                checkpoint(rows, levels, level_counts, meta)
                if time.time() < deadline:
                    nxt = next(it, None)
                    if nxt is not None:
                        futs[pool.submit(vb.worker,
                                         make_task(nxt, rows_by_design))] = nxt
                        submitted += 1
                else:
                    skipped = len(pending) - submitted
                break

    payload = checkpoint(rows, levels, level_counts, meta)
    s = payload["summary"]
    print(f"\n[r37] {s['orbits']} orbits in the record "
          f"({completed} solved here, {len(archived)} reused, "
          f"{skipped} never submitted)")
    print(f"[r37] highest complete level: "
          f"{s['highest_complete_level_orbits']} orbits "
          f"(level {s['highest_complete_level_per_design']} per design)")
    print(f"[r37] sign agreement {s['sign_agreement']}/{s['sign_scored']}; "
          f"predicted favors radial on {s['predicted_favors_radial']}")
    if s["sign_disagreements"]:
        print("[r37] sign disagreements:")
        for d in s["sign_disagreements"]:
            print(f"    {d['design']}{d['sobol_index']:03d} hp={d['hp_km']:.1f} "
                  f"pred={d['predicted']:.4g} meas={d['measured']:.4g}")
    c = s["calibration_ratio_fixed"]
    if c:
        print(f"[r37] calibration channel median {c['median']:.4f} "
              f"[{c['min']:.3f}, {c['max']:.3f}], "
              f"{len(c['outside_archived_band'])} outside the 0.94-1.04 band")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
