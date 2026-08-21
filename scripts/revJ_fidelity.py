"""Is the J-campaign harness the same instrument as the manuscript's?

The three JGCD campaigns are run by new code. New code invites a fair question:
when J1 reports that the radial policy wins the force metric on a second
gravity solution, is that the same measurement the manuscript made, or a
lookalike?

The question is settled rather than argued. The harness is pointed at the
*archived* population, the *archived* reference trajectories and the *archived*
accuracy parameters, and asked to recompute the deterministic force defect that
the frozen budget records already contain. Nothing is re-propagated and nothing
is re-calibrated: the same inputs go in, and the numbers that come out either
match the archive or they do not.

The defect is deterministic -- it is a difference of two kernel evaluations at
archived positions -- so the tolerance for a match is not statistical. Anything
above a rounding-level relative deviation means the two instruments differ, and
the campaigns would have to be read against the new instrument rather than
against the manuscript.

Usage:
    python revJ_fidelity.py --workers 11
"""

from __future__ import annotations

import argparse
import json
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import revJ_common as J

J.select_field("JGGRX_1800F")
J.install_field()

import rev12_atallah as at                                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "rJ_fidelity_check.json"
LOG = Path(__file__).resolve().parent / "rJ_fidelity.log"

PARETO = METRICS / "r14_budget_pareto.json"
SOURCES = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "raw": METRICS / "r11_raw" / "convergence",
          "truth_stem": "truth"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "raw": METRICS / "r11_raw" / "designB_convergence",
          "truth_stem": "truth"},
}
BETA = 1.00
TOLERANCE = 1.0e-9        # relative; the quantity is deterministic


def log(msg: str) -> None:
    J.log_line(LOG, f"fidelity {msg}")


def geometry(design: str) -> dict:
    payload = json.loads(SOURCES[design]["rows"].read_text(encoding="utf-8"))
    out = {}
    for row in payload["rows"]:
        g = row.get("design_point", row)
        out[int(row["sobol_index"])] = {
            "hp_km": float(g["hp_km"]), "ha_km": float(g["ha_km"]),
            "adopted_truth_degree": int(row["adopted_truth_degree"])}
    return out


def task(payload: dict) -> dict:
    design = payload["design"]
    index = int(payload["index"])
    try:
        g = payload["geom"]
        budget = payload["budget"]
        adopted = int(g["adopted_truth_degree"])
        model, args = J.model_for(adopted)
        raw = (SOURCES[design]["raw"] / f"sobolA_{index:03d}"
               / f"{SOURCES[design]['truth_stem']}_tighter.npz")
        if not raw.exists():
            return {"ok": False, "design": design, "index": index,
                    "error": f"missing archived reference {raw}"}
        Y, t = J.load_states(raw), J.load_times(raw)
        h_km = (np.linalg.norm(Y[:3], axis=0) - model.r_ref) / 1e3
        _, table = at.atallah_binned_schedule(
            model, J.atallah_g(adopted),
            float(budget["atallah"]["tol_accel_m_s2"]),
            g["hp_km"], g["ha_km"], floor=J.FLOOR, cap=adopted,
            bin_km=J.BIN_KM)
        table = {float(k): int(v) for k, v in table.items()}
        degrees = {"atallah": J.degrees_from_table(table, h_km),
                   "fixed": np.full(len(h_km),
                                    int(budget["fixed"]["degree"]), dtype=int)}
        mine = J.force_defects(t, Y[:3], degrees, adopted, args)
        rec = {"ok": True, "design": design, "index": index,
               "adopted_truth_degree": adopted, "n_epochs": int(len(t))}
        worst = 0.0
        for policy in ("atallah", "fixed"):
            archived = float(budget[policy]["defect"]["defect_rms_m_s2"])
            recomputed = mine[policy]["J_force_rms_m_s2"]
            dev = abs(recomputed / archived - 1.0) if archived else float("nan")
            worst = max(worst, dev)
            rec[policy] = {"archived_defect_rms_m_s2": archived,
                           "recomputed_defect_rms_m_s2": recomputed,
                           "relative_deviation": dev,
                           "J_force_mean_m_s2": mine[policy]["J_force_mean_m_s2"]}
        rec["max_relative_deviation"] = worst
        rec["matches"] = bool(worst <= TOLERANCE)
        rec["mean_defect_ratio_radial_over_constant"] = (
            mine["atallah"]["J_force_mean_m_s2"]
            / mine["fixed"]["J_force_mean_m_s2"])
        rec["archived_rms_ratio_constant_over_radial"] = float(
            budget["ratios"]["R_a_defect_rms"])
        return rec
    except Exception:
        return {"ok": False, "design": design, "index": index,
                "error": traceback.format_exc()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()

    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    tasks = []
    censored = []
    for design in ("A", "B"):
        geom = geometry(design)
        for row in pareto["designs"][design]["rows"]:
            index = int(row["sobol_index"])
            budget = row["budgets"].get(f"beta_{BETA:.2f}")
            if budget is None or budget.get("censored"):
                censored.append({"design": design, "index": index})
                continue
            if index not in geom:
                continue
            tasks.append({"design": design, "index": index,
                          "geom": geom[index], "budget": budget})
    log(f"{len(tasks)} archived cases to recompute, {len(censored)} censored")

    results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fut in as_completed([pool.submit(task, t) for t in tasks]):
            results.append(fut.result())
    good = [r for r in results if r.get("ok")]
    bad = [r for r in results if not r.get("ok")]
    matched = [r for r in good if r["matches"]]
    worst = max((r["max_relative_deviation"] for r in good), default=None)

    payload = {
        "schema": "rJ_fidelity_check_v1", "created_utc": J.utc_now(),
        "purpose": ("recompute the archived deterministic force defect with "
                    "the J-campaign harness, on archived inputs, so that the "
                    "new campaigns can be read against the manuscript's own "
                    "instrument rather than a lookalike"),
        "what_is_held_fixed": ["population", "reference trajectories",
                               "accuracy parameters", "comparator degrees",
                               "binning convention", "reference level"],
        "match_tolerance_relative": TOLERANCE,
        "counts": {"recomputed": len(good), "matched": len(matched),
                   "failed": len(bad), "censored_and_skipped": len(censored)},
        "max_relative_deviation": worst,
        "verdict": ("identical_instrument"
                    if good and len(matched) == len(good)
                    else "instruments_differ"),
        "complete": bool(good and not bad),
        "rows": sorted(good, key=lambda r: (r["design"], r["index"])),
        "failures": bad, "censored": censored,
        "provenance": J.provenance(),
    }
    J.atomic_json(OUT, payload)
    log(f"{payload['verdict']}: {len(matched)}/{len(good)} match, worst "
        f"relative deviation {worst:.3e}" if worst is not None
        else "no cases recomputed")
    return 0 if payload["verdict"] == "identical_instrument" else 1


if __name__ == "__main__":
    raise SystemExit(main())
