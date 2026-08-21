"""Band admissibility for R68: what the out-of-band cells are, and what they do.

The registration prescribes a 0.95-1.05 time-match band and says explicitly
that a cell still outside it after refinement "keeps its nearest integer match
and is flagged; misses are counted in the record and in the manuscript
sentence, not absorbed". Misses are therefore *not* an exclusion criterion,
and the primary counts include them. This file exists so that statement can be
checked rather than trusted: it lists every out-of-band cell, how far out and
in which direction it sits, what verdict it carries, and what the tally would
be without it.

The direction matters more than the count. A ratio above unity gives the
constant comparator more machine time than the member spent, so an excess of
high-side misses biases a comparison towards the constant degree and against
the variable member; a deficit does the reverse. Reporting the sensitivity
tally alone would hide which way the residual mismatch leans.

Usage:  python rev68_band_audit.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rev68_timing_full as r68            # noqa: E402

ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r68_band_admissibility.json"


def member_name(arm: str) -> str:
    return "radial" if arm == "endpoint" else "interior"


def tally(rows, member: str) -> dict:
    res = [r for r in rows if r["resolved"]]
    return {"resolved": len(res),
            f"{member}_wins": sum(1 for r in res if r["winner"] == member),
            "constant_wins": sum(1 for r in res if r["winner"] == "constant"),
            "unresolved": len(rows) - len(res)}


def main() -> int:
    prereg = json.loads((METRICS / "r68_preregistration.json"
                         ).read_text(encoding="utf-8"))
    lo, hi = r68.BAND_LO, r68.BAND_HI
    payload = {
        "schema": "r68_band_admissibility_v1",
        "created_utc": None,
        "registered_band": [lo, hi],
        "registered_rule": prereg["locked_choices"]["timing_match_miss"],
        "misses_are_an_exclusion_criterion": False,
        "how_to_read": (
            "the primary counts are the registered ones and include every "
            "scored cell, in band or not; the without-misses tally below is a "
            "sensitivity check and not the reported verdict. A ratio above "
            "unity means the constant comparator was given more machine time "
            "than the member spent, which favours the constant degree."),
        "arms": {},
    }
    for arm in sorted(r68.ARMS):
        p = r68.out_path(arm)
        if not p.exists():
            continue
        rec = json.loads(p.read_text(encoding="utf-8"))
        payload["created_utc"] = rec["created_utc"]
        member = member_name(arm)
        arm_out = {"member_k": rec["member_k"], "member": member,
                   "by_design": {}, "misses": []}
        for design in r68.DESIGNS:
            rows = [r for r in rec["rows"] if r["design"] == design]
            miss = [r for r in rows if r["timing_match_miss"]]
            keep = [r for r in rows if not r["timing_match_miss"]]
            ratios = np.array([r["achieved_time_ratio"] for r in rows])
            arm_out["by_design"][design] = {
                "cells": len(rows),
                "in_band": len(keep),
                "out_of_band": len(miss),
                "out_high": sum(1 for r in miss
                                if r["achieved_time_ratio"] > hi),
                "out_low": sum(1 for r in miss
                               if r["achieved_time_ratio"] < lo),
                "worst_ratio": (float(max(ratios, key=lambda x: abs(x - 1.0)))
                                if len(ratios) else None),
                "median_ratio": float(np.median(ratios)) if len(ratios) else None,
                "mean_abs_deviation_from_unity":
                    float(np.mean(np.abs(ratios - 1.0))) if len(ratios) else None,
                "primary_registered": tally(rows, member),
                "sensitivity_without_misses": tally(keep, member),
                "verdict_changes": (
                    tally(rows, member)[f"{member}_wins"] >
                    tally(rows, member)["constant_wins"]) != (
                    tally(keep, member)[f"{member}_wins"] >
                    tally(keep, member)["constant_wins"]),
            }
            for r in sorted(miss, key=lambda r: -abs(
                    r["achieved_time_ratio"] - 1.0)):
                arm_out["misses"].append({
                    "design": design, "sobol_index": r["sobol_index"],
                    "hp_km": r["hp_km"], "ha_km": r["ha_km"],
                    "comparator_degree": r["comparator_degree"],
                    "achieved_time_ratio": r["achieved_time_ratio"],
                    "side": ("high, favouring the constant degree"
                             if r["achieved_time_ratio"] > hi
                             else "low, favouring the member"),
                    "refinement_passes": r["refinement_passes"],
                    "resolved": r["resolved"], "winner": r["winner"]})
        payload["arms"][arm] = arm_out

    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    for arm, a in payload["arms"].items():
        print(f"\n{arm} (member {a['member']}, k={a['member_k']})")
        for design, d in a["by_design"].items():
            pr, se = d["primary_registered"], d["sensitivity_without_misses"]
            print(f"  {design}: {d['out_of_band']}/{d['cells']} out of band "
                  f"({d['out_high']} high, {d['out_low']} low), "
                  f"median ratio {d['median_ratio']:.3f}")
            print(f"     primary     {pr[f'{a['member']}_wins']}-"
                  f"{pr['constant_wins']} of {pr['resolved']} resolved")
            print(f"     sensitivity {se[f'{a['member']}_wins']}-"
                  f"{se['constant_wins']} of {se['resolved']} resolved"
                  f"   verdict changes: {d['verdict_changes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
