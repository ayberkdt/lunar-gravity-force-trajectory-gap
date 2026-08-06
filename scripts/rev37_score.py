"""R37 scoring: re-derive the enlarged panel's verdict from the stored rows.

The driver writes every orbit it solves, with that orbit's measured ratio, its
numerical envelopes and the resolution flag the R14 campaign computed. What it
does not do correctly is tally them: it scores sign agreement against every
measured ratio, including ratios the resolution rule leaves undecided. An
undecided measurement cannot confirm a prediction and cannot contradict one.

This script applies the manuscript's own convention instead --- sign agreement
among resolved comparisons only, unresolved ones counted separately as carrying
no verdict --- and reports both tallies. See metrics/r37_scoring_amendment.json
for why that is a correction rather than a choice.

Nothing is recomputed: this reads the archive the run has written so far, so it
can be run against a partial record while the run is still going.

Usage:  python rev37_score.py [--json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RECORD = METRICS / "r37_variational_extension.json"
OUT = METRICS / "r37_panel_verdict.json"

CAL_BAND = (0.94, 1.04)


def classify(r):
    m = r.get("measured") or {}
    pred = r.get("predicted_ratio_fixed_over_atallah")
    if pred is None or not m:
        return None
    meas = m["fixed_budget"] / m["atallah_budget"]
    cal = r.get("calibration_ratio_fixed")
    return {
        "design": r["design"], "sobol_index": r["sobol_index"],
        "hp_km": r["hp_km"], "level": r.get("level"),
        "reused_from_r14": bool(r.get("reused_from_r14")),
        "predicted": pred, "measured": meas,
        "resolved": bool(m.get("resolved")),
        "calibration_ratio_fixed": cal,
        "calibrated": (cal is not None and CAL_BAND[0] <= cal <= CAL_BAND[1]),
        "sign_agrees": (pred < 1.0) == (meas < 1.0),
        "predicted_favors_radial": pred > 1.0,
        "measured_favors_radial": meas > 1.0,
    }


def levels_complete(rec, cells):
    """Highest declared level whose every member is present in the record."""
    prereg = json.loads(
        (METRICS / "r37_preregistration.json").read_text(encoding="utf-8"))
    levels = prereg["selection_rule"]["levels_per_design"]
    pareto = json.loads(
        (METRICS / "r14_budget_pareto.json").read_text(encoding="utf-8"))
    have = {(c["design"], c["sobol_index"]) for c in cells}
    best, members_of = None, {}
    for n_per in levels:
        members = set()
        for d in ("A", "B"):
            pr = [r for r in pareto["designs"][d]["rows"]
                  if not r["budgets"]["beta_1.00"]["censored"]]
            pr.sort(key=lambda r: r["hp_km"])
            idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
            members |= {(d, int(pr[i]["sobol_index"])) for i in idx}
        members_of[n_per] = members
        if members <= have:
            best = n_per
    nxt = next((n for n in levels if n > (best or 0)), None)
    return {
        "highest_complete_per_design": best,
        "highest_complete_orbits": 2 * best if best else 0,
        "next_level_per_design": nxt,
        "next_level_orbits": 2 * nxt if nxt else None,
        "next_level_present": (len(members_of[nxt] & have) if nxt else None),
    }


def level_members(n_per):
    pareto = json.loads(
        (METRICS / "r14_budget_pareto.json").read_text(encoding="utf-8"))
    members = set()
    for d in ("A", "B"):
        pr = [r for r in pareto["designs"][d]["rows"]
              if not r["budgets"]["beta_1.00"]["censored"]]
        pr.sort(key=lambda r: r["hp_km"])
        idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
        members |= {(d, int(pr[i]["sobol_index"])) for i in idx}
    return members


def tally(cells):
    res = [c for c in cells if c["resolved"]]
    cal = [c["calibration_ratio_fixed"] for c in res
           if c["calibration_ratio_fixed"] is not None]
    return {
        "orbits": len(cells),
        "resolved": len(res),
        "unresolved": len(cells) - len(res),
        "sign_agreement": sum(1 for c in res if c["sign_agrees"]),
        "sign_disagreement": sum(1 for c in res if not c["sign_agrees"]),
        # the raw tally the driver produced before the scoring amendment, kept
        # so both readings are archived and neither has to be taken on trust
        "sign_agreement_raw": sum(1 for c in cells if c["sign_agrees"]),
        "sign_disagreement_raw": sum(1 for c in cells if not c["sign_agrees"]),
        # the instrument's own verdict does not depend on whether the
        # propagated comparison resolved, so it is counted over every orbit,
        # matching how the ladder's other rungs are counted
        "predicted_favors_radial_all": sum(c["predicted_favors_radial"]
                                           for c in cells),
        "predicted_favors_radial": sum(c["predicted_favors_radial"] for c in res),
        "measured_favors_radial": sum(c["measured_favors_radial"] for c in res),
        "calibration_median": float(np.median(cal)) if cal else None,
        "calibration_min": float(np.min(cal)) if cal else None,
        "calibration_max": float(np.max(cal)) if cal else None,
        "calibrated_within_band": sum(1 for c in res if c["calibrated"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    rec = json.loads(RECORD.read_text(encoding="utf-8"))
    cells = [c for c in (classify(r) for r in rec["rows"]) if c]
    res = [c for c in cells if c["resolved"]]
    unres = [c for c in cells if not c["resolved"]]
    agree = [c for c in res if c["sign_agrees"]]
    disagree = [c for c in res if not c["sign_agrees"]]
    uncal = [c for c in cells if c["calibration_ratio_fixed"] is not None
             and not c["calibrated"]]
    cal_res = [c["calibration_ratio_fixed"] for c in res
               if c["calibration_ratio_fixed"] is not None]

    lv = levels_complete(rec, cells)
    panel_keys = (level_members(lv["highest_complete_per_design"])
                  if lv["highest_complete_per_design"] else set())
    for c in cells:
        c["in_panel"] = (c["design"], c["sobol_index"]) in panel_keys
    panel_cells = ([c for c in cells
                    if (c["design"], c["sobol_index"])
                    in level_members(lv["highest_complete_per_design"])]
                   if lv["highest_complete_per_design"] else [])

    out = {
        "schema": "r37_panel_verdict_v1",
        "scored_under": "r37_scoring_amendment.json",
        "orbits_in_record": len(cells),
        "levels": lv,
        "panel": {
            "definition": ("the highest fully completed level of the "
                           "registered chain; orbits beyond it are in "
                           "'record' but not in the panel"),
            "orbits": lv["highest_complete_orbits"],
            **tally(panel_cells)},
        "record": {"definition": "every orbit solved, complete level or not",
                   **tally(cells)},
        "resolved": {
            "n": len(res),
            "sign_agreement": len(agree),
            "sign_disagreement": len(disagree),
            "predicted_favors_radial": sum(c["predicted_favors_radial"]
                                           for c in res),
            "measured_favors_radial": sum(c["measured_favors_radial"]
                                          for c in res),
            "calibration_median": (float(np.median(cal_res))
                                   if cal_res else None),
            "calibration_min": float(np.min(cal_res)) if cal_res else None,
            "calibration_max": float(np.max(cal_res)) if cal_res else None,
        },
        "unresolved": {
            "n": len(unres),
            "note": ("carry no verdict: the propagated comparison lies inside "
                     "the numerical envelope, so it can neither confirm nor "
                     "contradict a prediction"),
            "orbits": [{k: c[k] for k in
                        ("design", "sobol_index", "hp_km", "predicted",
                         "measured", "calibration_ratio_fixed", "in_panel")}
                       for c in sorted(unres, key=lambda c: c["hp_km"])],
        },
        "sign_disagreements_resolved": [
            {k: c[k] for k in ("design", "sobol_index", "hp_km", "predicted",
                               "measured", "calibration_ratio_fixed",
                               "calibrated")}
            for c in sorted(disagree, key=lambda c: c["hp_km"])],
        "outside_calibration_band": [
            {k: c[k] for k in ("design", "sobol_index", "hp_km",
                               "calibration_ratio_fixed", "resolved",
                               "in_panel")}
            for c in sorted(uncal, key=lambda c: c["hp_km"])],
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    if a.json:
        print(json.dumps(out, indent=2))
        return 0

    L = out["levels"]
    r = out["resolved"]
    print(f"orbits in record: {out['orbits_in_record']}")
    print(f"highest complete level: {L['highest_complete_orbits']} orbits "
          f"(level {L['highest_complete_per_design']} per design)")
    if L["next_level_orbits"]:
        print(f"  next level {L['next_level_orbits']} orbits: "
              f"{L['next_level_present']} present")
    print()
    print(f"resolved comparisons: {r['n']}")
    print(f"  sign agreement:    {r['sign_agreement']}/{r['n']}")
    print(f"  sign disagreement: {r['sign_disagreement']}")
    print(f"  predicted favors radial {r['predicted_favors_radial']}, "
          f"measured favors radial {r['measured_favors_radial']}")
    if r["calibration_median"] is not None:
        print(f"  calibration channel median {r['calibration_median']:.4f} "
              f"[{r['calibration_min']:.3f}, {r['calibration_max']:.3f}]")
    print(f"unresolved (no verdict): {out['unresolved']['n']}")
    for c in out["unresolved"]["orbits"]:
        print(f"    {c['design']}{c['sobol_index']:03d} hp={c['hp_km']:6.1f} "
              f"pred={c['predicted']:8.4f} meas={c['measured']:8.4f}")
    if out["sign_disagreements_resolved"]:
        print("SIGN DISAGREEMENTS on resolved comparisons:")
        for c in out["sign_disagreements_resolved"]:
            print(f"    {c['design']}{c['sobol_index']:03d} "
                  f"hp={c['hp_km']:6.1f} pred={c['predicted']:8.4f} "
                  f"meas={c['measured']:8.4f} "
                  f"calib={c['calibration_ratio_fixed']:.3f} "
                  f"{'calibrated' if c['calibrated'] else 'NOT calibrated'}")
    if out["outside_calibration_band"]:
        print(f"outside the {CAL_BAND[0]}-{CAL_BAND[1]} calibration band: "
              f"{len(out['outside_calibration_band'])}")
        for c in out["outside_calibration_band"]:
            print(f"    {c['design']}{c['sobol_index']:03d} "
                  f"hp={c['hp_km']:6.1f} "
                  f"calib={c['calibration_ratio_fixed']:.3f} "
                  f"resolved={c['resolved']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
