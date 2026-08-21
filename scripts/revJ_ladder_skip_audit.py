"""Complete the ladder's skip accounting and cross-check its censoring.

rJ_ladder_primary.json stores counts = {scored, skipped} but only the first 50
skip entries, so 56 of the 106 skips carry no recorded reason. This script
re-derives the reason for every skipped (design, orbit, budget) pair without
integrating anything: the incomplete-archive reason is a filesystem check and
the censor reason is the same degree-schedule computation the ladder ran, which
touches no trajectory.

It also does what the censor itself cannot: for every censored pair it looks up
the archived R14 verdict at that budget, so the record states how many resolved
comparisons -- and how many radial-favoring ones -- the censor removed. That
number is what the manuscript must disclose next to the ladder's tallies.

Output: metrics/rJ_ladder_skip_audit.json  (amendment; the ladder record is
not rewritten).

Usage:  python revJ_ladder_skip_audit.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import revJ_common as J

J.select_field("JGGRX_1800F")
J.install_field()

import rev12_atallah as at                                        # noqa: E402
import revJ_ladder_primary as lad                                 # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "rJ_ladder_skip_audit.json"


def r14_verdict(design: str, beta: float, index: int) -> dict | None:
    """The archived R14 comparison for this pair, at the archived error_tight
    level and in the archived fixed/radial ratio convention."""
    p = METRICS / f"r14_trajectory_{design}_beta_{beta:.2f}.json"
    if not p.exists():
        return None
    rec = json.loads(p.read_text(encoding="utf-8"))
    for row in rec["rows"]:
        if int(row["sobol_index"]) == index:
            c = row["comparison"]
            return {"resolved": c["resolved"],
                    "resolved_winner": c["resolved_winner"],
                    "raw_winner": c["raw_winner"],
                    "rho_budget_fixed_over_radial": c["rho_budget"]}
    return None


def classify(payload: dict) -> tuple[str, dict]:
    """Re-run the ladder's own skip logic for one pair, integration-free."""
    design, index, beta = payload["design"], payload["index"], payload["beta"]
    g, budget = payload["geom"], payload["budget"]
    adopted = int(g["adopted_truth_degree"])
    paths = {(p, lv): lad._raw(design, beta, index, p, lv)
             for p in ("reference", "constant", "radial")
             for lv in ("tight", "tighter")}
    if any(v is None for v in paths.values()):
        return "incomplete archived record at this budget", {}
    model, _ = J.model_for(adopted)
    Y = J.load_states(paths[("reference", "tighter")])
    R = Y[:3]
    h_km = (np.linalg.norm(R, axis=0) - model.r_ref) / 1e3
    _, table = at.atallah_binned_schedule(
        model, J.atallah_g(adopted),
        float(budget["atallah"]["tol_accel_m_s2"]),
        g["hp_km"], g["ha_km"], floor=J.FLOOR, cap=adopted, bin_km=J.BIN_KM)
    table = {float(k): int(v) for k, v in table.items()}
    n_f = int(budget["fixed"]["degree"])
    deg_r = J.degrees_from_table(table, h_km)
    detail = {"fixed_degree": n_f, "radial_max_degree": int(deg_r.max()),
              "adopted_truth_degree": adopted,
              "constant_triggers": bool(n_f >= adopted),
              "radial_triggers": bool(int(deg_r.max()) >= adopted)}
    if n_f >= adopted or int(deg_r.max()) >= adopted:
        return "policy reaches the reference degree; censored", detail
    return "scored", detail


def main() -> int:
    ladder = json.loads((METRICS / "rJ_ladder_primary.json").read_text(
        encoding="utf-8"))
    scored = {(r["design"], r["beta"], r["index"]) for r in ladder["rows"]}

    pareto = json.loads(lad.PARETO.read_text(encoding="utf-8"))
    pairs = []
    for design in ("A", "B"):
        src = json.loads(lad.ARCHIVED[design]["rows"].read_text(
            encoding="utf-8"))
        crit = {int(r["sobol_index"]): r for r in src["rows"]}
        par = {int(r["sobol_index"]): r
               for r in pareto["designs"][design]["rows"]}
        for index, row in crit.items():
            g = row.get("design_point", row)
            geom = {"hp_km": float(g["hp_km"]), "ha_km": float(g["ha_km"]),
                    "adopted_truth_degree": int(row["adopted_truth_degree"])}
            for beta in lad.BUDGETS:
                budget = par[index]["budgets"].get(f"beta_{beta:.2f}")
                if budget is None or budget.get("censored"):
                    continue
                pairs.append({"design": design, "index": index, "beta": beta,
                              "geom": geom, "budget": budget})

    skips, censor_removed = [], {}
    for p in pairs:
        key = (p["design"], p["beta"], p["index"])
        if key in scored:
            continue
        reason, detail = classify(p)
        if reason == "scored":
            reason = ("classified as scorable on re-check but absent from the "
                      "ladder record")
        entry = {"design": p["design"], "beta": p["beta"],
                 "index": p["index"], "reason": reason, **detail}
        if reason.startswith("policy reaches"):
            v = r14_verdict(p["design"], p["beta"], p["index"])
            entry["archived_r14_verdict_error_tight"] = v
            b = censor_removed.setdefault(
                f"beta_{p['beta']:.2f}",
                {"censored": 0, "resolved_in_r14": 0,
                 "radial_wins_removed": 0, "fixed_wins_removed": 0})
            b["censored"] += 1
            if v and v["resolved"]:
                b["resolved_in_r14"] += 1
                if v["resolved_winner"] == "atallah":
                    b["radial_wins_removed"] += 1
                else:
                    b["fixed_wins_removed"] += 1
        skips.append(entry)

    reasons = {}
    for s in skips:
        reasons[s["reason"]] = reasons.get(s["reason"], 0) + 1

    # the ladder's stored (truncated) skip list, for cross-checking
    stored = {(s["design"], s["beta"], s["index"]): s["reason"]
              for s in ladder.get("skipped", [])}
    agree = sum(1 for s in skips
                if stored.get((s["design"], s["beta"], s["index"]))
                == s["reason"])

    payload = {
        "schema": "rJ_ladder_skip_audit_v1", "created_utc": J.utc_now(),
        "amends": {"file": "rJ_ladder_primary.json",
                   "counts_in_record": ladder["counts"],
                   "skip_entries_in_record": len(ladder.get("skipped", [])),
                   "note": ("the record stored only the first 50 skip entries; "
                            "this audit derives all of them, integration-free, "
                            "and does not rewrite the record")},
        "censor_rule": ("skip when fixed_degree >= adopted reference degree OR "
                        "max radial scheduled degree >= adopted reference "
                        "degree. In this archive the constant branch never "
                        "fires (fixed_degree <= 132 against reference 300); "
                        "every censor is triggered by the radial schedule "
                        "reaching the reference degree"),
        "totals": {"pairs": len(pairs), "scored": len(scored),
                   "skipped": len(skips)},
        "reasons": reasons,
        "censor_removed_by_budget": censor_removed,
        "stored_skip_entries_agreeing": agree,
        "skips": sorted(skips, key=lambda s: (s["beta"], s["design"],
                                              s["index"])),
        "provenance": J.provenance()}
    J.atomic_json(OUT, payload)

    print(f"[written] {OUT.name}")
    print(f"  pairs {len(pairs)} = scored {len(scored)} + skipped {len(skips)}")
    for k, v in sorted(reasons.items()):
        print(f"  reason: {k} -> {v}")
    ct = any(s.get("constant_triggers") for s in skips)
    print(f"  constant branch ever triggers censor: {ct}")
    for b, v in sorted(censor_removed.items()):
        print(f"  {b}: censored {v['censored']}, resolved in R14 "
              f"{v['resolved_in_r14']}, radial wins removed "
              f"{v['radial_wins_removed']}, fixed wins removed "
              f"{v['fixed_wins_removed']}")
    print(f"  stored 50 entries re-derived identically: {agree}/50")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
