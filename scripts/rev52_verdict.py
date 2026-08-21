#!/usr/bin/env python3
"""R52: score block B's cap-lifted ladder, against its own parent and block A.

Written while the control's base is still propagating and before any block-B
ladder comparison exists, for the reason rev50_verdict.py and rev51_verdict.py
were: what counts as "the blocks replicate" has to be arithmetic fixed in
advance rather than a judgement made once the numbers are in view.

Two comparisons, in this order:

  1. block B uncapped against block B capped, level by level. This is the same
     reading rev51_verdict.py makes of block A, and it is what says whether the
     ceiling was doing work on this block.

  2. block B uncapped against block A uncapped, level by level. This is the
     replication question, and it is the one the registered outcome classes
     P/Q/R/S are about.

The registered outcome classes, from r52_preregistration.json:

  P_replicates              same per-level verdicts as block A, and the score
                            still orders in radial span at every budget that
                            decides
  Q_replicates_ordering_only  verdicts agree with block A but the wide-level
                            magnitudes do not
  R_blocks_disagree         a level changes hands in block B where it did not
                            in block A, or the two blocks place the turnover
                            between different level pairs at the same budget
  S_undecided               too few comparisons resolve per level to order
                            block B's levels

The registration names those classes; it does not fix the numeric threshold
that separates P from Q, exactly as R51's registration left MAGNITUDE_DROP to
its scoring script. That threshold is BLOCK_MAGNITUDE_FACTOR below, it is set
here before any block-B ladder exists, and it is published in the verdict
record so a reader can see what it was rather than infer it.

Usage:
    python rev52_verdict.py
    python rev52_verdict.py --beta 0.75
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

CAPPED_KEY, UNCAPPED_KEY = "RS2", "RS2U"
CAPPED_DESIGN = "r50_span_ladder_b_design_frozen.json"
SIBLING_VERDICT = METRICS / "r51_verdict.json"
BETAS = ["1.00", "0.75", "0.62", "0.50"]

# A factor of ten between the two blocks' median ratios at the same wide cell
# is a magnitude the blocks do not share. Ten mirrors the order-of-magnitude
# criterion R51 used for its own X/Y split, so the two campaigns are not scored
# on different scales.
BLOCK_MAGNITUDE_FACTOR = 10.0


def levels() -> dict:
    d = json.loads((METRICS / CAPPED_DESIGN).read_text(encoding="utf-8"))
    return {o["sobol_index"]: o["apolune_level_km"] for o in d["orbits"]}


def cells(key: str, beta: str, index: dict) -> dict | None:
    """Level cells for one design at one budget, read as rev50_verdict reads."""
    p = METRICS / f"r14_trajectory_{key}_beta_{beta}.json"
    if not p.exists():
        return None
    out: dict[float, dict] = {}
    for r in json.loads(p.read_text(encoding="utf-8"))["rows"]:
        c = r["comparison"]
        lvl = index[r["sobol_index"]]
        e = out.setdefault(lvl, {"orbits": 0, "resolved": 0, "radial": 0,
                                 "fixed": 0, "rho": []})
        e["orbits"] += 1
        if c.get("rho_budget") is not None:
            e["rho"].append(float(c["rho_budget"]))
        if c.get("resolved"):
            e["resolved"] += 1
            won = (c.get("resolved_winner") or c.get("raw_winner")) == "atallah"
            e["radial" if won else "fixed"] += 1
    for e in out.values():
        e["median_rho"] = median(e["rho"]) if e["rho"] else None
        e.pop("rho")
        e["score"] = ((e["radial"] - e["fixed"]) / e["resolved"]
                      if e["resolved"] else None)
        e["verdict"] = ("undecided" if not e["resolved"]
                        else "radial" if e["radial"] > e["fixed"]
                        else "constant" if e["fixed"] > e["radial"]
                        else "split")
    return dict(sorted(out.items()))


def kendall_tau(pairs) -> float | None:
    if len(pairs) < 3:
        return None
    con = dis = 0
    for (x1, y1), (x2, y2) in combinations(pairs, 2):
        s = (x2 - x1) * (y2 - y1)
        con += s > 0
        dis += s < 0
    n = con + dis
    return (con - dis) / n if n else 0.0


def turnover(verdicts: dict) -> tuple | None:
    """The level pair the verdict first flips across, going up in span."""
    lv = sorted(verdicts)
    for a, b in zip(lv, lv[1:]):
        if verdicts[a] != verdicts[b] and "undecided" not in (verdicts[a],
                                                              verdicts[b]):
            return (a, b)
    return None


def sibling() -> dict:
    """Block A's ceiling-free cells, keyed (beta, level)."""
    if not SIBLING_VERDICT.exists():
        return {}
    d = json.loads(SIBLING_VERDICT.read_text(encoding="utf-8"))
    out = {}
    for beta, blk in d["by_budget"].items():
        for lvl, cell in blk["levels"].items():
            out[(beta, float(lvl))] = cell["uncapped"]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", default=None)
    a = ap.parse_args()

    prereg_p = METRICS / "r52_preregistration.json"
    if not prereg_p.exists():
        print("r52_preregistration.json missing; the control is not registered")
        return 0
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))
    index = levels()
    block_a = sibling()
    betas = [a.beta] if a.beta else BETAS

    report, changed_hands, taus, drops = {}, [], [], []
    disagree, magnitude_gap, turnover_gap = [], [], []

    for beta in betas:
        cap = cells(CAPPED_KEY, beta, index)
        unc = cells(UNCAPPED_KEY, beta, index)
        if cap is None or unc is None:
            continue
        rows = {}
        for lvl in sorted(unc):
            c, u = cap.get(lvl), unc[lvl]
            row = {"capped": c, "uncapped": u,
                   "verdict_changed": bool(c and c["verdict"] != u["verdict"])}
            if row["verdict_changed"]:
                changed_hands.append({"beta": beta, "level_km": lvl,
                                      "from": c["verdict"], "to": u["verdict"]})
            if c and c["median_rho"] and u["median_rho"]:
                row["rho_factor"] = c["median_rho"] / u["median_rho"]
                if lvl >= 1200.0:
                    drops.append(row["rho_factor"])

            # the replication comparison
            peer = block_a.get((beta, lvl))
            if peer:
                row["block_a_uncapped"] = peer
                row["agrees_with_block_a"] = peer["verdict"] == u["verdict"]
                if not row["agrees_with_block_a"]:
                    disagree.append({"beta": beta, "level_km": lvl,
                                     "block_a": peer["verdict"],
                                     "block_b": u["verdict"]})
                if (lvl >= 1200.0 and peer["median_rho"] and u["median_rho"]):
                    f = max(peer["median_rho"], u["median_rho"]) / min(
                        peer["median_rho"], u["median_rho"])
                    row["block_ratio_factor"] = f
                    if f >= BLOCK_MAGNITUDE_FACTOR:
                        magnitude_gap.append({"beta": beta, "level_km": lvl,
                                              "factor": f})
            rows[lvl] = row

        deciding = [(lvl, r["uncapped"]["score"]) for lvl, r in rows.items()
                    if r["uncapped"]["resolved"]]
        tau = kendall_tau(deciding) if len(deciding) >= 3 else None
        if tau is not None:
            taus.append(tau)

        t_b = turnover({lvl: r["uncapped"]["verdict"] for lvl, r in rows.items()})
        t_a = None
        if block_a:
            av = {lvl: block_a[(beta, lvl)]["verdict"] for lvl in rows
                  if (beta, lvl) in block_a}
            t_a = turnover(av) if av else None
        if t_a != t_b:
            turnover_gap.append({"beta": beta, "block_a": t_a, "block_b": t_b})

        report[beta] = {"levels": rows, "kendall_tau_uncapped": tau,
                        "levels_deciding": len(deciding),
                        "turnover_block_b": t_b, "turnover_block_a": t_a}

    if not report:
        print("no block-B uncapped ladder on disk yet; nothing to read")
        return 0

    # P carries two clauses, not one: the verdicts must match block A *and* the
    # score must order in radial span at every budget that decides. An earlier
    # version of this file tested only the first and would have returned P with
    # the second failing, which is the defect R51's registration was faulted
    # for. Both are tested here and the second is reported separately, because
    # a clause that fails silently is worse than one that fails loudly.
    ordering_clause = [b for b, blk in report.items()
                       if blk["kendall_tau_uncapped"] is not None
                       and abs(blk["kendall_tau_uncapped"]) != 1.0]

    if not taus:
        outcome = "S_undecided"
    elif disagree or turnover_gap:
        outcome = "R_blocks_disagree"
    elif magnitude_gap:
        outcome = "Q_replicates_ordering_only"
    else:
        outcome = "P_replicates"

    payload = {
        "schema": "r52_verdict_v1",
        "preregistration_sha256": prereg["preregistration_sha256"],
        "rule": prereg["verdict_rule"],
        "budgets_read": sorted(report),
        "by_budget": report,
        "levels_changing_hands_vs_capped_parent": changed_hands,
        "levels_disagreeing_with_block_a": disagree,
        "turnover_disagreements": turnover_gap,
        "wide_cells_beyond_block_magnitude_factor": magnitude_gap,
        "block_magnitude_factor": BLOCK_MAGNITUDE_FACTOR,
        "kendall_tau_uncapped": taus,
        "median_rho_factor_wide_levels": median(drops) if drops else None,
        "outcome": outcome,
        "outcome_text": prereg["outcomes"][outcome],
        "p_ordering_clause_holds": not ordering_clause,
        "budgets_failing_p_ordering_clause": ordering_clause,
        "registered_classes_cover_result": not (outcome == "P_replicates"
                                                and ordering_clause),
    }
    (METRICS / "r52_verdict.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    for beta, blk in report.items():
        print(f"\nbeta {beta}   (block B: capped -> uncapped | block A uncapped)")
        for lvl, r in blk["levels"].items():
            u = r["uncapped"]
            peer = r.get("block_a_uncapped")
            peer_txt = ("--" if not peer else
                        f"{peer['verdict']:<8} {peer['radial']:2d}-"
                        f"{peer['fixed']:<2d}")
            flag = "" if r.get("agrees_with_block_a", True) else "   <-- differs"
            print(f"  ha {lvl:6.0f} km  B {u['verdict']:<8} "
                  f"{u['radial']:2d}-{u['fixed']:<2d} of {u['resolved']:2d}"
                  f"   | A {peer_txt}{flag}")
        print(f"  tau(B uncapped) = {blk['kendall_tau_uncapped']}   "
              f"turnover A {blk['turnover_block_a']} vs B "
              f"{blk['turnover_block_b']}")
    print(f"\nlevels disagreeing with block A: {disagree or 'none'}")
    print(f"turnover disagreements: {turnover_gap or 'none'}")
    print(f"OUTCOME: {outcome}")
    print(f"  {payload['outcome_text']}")
    if outcome == "P_replicates" and ordering_clause:
        print("\n  *** DECLARED DEPARTURE ***")
        print("  P has two clauses. The verdicts clause holds at every cell.")
        print(f"  The span-ordering clause does not, at beta "
              f"{', '.join(ordering_clause)}.")
        print("  No registered class covers a result that replicates the")
        print("  verdicts, the turnover and the magnitudes while failing the")
        print("  ordering statistic, so the registration does not partition")
        print("  this outcome space and the gap is recorded rather than")
        print("  resolved by choosing a class after the fact.")
    print("[written] r52_verdict.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
