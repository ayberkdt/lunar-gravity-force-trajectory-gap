"""R50: read the paired span ladder level by level, under a rule fixed first.

Written while the first ladder was still propagating and before any comparison
of this population existed, for the reason rev30_verdict exists: the sentence
that reaches the manuscript should be the one the tallies support, not the one
that reads well. What counts as a dependence is therefore arithmetic here rather
than judgement later.

Three readouts, all per apolune level and per budget:

  primary    the budget-calibrated radial endpoint against its equal-budget
             constant degree, from the R14 record. This is the comparison the
             discussion's geometry sentence is about.
  secondary  the interior member of the span family against its work-matched
             constant degree, from the R19 record. This is the comparison the
             budget-axis sentence is about.
  paired     within each identity, the levels at which each verdict changes,
             so the dependence is read inside identities and not only across
             level medians.

The score of a level is s = (wins - losses) / resolved, in [-1, 1], positive
when the varying-degree policy wins. A level decides when it resolves at least
one comparison.

Outcome rule, applied to the primary readout:

  W_undecided        fewer than three levels decide at every budget
  V_no_dependence    no budget changes its verdict across levels
  T_span_dependence  some budget changes verdict, and at every budget that
                     decides three or more levels the score is perfectly
                     ordered in radial span (Kendall tau = +-1)
  U_threshold        some budget changes verdict, but the ordering is not
                     perfect: a flip rather than a gradient

Blocks A and B are pooled level by level, which the registration allows because
they are the same design at the same levels, and each block is also reported on
its own so the pooling can be checked rather than trusted.

Usage:
    python rev50_verdict.py
    python rev50_verdict.py --beta 0.75
    python rev50_verdict.py --all
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

import population_registry as registry                        # noqa: E402

REG = "r50"
BETAS = [1.00, 0.75, 0.62, 0.50, 1.25, 1.50]


def levels_of(name: str, spec: dict) -> dict:
    """sobol_index -> (identity, apolune level, radial span) from the design."""
    d = json.loads((METRICS / spec["file"]).read_text(encoding="utf-8"))
    if d["design_sha256"] != spec["design_sha256"]:
        raise SystemExit(f"{spec['file']} does not match the registered hash")
    return {o["sobol_index"]: {"identity": o["identity_index"],
                               "level_km": o["apolune_level_km"],
                               "span_km": o["radial_span_km"],
                               "in_box": o["apolune_level_inside_factor_box"]}
            for o in d["orbits"]}


def primary_rows(key: str, beta: float) -> list[dict] | None:
    """Radial endpoint against equal-budget constant, from the R14 record."""
    p = METRICS / f"r14_trajectory_{key}_beta_{beta:.2f}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    out = []
    for r in d["rows"]:
        c = r["comparison"]
        winner = c.get("resolved_winner") or c.get("raw_winner")
        out.append({"sobol_index": r["sobol_index"],
                    "resolved": bool(c["resolved"]),
                    "varying_wins": winner == "atallah",
                    "rho": c.get("rho_budget")})
    return out


def secondary_rows(key: str, beta: float) -> list[dict] | None:
    """Interior member against work-matched constant, from the R19 record."""
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    if "summary" not in d:
        return None
    return [{"sobol_index": r["sobol_index"],
             "resolved": bool(r["resolved"]),
             "varying_wins": r["winner"] == "interior",
             "rho": r.get("rho_workmatched")}
            for r in d["rows"]]


def group_by_level(rows: list[dict], index: dict) -> dict:
    cells: dict[float, dict] = {}
    for r in rows:
        meta = index[r["sobol_index"]]
        c = cells.setdefault(meta["level_km"], {
            "level_km": meta["level_km"], "inside_factor_box": meta["in_box"],
            "orbits": 0, "resolved": 0, "varying_wins": 0, "fixed_wins": 0,
            "rho": []})
        c["orbits"] += 1
        if r["rho"] is not None:
            c["rho"].append(float(r["rho"]))
        if not r["resolved"]:
            continue
        c["resolved"] += 1
        if r["varying_wins"]:
            c["varying_wins"] += 1
        else:
            c["fixed_wins"] += 1
    for c in cells.values():
        c["unresolved"] = c["orbits"] - c["resolved"]
        c["median_rho"] = median(c["rho"]) if c["rho"] else None
        c.pop("rho")
        c["score"] = ((c["varying_wins"] - c["fixed_wins"]) / c["resolved"]
                      if c["resolved"] else None)
        c["verdict"] = ("undecided" if not c["resolved"]
                        else "varying" if c["varying_wins"] > c["fixed_wins"]
                        else "constant" if c["fixed_wins"] > c["varying_wins"]
                        else "split")
    return dict(sorted(cells.items()))


def kendall_tau(pairs: list[tuple[float, float]]) -> float | None:
    """Rank correlation of score against radial span, on four points."""
    if len(pairs) < 3:
        return None
    con = dis = 0
    for (x1, y1), (x2, y2) in combinations(pairs, 2):
        s = (x2 - x1) * (y2 - y1)
        if s > 0:
            con += 1
        elif s < 0:
            dis += 1
    n = con + dis
    return (con - dis) / n if n else 0.0


def paired_turning_points(rows: list[dict], index: dict) -> dict:
    """Within an identity, the lowest level whose comparison the varying
    degree wins, and whether the identity is ordered once it starts winning."""
    by_identity: dict[int, list[tuple[float, dict]]] = {}
    for r in rows:
        meta = index[r["sobol_index"]]
        by_identity.setdefault(meta["identity"], []).append((meta["level_km"], r))
    turning, ordered, undecided = {}, 0, 0
    for identity, members in by_identity.items():
        members.sort(key=lambda m: m[0])
        seq = [(lvl, r) for lvl, r in members if r["resolved"]]
        if not seq:
            undecided += 1
            turning[identity] = None
            continue
        wins = [lvl for lvl, r in seq if r["varying_wins"]]
        turning[identity] = min(wins) if wins else None
        # ordered means: once the varying degree starts winning it keeps
        # winning at every wider level this identity resolves.
        flags = [r["varying_wins"] for _, r in seq]
        first = flags.index(True) if True in flags else len(flags)
        if all(flags[first:]):
            ordered += 1
    counted = len(by_identity) - undecided
    return {"turning_level_km": turning,
            "identities": len(by_identity),
            "identities_with_no_resolved_member": undecided,
            "identities_ordered_after_first_win": ordered,
            "identities_counted": counted,
            "fraction_ordered": (ordered / counted) if counted else None}


def outcome_of(per_beta: dict) -> tuple[str, dict]:
    deciding = {b: v for b, v in per_beta.items()
                if sum(1 for c in v["levels"].values() if c["resolved"]) >= 3}
    evidence = {
        "budgets_present": sorted(per_beta),
        "budgets_deciding_three_or_more_levels": sorted(deciding),
        "verdict_changes_across_levels": {},
        "kendall_tau_score_vs_span": {},
    }
    if not deciding:
        return "W_undecided", evidence

    changed, taus = False, []
    for b, v in deciding.items():
        verdicts = {c["verdict"] for c in v["levels"].values() if c["resolved"]}
        changes = len(verdicts - {"undecided"}) > 1
        # the budget keys of this report are the strings the tables are keyed
        # by, not floats; formatting them as floats raised here on the first
        # complete grid.
        evidence["verdict_changes_across_levels"][str(b)] = changes
        changed = changed or changes
        tau = kendall_tau([(c["level_km"], c["score"])
                           for c in v["levels"].values() if c["resolved"]])
        evidence["kendall_tau_score_vs_span"][str(b)] = tau
        if tau is not None:
            taus.append(tau)
    if not changed:
        return "V_no_dependence", evidence
    if taus and all(abs(t) == 1.0 for t in taus):
        return "T_span_dependence", evidence
    return "U_threshold", evidence


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=1.00)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    prereg = registry.registration(REG)
    pops = registry.populations(REG)
    index = {name: levels_of(name, spec) for name, spec in pops.items()}
    betas = BETAS if a.all else [a.beta]

    report: dict = {"primary": {}, "secondary": {}, "blocks": {}}
    for beta in betas:
        for readout, reader in (("primary", primary_rows),
                                ("secondary", secondary_rows)):
            pooled, blocks = [], {}
            for name, spec in pops.items():
                rows = reader(spec["design_key"], beta)
                if rows is None:
                    continue
                blocks[name] = {
                    "levels": group_by_level(rows, index[name]),
                    "paired": paired_turning_points(rows, index[name])}
                pooled += [dict(r, sobol_index=r["sobol_index"],
                                _block=name) for r in rows]
            if not blocks:
                continue
            merged: dict[int, dict] = {}
            merged_index: dict[int, dict] = {}
            for n, (name, _) in enumerate(pops.items()):
                if name not in blocks:
                    continue
                for r in (reader(pops[name]["design_key"], beta) or []):
                    uid = n * 1000 + r["sobol_index"]
                    merged[uid] = dict(r, sobol_index=uid)
                    m = index[name][r["sobol_index"]]
                    merged_index[uid] = dict(m, identity=n * 1000 + m["identity"])
            report[readout][f"{beta:.2f}"] = {
                "levels": group_by_level(list(merged.values()), merged_index),
                "paired": paired_turning_points(list(merged.values()),
                                                merged_index),
                "blocks_pooled": sorted(blocks)}
            report["blocks"].setdefault(f"{beta:.2f}", {})[readout] = blocks

    if not report["primary"]:
        print("no R14 record for this population yet; nothing to read")
        return 0

    outcome, evidence = outcome_of(report["primary"])
    payload = {
        "schema": "r50_verdict_v1",
        "registry": REG,
        "preregistration_sha256": prereg["preregistration_sha256"],
        "amendment": ("r50_budget_extension_amendment.json"
                      if (METRICS
                          / "r50_budget_extension_amendment.json").exists()
                      else None),
        "rule": prereg["verdict_rule"],
        "readouts": prereg["readouts"],
        "pooling": prereg["pooling"],
        "primary_by_budget": report["primary"],
        "secondary_by_budget": report["secondary"],
        "by_block": report["blocks"],
        "outcome": outcome,
        "outcome_text": prereg["outcomes"][outcome],
        "outcome_evidence": evidence,
    }
    out = METRICS / "r50_verdict.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    for readout in ("primary", "secondary"):
        if not report[readout]:
            continue
        label = ("radial endpoint vs equal-budget constant"
                 if readout == "primary"
                 else "interior member vs work-matched constant")
        print(f"\n{readout.upper()}  ({label})")
        for beta, v in report[readout].items():
            print(f"  beta {beta}   pooled over {', '.join(v['blocks_pooled'])}")
            for c in v["levels"].values():
                rho = "     -" if c["median_rho"] is None else f"{c['median_rho']:6.3f}"
                box = "in " if c["inside_factor_box"] else "out"
                print(f"    ha {c['level_km']:>6.0f} km [{box}]  "
                      f"{c['verdict']:<9} {c['varying_wins']:>2}-"
                      f"{c['fixed_wins']:<2} of {c['resolved']:>2} resolved / "
                      f"{c['orbits']:>2} orbits   rho {rho}")
            p = v["paired"]
            print(f"    paired: {p['identities_ordered_after_first_win']}/"
                  f"{p['identities_counted']} identities keep winning at every "
                  f"wider level once they start")
    print(f"\nOUTCOME: {outcome}")
    print(f"  {payload['outcome_text']}")
    print(f"  evidence: {json.dumps(evidence['kendall_tau_score_vs_span'])}")
    print(f"[written] {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
