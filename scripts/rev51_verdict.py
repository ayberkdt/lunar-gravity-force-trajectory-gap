"""R51: read the cap-lifted ladder against its capped parent, level by level.

Written while the control's base was still propagating and before any of its
comparisons existed, for the same reason rev50_verdict.py was: the sentence
that reaches the manuscript should be the one the tallies support, and what
counts as "the ordering survived" has to be arithmetic fixed in advance rather
than a judgement made once the numbers are in view.

The comparison is paired twice over. Each orbit of this population is the same
orbit as in the parent, with the same identity and the same apolune level, so a
level cell here faces exactly the level cell there; and within a level the
members are the same identities. Nothing is pooled across levels.

The registered outcome classes, from r51_preregistration.json:

  X_ordering_survives                 same policy carries the same levels, and
                                      the score still orders in radial span at
                                      every budget that decides
  Y_ordering_survives_magnitude_falls ordering unchanged, but the median ratio
                                      at the wide levels drops by an order of
                                      magnitude or more
  Z_ordering_breaks                   a level changes hands, or the score stops
                                      ordering in span
  W_undecided                         too few comparisons resolve per level to
                                      order the levels

Y is a refinement of X rather than a rival: the ordering test decides first,
and the magnitude test then chooses between X and Y. That is stated here
because a reader can otherwise not tell which of the two a run returned.

Usage:
    python rev51_verdict.py
    python rev51_verdict.py --beta 0.75
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

CAPPED_KEY, UNCAPPED_KEY = "RS1", "RS1U"
CAPPED_DESIGN = "r50_span_ladder_a_design_frozen.json"
BETAS = ["1.00", "0.75", "0.62", "0.50"]
MAGNITUDE_DROP = 10.0          # "an order of magnitude or more"


def levels() -> dict:
    d = json.loads((METRICS / CAPPED_DESIGN).read_text(encoding="utf-8"))
    return {o["sobol_index"]: o["apolune_level_km"] for o in d["orbits"]}


def cells(key: str, beta: str, index: dict) -> dict | None:
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", default=None)
    a = ap.parse_args()

    prereg_p = METRICS / "r51_preregistration.json"
    if not prereg_p.exists():
        print("r51_preregistration.json missing; the control is not registered")
        return 0
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))
    index = levels()
    betas = [a.beta] if a.beta else BETAS

    report, changed_hands, taus, drops = {}, [], [], []
    for beta in betas:
        cap = cells(CAPPED_KEY, beta, index)
        unc = cells(UNCAPPED_KEY, beta, index)
        if cap is None or unc is None:
            continue
        rows = {}
        for lvl in sorted(unc):
            c, u = cap.get(lvl), unc[lvl]
            rows[lvl] = {"capped": c, "uncapped": u,
                         "verdict_changed": bool(c and c["verdict"]
                                                 != u["verdict"])}
            if rows[lvl]["verdict_changed"]:
                changed_hands.append({"beta": beta, "level_km": lvl,
                                      "from": c["verdict"],
                                      "to": u["verdict"]})
            if c and c["median_rho"] and u["median_rho"]:
                rows[lvl]["rho_factor"] = c["median_rho"] / u["median_rho"]
                if lvl >= 1200.0:
                    drops.append(rows[lvl]["rho_factor"])
        deciding = [(lvl, r["uncapped"]["score"]) for lvl, r in rows.items()
                    if r["uncapped"]["resolved"]]
        tau = kendall_tau(deciding) if len(deciding) >= 3 else None
        if tau is not None:
            taus.append(tau)
        report[beta] = {"levels": rows, "kendall_tau_uncapped": tau,
                        "levels_deciding": len(deciding)}

    if not report:
        print("no uncapped ladder on disk yet; nothing to read")
        return 0

    if not taus:
        outcome = "W_undecided"
    elif changed_hands or not all(abs(t) == 1.0 for t in taus):
        outcome = "Z_ordering_breaks"
    elif drops and median(drops) >= MAGNITUDE_DROP:
        outcome = "Y_ordering_survives_magnitude_falls"
    else:
        outcome = "X_ordering_survives"

    payload = {
        "schema": "r51_verdict_v1",
        "preregistration_sha256": prereg["preregistration_sha256"],
        "rule": prereg["verdict_rule"],
        "budgets_read": sorted(report),
        "by_budget": report,
        "levels_changing_hands": changed_hands,
        "kendall_tau_uncapped": taus,
        "median_rho_factor_wide_levels": median(drops) if drops else None,
        "magnitude_drop_threshold": MAGNITUDE_DROP,
        "outcome": outcome,
        "outcome_text": prereg["outcomes"][outcome],
    }
    (METRICS / "r51_verdict.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    for beta, blk in report.items():
        print(f"\nbeta {beta}   (capped -> uncapped)")
        for lvl, r in blk["levels"].items():
            c, u = r["capped"], r["uncapped"]
            cap_txt = ("--" if not c else
                       f"{c['verdict']:<8} {c['radial']:2d}-{c['fixed']:<2d} "
                       f"rho {c['median_rho']:.4g}")
            print(f"  ha {lvl:6.0f} km  {cap_txt}")
            print(f"                 {u['verdict']:<8} "
                  f"{u['radial']:2d}-{u['fixed']:<2d} of {u['resolved']:2d}  "
                  f"rho {u['median_rho']:.4g}"
                  + (f"   ratio falls x{r['rho_factor']:.3g}"
                     if "rho_factor" in r else ""))
        print(f"  tau(uncapped) = {blk['kendall_tau_uncapped']}")
    print(f"\nlevels changing hands: {changed_hands or 'none'}")
    print(f"OUTCOME: {outcome}")
    print(f"  {payload['outcome_text']}")
    print("[written] r51_verdict.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
