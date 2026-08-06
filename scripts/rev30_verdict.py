"""R30: apply the pre-registered stratum rule to the geometry panel.

Same shape as rev29_verdict.py and the same reason for existing: the rule is
read out of the registration and applied by a script, so that the sentence that
ends up in the manuscript is the one the tallies support rather than the one
that reads well.

Per stratum the R19 realized-work tally decides: interior, constant, split, or
undecided when nothing resolves. The panel returns K only when every propagated
stratum returns interior; L when they disagree; M when a stratum cannot decide.

Strata are never pooled, here or anywhere else. Their sub-boxes overlap each
other and the coverage box, so a pooled count would weight the geometry that is
oversampled twice. The coverage designs are printed alongside for reference and
take no part in the panel verdict.

Usage:
    python rev30_verdict.py
    python rev30_verdict.py --beta 0.75
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

import population_registry as registry


def tally(path: Path) -> dict:
    if not path.exists():
        return {"verdict": "not propagated", "record": None}
    d = json.loads(path.read_text(encoding="utf-8"))
    s = d["summary"]
    i, f = int(s["resolved_interior_wins"]), int(s["resolved_fixed_wins"])
    if i + f == 0:
        v = "undecided"
    elif i > f:
        v = "interior"
    elif f > i:
        v = "constant"
    else:
        v = "split"
    return {"verdict": v, "record": path.name, "orbits": s["orbits"],
            "resolved": s["resolved"], "interior_wins": i, "fixed_wins": f,
            "unresolved": int(s["unresolved"]), "median_rho": s["median_rho"],
            "achieved_work_ratio_median": s["achieved_work_ratio"]["median"]}


def coverage_path(design: str, beta: float) -> Path:
    if design == "C":
        return METRICS / f"r19_equal_total_work_C_beta_{beta:.2f}.json"
    suffix = "" if abs(beta - 1.0) < 1e-12 else f"_beta_{beta:.2f}"
    return METRICS / f"r19_equal_total_work_{design}{suffix}.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=1.00)
    ap.add_argument("--registry", default="r30")
    a = ap.parse_args()

    prereg = registry.registration(a.registry)
    output = METRICS / f"{a.registry}_verdict.json"
    strata = {}
    for name, spec in registry.populations(a.registry).items():
        p = (METRICS / f"r19_equal_total_work_{spec['design_key']}"
                       f"_beta_{a.beta:.2f}.json")
        strata[name] = {**tally(p), "design_key": spec["design_key"],
                        "sub_box": (prereg.get("sub_boxes", {}).get(name)
                                    or prereg.get("design", {}).get(
                                        "factor_box"))}

    live = {k: v for k, v in strata.items() if v["verdict"] != "not propagated"}
    verdicts = {v["verdict"] for v in live.values()}
    keys = sorted(prereg["outcomes"])          # [interior-holds, disagrees, undecided]
    holds, disagrees, undecided = keys[0], keys[1], keys[2]
    # A panel verdict is a statement about every declared population. With some
    # of them unpropagated it is not one, and reporting the propagated subset
    # under a pre-registered outcome would read as a result the campaign never
    # reached. Same distinction rev29_verdict.py draws for an unfinished ladder.
    missing = [k for k, v in strata.items() if v["verdict"] == "not propagated"]
    if missing:
        outcome = "incomplete_not_an_outcome"
    elif "undecided" in verdicts:
        outcome = undecided
    elif verdicts == {"interior"}:
        outcome = holds
    else:
        outcome = disagrees

    coverage = {d: tally(coverage_path(d, a.beta)) for d in ("A", "B", "C")}
    payload = {
        "schema": f"{a.registry}_verdict_v1",
        "registry": a.registry,
        "beta": a.beta,
        "preregistration_sha256": prereg["preregistration_sha256"],
        "rule": prereg["verdict_rule"],
        "pooling": prereg.get("pooling") or prereg.get("why_separate_from_R30"),
        "strata": strata,
        "coverage_designs_for_reference": coverage,
        "strata_propagated": sorted(live),
        "outcome": outcome,
        "outcome_text": (prereg["outcomes"].get(outcome) or
                         ("declared populations not propagated: "
                          + ", ".join(missing)
                          + ". A panel verdict is a statement about all of "
                            "them, so this is a state of the campaign and not "
                            "a result of it; the propagated populations are "
                            "reported individually.")),
        "declared_but_not_propagated": missing,
    }
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"budget beta = {a.beta:.2f}   registration {a.registry}\n")
    print("populations:")
    for name, v in strata.items():
        if v["verdict"] == "not propagated":
            print(f"  {name:<14} not propagated")
            continue
        print(f"  {name:<14} {v['verdict']:<9} "
              f"{v['interior_wins']:>2}-{v['fixed_wins']:<2} of "
              f"{v['resolved']:>2} resolved / {v['orbits']} orbits, "
              f"rho {v['median_rho']:.3f}")
    print("\ncoverage designs (reference, not part of the panel verdict):")
    for d, v in coverage.items():
        if v["verdict"] == "not propagated":
            print(f"  design {d}      not propagated")
            continue
        print(f"  design {d}      {v['verdict']:<9} "
              f"{v['interior_wins']:>2}-{v['fixed_wins']:<2} of "
              f"{v['resolved']:>2} resolved / {v['orbits']} orbits, "
              f"rho {v['median_rho']:.3f}")
    print(f"\nOUTCOME: {outcome}")
    print(f"  {payload['outcome_text']}")
    print(f"[written] {output.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
