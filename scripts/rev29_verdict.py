"""R29: apply the pre-registered verdict rule to the design-C ladder.

The rule is not restated here in prose and then implemented separately; it is
read out of r29_preregistration.json and applied. What the script decides:

  per budget   interior if resolved_interior_wins > resolved_fixed_wins,
               constant if the reverse, split if equal, unresolved if the
               record is missing or has no resolved comparisons;
  per design   the crossing bracket is (largest grid budget with a constant
               verdict, smallest grid budget with an interior verdict];
  campaign     H if design C brackets the crossing where A and B do,
               I if it brackets it somewhere else, J if it cannot place it.

Three things this deliberately does not do. It does not average designs. It does
not let the median error ratio decide anything -- that number is recorded as a
secondary signature because it moves continuously where the tally jumps, and the
campaign has decided by tally since R19. And it does not read the post-hoc
budget into the bracket: beta = 0.62 is reported beside the grid, never inside
it, on every design.

Usage:
    python rev29_verdict.py
    python rev29_verdict.py --json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PREREG = METRICS / "r29_preregistration.json"
R26_PREREG = METRICS / "r26_preregistration.json"
OUTPUT = METRICS / "r29_verdict.json"

POST_HOC = 0.62


def record_path(design: str, beta: float) -> Path:
    """Design C always carries a budget suffix; A and B leave beta = 1 bare."""
    if design == "C":
        return METRICS / f"r19_equal_total_work_C_beta_{beta:.2f}.json"
    suffix = "" if abs(beta - 1.0) < 1e-12 else f"_beta_{beta:.2f}"
    return METRICS / f"r19_equal_total_work_{design}{suffix}.json"


def verdict_at(design: str, beta: float) -> dict:
    p = record_path(design, beta)
    post_hoc = abs(beta - POST_HOC) < 1e-12
    if not p.exists():
        return {"beta": beta, "verdict": "not propagated", "record": None,
                "post_hoc": post_hoc}
    d = json.loads(p.read_text(encoding="utf-8"))
    s = d["summary"]
    i, f = int(s["resolved_interior_wins"]), int(s["resolved_fixed_wins"])
    if i + f == 0:
        v = "unresolved"
    elif i > f:
        v = "interior"
    elif f > i:
        v = "constant"
    else:
        v = "split"
    return {"beta": beta, "verdict": v, "record": p.name,
            "orbits": s["orbits"], "resolved": s["resolved"],
            "interior_wins": i, "fixed_wins": f,
            "unresolved": int(s["unresolved"]),
            "median_rho": s["median_rho"],
            "achieved_work_ratio_median": s["achieved_work_ratio"]["median"],
            "post_hoc": post_hoc}


def bracket(entries: list[dict]) -> dict:
    """(largest constant budget, smallest interior budget], grid budgets only."""
    grid = [e for e in entries if not e["post_hoc"]
            and e["verdict"] in ("interior", "constant", "split")]
    const = [e["beta"] for e in grid if e["verdict"] == "constant"]
    inter = [e["beta"] for e in grid if e["verdict"] == "interior"]
    if not const or not inter:
        return {"placed": False,
                "why": ("no constant verdict below an interior one among the "
                        "propagated grid budgets")}
    lo, hi = max(const), min(inter)
    if lo >= hi:
        return {"placed": False,
                "why": (f"verdicts are not monotone: constant at {lo:.2f} is "
                        f"not below interior at {hi:.2f}")}
    return {"placed": True, "lower": lo, "upper": hi,
            "text": f"({lo:.2f}, {hi:.2f}]"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    r26 = json.loads(R26_PREREG.read_text(encoding="utf-8"))
    grid = list(prereg["stages"]["calibration"]["grid"])
    betas = sorted(set(grid + [POST_HOC]))

    designs = {}
    for d in ("A", "B", "C"):
        entries = [verdict_at(d, b) for b in betas]
        designs[d] = {"entries": entries, "bracket": bracket(entries)}

    ref = prereg["verdict_rule"]["reference_bracket"]
    ref_text = f"({ref['A'][0]:.2f}, {ref['A'][1]:.2f}]"
    c = designs["C"]["bracket"]
    propagated = [e for e in designs["C"]["entries"]
                  if e["verdict"] != "not propagated"]
    # "cannot place the crossing" and "has not yet propagated the budgets that
    # would place it" are different states, and only the first is an outcome.
    # A campaign still running would otherwise read as a pre-registered result.
    grid_done = [e for e in designs["C"]["entries"]
                 if not e["post_hoc"] and e["verdict"] not in
                 ("not propagated", "unresolved")]
    incomplete = len(grid_done) < 2
    if incomplete:
        outcome = "incomplete_not_an_outcome"
    elif not c["placed"]:
        outcome = "J_designC_unresolved"
    elif c["text"] == ref_text:
        outcome = "H_designC_agrees"
    else:
        outcome = "I_designC_disagrees"

    payload = {
        "schema": "r29_verdict_v1",
        "preregistration_sha256": prereg["preregistration_sha256"],
        "parent_preregistration_sha256": r26["preregistration_sha256"],
        "rule": prereg["verdict_rule"],
        "designs": designs,
        "reference_bracket": ref_text,
        "designC_bracket": c.get("text"),
        "outcome": outcome,
        "outcome_text": (r26["outcomes"].get(outcome) or
                         (f"the ladder has propagated {len(grid_done)} grid "
                          f"budget(s); at least two are needed before the "
                          f"pre-registered outcomes apply. This is a state of "
                          f"the campaign, not a result of it.")),
        "grid_budgets_with_a_verdict": [e["beta"] for e in grid_done],
        "designC_budgets_propagated": [e["beta"] for e in propagated],
        "completion": {
            e["beta"]: (f"{e.get('orbits')} orbits" if e["record"] else "none")
            for e in designs["C"]["entries"]},
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if a.json:
        print(json.dumps(payload, indent=2))
        return 0

    print(f"pre-registered reference bracket (designs A, B): {ref_text}\n")
    for d in ("A", "B", "C"):
        print(f"design {d}:")
        for e in designs[d]["entries"]:
            if e["verdict"] == "not propagated":
                continue
            flag = "  [post hoc, outside the bracket rule]" if e["post_hoc"] else ""
            print(f"  beta {e['beta']:.2f}  {e['verdict']:<9} "
                  f"{e['interior_wins']:>2}-{e['fixed_wins']:<2} of "
                  f"{e['resolved']:>2} resolved / {e['orbits']} orbits, "
                  f"rho {e['median_rho']:.3f}{flag}")
        b = designs[d]["bracket"]
        print(f"  bracket: {b['text'] if b['placed'] else 'not placed - ' + b['why']}\n")
    print(f"OUTCOME: {outcome}")
    print(f"  {payload['outcome_text']}")
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
