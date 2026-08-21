"""R53 verdict: does the post-hoc budget sit between its neighbours?

The registration asks one question per cell. Each of the seven populations
already had a score at 0.50 and a score at 0.75, and the budget ordering
printed for it was read off that pair. A point at 0.62 either falls inside that
interval, in which case the ordering is monotone through it and the crossing
already reported is located more finely, or it does not, in which case the
population's bracket sentence has to say so.

Two comparisons are scored, not one, because the regime map draws two panels
and a population appears in both: the radial endpoint against its equal-budget
constant degree, which is where the geometry-axis result lives, and the interior
member against the constant degree matched on realized total quadratic work,
which is where the crossing bracket lives. Neither is nominated as the primary
after the fact. A cell is bracketed only if it is bracketed on both, so the
outcome cannot be bought by choosing the panel that cooperates.

The score is the one used everywhere else, (wins - losses) / resolved, and it is
undefined when too few comparisons resolve, which is the registration's outcome
Y rather than a number to be rounded into one.

Usage:  python rev53_verdict.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

REG = "r53"
BETA = "0.62"
NEIGHBOURS = ("0.50", "0.75")

# the resolved count below which the regime map leaves a cell unfilled; the same
# floor decides outcome Y here so that a cell the figure will not draw is not
# scored in prose either.
MIN_RESOLVED = 6


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def endpoint_panel(key: str, beta: str) -> dict | None:
    """The radial endpoint against its equal-budget constant degree."""
    p = METRICS / f"r14_trajectory_{key}_beta_{beta}.json"
    if not p.exists():
        return None
    s = json.loads(p.read_text(encoding="utf-8"))["summary"]
    return {"wins": s["resolved_atallah_wins"],
            "losses": s["resolved_fixed_wins"],
            "resolved": s["resolved_atallah_wins"] + s["resolved_fixed_wins"],
            "median_rho": s["rho_budget"]["median"]}


def interior_panel(key: str, beta: str) -> dict | None:
    """The interior member against the work-matched constant degree."""
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta}.json"
    if not p.exists():
        return None
    s = json.loads(p.read_text(encoding="utf-8"))["summary"]
    return {"wins": s["resolved_interior_wins"],
            "losses": s["resolved_fixed_wins"],
            "resolved": s["resolved_interior_wins"] + s["resolved_fixed_wins"],
            "median_rho": s["median_rho"]}


def score(cell: dict | None) -> float | None:
    if cell is None or cell["resolved"] < MIN_RESOLVED:
        return None
    return (cell["wins"] - cell["losses"]) / cell["resolved"]


def bracketed(lo: float | None, mid: float | None,
              hi: float | None) -> bool | None:
    """Is the middle score inside the interval its neighbours span?

    Direction is taken from the neighbours rather than assumed: a population
    whose score falls with budget brackets in the other order, and nothing in
    the registration says which way a given population runs.
    """
    if lo is None or mid is None or hi is None:
        return None
    return min(lo, hi) <= mid <= max(lo, hi)


def informative(lo: float | None, hi: float | None) -> bool | None:
    """Does the neighbouring pair leave an interval to fall inside?

    On the wide-elliptic population the interior member beats its work-matched
    comparator on every resolved orbit at 0.50, 0.62 and 0.75 alike, so that
    panel's score is pinned at +1 and the bracket test passes without being
    tested. It is the same saturation the ordering statistic runs into on the
    span ladder, and it is recorded rather than counted: a degenerate interval
    is not evidence of monotonicity, and a cell whose panels are both
    degenerate has not been checked by this campaign at all.
    """
    if lo is None or hi is None:
        return None
    return lo != hi


def main() -> int:
    prereg_p = METRICS / f"{REG}_preregistration.json"
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))

    report, unrun, unscored, outside, degenerate = {}, [], [], [], []
    for spec in sorted(prereg["cells"], key=lambda c: c["order"]):
        key = spec["design_key"]
        panels = {}
        for name, fn in (("endpoint", endpoint_panel),
                         ("interior", interior_panel)):
            cells = {b: fn(key, b) for b in (*NEIGHBOURS, BETA)}
            scores = {b: score(cells[b]) for b in cells}
            panels[name] = {
                "cells": cells,
                "scores": scores,
                "bracketed": bracketed(scores[NEIGHBOURS[0]], scores[BETA],
                                       scores[NEIGHBOURS[1]]),
                "informative": informative(scores[NEIGHBOURS[0]],
                                           scores[NEIGHBOURS[1]]),
            }
        present = [panels[n]["cells"][BETA] is not None for n in panels]
        ran = all(present)
        # a cell that wrote one stage and not the other is not a cell that was
        # never started, and the two are not reported as the same thing
        state = "complete" if ran else ("partial" if any(present)
                                        else "not run")
        verdicts = [panels[n]["bracketed"] for n in panels]
        live = [n for n in panels if panels[n]["informative"] is True]
        report[key] = {
            "population": spec["population"],
            "order": spec["order"],
            "ran": ran,
            "state": state,
            "panels": panels,
            "informative_panels": live,
            "bracketed_on_both": (all(v is True for v in verdicts)
                                  if ran else None),
        }
        if ran and not live:
            degenerate.append(key)
        if not ran:
            unrun.append(key)
        elif any(v is None for v in verdicts):
            unscored.append(key)
        elif not all(verdicts):
            outside.append(key)

    # W over an empty set is true and says nothing. A campaign that reached no
    # cell has an outcome, and it is not the one about monotonicity.
    if not any(v["ran"] for v in report.values()):
        outcome, text = "Z_none_run", (
            "no cell of the column was completed, so nothing is scored and no "
            "statement about the budget ordering follows from this campaign.")
    elif outside:
        outcome, text = "X_not_monotone", (
            "at least one completed cell falls outside the interval its own "
            "0.50 and 0.75 scores span. The non-monotonicity is reported per "
            "population and the affected bracket sentence is qualified.")
    elif unscored:
        outcome, text = "Y_underresolved", (
            "every completed cell that scores is bracketed, and at least one "
            "resolves too few comparisons to score. It is drawn as unresolved "
            "and no verdict is read from it.")
    else:
        outcome, text = "W_monotone", (
            "every completed cell scores inside the interval its own 0.50 and "
            "0.75 scores span, on both panels. The budget ordering already "
            "reported for each population is monotone through the new point, "
            "and the crossing intervals already printed stand unchanged and "
            "better located.")

    payload = {
        "schema": f"{REG}_verdict_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "beta": float(BETA),
        "neighbours": list(NEIGHBOURS),
        "preregistration_sha256": sha(prereg_p),
        "rule": ("per population and per panel, the score (wins - losses) / "
                 "resolved at 0.62 against the interval its 0.50 and 0.75 "
                 "scores span; a cell counts as bracketed only if it is "
                 "bracketed on both panels"),
        "min_resolved": MIN_RESOLVED,
        "cells_run": [k for k, v in report.items() if v["ran"]],
        "cells_declared_and_not_run": unrun,
        "cells_underresolved": unscored,
        "cells_with_no_informative_panel": degenerate,
        "cells_outside_their_interval": outside,
        "outcome": outcome,
        "outcome_text": text,
        "outcome_note": ("outcome Z of the registration, a cell the clock did "
                         "not reach, is not exclusive of the others and is "
                         "carried by cells_declared_and_not_run"),
        "report": report,
    }
    out = METRICS / f"{REG}_verdict.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    print(f"[{REG}] {outcome}")
    for key, v in sorted(report.items(), key=lambda kv: kv[1]["order"]):
        if not v["ran"]:
            print(f"  {key:<4} {v['state']}")
            continue
        bits = []
        for name in ("endpoint", "interior"):
            p = v["panels"][name]
            c = p["cells"][BETA]
            s = p["scores"]
            bits.append(f"{name} {c['wins']}-{c['losses']}"
                        f" s={s[BETA]:+.2f} in "
                        f"[{s[NEIGHBOURS[0]]:+.2f},{s[NEIGHBOURS[1]]:+.2f}]"
                        if s[BETA] is not None
                        else f"{name} underresolved")
        mark = {True: "ok", False: "OUTSIDE", None: "?"}[v["bracketed_on_both"]]
        if not v["informative_panels"]:
            mark = "degenerate"
        print(f"  {key:<4} [{mark}] " + " | ".join(bits))
    print(f"[written] {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
