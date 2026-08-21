"""R55 registration: complete the budget grid at beta = 1.25 and 1.50.

Seven rows of the regime map stop at beta = 1.00 while the two confirmatory
designs and the four ladder levels carry 1.25 and 1.50. The grid is therefore
ragged on its compute-rich side, and every statement the paper makes about what
happens as compute becomes plentiful rests on the two designs alone.

This is the same shape as R53 and is registered the same way. It adds a column
to a grid rather than a degree of freedom: no population, orbit, parameter or
reference degree is added, the budgets are two values of the grid registered in
r14_preregistration.json, and every population's calibration record already
carries both of them among its computed budgets, which is checked here and
recorded with the file digest that proves it. What runs is the ladder those
populations have already run at 0.50, 0.62, 0.75 and 1.00, at two further
budgets.

The cell order is budget-major on purpose. A run that stops early then leaves a
complete beta = 1.25 column across all seven populations rather than a ragged
edge of whichever populations happened to be reached, and a complete column is
readable where a ragged one is not.

Outcome classes are fixed here, before anything propagates:

  W  every new cell continues the budget ordering its own 0.75 and 1.00 scores
     establish, on both panels.
  X  at least one population's tally reverses back across 1.00, which would
     make the compute-rich end non-monotone and is reported as such.
  Y  a cell resolves fewer than `min_resolved` comparisons on both panels, so
     it carries no score rather than a rounded one.
  Z  a cell is declared and not run. Its trigger is the clock **or** a resource
     floor. R53 registered Z by the clock alone and then stopped on a disk
     floor inside its window, so its classes did not partition the ways it
     could stop; this registration fixes that defect rather than inheriting it.

Prohibited: adding, dropping or reordering a cell once this file is written;
recomputing any earlier verdict; quoting a crossing more precisely than the
sampled grid supports; and reading any cell of this campaign as enlarging the
confirmatory population.

Usage:  python rev55_preregister.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

REG = "r55"
BUDGETS = (1.25, 1.50)
MIN_RESOLVED = 6

# population -> (design key, parent registry, ladder driver, calibration record)
POPULATIONS = [
    ("designC", "C", "r29", "rev29_designC_ladder.py",
     "r29_budget_pareto_designC.json"),
    ("polar", "SP", "r30", "rev30_stratum_ops.py",
     "r30_budget_pareto_polar.json"),
    ("equatorial", "SE", "r30", "rev30_stratum_ops.py",
     "r30_budget_pareto_equatorial.json"),
    ("frozen_like", "SF", "r30", "rev30_stratum_ops.py",
     "r30_budget_pareto_frozen_like.json"),
    ("high_apolune", "SH", "r30", "rev30_stratum_ops.py",
     "r30_budget_pareto_high_apolune.json"),
    ("operational_elliptical", "OE", "r31", "rev30_stratum_ops.py",
     "r31_budget_pareto_operational_elliptical.json"),
    ("operational_elliptical_uncapped", "OEU", "r38", "rev30_stratum_ops.py",
     "r38_budget_pareto_operational_elliptical_uncapped.json"),
]

# Measured, not guessed: the frozen-like cell ran its whole ladder in 119.0 min
# at four workers on 15 August, and design C's took 123.2 min at eleven. The
# prior is the larger of those with a margin, and the campaign refuses to start
# a cell it cannot finish inside it.
PRIOR_MIN = 150.0


def sha(name: str) -> str:
    return hashlib.sha256((METRICS / name).read_bytes()).hexdigest()


def main() -> int:
    out = METRICS / f"{REG}_preregistration.json"
    if out.exists():
        print(f"[abort] {out.name} exists; a registration is written once, "
              f"before the first propagation, and never rewritten")
        return 2

    calibrations, missing = {}, []
    for pop, key, registry, driver, cal in POPULATIONS:
        p = METRICS / cal
        if not p.exists():
            missing.append(f"{key}: {cal} is not on disk")
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        have = d.get("budgets_computed") or []
        for b in BUDGETS:
            if b not in have:
                missing.append(f"{key}: calibration does not carry beta={b}")
        calibrations[key] = {"record": cal, "sha256": sha(cal),
                             "budgets_computed": have}
    if missing:
        print("[abort] this campaign only adds budgets that are already "
              "calibrated, and these are not:")
        for m in missing:
            print("  " + m)
        return 2

    cells = []
    order = 0
    for beta in BUDGETS:                      # budget-major: see the docstring
        for pop, key, registry, driver, cal in POPULATIONS:
            order += 1
            cells.append({"order": order, "budget": beta, "population": pop,
                          "design_key": key, "registry": registry,
                          "driver": driver, "prior_min": PRIOR_MIN,
                          "measured_prior":
                          "r53 frozen-like ladder 119.0 min at 4 workers"})

    payload = {
        "schema": f"{REG}_preregistration_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "campaign": REG,
        "budgets": list(BUDGETS),
        "what_this_is": (
            "the compute-rich half of the registered budget grid, run on the "
            "seven populations whose regime-map rows stop at 1.00. A column "
            "added to a grid, not a new degree of freedom."),
        "what_is_not_added": (
            "no population, orbit, initial state, parameter, tolerance or "
            "reference degree. Both budgets are values of the grid registered "
            "in r14_preregistration.json and both are already present in every "
            "population's calibration record, whose digests are indexed here."),
        "cell_order_rule": (
            "budget-major, so that a run stopping early leaves a complete "
            "beta=1.25 column across all seven populations rather than a "
            "ragged edge."),
        "cells": cells,
        "calibrations_read": calibrations,
        "min_resolved": MIN_RESOLVED,
        "scoring_rule": (
            "per population and per panel, the score (wins - losses) / "
            "resolved, read against the ordering that population's own 0.75 "
            "and 1.00 scores establish. Both panels are scored and neither is "
            "nominated as primary after the fact."),
        "outcomes": {
            "W": "every new cell continues its population's own 0.75-to-1.00 "
                 "budget ordering, on both panels",
            "X": "at least one population's tally reverses back across 1.00, "
                 "making the compute-rich end non-monotone",
            "Y": "a cell resolves fewer than min_resolved comparisons on both "
                 "panels and carries no score",
            "Z": "a cell is declared and not run, its trigger being the clock "
                 "or a resource floor. R53 registered Z by the clock alone and "
                 "stopped on a disk floor inside its window; that defect is "
                 "corrected here rather than inherited",
        },
        "prohibited": [
            "adding, dropping or reordering a cell after this file is written",
            "recomputing any earlier campaign's verdict",
            "quoting a crossing more precisely than the sampled grid supports",
            "reading any cell here as enlarging the confirmatory population",
            "pooling any population with another",
        ],
        "expected_weakness": (
            "the resolution rule decides fewer comparisons as the budget "
            "grows, because both policies descend toward the numerical floor. "
            "Cells at these budgets are therefore expected to be thinner than "
            "those below them, and outcome Y is a real possibility rather than "
            "a formality. A thin cell is reported as thin."),
    }
    payload["preregistration_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] {out.name}: {len(cells)} cells, "
          f"budgets {', '.join(str(b) for b in BUDGETS)}")
    print(f"  order: {', '.join(c['design_key'] + '@' + str(c['budget']) for c in cells)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
