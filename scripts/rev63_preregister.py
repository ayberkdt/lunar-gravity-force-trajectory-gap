"""Freeze the R63 (O55) cells, statistics and outcome classes before any
propagation.

One lesson from (O54) is written into this registration rather than repeated.
(O54)'s outcome classes were keyed to a single budget, beta = 0.50, and that
budget turned out to sit below the crossing at both apolune levels on one
block under the very accounting the campaign was testing. The class it
returned was correct under its own rule and understated what the campaign
measured. Here the classes are keyed to the shape of the margin across the
apolune levels, evaluated at every propagated budget, and the class is
declared in advance to be a summary of the beta = 1.00 row that must be read
beside the full per-budget table, never instead of it.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r63_preregistration.json"

BLOCKS = [("span_ladder_a_uncapped", "RS1U", "r51"),
          ("span_ladder_b_uncapped", "RS2U", "r52")]
BETAS = (0.50, 0.75, 1.00)
APOLUNE_KM = [300.0, 600.0, 1200.0, 2400.0]
# Descending in beta, so the cells array *is* the attempt order the prose
# below declares. The supervisor walks this array; keeping the two in one
# place is what stops them disagreeing.
CELLS = [{"population": p, "design_key": k, "registry": r, "beta": b,
          "apolune_km": APOLUNE_KM, "identities": 64}
         for b in sorted(BETAS, reverse=True) for p, k, r in BLOCKS]


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r63-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for name in ("rev63_ladder_uncapped_rematch.py", "rev63_campaign.py",
                 "rev44_equal_work_tighter.py", "rev18_span_sweep.py",
                 "population_registry.py"):
        inputs[f"python_codes/{name}"] = sha(HERE / name)
    for c in CELLS:
        tag = f"beta_{c['beta']:.2f}"
        for pre in ("r18_span_sweep", "r19_equal_total_work"):
            f = METRICS / f"{pre}_{c['design_key']}_{tag}.json"
            inputs[f"metrics/{f.name}"] = sha(f)
    for pop, _, pre in BLOCKS:
        f = METRICS / f"{pre}_{pop}_rows.json"
        inputs[f"metrics/{f.name}"] = sha(f)

    payload = {
        "schema": "r63_preregistration_v1",
        "campaign": "R63 (O55): the ceiling-free apolune ladder, interior "
                    "panel, re-matched at the scoring tolerance",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O54) carried the level-consistent match to the controlled "
            "apolune ladder but stopped at 300 and 600 km, because on the "
            "capped blocks the reference-degree ceiling binds above that and "
            "would confound the accounting change with the ceiling. The "
            "ceiling-free blocks of (O50) and (O51) remove that confound, so "
            "the two widest levels can be run on them.",
        "two_questions":
            "at 300 and 600 km this replicates (O54) on independent "
            "ceiling-free identities; at 1200 and 2400 km it asks what (O54) "
            "could not, whether the interior member's advantage keeps "
            "widening with apolune once both the ceiling and the accounting "
            "objection are removed.",
        "cells": CELLS,
        "attempt_order": "beta 1.00 on both blocks first, then 0.75, then "
                         "0.50: the declared budget carries the summary class "
                         "and the widest levels are decidable there.",
        "frozen_rules": {
            "construction": "(O42)/(O53)/(O54) unchanged: the member is "
                            "k = 0.50 of the archived ladder sweep, reused, "
                            "and only the matched constant-degree comparator "
                            "is propagated, at both tolerance levels",
            "primary_statistic": "per block, budget and apolune level, the "
                                 "resolved interior--fixed tally and the "
                                 "median error ratio, never pooled across "
                                 "levels or across blocks",
            "robustness": "every level re-tallied at resolution cuts 0.5, 1 "
                          "and 2, nothing repropagated; a level whose leading "
                          "side moves across the cuts is reported as "
                          "cut-sensitive and never as a verdict",
            "censoring": "a comparator at or above the adopted reference "
                         "degree is censored and reported, never clamped. The "
                         "plan subcommand measured 0 censored in all six "
                         "cells, which is what ceiling-free means here.",
            "no_outcome_dependent_stopping":
                "all six cells run to completion or to the window deadline; "
                "the attempt order is fixed here and does not depend on any "
                "result this campaign produces.",
        },
        "class_is_a_summary_not_a_verdict":
            "the classes below summarise the beta = 1.00 row. They are "
            "declared to be read beside the full per-budget, per-level table "
            "and never instead of it. (O54) showed why: a class keyed to one "
            "budget can be correct under its own rule and still understate "
            "the campaign, because the accounting under test is itself what "
            "moves the crossing along the budget axis.",
        "declared_outcomes": {
            "A_widening": "at beta = 1.00, in both blocks, the resolved "
                          "interior-minus-fixed margin is non-decreasing "
                          "across 300, 600, 1200 and 2400 km and the interior "
                          "member leads at 600 km and above. Reading: the "
                          "radial-span dependence is monotone in apolune "
                          "under level-consistent matching and free of the "
                          "ceiling, which is the strongest form of the "
                          "controlled claim.",
            "B_saturating": "at beta = 1.00 the interior member leads at "
                            "600 km and above in both blocks, but the margin "
                            "is non-monotone above 600 km in at least one "
                            "block. Reading: the dependence is real but "
                            "saturates, and the main text says the advantage "
                            "widens up to a level rather than with span "
                            "generally.",
            "C_reversal": "at beta = 1.00 the leading side returns to the "
                          "constant degree at 1200 or 2400 km in at least one "
                          "block. Reading: the advantage is not monotone in "
                          "apolune, the widest arcs are their own regime, and "
                          "the geometry sentence is bounded to the levels "
                          "where it holds.",
            "note": "all three are publishable and all three change the "
                    "ladder discussion and the geometry paragraph of Section "
                    "IX.B.",
        },
        "plan_validation":
            "the plan subcommand was run on all six cells before this "
            "registration was written: 64 identities each, 0 missing "
            "telemetry, 0 censored, 0 degenerate.",
        "driver": "python_codes/rev63_ladder_uncapped_rematch.py",
        "supervisor": "python_codes/rev63_campaign.py",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r63-prereg] wrote {OUT.name}: {len(CELLS)} cells, "
          f"{len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
