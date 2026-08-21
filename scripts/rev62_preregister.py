"""Freeze the R62 (O54) cells, statistics and outcome classes before any
propagation.

The outcome classes are pinned to the sharpest contrast the tight-level
baseline already contains, so that no class can be chosen after the fact: at
beta = 0.50 the 300 km subset leads the constant degree and the 600 km subset
leads the interior member, in both blocks (2-14 -> 11-1 on block A, 1-15 ->
13-1 on block B). Whether that controlled apolune step survives the
scoring-tolerance match is the whole question, and the three classes below
partition the possible answers to it.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r62_preregistration.json"

BLOCKS = [("span_ladder_a", "RS1"), ("span_ladder_b", "RS2")]
BETAS = (0.50, 0.75, 1.00)
APOLUNE_KM = [300.0, 600.0]
CELLS = [{"population": p, "design_key": k, "beta": b,
          "apolune_km": APOLUNE_KM, "identities": 32}
         for b in BETAS for p, k in BLOCKS]

# The tight-level baseline these cells are read against, computed from the
# sealed R19 ladder records and pinned here so the comparison cannot drift.
BASELINE = {
    "RS1": {"0.50": {"300": [2, 14], "600": [11, 1]},
            "0.75": {"300": [7, 4], "600": [11, 0]},
            "1.00": {"300": [3, 3], "600": [15, 0]}},
    "RS2": {"0.50": {"300": [1, 15], "600": [13, 1]},
            "0.75": {"300": [1, 5], "600": [12, 1]},
            "1.00": {"300": [3, 3], "600": [11, 0]}},
}


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r62-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for name in ("rev62_ladder_interior_rematch.py", "rev62_campaign.py",
                 "rev44_equal_work_tighter.py", "rev18_span_sweep.py",
                 "population_registry.py"):
        inputs[f"python_codes/{name}"] = sha(HERE / name)
    for c in CELLS:
        tag = f"beta_{c['beta']:.2f}"
        for pre in ("r18_span_sweep", "r19_equal_total_work"):
            f = METRICS / f"{pre}_{c['design_key']}_{tag}.json"
            inputs[f"metrics/{f.name}"] = sha(f)
    for pop, _ in BLOCKS:
        inputs[f"metrics/r50_{pop}_rows.json"] = sha(
            METRICS / f"r50_{pop}_rows.json")

    payload = {
        "schema": "r62_preregistration_v1",
        "campaign": "R62 (O54): the controlled apolune ladder's interior "
                    "panel, re-matched at the scoring tolerance",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O53) showed that the tolerance the realized-work match is made "
            "at changes which populations share the crossing bracket, and "
            "that two of those changes survive every resolution cut. The "
            "geometry strata therefore now carry the level-consistent match "
            "while the controlled test of the radial-span direction, the "
            "paired apolune ladder of (O49), still carries the older "
            "tight-level convention on its interior panel. The uncontrolled "
            "evidence is held to a stricter accounting than the controlled "
            "evidence. This campaign removes that asymmetry.",
        "question":
            "Does a controlled apolune change, perilune and angular geometry "
            "held fixed, still move the interior crossing once the "
            "comparator is matched at the tolerance the errors are scored "
            "at?",
        "cells": CELLS,
        "scope_and_why_it_stops_here":
            "Only the 300 and 600 km levels are run. (O49) reports the "
            "reference-degree ceiling beginning to bind at 1200 km and "
            "binding on every orbit at 2400 km, so a rematch there would "
            "confound the accounting change with the ceiling. The 300 to 600 "
            "step alone answers the question, and the plan subcommand "
            "confirmed 0 censored comparators at both levels in all six "
            "cells, which is the same statement measured rather than "
            "asserted.",
        "tight_level_baseline_interior_minus_fixed": BASELINE,
        "frozen_rules": {
            "construction": "(O42)/(O53) unchanged: the member is k = 0.50 of "
                            "the archived R18 ladder sweep, reused, and only "
                            "the matched constant-degree comparator is "
                            "propagated, at both tolerance levels",
            "match_target": "W_k^tighter from the archived tighter-level "
                            "member telemetry; N* = round(N_0 sqrt(W_k/W_0)) "
                            "with both works at that level",
            "primary_statistic": "per block and budget, the resolved "
                                 "interior--fixed tally at 300 km and at "
                                 "600 km separately, never pooled across "
                                 "levels, with the median error ratio at "
                                 "each level",
            "robustness": "every cell re-tallied at resolution cuts 0.5, 1 "
                          "and 2, nothing repropagated; a level whose leading "
                          "side moves across the cuts is reported as "
                          "cut-sensitive and never as a verdict",
            "censoring": "a comparator at or above the adopted reference "
                         "degree is censored and reported, never clamped",
            "no_outcome_dependent_stopping":
                "all six cells run to completion or to the window deadline; "
                "the attempt order is fixed here and does not depend on any "
                "result this campaign produces",
        },
        "attempt_order": "beta 0.50 on both blocks first, then 0.75, then "
                         "1.00: the 0.50 row carries the outcome classes, so "
                         "a window that ends early still decides the "
                         "campaign.",
        "declared_outcomes": {
            "A_shift_preserved":
                "in both blocks at beta = 0.50 the 300 km subset leads the "
                "constant degree and the 600 km subset leads the interior "
                "member, and both leading sides hold at all three resolution "
                "cuts. Reading: the controlled leg of the radial-span "
                "interpretation survives the accounting change, and the "
                "geometry claim is stated without an accounting caveat on "
                "its controlled evidence.",
            "B_shift_attenuated":
                "the 300 to 600 step still moves the tally toward the "
                "interior member in both blocks at beta = 0.50, but in at "
                "least one block the leading side no longer differs between "
                "the levels, or the difference is cut-sensitive. Reading: "
                "the radial-span interpretation is partly a property of the "
                "cost accounting, and the main text's geometry sentences are "
                "narrowed to say so.",
            "C_shift_absent":
                "at beta = 0.50 the leading side is the same at both levels "
                "in both blocks. Reading: for the interior policy the "
                "radial-span causal reading does not survive level-consistent "
                "matching; the claim is withdrawn for the interior member and "
                "retained only for the endpoint ladder, which this campaign "
                "does not touch.",
            "note": "all three are publishable and all three change the same "
                    "places: the geometry paragraph of Section IX.B, the "
                    "ladder discussion in the supplement, and (O49)'s entry "
                    "in the experiment contract.",
        },
        "plan_validation":
            "the plan subcommand was run on all six cells before this "
            "registration was written: 32 identities each, 0 missing "
            "telemetry, 0 censored, 0 degenerate; about 11 min per cell at "
            "10 workers.",
        "driver": "python_codes/rev62_ladder_interior_rematch.py",
        "supervisor": "python_codes/rev62_campaign.py",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r62-prereg] wrote {OUT.name}: {len(CELLS)} cells, "
          f"{len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
