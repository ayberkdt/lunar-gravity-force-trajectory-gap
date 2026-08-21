"""Freeze the R61 (O42-ext) cell list and outcome classes before propagation.

Mirrors rev44's registration: the cells, the order they are attempted in, the
frozen construction rules and the outcome classes are fixed here, and every
input the campaign reads is hashed as read, so a later edit to any of them
would surface rather than pass.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r61_preregistration.json"

POPS = [("C", "C"), ("low_perilune", "SL"), ("polar", "SP"),
        ("equatorial", "SE"), ("frozen_like", "SF"), ("high_apolune", "SH")]

# Bracket pass first, declared-budget pass second: every population gets the
# pair that locates its crossing before any population gets a second budget.
ORDER = ([{"population": p, "design_key": k, "beta": b, "pass": "bracket"}
          for b in (0.50, 0.75) for p, k in POPS]
         + [{"population": p, "design_key": k, "beta": 1.00,
             "pass": "declared"} for p, k in POPS])


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r61-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for name in ("rev61_equal_work_tighter_ext.py",
                 "rev44_equal_work_tighter.py", "rev18_span_sweep.py",
                 "rev14_budget_trajectory.py", "population_registry.py"):
        inputs[f"python_codes/{name}"] = sha(HERE / name)
    for cell in ORDER:
        tag = f"beta_{cell['beta']:.2f}"
        f = METRICS / f"r18_span_sweep_{cell['design_key']}_{tag}.json"
        inputs[f"metrics/{f.name}"] = sha(f)
    for pop, _ in POPS:
        f = (METRICS / "r26_designC_rows.json" if pop == "C"
             else METRICS / f"r30_{pop}_rows.json")
        inputs[f"metrics/{f.name}"] = sha(f)

    payload = {
        "schema": "r61_preregistration_v1",
        "campaign": "R61 (O42-ext): the scoring-tolerance rematch on the "
                    "third coverage design and the five geometry strata",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "R44 (O42) re-established the equal-realized-work match at the "
            "tolerance the errors are scored at, on designs A and B only, and "
            "one of those two moved its tally crossing above 0.75. Section "
            "IX.B nevertheless reads the crossing bracket across seven "
            "populations and attributes it to the matching convention as much "
            "as to the dynamics, on that two-population base; council round "
            "22 item 3(c) named the same gap from the other side. This "
            "campaign supplies the six missing populations so the convention "
            "statement is measured on the populations it generalises over "
            "rather than inferred from two.",
        "relationship_to_r44":
            "method identical, population set disjoint. R44's eight A/B cells "
            "and its sealed manifest are untouched; this campaign writes r61_* "
            "records and its own case and raw trees, and reads r18_cases "
            "read-only.",
        "cells": ORDER,
        "attempt_order":
            "as listed. The bracket pass (beta 0.50 and 0.75) runs for all six "
            "populations before the declared-budget pass (beta 1.00), so a "
            "window that ends early still leaves every population with the "
            "pair that locates its crossing. Populations are ordered C, low "
            "perilune, polar, equatorial, frozen-like, high apolune: the "
            "seven-population bracket statement covers the first five, and "
            "high apolune is the declared exception to it.",
        "frozen_rules": {
            "member": "k = 0.50 of the archived R18 span sweep, reused "
                      "unchanged",
            "match_target": "W_k^tighter = mean(N_k^2)_calls x n_RHS(k) from "
                            "the archived R18 tighter-level telemetry",
            "comparator_estimate":
                "N* = round(N_0 x sqrt(W_k^tighter / W_0^tighter)), both works "
                "at the tighter level; where the beta-specific fixed sidecar "
                "is not addressable the orbit's fixed_critical run supplies "
                "n_RHS and the source is recorded per case. This affects the "
                "first estimate only; the achieved ratio is measured from the "
                "propagated runs at both levels.",
            "propagation": "each comparator at both tolerance levels, the "
                           "R18/R19 contract unchanged",
            "censoring": "a work-matched degree at or above the orbit's "
                         "adopted reference degree is censored and reported, "
                         "never clamped",
            "scoring": "the existing reference-inclusive envelope rule, errors "
                       "at the tighter level against archived truths; no new "
                       "decision logic",
            "no_outcome_dependent_stopping":
                "cells are attempted in the registered order and each runs to "
                "completion or to the window deadline; a cell cut by the "
                "deadline is reported as partial and never as a verdict. The "
                "order does not depend on any result this campaign produces.",
        },
        "declared_outcomes": {
            "A": "no additional population moves its crossing out of "
                 "(0.50, 0.75] under level-consistent matching: the convention "
                 "caveat narrows to design B, and Section IX.B says so with "
                 "the count of populations tested.",
            "B": "one or more further populations move the crossing above "
                 "0.75: the convention reading is corroborated and reported "
                 "with the per-population brackets, and the bracket sentence "
                 "carries the convention it holds under.",
            "C": "cells are unresolved-dominated so that no bracket can be "
                 "read: the populations concerned are reported as undecidable "
                 "under this match with their resolved counts, and the "
                 "two-population caveat stands unchanged.",
            "note": "all three are publishable and all three change the same "
                    "three places: Section IX.B's bracket sentence, Section "
                    "VIII.C, and supp_strata.tex's 'exists for designs A and B "
                    "only' clause.",
        },
        "plan_validation":
            "the plan subcommand was run on all eighteen cells on 2026-08-19 "
            "before this registration was written: 64 orbits each, 0 missing "
            "telemetry, 0 censored at estimate time, 0 degenerate; estimated "
            "3.6 h single-thread per cell, about 21 min at 10 workers.",
        "driver": "python_codes/rev61_equal_work_tighter_ext.py",
        "supervisor": "python_codes/rev61_campaign.py",
        "window": "one night, supervisor stop-at 2026-08-19T11:00 local, hard "
                  "cap 12:00; cells not reached are left unrun and listed as "
                  "such.",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r61-prereg] wrote {OUT.name}: {len(ORDER)} cells, "
          f"{len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
