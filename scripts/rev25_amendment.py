"""Amendment to the R25 registration: one bisection step at beta = 0.62.

Written after the beta = 0.75 result was known, and it says so. The original
registration asked where the sign change lies and fixed the wording for four
outcomes; design A answered outcome A, putting the change between beta = 0.50
and beta = 0.75. This amendment adds the midpoint of that bracket to halve it.

The midpoint is taken as 0.62 rather than 0.625 for a naming reason, and the
reason is worth stating because it is a trap. Every driver tags its records
with f"beta_{beta:.2f}", so a run at 0.625 would write files named beta_0.62
while carrying 0.625 in the record: the file name would contradict its own
contents and would collide with any later run at 0.62. The two values split
the bracket into 0.12 and 0.13 rather than 0.125 and 0.125, which is not a
distinction any measurement here can see.

Why this is a legitimate extension and not a search. The quantity being
measured is a location, not a verdict: the run reports which half of the
bracket contains the crossing, and both halves are equally reportable. What it
cannot do is change the direction of any claim, because the claim's direction
at beta = 0.50, 0.75 and 1.00 is already fixed by records that this run does
not touch.

What it can change is how favorable the interval looks, which is why both
outcomes are written down here before the run. A crossing in the upper half
narrows the interval on which the constructive claim holds and is the less
convenient result; a crossing in the lower half widens it.

Usage:  python rev25_amendment.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "r25_preregistration_amendment.json"
PARENT = METRICS / "r25_preregistration.json"

OUTCOMES = {
    "E_crossing_in_the_upper_half": (
        "The constant degree wins on resolved realized-work comparisons at "
        "beta = 0.62. The sign change then lies between 0.62 and 0.75, the "
        "interval on which the constructive claim holds is narrower than the "
        "beta = 0.75 result alone suggested, and the manuscript states the "
        "bracket as (0.62, 0.75]."),
    "F_crossing_in_the_lower_half": (
        "The interior member wins on resolved realized-work comparisons at "
        "beta = 0.62. The sign change then lies between 0.50 and 0.62, the "
        "claim holds on a wider interval than the two archived budgets could "
        "show, and the manuscript states the bracket as (0.50, 0.62]."),
    "G_unresolved_or_split": (
        "Too few comparisons resolve, or they divide without a majority. The "
        "bracket is then reported as unchanged at (0.50, 0.75) with the "
        "midpoint counts given, and the sampled grid is stated as unable to "
        "localize further."),
}


def main() -> int:
    parent = (json.loads(PARENT.read_text(encoding="utf-8"))
              if PARENT.exists() else {})
    payload = {
        "schema": "r25_preregistration_amendment_v1",
        "created_utc": base.utc_now(),
        "amends": {
            "file": PARENT.name,
            "preregistration_sha256": parent.get("preregistration_sha256"),
            "unchanged": ("the parent registration's protocol, panel and "
                          "stopping rule apply verbatim to this budget; only "
                          "the value of beta is new"),
        },
        "campaign": "R25 (O33) bisection step: beta = 0.62",
        "written_after_seeing": (
            "the beta = 0.75 design-A result (27 interior against 15 fixed on "
            "42 resolved comparisons, median ratio 1.62), which placed the "
            "sign change between beta = 0.50 and beta = 0.75. This amendment "
            "is therefore not blind and is not presented as blind."),
        "why_this_is_localization_not_search": (
            "the reported quantity is which half of an existing bracket holds "
            "the crossing. Both halves are reportable, and neither changes the "
            "direction of any verdict already recorded at 0.50, 0.75 or 1.00, "
            "none of which this run touches. The stake is the width of the "
            "interval, and the wording for both widths is fixed below before "
            "any number exists."),
        "protocol": (
            "the full chain, since no design has an archived budget record at "
            "this value: rev14_budget_trajectory, then rev18_span_sweep, then "
            "rev19_equal_total_work, each called with --beta 0.62 and "
            "otherwise unchanged. Design A first; design B only if the window "
            "allows a complete chain."),
        "outcomes": OUTCOMES,
        "reporting_commitment": (
            "the bracket is reported at whatever width the run leaves it, in "
            "the same sentence of the main text that currently states it, and "
            "a truncated chain is reported with its completion fraction rather "
            "than as a located crossing."),
    }
    payload["amendment_sha256"] = base.object_hash(payload)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r25 amend] {OUT.name} sha256={payload['amendment_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
