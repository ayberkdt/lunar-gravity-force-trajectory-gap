"""Freeze the R25 (O33) registration before the crossover run starts.

Written and hashed before stage 1 is launched. It is not blind: the direction is
partly guessable, since the interior member wins at beta = 1 and loses at
beta = 0.5, so the midpoint is expected to fall somewhere between. What the
registration fixes is the sentence that goes into the manuscript for each
outcome, including the two that narrow the claim.

Usage:  python rev25_preregister.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "r25_preregistration.json"

OUTCOMES = {
    "A_interior_wins_at_0.75": (
        "The interior member holds a resolved advantage at beta = 0.75. The "
        "sign change then lies between beta = 0.50 and beta = 0.75, and the "
        "constructive claim is stated on that interval rather than at the "
        "anchor budget alone."),
    "B_constant_wins_at_0.75": (
        "The constant degree wins on resolved comparisons at beta = 0.75. The "
        "sign change then lies between beta = 0.75 and beta = 1.00, the "
        "claim's interval is narrower than the paper currently implies, and "
        "the abstract and conclusion are tightened to say so."),
    "C_mostly_unresolved": (
        "Too few comparisons clear the summed envelope to read a direction. "
        "Reported as a limit of the measurement at this budget, not as "
        "support either way, with the resolved and undecided counts given."),
    "D_split_without_a_majority": (
        "Resolved comparisons divide without a clear majority. Reported as "
        "the crossing lying at or near beta = 0.75, with the counts, and "
        "without claiming a located optimum."),
}


def main() -> int:
    payload = {
        "schema": "r25_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R25 (O33): the realized-work comparison at the midpoint "
                     "budget beta = 0.75 on design A"),
        "question": (
            "The constructive result holds at beta = 1 and reverses at "
            "beta = 0.5. Nothing in the paper says where between them it turns "
            "over, so the claim is currently stated on an interval whose width "
            "is unmeasured. This measures the midpoint."),
        "protocol": {
            "stage_1": ("the R18 span sweep run at beta = 0.75 exactly as it "
                        "was run at every other budget, same family, same "
                        "work calibration, same two tolerance levels; it "
                        "produces the interior member and its realized work"),
            "stage_2": ("the R19 realized-work comparison run verbatim at "
                        "beta = 0.75: the comparator degree is estimated as "
                        "round[N_0 sqrt(W_k/W_0)], propagated fresh at both "
                        "levels rather than assumed, and its achieved work "
                        "ratio is measured"),
            "nothing_new": ("no new machinery is introduced. Both drivers are "
                            "the archived ones, called with a budget argument "
                            "they already accept."),
        },
        "panel": {
            "design": "A",
            "orbits": 64,
            "censoring": ("R19's existing rule: an orbit whose work-matched "
                          "degree lands at or above the adopted truth degree "
                          "is censored and listed, never clamped"),
            "design_B_not_run": (
                "design B has no archived R14 trajectory record at "
                "beta = 0.75. Building one would be a new budget campaign "
                "rather than a control, so this is a single-design result and "
                "is reported as one."),
        },
        "stopping_rule": (
            "a wall clock stop is set in advance. If it truncates either "
            "stage, the aggregate covers the completed orbits only, the "
            "completion fraction is quoted wherever the result is, and the "
            "incomplete orbits are named. The panel is walked in archived "
            "Sobol index order, which does not correlate with cost, so an "
            "early stop shortens it without tilting it."),
        "outcomes": OUTCOMES,
        "not_blind": (
            "The interior member wins at beta = 1 and loses at beta = 0.5, so "
            "some crossing between them is expected; the registration's value "
            "is that it fixes the manuscript sentence for each outcome, "
            "including B, which narrows the claim, before any number exists."),
        "reporting_commitment": (
            "the result is reported whichever way it comes out, in the main "
            "text rather than the supplement, and the abstract's statement of "
            "the budget interval is updated to match."),
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r25 prereg] {OUT.name} "
          f"sha256={payload['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
