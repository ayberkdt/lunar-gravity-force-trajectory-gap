"""Pre-registration of the R23 constructive-claim controls.

Written and frozen BEFORE the beta = 0.5 propagation is released and before any
R23 aggregate is inspected. It does not alter the R14 protocol or its R15
amendment, both of which stand as issued and as hashed; it adds the rules the
two R23 items need.

Why this registration exists
----------------------------
The R14/R15 protocols were written around the paper's *negative* claim -- that a
budget-calibrated radial degree history loses to a constant degree at equal
budget -- and they are strict about it. Two controls were invented for that
claim and applied to it: realized-total-work matching (R19), which removes the
per-call abstraction, and the post-hoc fixed oracle (R15-A), which represents the
constant family by its best member rather than by one nominated degree.

Neither control was ever applied to the paper's *constructive* claim, that at a
fixed budget an interior member of the span family beats both endpoints. That
claim is currently supported at beta = 0.5 under nominal per-call accounting
only, and the interior member is known from the R18 records to consume between
1.05 and 2.68 times the constant endpoint's realized work at that scale, median
1.55. A control the authors considered decisive against one claim cannot be
withheld from the other.

Disclosure -- this is not a blind registration
----------------------------------------------
An archive-based reconstruction carried out during review already indicates
outcome D below at beta = 0.5 for design A: with realized total work equalized,
the constant comparator was more accurate on 44 of 64 orbits, with a median
error ratio near 0.49 and a work ratio near 0.98. That reconstruction did not
propagate a fresh comparator, which is why R23-A is being run. This document is
therefore not a blind pre-registration and is not presented as one. Its purpose
is narrower and still worth having: it fixes, before the fresh propagation
exists, what the manuscript will say for every outcome, so that the fresh
numbers cannot be re-framed after they land.

Usage:  python rev23_preregister.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r23_preregistration.json"

PAYLOAD = {
    "schema": "r23_preregistration_v1",
    "campaign": "R23 -- controls applied to the constructive interior-optimum claim",
    "amends": ["r14_preregistration.json", "r15_preregistration_amendment.json"],
    "status": (
        "additive registration; the R14 protocol and its R15 amendment are "
        "unchanged and their hashes remain valid. Written before the beta = 0.5 "
        "propagation was released and before any R23 aggregate was inspected."),

    "blindness_disclosure": {
        "is_blind": False,
        "what_is_already_known": (
            "an archive-based reconstruction during review, without fresh "
            "propagation, indicated the constant comparator wins at beta = 0.5 "
            "on design A: 44 of 64 orbits more accurate, median error ratio "
            "about 0.49, realized work ratio about 0.98"),
        "why_the_run_still_happens": (
            "the reconstruction scored a comparator that was never propagated. "
            "The paper's own resolution rule requires a propagated trajectory "
            "against the archived truth at two tolerance levels; an inferred "
            "error is not admissible evidence under that rule, for this claim "
            "or against it"),
        "what_this_document_buys": (
            "the mapping from outcome to manuscript wording is fixed before the "
            "propagated numbers exist, so that a result already expected to be "
            "unfavorable cannot be narrowed, re-scoped, or moved to the "
            "supplement after it lands"),
    },

    "R23-A: realized-work control at the constructive budget": {
        "frozen_before": "any beta = 0.5 comparator propagation",
        "question": (
            "At beta = 0.5, does the interior span member k = 0.5 still beat a "
            "constant degree once the comparison is equalized on realized total "
            "quadratic work rather than on nominal per-call work?"),
        "protocol": (
            "the archived R19 protocol, unchanged: for each orbit read the "
            "interior member's realized total quadratic work at the tight level "
            "from the R18 record, take N* = round(N_0 * sqrt(W_k / W_0)) as the "
            "first estimate, propagate that constant degree fresh at both "
            "tolerance levels, and report the achieved work ratio as measured "
            "rather than assumed"),
        "code_identity": (
            "rev19_equal_total_work.py with a --beta argument added. The beta = 1 "
            "record paths are unchanged and the beta = 1 config hashes are "
            "unaffected, so the archived campaign is not disturbed"),
        "populations": "design A (required) and design B (run if the window allows)",
        "censoring_rule": (
            "inherited from R19: an orbit whose work-matched degree lands at or "
            "above its adopted truth degree is censored and reported, never "
            "clamped, because clamping would hand the comparison to the fixed "
            "side. At beta = 0.5 the preview shows zero orbits censored on "
            "either design"),
        "resolution_rule": (
            "inherited unchanged: |E_int - E_fix| > (env_int) + (env_fix), with "
            "unresolved comparisons reported as unresolved and never counted as "
            "a tie or a win for either side"),
        "primary_statistics_frozen": [
            "orbits, resolved, resolved interior wins, resolved fixed wins, unresolved",
            "median of the per-orbit error ratio rho = E_fix / E_int",
            "achieved total-work ratio: median, min, max",
        ],
        "statistic_convention": (
            "rho is reported as the median of per-orbit ratios. Where the "
            "manuscript compares this against the per-call accounting, the "
            "per-call number must be recomputed as a median of ratios too; the "
            "existing 'sixfold' figure is a ratio of medians and the two are not "
            "interchangeable"),
    },

    "R23-B: oracle control on the constructive claim": {
        "frozen_before": "any interior-versus-oracle aggregate being computed",
        "question": (
            "Does the interior member beat the constant-degree *family*, or only "
            "the single nominated budget-saturating degree it was compared with?"),
        "protocol": (
            "reuse the R15-A fixed-oracle archive on its 16-orbit panel at "
            "beta = 1 and score the interior member against N_oracle on the same "
            "orbits, with the same truth, the same tolerance pair and the same "
            "resolution rule"),
        "ladder_is_pre_declared": (
            "the degree ladder is the R15-A ladder, offsets [0, 1, 2, 3, 4, 6, 8, "
            "12, 16, 24] below the budget-saturating degree, frozen when R15 was "
            "written and reused here unchanged. It was not chosen after seeing "
            "any interior result"),
        "no_new_propagation_expected": (
            "the ladder errors are already archived per orbit. If and only if a "
            "convention check shows the R15 and R18 errors are not computed "
            "against the same truth on the same grid with the same statistic, "
            "the comparison is abandoned rather than patched, and the item is "
            "answered by narrowing the claim instead"),
        "symmetry_statement_required": (
            "the oracle is post-hoc for whichever side it is applied to. The "
            "manuscript must describe it identically in both places: a lower "
            "envelope over the constant family, not a selectable policy"),
    },

    "decision_logic": {
        "A_interior_wins_under_realized_work": (
            "if the interior member retains a resolved majority at beta = 0.5 "
            "under realized-work matching, the constructive claim is stated "
            "without the per-call qualifier, and Section 8.3 reports both "
            "accountings side by side"),
        "B_interior_wins_but_only_at_some_scales": (
            "state the claim with its budget range attached, naming the scales "
            "tested and the scales at which it does not hold. The abstract "
            "carries the range, not the headline"),
        "C_unresolved_dominates": (
            "if the comparison is mostly unresolved, report it as unresolved and "
            "withdraw the constructive claim to a statement about the per-call "
            "accounting only, saying so explicitly"),
        "D_constant_degree_wins_under_realized_work": (
            "report it plainly and in the main text, in the same voice the paper "
            "uses for its negative result: at beta = 0.5 the interior member's "
            "apparent advantage over the constant degree is an artifact of "
            "per-call accounting, and under realized total work the constant "
            "degree is the more accurate spend. The constructive claim is then "
            "restricted to the budget scales where it survives this control, the "
            "restriction is carried into the abstract and Section 10, and the "
            "sentence claiming the interior optimum is 'best evidenced where it "
            "matters most' is removed rather than softened"),
        "E_oracle_beats_interior": (
            "if the interior member does not beat N_oracle on the panel, the "
            "claim is written as a claim against the nominated budget-saturating "
            "constant degree, not against the constant family, and the asymmetry "
            "the paper currently relies on is removed from the radial discussion "
            "as well, so that both claims are stated at the same strength"),
    },

    "prohibited": [
        "reporting beta = 0.5 only if it is favorable",
        "moving an unfavorable realized-work result to the supplement",
        "substituting the archived constant endpoint for the work-matched "
        "comparator on any orbit whose required degree does not round back to it",
        "counting an unresolved comparison as a tie or as a win",
        "quoting a median of ratios against a ratio of medians",
        "describing the fixed oracle as operational or selectable in advance",
        "applying a control to the radial endpoint that is withheld from the "
        "interior member, or the reverse",
        "re-running either campaign with different tolerances after seeing the result",
    ],

    "numerical_contract": {
        "tight": {"rtol": 1.0e-12, "atol_position_m": 1.0e-5,
                  "atol_velocity_m_s": 1.0e-8},
        "tighter": {"rtol": 1.0e-13, "atol_position_m": 1.0e-6,
                    "atol_velocity_m_s": 1.0e-9},
        "atol_kind": "vector", "max_step_s": 60.0,
        "duration_s": 604800.0, "output_step_s": 120.0,
        "timing_comparable": False,
        "note": ("R23-A is scored on position error only. The propagation runs "
                 "under an idle-gated queue for the sake of the unrelated job "
                 "sharing the machine, not because any R23 statistic depends on "
                 "wall clock"),
    },
}


def main() -> int:
    METRICS.mkdir(parents=True, exist_ok=True)
    payload = dict(PAYLOAD)
    payload["created_utc"] = base.utc_now()
    payload["frozen_sha256"] = base.object_hash(PAYLOAD)
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r23] wrote {OUTPUT.name}")
    print(f"[r23] frozen hash {payload['frozen_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
