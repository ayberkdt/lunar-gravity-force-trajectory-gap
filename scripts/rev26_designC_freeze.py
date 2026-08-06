"""Freeze a third scrambled-Sobol design, and the registration that governs it.

R10 froze two seeds, 20260723 for design A and 20260724 for design B, and
declared B a "frozen unpropagated future replicate" before any of its orbits
were run. That declaration is what makes design B an independent replication
rather than a second look. No third seed was declared, so a design C chosen
now is chosen after seeing that A and B agree, and the choice has to be made in
a way that cannot be steered.

The rule is therefore arithmetic and public: design C takes the next integer in
the sequence the archive already contains, 20260725. It is not selected by
inspecting anything. Anyone can check it against the two seeds R10 froze.

The commitment that matters more than the rule is below in the registration:
this is the only design C that will be drawn. If its numbers disagree with
designs A and B, the disagreement is the result and it is reported as the
result. A second seed will not be tried.

This script propagates nothing. It writes the frozen design and the
registration, both hashed, so that everything downstream can be checked against
a design that existed before any of it ran.

Usage:  python rev26_designC_freeze.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SEED_A = 20260723
SEED_B = 20260724
SEED_C = SEED_B + 1          # the rule, stated as arithmetic and not as choice

DESIGN_OUT = METRICS / "r26_sobolC_design_frozen.json"
PREREG_OUT = METRICS / "r26_preregistration.json"

OUTCOMES = {
    "H_designC_agrees": (
        "Design C places the sign change in the same bracket as designs A and "
        "B. The crossing is then reported as replicated on three independent "
        "scrambled draws, and the manuscript says three rather than two."),
    "I_designC_disagrees": (
        "Design C places the sign change in a different bracket. The crossing "
        "is then reported as design-dependent, the interval quoted in the "
        "abstract is widened to cover all three designs, and the disagreement "
        "is described rather than averaged away."),
    "J_designC_unresolved": (
        "Design C resolves too few comparisons to place the crossing. "
        "Reported as a design on which the measurement does not decide, with "
        "its resolved and undecided counts, and the two-design result stands "
        "unchanged."),
}


def main() -> int:
    protocol_sha = base.object_hash(
        {"rule": "next integer after the two seeds frozen in R10",
         "seed_a": SEED_A, "seed_b": SEED_B, "seed_c": SEED_C})

    model = base.load_model(base.TRUTH_DEGREE_DEFAULT) \
        if hasattr(base, "TRUTH_DEGREE_DEFAULT") else base.load_model(300)
    design = base.make_design(SEED_C, "sobolC", "third_independent_replicate",
                              model, protocol_sha)
    design["seed_rule"] = (
        f"SEED_C = SEED_B + 1 = {SEED_C}; the next integer after the two seeds "
        f"R10 froze ({SEED_A}, {SEED_B}). Arithmetic, not selected.")
    design["propagation_status"] = "frozen_pending_base_generation"
    design["design_sha256"] = base.object_hash(
        {k: v for k, v in design.items() if k != "design_sha256"})
    DESIGN_OUT.write_text(json.dumps(design, indent=2), encoding="utf-8")

    prereg = {
        "schema": "r26_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R26: a third independent scrambled-Sobol design, drawn "
                     "to test whether the budget crossing replicates beyond "
                     "the two designs the paper already carries"),
        "seed_rule": design["seed_rule"],
        "not_blind": (
            "R10 declared two seeds and no third, so this design is drawn "
            "after seeing that designs A and B agree. That is stated rather "
            "than hidden. The protection is not blindness but the "
            "single-draw commitment below."),
        "single_draw_commitment": (
            "this is the only design C that will be drawn. If its numbers "
            "disagree with A and B, the disagreement is the result and is "
            "reported as such. No second seed will be tried, and no orbit of "
            "this design will be dropped after its numbers are known."),
        "design": {
            "seed": SEED_C,
            "family": "sobolC",
            "orbits": len(design["orbits"]),
            "generator": design["generator"],
            "design_sha256": design["design_sha256"],
            "frozen_before": ("any propagation of this design; the base "
                              "generation reads this file and does not "
                              "regenerate the points"),
        },
        "protocol": (
            "identical to designs A and B: the same truth-degree rule, the "
            "same two vector-tolerance levels, the same seven-day arc, the "
            "same output grid, and the same truth-inclusive resolution rule. "
            "The base is the truth trajectories and the critical-degree "
            "comparator; the budget ladder then runs the archived R14, R18 "
            "and R19 drivers with --design C."),
        "scope": (
            "a third design strengthens or contradicts the population-level "
            "statements only. It does not extend any claim to a new body, "
            "arc length, or force model."),
        "outcomes": OUTCOMES,
        "stopping_rule": (
            "the base is generated whenever the machine would otherwise idle. "
            "A partially generated base is reported as such and no budget "
            "result is quoted from a design whose base is incomplete."),
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")

    incl = design["realized_inclination_range_deg"]
    print(f"[r26] {DESIGN_OUT.name}  seed={SEED_C}  "
          f"design_sha256={design['design_sha256'][:16]}")
    print(f"      {len(design['orbits'])} orbits, inclination "
          f"{incl[0]:.1f}-{incl[1]:.1f} deg, "
          f"regimes {design['inclination_regime_counts']}")
    print(f"[r26] {PREREG_OUT.name}  "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
