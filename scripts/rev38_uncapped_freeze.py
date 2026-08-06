"""R38: freeze the cap-lifted control on the operational elliptical population.

Round 6 of the review found that the one population where the paper's negative
result reverses -- the operational elliptical geometry, where the radial rule
wins 49-10 at beta = 1 with a median error ratio of 81 -- is confounded with the
degree ceiling. Thirty-eight of its sixty-four orbits reach the adopted
reference degree somewhere on the arc, and where a policy's degree reaches the
reference degree its force model *is* the reference model, so its defect
vanishes by construction. All thirty-seven resolved capped comparisons fall to
the radial rule. On the twenty-six orbits that never touch the ceiling the same
tally is 12-10 at a median ratio of 1.46, which decides nothing.

The paper already audits the ceiling wherever the radial rule loses. It has
never audited it where the radial rule wins. This registration runs that audit
the expensive way, by removing the ceiling rather than by subsetting around it.

What changes and what does not
------------------------------
Nothing about the population changes. The same sixty-four frozen orbits, the
same initial states bit for bit, the same seven-day arc, the same tolerance
levels, the same resolution rule, the same two-policy base scope, and -- this is
the point of the design -- ``n_critical`` and ``n_work`` copied unchanged from
the R31 prepass rather than recomputed. The comparator is therefore the same
comparator, not a re-derived one, and the two campaigns are paired orbit by
orbit.

The single change is the adopted reference degree, 300 -> 600. That one number
does two things at once, and it has to do both or the control is not a control:
it raises the truth the error is measured against, and it lifts the ceiling the
calibrated radial schedule is clamped to. A ceiling below the reference degree
would be arbitrary; a reference degree below the ceiling would be incoherent.

Why 600 and not 900
-------------------
Checked before propagation, on the archived altitude histories, with no
propagation of any kind: at cap 600 the equal-work radial schedule requests at
most 522 degrees, on sixty-four orbits of sixty-four, at every one of the three
budgets. ``fraction_at_cap`` is 0.000 everywhere. Six hundred clears the demand
with margin and costs four times a truth at 300; nine hundred would cost nine
times and buy nothing. The probe also confirms the confound is real rather than
formal: the orbits now asking for 418 to 522 degrees near perilune are exactly
the ones the old ceiling was clamping to 300.

Usage:
    python rev38_uncapped_freeze.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PARENT_DESIGN = METRICS / "r31_operational_elliptical_design_frozen.json"
PARENT_ROWS = METRICS / "r31_operational_elliptical_rows.json"
PARENT_PREREG = METRICS / "r31_preregistration.json"

NAME = "operational_elliptical_uncapped"
DESIGN_KEY = "OEU"
PARENT_KEY = "OE"
NEW_DEGREE = 600
PARENT_DEGREE = 300

DESIGN_OUT = METRICS / f"r38_{NAME}_design_frozen.json"
ROWS_OUT = METRICS / f"r38_{NAME}_rows.json"
PREREG_OUT = METRICS / "r38_preregistration.json"

# Measured by the pre-propagation probe described in the module docstring.
PROBE = {
    "method": ("the calibration bisection of rev14_budget_pareto re-run at "
               "cap 600 on the archived R31 truth altitude histories; no "
               "propagation, no new trajectory, nothing written to the "
               "archive"),
    "orbits_probed": 64,
    "budget_probed": 1.00,
    "max_requested_degree": 522,
    "orbits_touching_cap_600": 0,
    "all_budgets_spot_checked": [1.00, 0.75, 0.50],
    "max_requested_degree_six_lowest_perilune_all_budgets": 522,
    "conclusion": ("600 clears the demand on every orbit at every budget, so "
                   "the control removes the ceiling rather than raising it"),
}

OUTCOMES = {
    "P_ceiling_carried_it": (
        "The uncapped tally collapses towards the cap-free subset, that is, it "
        "no longer resolves in favour of the radial rule. The 49-10 result and "
        "the median ratio of 81 are then reported as artefacts of the ceiling. "
        "The uncapped tally becomes the operational-elliptical result "
        "everywhere it is quoted -- Section 8.3, the recommendation table, the "
        "regime-map figure and the conclusion -- and the capped numbers stay "
        "only as the audit that retired them."),
    "Q_survives_intact": (
        "The uncapped tally still favours the radial rule by a margin of the "
        "same order. The ceiling is then reported as having flattered the "
        "magnitude without creating the direction, both tallies are printed "
        "side by side, and the reversal stands as the paper's bound on its own "
        "negative claim."),
    "R_direction_holds_magnitude_falls": (
        "The uncapped tally favours the radial rule but the median ratio drops "
        "by an order of magnitude or more. The uncapped numbers are the ones "
        "quoted; the sentence that carries them says the magnitude was a "
        "ceiling effect and the direction was not."),
    "S_reversal": (
        "The uncapped tally favours the constant degree. The paper's one "
        "positive geometry result is withdrawn, the recommendation-table row "
        "is rewritten against it, and the negative claim is reported as "
        "holding on this geometry too."),
}

PROHIBITED = [
    "choosing between the capped and the uncapped tally after seeing both",
    "quoting the capped 49-10 or the median ratio 81 without the uncapped "
    "tally in the same sentence, whatever the outcome",
    "dropping an orbit, changing a budget or re-deriving n_critical or n_work "
    "after propagation starts",
    "pooling this population with the strata, the coverage designs or its own "
    "capped parent",
    "treating the cap-free subset of the capped run as a substitute for this "
    "control; it is a subset of a confounded comparison, not a clean one",
]


def main() -> int:
    parent = json.loads(PARENT_DESIGN.read_text(encoding="utf-8"))
    prereg31 = json.loads(PARENT_PREREG.read_text(encoding="utf-8"))
    if parent["design_sha256"] != prereg31["design"]["design_sha256"]:
        raise SystemExit("parent design does not match its own registration")

    orbits = []
    for o in parent["orbits"]:
        q = dict(o)
        if int(q["truth_degree"]) != PARENT_DEGREE:
            raise SystemExit(f"orbit {q['sobol_index']} is not a "
                             f"{PARENT_DEGREE}-degree orbit; the control "
                             f"assumes a uniform parent reference degree")
        q["truth_degree"] = NEW_DEGREE
        q["parent_truth_degree"] = PARENT_DEGREE
        orbits.append(q)

    design = dict(parent)
    design.pop("design_sha256", None)
    design.update({
        "schema": "r38_uncapped_control_design_v1",
        "population": NAME,
        "design_key": DESIGN_KEY,
        "role": ("the R31 operational elliptical population re-referenced at "
                 "degree 600 so that the calibrated radial schedule is never "
                 "clamped; a paired control, not a new draw"),
        "derived_from": {
            "population": prereg31["population"],
            "design_key": PARENT_KEY,
            "file": PARENT_DESIGN.name,
            "design_sha256": parent["design_sha256"],
            "orbits_identical": True,
            "initial_states_identical": True,
            "what_differs": ("the adopted reference degree only: "
                             f"{PARENT_DEGREE} -> {NEW_DEGREE}"),
        },
        "seed_rule": ("no draw. The seed field is the parent's and is carried "
                      "for provenance; this population samples nothing."),
        "propagation_status": "frozen_pending_base_generation",
        "ceiling_probe": PROBE,
        "orbits": orbits,
    })
    design["design_sha256"] = base.object_hash(design)
    DESIGN_OUT.write_text(json.dumps(design, indent=2), encoding="utf-8")

    # ---- rows: the parent prepass, with the adopted reference degree raised.
    parent_rows = json.loads(PARENT_ROWS.read_text(encoding="utf-8"))
    rows = []
    for r in parent_rows["rows"]:
        q = dict(r)
        if int(q["adopted_truth_degree"]) != PARENT_DEGREE:
            raise SystemExit(f"row {q['sobol_index']} does not carry the "
                             f"parent reference degree")
        q["adopted_truth_degree"] = NEW_DEGREE
        q["parent_adopted_truth_degree"] = PARENT_DEGREE
        rows.append(q)

    out = dict(parent_rows)
    out.update({
        "schema": "r11_designB_rows_v1",
        "created_utc": base.utc_now(),
        "design_frozen_sha256": base.file_hash(DESIGN_OUT),
        "script_sha256": base.file_hash(Path(__file__)),
        "truth_degree_rule": (
            f"overridden for this control: adopted truth {NEW_DEGREE} on every "
            f"orbit, in place of the archived rule's {PARENT_DEGREE}. The "
            f"schedule basis (original_truth_degree) is untouched, and neither "
            f"n_critical nor n_work is recomputed: both are the R31 prepass "
            f"values, copied, so the comparator is the same comparator."),
        "derived_from_rows": {
            "file": PARENT_ROWS.name,
            "sha256": base.file_hash(PARENT_ROWS),
            "fields_changed": ["adopted_truth_degree"],
            "fields_copied_unchanged": ["n_critical", "n_work",
                                        "original_truth_degree",
                                        "design_point", "empirical_prepass"],
        },
        "prepass_rerun": False,
        "rows": rows,
    })
    ROWS_OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")

    prereg = {
        "schema": "r38_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R38: the operational elliptical population re-run with "
                     "its reference degree raised from 300 to 600, so that the "
                     "calibrated radial schedule is never clamped to the "
                     "reference model"),
        "question": ("the radial rule wins 49-10 on this population at the "
                     "declared budget, at a median error ratio of 81, and it "
                     "is the only population where the paper's negative result "
                     "reverses. Thirty-eight of the sixty-four orbits reach "
                     "the reference degree, where the radial force model is "
                     "the reference model and its defect vanishes by "
                     "construction, and all thirty-seven resolved capped "
                     "comparisons fall to the radial rule. Does the reversal "
                     "survive when the ceiling is removed?"),
        "why_not_the_subset": (
            "the cap-free subset of the capped run answers a weaker question. "
            "It is twenty-six orbits, it resolves 12-10, and it is drawn from "
            "a comparison whose calibration was itself computed against the "
            "ceiling. Removing the ceiling recalibrates the schedule on all "
            "sixty-four orbits at equal realized work, which is the comparison "
            "the paper claims to be making."),
        "why_a_paired_control": (
            "same orbits, same initial states, same comparator degrees, same "
            "tolerances, same arc. One number differs. Anything that moves is "
            "attributable to the ceiling and to nothing else."),
        "not_blind": (
            "the capped result is known and so is the cap-free subset of it. "
            "The protections are that no orbit is drawn, no parameter is "
            "chosen, the reference degree was fixed by a demand probe rather "
            "than by preference, and the four outcomes below are written "
            "before the first trajectory."),
        "ceiling_probe": PROBE,
        "protocol": ("identical to R31 in every respect except the adopted "
                     "reference degree: seven-day arc, the same two tolerance "
                     "levels, the same output grid, the same resolution rule, "
                     "the same two-policy base scope inherited from "
                     "r30_preregistration.json, the same R12-rule operating "
                     "point regeneration, the same calibration and the same "
                     "R19 ladder"),
        "budget": ("beta = 1.00 first, because it carries the published 49-10 "
                   "and the median ratio of 81; then 0.75, then 0.50, in that "
                   "order, and only in that order"),
        "population": NAME,
        "design": {"seed": design["seed"], "design_key": DESIGN_KEY,
                   "design_sha256": design["design_sha256"],
                   "file": DESIGN_OUT.name,
                   "factor_box": design["factor_box"],
                   "derived_from": design["derived_from"]},
        "outcomes": OUTCOMES,
        "verdict_rule": ("the R19 realized-work tally, as everywhere else: "
                         "radial if resolved_atallah_wins exceeds "
                         "resolved_fixed_wins, constant if the reverse, split "
                         "if equal, undecided if nothing resolves. The "
                         "resolution rule is the archive's and is not "
                         "recomputed here."),
        "reporting_commitment": (
            "whatever comes out is printed next to the capped numbers it "
            "controls, in the main text and not only in the supplement, "
            "including the case where the reversal survives untouched. The "
            "out-of-box qualifier that R31 attached to this population travels "
            "with the control as well."),
        "prohibited": PROHIBITED,
        "parent_registration": {
            "file": PARENT_PREREG.name,
            "sha256": prereg31["preregistration_sha256"],
        },
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")

    print(f"[r38] {DESIGN_OUT.name}  key={DESIGN_KEY}  "
          f"sha={design['design_sha256'][:16]}  "
          f"orbits={len(orbits)}  reference degree {NEW_DEGREE}")
    print(f"[r38] {ROWS_OUT.name}  rows={len(rows)}  prepass reused from "
          f"{PARENT_ROWS.name}")
    print(f"[r38] {PREREG_OUT.name}  "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
