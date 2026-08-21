"""R51: freeze the cap-lifted control on the paired radial-span ladder.

R50 produces the largest error ratios in the campaign at its widest apolune
levels, and the degree ceiling is the first thing a referee will suspect of
producing them. The suspicion is not idle: round 6 of the review raised exactly
this against the operational elliptical population, and R38 found that lifting
the ceiling there cut the median ratio by a factor of 14.6 while leaving the
direction standing. A span ladder whose wide end is quoted at a ratio of order
ten thousand cannot go to a referee without the same audit.

What changes and what does not
------------------------------
Nothing about the population changes. The same sixty-four frozen orbits of block
A, the same initial states bit for bit, the same identities and apolune levels,
the same seven-day arc, the same tolerance levels, the same resolution rule, the
same two-policy base scope, and -- the point of the design -- ``n_critical`` and
``n_work`` copied from the R50 prepass rather than recomputed, so the comparator
is the same comparator and the two campaigns are paired orbit by orbit.

The single change is the adopted reference degree, 300 -> 600. That one number
raises the truth the error is measured against and lifts the ceiling the
calibrated radial schedule is clamped to. In this code the two are the same
number (``cap`` is ``adopted_truth_degree`` in rev14_budget_pareto.worker), so
they cannot be varied separately, and a control that varied only one would not
be coherent.

Why the probe has to have run
-----------------------------
The reference degree of a control is fixed by measured demand, not by
preference. rev51_ceiling_probe.py measures, without propagating anything, how
much of the calibrated schedule is clamped today and what the schedule asks for
once the cap is 600. This script refuses to register a control without that
record on disk, and refuses to register one at 600 if the probe says 600 does
not clear the demand or that the cap never binds. In the second case the control
would be answering a question the data says does not arise.

The probe is read per apolune level, because R50's levels differ in exactly the
way that matters here: if the cap binds only at the wide levels, then the
ceiling is confounded with the span axis itself, and that is the strongest
possible reason to run this control rather than a weaker one.

Usage:  python rev51_uncapped_freeze.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PARENT_BLOCK = "span_ladder_a"
PARENT_KEY = "RS1"
PARENT_DESIGN = METRICS / f"r50_{PARENT_BLOCK}_design_frozen.json"
PARENT_ROWS = METRICS / f"r50_{PARENT_BLOCK}_rows.json"
PARENT_PREREG = METRICS / "r50_preregistration.json"
PROBE_RECORD = METRICS / "r51_ceiling_probe.json"

NAME = "span_ladder_a_uncapped"
DESIGN_KEY = "RS1U"
NEW_DEGREE = 600
PARENT_DEGREE = 300
BUDGETS = [1.00, 0.75, 0.62, 0.50]

DESIGN_OUT = METRICS / f"r51_{NAME}_design_frozen.json"
ROWS_OUT = METRICS / f"r51_{NAME}_rows.json"
PREREG_OUT = METRICS / "r51_preregistration.json"

OUTCOMES = {
    "X_ordering_survives": (
        "The level ordering of the span ladder is unchanged: the same levels "
        "favour the same policy and the score is ordered in radial span at "
        "every budget that decides. The ceiling is then reported as having "
        "flattered the magnitudes without creating the ordering, both sets of "
        "numbers are printed side by side, and the span dependence is quoted "
        "from the uncapped run."),
    "Y_ordering_survives_magnitude_falls": (
        "The ordering is unchanged but the median ratios at the wide levels "
        "drop by an order of magnitude or more, as they did on the operational "
        "elliptical population. The uncapped numbers become the quoted ones and "
        "the sentence carrying them says the magnitude was a ceiling effect and "
        "the ordering was not."),
    "Z_ordering_breaks": (
        "The level ordering changes: a level that favoured the radial endpoint "
        "under the cap no longer does, or the score stops ordering in span. The "
        "span dependence is then reported as confounded with the ceiling, the "
        "capped ladder is retired from every place it is quoted, and the "
        "discussion says that the geometry axis could not be separated from the "
        "reference-degree ceiling by this campaign."),
    "W_undecided": (
        "Too few comparisons resolve per level under the raised reference to "
        "order the levels. Reported with its resolved and unresolved counts, "
        "and read as leaving the capped result unaudited rather than as "
        "confirming it."),
}

PROHIBITED = [
    "choosing between the capped and the uncapped ladder after seeing both",
    "quoting a capped median ratio from a wide level without the uncapped one "
    "in the same sentence, whatever the outcome",
    "dropping an orbit, changing a level, or re-deriving n_critical or n_work "
    "after propagation starts",
    "pooling this population with its capped parent, with the strata, with the "
    "coverage designs or with R31",
    "treating the subset of orbits that never touch the ceiling in the capped "
    "run as a substitute for this control",
]


def probe() -> dict:
    if not PROBE_RECORD.exists():
        raise SystemExit(
            f"{PROBE_RECORD.name} is missing. The reference degree of a "
            f"control is fixed by measured demand, so run "
            f"rev51_ceiling_probe.py first.")
    p = json.loads(PROBE_RECORD.read_text(encoding="utf-8"))
    if p["raised_reference_degree"] != NEW_DEGREE:
        raise SystemExit(f"the probe measured {p['raised_reference_degree']}, "
                         f"this control is written for {NEW_DEGREE}")
    if p["max_requested_degree_at_raised_cap"] >= NEW_DEGREE:
        raise SystemExit(
            f"the probe requests {p['max_requested_degree_at_raised_cap']} at "
            f"cap {NEW_DEGREE}: the control would still be capped, so it is "
            f"not a cap-lifted control")
    clamped = p["orbits_clamped_at_parent_cap_by_level"]
    if not any(clamped.values()):
        raise SystemExit(
            "the probe finds no orbit clamped at the parent cap, so the "
            "ceiling cannot be what produces the capped ratios and this "
            "control buys nothing. Record that finding instead of running it.")
    return p


def main() -> int:
    p = probe()
    parent = json.loads(PARENT_DESIGN.read_text(encoding="utf-8"))
    prereg50 = json.loads(PARENT_PREREG.read_text(encoding="utf-8"))
    registered = prereg50["strata"][PARENT_BLOCK]
    if parent["design_sha256"] != registered["design_sha256"]:
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
        "schema": "r51_uncapped_control_design_v1",
        "population": NAME,
        "design_key": DESIGN_KEY,
        "role": ("the R50 block A span ladder re-referenced at degree 600 so "
                 "that the calibrated radial schedule is never clamped; a "
                 "paired control, not a new draw"),
        "derived_from": {
            "population": PARENT_BLOCK,
            "design_key": PARENT_KEY,
            "file": PARENT_DESIGN.name,
            "design_sha256": parent["design_sha256"],
            "orbits_identical": True,
            "initial_states_identical": True,
            "levels_identical": parent["apolune_levels_km"],
            "what_differs": ("the adopted reference degree only: "
                             f"{PARENT_DEGREE} -> {NEW_DEGREE}"),
        },
        "seed_rule": ("no draw. The seed field is the parent's and is carried "
                      "for provenance; this population samples nothing."),
        "propagation_status": "frozen_pending_base_generation",
        "ceiling_probe": p,
        "orbits": orbits,
    })
    design["design_sha256"] = base.object_hash(design)
    DESIGN_OUT.write_text(json.dumps(design, indent=2), encoding="utf-8")

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
            f"n_critical nor n_work is recomputed: both are the R50 prepass "
            f"values, copied, so the comparator is the same comparator."),
        "derived_from_rows": {
            "file": PARENT_ROWS.name,
            "sha256": base.file_hash(PARENT_ROWS),
            "fields_changed": ["adopted_truth_degree"],
            "fields_copied_unchanged": ["n_critical", "n_work",
                                        "original_truth_degree"],
        },
        "prepass_rerun": False,
        "rows": rows,
    })
    ROWS_OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")

    prereg = {
        "schema": "r51_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R51: block A of the paired radial-span ladder re-run "
                     "with its reference degree raised from 300 to 600, so "
                     "that the calibrated radial schedule is never clamped to "
                     "the reference model"),
        "question": ("R50 finds the budget crossing moving with radial span, "
                     "and the widest level carries the largest error ratios in "
                     "the campaign. Where a policy's degree reaches the "
                     "reference degree its force model is the reference model "
                     "and its defect vanishes by construction. Does the level "
                     "ordering survive when the ceiling is removed?"),
        "why_a_paired_control": (
            "same orbits, same identities, same apolune levels, same initial "
            "states, same comparator degrees, same tolerances, same arc. One "
            "number differs, so anything that moves is attributable to the "
            "ceiling and to nothing else."),
        "why_block_A_only": (
            "the ordering is a property of the levels, and block A carries all "
            "four of them with sixty-four orbits. Auditing block B as well "
            "would double the cost to sharpen a count rather than to test the "
            "claim. If the ordering breaks here, block B's capped numbers fall "
            "with it and are reported as falling."),
        "not_blind": (
            "the capped ladder is known. The protections are that no orbit is "
            "drawn, no parameter is chosen, the reference degree was fixed by "
            "the demand probe rather than by preference, and the outcomes "
            "below are written before the first trajectory of this control."),
        "ceiling_probe": p,
        "protocol": ("identical to R50 in every respect except the adopted "
                     "reference degree: the same seven-day arc, tolerance "
                     "levels, output grid, resolution rule, two-policy base "
                     "scope inherited from r30_preregistration.json, R12-rule "
                     "operating point regeneration, calibration and ladder"),
        "budgets": BUDGETS,
        "budget_order": ("1.00 first, because it is the budget the geometry "
                         "sentence is written at, then 0.75, 0.62, 0.50. The "
                         "amendment budgets of R50 are not run here unless the "
                         "registered four are carried first, and a budget the "
                         "clock does not reach is reported as not run."),
        "population": NAME,
        "design": {"design_key": DESIGN_KEY,
                   "design_sha256": design["design_sha256"],
                   "file": DESIGN_OUT.name},
        "strata": {NAME: {"file": DESIGN_OUT.name,
                          "design_sha256": design["design_sha256"],
                          "seed": design.get("seed"),
                          "design_key": DESIGN_KEY}},
        "sub_boxes": {NAME: ("the R50 block A identities at the same four "
                             "apolune levels, referenced at degree 600")},
        "outcomes": OUTCOMES,
        "verdict_rule": ("per level, the R19 realized-work tally and the R14 "
                         "endpoint tally, read exactly as rev50_verdict.py "
                         "reads them, and compared level by level against the "
                         "capped parent"),
        "reporting_commitment": (
            "the capped and uncapped ladders are printed together wherever a "
            "wide-level ratio is quoted. If the ordering breaks, the capped "
            "ladder is retired rather than kept beside its own audit."),
        "prohibited": PROHIBITED,
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")

    print(f"[r51] {DESIGN_OUT.name}  sha={design['design_sha256'][:16]}")
    print(f"      {len(orbits)} orbits, reference degree "
          f"{PARENT_DEGREE} -> {NEW_DEGREE}, levels "
          f"{design['apolune_levels_km']}")
    print(f"[r51] {ROWS_OUT.name}  rows={len(rows)}, prepass not re-run")
    print(f"[r51] {PREREG_OUT.name} "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    print(f"      probe: max requested {p['max_requested_degree_at_raised_cap']}"
          f" at cap {NEW_DEGREE}; clamped today by level "
          f"{p['orbits_clamped_at_parent_cap_by_level']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
