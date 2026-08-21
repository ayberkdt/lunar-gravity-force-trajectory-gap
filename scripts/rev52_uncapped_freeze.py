#!/usr/bin/env python3
"""Freeze the ceiling-free control on block B of the paired radial-span ladder.

R51 lifted the reference-degree ceiling on block A only, and said why: the
ordering is a property of the levels and one block carries all four of them.
That left the pooled ladder of R50 -- which is what the primary tables and the
regime map draw -- half audited, and the manuscript has to say so in four
places. This campaign closes that gap by running the same control on the second
identity block.

It is a replication, not a repeat. Block B is an independent draw of sixteen
identities at the same four apolune levels, so it can disagree with block A,
and the outcome classes below are written around that possibility rather than
around confirmation. They are frozen before the first trajectory of this
control exists; the script refuses to run if any block-B ceiling-free record is
already on disk.

Nothing under metrics/r50_* or metrics/r51_* is written. The parent design,
initial states, comparator degrees and prepass values are copied from R50 block
B rather than re-derived, exactly as R51 copied them from block A, so the only
number that differs from the capped parent is the adopted reference degree.

Usage:  python rev52_uncapped_freeze.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev10_sobol_confirmatory as base                       # noqa: E402

PARENT_BLOCK = "span_ladder_b"
PARENT_KEY = "RS2"
PARENT_DESIGN = METRICS / f"r50_{PARENT_BLOCK}_design_frozen.json"
PARENT_ROWS = METRICS / f"r50_{PARENT_BLOCK}_rows.json"
PARENT_PREREG = METRICS / "r50_preregistration.json"
PROBE_RECORD = METRICS / "r52_ceiling_probe.json"

SIBLING_PREREG = METRICS / "r51_preregistration.json"
SIBLING_VERDICT = METRICS / "r51_verdict.json"

NAME = "span_ladder_b_uncapped"
DESIGN_KEY = "RS2U"
NEW_DEGREE = 600
PARENT_DEGREE = 300
BUDGETS = [1.00, 0.75, 0.62, 0.50]

DESIGN_OUT = METRICS / f"r52_{NAME}_design_frozen.json"
ROWS_OUT = METRICS / f"r52_{NAME}_rows.json"
PREREG_OUT = METRICS / "r52_preregistration.json"

GUARD = [DESIGN_OUT, ROWS_OUT, PREREG_OUT,
         METRICS / f"r52_{NAME}_convergence.json",
         METRICS / f"r52_{NAME}_operating_point.json",
         METRICS / f"r52_budget_pareto_{NAME}.json",
         METRICS / "r52_verdict.json"]

OUTCOMES = {
    "P_replicates": (
        "Block B returns the same per-level verdicts under the raised "
        "reference as block A did, and its score orders in radial span at "
        "every budget that decides. The ceiling-free ladder is then reported "
        "pooled over both blocks, the c marks come off the regime map, and the "
        "manuscript stops scoping the control to one block."),
    "Q_replicates_ordering_only": (
        "The per-level verdicts agree with block A but the magnitudes do not: "
        "the wide-level median ratios move by a factor the two blocks do not "
        "share. The ordering is then pooled and the magnitudes are reported "
        "per block, with the disagreement named rather than averaged away."),
    "R_blocks_disagree": (
        "A level changes hands in block B where it did not in block A, or the "
        "two blocks put the turnover between different level pairs at the same "
        "budget. The ceiling-free ladder is then reported unpooled, both "
        "blocks are printed, and the span dependence is quoted as replicated "
        "in direction but not in the cell that differs. It is not resolved by "
        "running a third block."),
    "S_undecided": (
        "Too few comparisons resolve per level under the raised reference to "
        "order block B's levels. Reported with its resolved and unresolved "
        "counts, and read as leaving block B unaudited rather than as "
        "agreeing with block A."),
}

PROHIBITED = [
    "pooling block B's ceiling-free cells with block A's before the per-block "
    "verdicts are read and recorded",
    "quoting a pooled ceiling-free ratio if the blocks return outcome R",
    "dropping an orbit, changing a level, or re-deriving n_critical or n_work "
    "after propagation starts",
    "reporting a budget the clock does not reach as anything other than "
    "declared and not run",
    "revising R51's block-A numbers in the light of block B",
]


def guard() -> None:
    existing = [p.name for p in GUARD if p.exists()]
    if existing:
        raise SystemExit(
            "a block-B ceiling-free record already exists, so this "
            "registration would not be written before the data: "
            + ", ".join(existing))


def probe() -> dict:
    if not PROBE_RECORD.exists():
        raise SystemExit(
            f"{PROBE_RECORD.name} is missing. The reference degree of a "
            f"control is fixed by measured demand on the block it controls, "
            f"so run rev51_ceiling_probe.py --block b first.")
    p = json.loads(PROBE_RECORD.read_text(encoding="utf-8"))
    if p["population"] != PARENT_BLOCK:
        raise SystemExit(f"the probe on disk is for {p['population']}, not "
                         f"{PARENT_BLOCK}")
    if p["raised_reference_degree"] != NEW_DEGREE:
        raise SystemExit(f"the probe measured {p['raised_reference_degree']}, "
                         f"this control is written for {NEW_DEGREE}")
    if p["max_requested_degree_at_raised_cap"] >= NEW_DEGREE:
        raise SystemExit(
            f"the probe requests {p['max_requested_degree_at_raised_cap']} at "
            f"cap {NEW_DEGREE}: the control would still be capped on block B, "
            f"so it is not a cap-lifted control")
    if not any(p["orbits_clamped_at_parent_cap_by_level"].values()):
        raise SystemExit(
            "the probe finds no block-B orbit clamped at the parent cap, so "
            "the ceiling cannot be what produces block B's capped ratios. "
            "Record that finding instead of running this control.")
    return p


def main() -> int:
    guard()
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
        "schema": "r52_uncapped_control_design_v1",
        "population": NAME,
        "design_key": DESIGN_KEY,
        "role": ("the R50 block B span ladder re-referenced at degree 600 so "
                 "that the calibrated radial schedule is never clamped; the "
                 "replication partner of R51, not a new draw"),
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
        "replicates": {
            "campaign": "r51",
            "population": "span_ladder_a_uncapped",
            "preregistration_sha256": json.loads(
                SIBLING_PREREG.read_text(encoding="utf-8")
            )["preregistration_sha256"],
            "verdict_sha256": (base.file_hash(SIBLING_VERDICT)
                               if SIBLING_VERDICT.exists() else None),
            "note": ("block A's result is known and is pinned here by digest "
                     "so that this registration cannot later be read as "
                     "having been written against a different one"),
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
        "schema": "r52_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R52: block B of the paired radial-span ladder re-run "
                     "with its reference degree raised from 300 to 600, the "
                     "replication partner of R51"),
        "question": ("R51 removed the reference-degree ceiling on block A and "
                     "found the level ordering standing while the wide-level "
                     "magnitudes fell. The pooled ladder that the primary "
                     "tables and the regime map draw covers both blocks, so "
                     "that control audits half of what is printed. Does block "
                     "B, an independent draw of identities at the same four "
                     "levels, return the same ceiling-free result?"),
        "why_this_completes_r51": (
            "R51's own registration deferred block B on the grounds that "
            "auditing it would double the cost to sharpen a count rather than "
            "to test the claim. That was true of the claim R51 tested. It is "
            "not true of the pooled numbers the manuscript prints, which is "
            "why this campaign is run and why its scope is the pooled ladder "
            "rather than the ordering alone."),
        "why_a_paired_control": (
            "same orbits, same identities, same apolune levels, same initial "
            "states, same comparator degrees, same tolerances, same arc. One "
            "number differs from the capped parent. The indirect effect "
            "recorded for block A also applies here: raising the reference "
            "moves the calibration bisection's search, so a ceiling-free level "
            "can still move by tens of per cent, and that is reported as the "
            "resolution of the comparison rather than as a ceiling effect."),
        "not_blind": (
            "block A's ceiling-free result is known and is pinned by digest in "
            "the design record. The protections are that no orbit is drawn, no "
            "parameter is chosen, the reference degree was fixed by a demand "
            "probe run on block B rather than inherited from block A, and the "
            "outcomes below are written before the first trajectory of this "
            "control."),
        "ceiling_probe": p,
        "protocol": ("identical to R50 and R51 in every respect except the "
                     "adopted reference degree: the same seven-day arc, "
                     "tolerance levels, output grid, resolution rule, "
                     "two-policy base scope inherited from "
                     "r30_preregistration.json, R12-rule operating point "
                     "regeneration, calibration and ladder"),
        "budgets": BUDGETS,
        "budget_order": ("1.00 first, because it is the budget the geometry "
                         "sentence is written at, then 0.75, 0.62, 0.50. The "
                         "amendment budgets of R50 are not run here: the "
                         "six-budget probe puts the largest request at 599 "
                         "against a ceiling of 600, so a ceiling-free control "
                         "at those budgets would need a higher reference than "
                         "this one and is a separate campaign. A budget the "
                         "clock does not reach is reported as declared and "
                         "not run."),
        "population": NAME,
        "design": {"design_key": DESIGN_KEY,
                   "design_sha256": design["design_sha256"],
                   "file": DESIGN_OUT.name},
        "strata": {NAME: {"file": DESIGN_OUT.name,
                          "design_sha256": design["design_sha256"],
                          "seed": design.get("seed"),
                          "design_key": DESIGN_KEY}},
        "sub_boxes": {NAME: ("the R50 block B identities at the same four "
                             "apolune levels, referenced at degree 600")},
        "outcomes": OUTCOMES,
        "verdict_rule": ("per level, the R19 realized-work tally and the R14 "
                         "endpoint tally, read exactly as rev50_verdict.py "
                         "reads them, compared level by level against the "
                         "capped block-B parent and then against block A's "
                         "ceiling-free cells"),
        "reporting_commitment": (
            "the per-block ceiling-free verdicts are recorded before any "
            "pooled ceiling-free number is quoted. If the blocks disagree at a "
            "level, the ceiling-free ladder is reported unpooled and the "
            "disagreeing cell is named wherever that level is discussed."),
        "prohibited": PROHIBITED,
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")

    print(f"[r52] {DESIGN_OUT.name}  sha={design['design_sha256'][:16]}")
    print(f"      {len(orbits)} orbits, reference degree "
          f"{PARENT_DEGREE} -> {NEW_DEGREE}, levels "
          f"{design['apolune_levels_km']}")
    print(f"[r52] {ROWS_OUT.name}  rows={len(rows)}, prepass not re-run")
    print(f"[r52] {PREREG_OUT.name} "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    print(f"      probe: max requested {p['max_requested_degree_at_raised_cap']}"
          f" at cap {NEW_DEGREE}; clamped today by level "
          f"{p['orbits_clamped_at_parent_cap_by_level']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
