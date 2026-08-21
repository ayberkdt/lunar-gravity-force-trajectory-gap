#!/usr/bin/env python3
"""Write the hashed record of R51's departure from its own registration.

The registered rule scored the ceiling-free control as outcome Z, whose
consequence is that the capped ladder be withdrawn everywhere it is quoted and
the geometry axis reported as inseparable from the reference-degree ceiling.
That consequence was not applied. The supplement says so in prose, but every
other declared deviation in this paper carries a machine-readable record whose
digest a manifest seals (R15, R23, R25, R28, R37, R50), and a departure argued
only in prose is the one a reader cannot check.

This file records what fired, what was done instead, and the numbers the
decision rests on, and it pins the registration and verdict it refers to by
digest so it cannot drift away from them. It asserts nothing new: every count
here is read back out of r51_verdict.json rather than typed in.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev10_sobol_confirmatory as base                       # noqa: E402

PREREG = METRICS / "r51_preregistration.json"
VERDICT = METRICS / "r51_verdict.json"
OUT = METRICS / "r51_registration_departure.json"


def main() -> int:
    for p in (PREREG, VERDICT):
        if not p.exists():
            raise SystemExit(f"{p.name} is missing; nothing to record")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    verdict = json.loads(VERDICT.read_text(encoding="utf-8"))

    # The verdict carries the registration's self-declared content digest, the
    # one the freezing script computed over the protocol before the file
    # existed; that is a different quantity from the digest of the file on
    # disk, and both are recorded below so neither can be read for the other.
    if verdict["preregistration_sha256"] != prereg["preregistration_sha256"]:
        raise SystemExit("the verdict was scored against a different "
                         "registration than the one on disk")

    cells = [(b, lvl, c)
             for b, blk in verdict["by_budget"].items()
             for lvl, c in blk["levels"].items()]
    changed = [(b, lvl, c) for b, lvl, c in cells if c["verdict_changed"]]
    tally_moved = [(b, lvl) for b, lvl, c in cells
                   if (c["capped"]["radial"], c["capped"]["fixed"])
                   != (c["uncapped"]["radial"], c["uncapped"]["fixed"])]

    record = {
        "schema": "r51_registration_departure_v1",
        "created_utc": base.utc_now(),
        "departs_from": PREREG.name,
        "departs_from_declared_sha256": prereg["preregistration_sha256"],
        "departs_from_file_sha256": base.file_hash(PREREG),
        "scored_by": VERDICT.name,
        "scored_by_sha256": base.file_hash(VERDICT),
        "outcome_returned": verdict["outcome"],
        "trigger_that_fired": (
            "a level favouring the radial endpoint under the ceiling no longer "
            "does: " + "; ".join(
                f"beta={b}, {float(lvl):.0f} km, {c['capped']['verdict']} -> "
                f"{c['uncapped']['verdict']}" for b, lvl, c in changed)),
        "trigger_that_did_not_fire": (
            "the score stops ordering in span; Kendall tau of the ceiling-free "
            "level scores is "
            + ", ".join(str(t) for t in verdict["kendall_tau_uncapped"])
            + " at the budgets read"),
        "registered_consequence": verdict["outcome_text"],
        "consequence_applied": False,
        "what_was_done_instead": (
            "both runs are printed together; the ceiling-free values are the "
            "ones interpreted quantitatively in prose; the cell that changed "
            "hands is named wherever the lowest budget is discussed; and the "
            "scope of the control is stated wherever its numbers are used."),
        "superseded_by": "r52_registration_departure.json",
        "superseded_note": (
            "written while the control covered block A only. R52 then ran the "
            "same control on block B, and the sentence above originally read "
            "that the pooled tables and the regime map keep the capped values. "
            "That is no longer what happens: the regime map draws the ladder "
            "rows from the ceiling-free records at the four registered budgets "
            "and marks the two amended budgets, which have no ceiling-free "
            "run, cell by cell. The ladder tables stay capped because they "
            "carry all six budgets. This record is kept as it was written "
            "apart from that one correction, because a departure record edited "
            "to match later evidence stops being a record of what was decided "
            "when."),
        "why": (
            "the registration attached one consequence to two triggers of "
            "different severity. Withdrawing the ladder is the response to the "
            "ordering failing, which did not happen. Applying it to a single "
            "cell changing hands would retire a result the control itself "
            "leaves standing."),
        "evidence_the_decision_rests_on": {
            "cells_covered": len(cells),
            "verdicts_unchanged": len(cells) - len(changed),
            "verdicts_changed": len(changed),
            "tallies_moved": len(tally_moved),
            "kendall_tau_uncapped": verdict["kendall_tau_uncapped"],
            "median_rho_factor_wide_levels":
                verdict["median_rho_factor_wide_levels"],
        },
        "what_the_departure_does_not_cover": (
            "block B and the two amended budgets have no ceiling-free "
            "counterpart, so the departure is scoped to block A at the "
            "budgets in r51_verdict.json and is not a statement about the "
            "cells it does not reach."),
        "reported_in": "supplement, span-ladder section, ceiling-removed "
                       "subsection",
    }
    OUT.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    print(f"  outcome {record['outcome_returned']}, "
          f"{record['evidence_the_decision_rests_on']['verdicts_changed']} of "
          f"{len(cells)} verdicts changed, "
          f"{len(tally_moved)} tallies moved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
