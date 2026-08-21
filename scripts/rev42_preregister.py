"""R42 registration: complete the R37 level chain instead of restarting it.

R37 fixed a nested level chain (8/14/26/50/80/112/128) and a wall clock, and the
clock stopped it inside level 56 with 102 orbits in the record and 80 of them
forming the highest completed level. The two remaining levels are not a new
experiment and this registration does not pretend they are: the selection rule,
the level chain, the scoring and the worker are R37's, pinned here by digest.
What is registered is the continuation itself -- which orbits are carried rather
than recomputed, what must hold before they may be carried, and what a second
truncated level would be allowed to be called.

The honest part of this file is that it is not blind. Levels up to 40 per design
are reported: 63/63 sign agreement on the 80-orbit panel, 96/102 over the whole
record. The 26 orbits this run can reach have never been solved, and the rule
that selects them was fixed in R37 before any of them existed, but nobody should
read this registration as though the outcome were unknown.

Usage:
    python rev42_preregister.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
R37_PREREG = METRICS / "r37_preregistration.json"
R37_MANIFEST = METRICS / "r37_final_experiment_manifest.json"
R37_RECORD = METRICS / "r37_variational_extension.json"
OUT = METRICS / "r42_preregistration.json"


def main() -> int:
    r37 = json.loads(R37_PREREG.read_text(encoding="utf-8"))
    man = json.loads(R37_MANIFEST.read_text(encoding="utf-8"))
    rec = json.loads(R37_RECORD.read_text(encoding="utf-8"))
    sealed = man["result_json"]["r37_variational_extension.json"]["sha256"]
    on_disk = base.file_hash(R37_RECORD)
    if sealed != on_disk:
        print(f"[r42 prereg] ABORT: {R37_RECORD.name} does not match the R37 "
              f"manifest ({on_disk[:16]} vs {sealed[:16]}); the record this "
              f"continuation would carry is not the sealed one.")
        return 1

    payload = {
        "schema": "r42_preregistration_v1",
        "campaign": ("R42 -- completion of the R37 forced-variational level "
                     "chain to levels 56 and 64 per design (112 and 128 "
                     "orbits)"),
        "question": ("R37's wall clock stopped inside level 56. Does the "
                     "mechanism rung hold on the 26 orbits the clock never "
                     "reached, and does the panel therefore reach both full "
                     "designs?"),
        "written_before": ("any of the 26 unsolved orbits was submitted. The "
                           "26 are exactly the orbits of levels 56 and 64 that "
                           "carry no row in the sealed R37 record."),
        "not_blind": (
            "levels 4 through 40 per design are already reported: sign "
            "agreement 63/63 on the 80-orbit panel and 96/102 over the whole "
            "102-orbit record, with six disagreements all on design A. The "
            "outcome on the 26 orbits reachable here is unknown, and the rule "
            "that selects them was fixed in R37 before any orbit outside the "
            "archived eight existed."),
        "inherited_from_r37": {
            "note": ("the selection rule, the nested level chain, the stopping "
                     "rule and the three outcomes are R37's, not restated here "
                     "in weakened form. They are pinned by digest and read at "
                     "run time."),
            "r37_preregistration.json": {
                "sha256": r37["preregistration_sha256"]},
            "selection_rule": r37["selection_rule"],
            "outcomes": r37["outcomes"],
        },
        "carried_rows": {
            "statement": ("the 102 rows of the sealed R37 record are carried "
                          "forward unchanged rather than recomputed, and each "
                          "carries the provenance flag it arrived with. No "
                          "carried row is re-scored, re-selected or dropped."),
            "r37_variational_extension.json": {"sha256": sealed,
                                               "rows": len(rec["rows"])},
            "why_not_recomputed": (
                "recomputing 94 solved orbits costs about a hundred core-hours "
                "and can only reproduce them. What that cost would buy is "
                "bought instead by the two admissibility checks below."),
        },
        "admissibility_self_check": {
            "digest": ("the carried record must match the digest sealed in "
                       "r37_final_experiment_manifest.json byte for byte, "
                       "checked before the first orbit is submitted and again "
                       "at the end. A mismatch aborts without writing."),
            "recomputed_orbit": "B005",
            "archived_predicted_ratio": 2.3193085489592344,
            "abort_threshold_rel": 0.001,
            "statement": ("R37's own admissibility orbit is recomputed under "
                          "the same threshold before any new orbit is "
                          "accepted. Together the two checks say that the "
                          "carried file is the sealed one and that the current "
                          "source still computes what produced it."),
        },
        "stopping_rule": {
            "statement": ("a wall clock is fixed before the run. Orbits are "
                          "submitted in level order, so level 56 is completed "
                          "before level 64 is begun. Whatever the clock leaves "
                          "incomplete is reported as incomplete."),
            "reporting": ("the highest fully completed level is the panel. If "
                          "the clock truncates level 64 as it truncated level "
                          "56, the panel is 112 orbits and the orbits beyond "
                          "it are reported separately with the completion "
                          "fraction of the truncated level. A truncated level "
                          "is never reported as though it were the panel."),
            "concurrency_note": (
                "this run shares its machine with the R30 stratum ladders, so "
                "its worker count is smaller than the 8 R37 used and its wall "
                "clock buys fewer orbits per hour. That is a scheduling fact "
                "and not a change to the protocol: which orbits are attempted, "
                "and in what order, does not depend on it."),
        },
        "writes": ["metrics/r42_variational_completion.json"],
        "prohibited": [
            "writing metrics/r37_variational_extension.json, which is sealed "
            "under the R37 manifest and is read, never touched",
            "changing the level chain, the selection rule or the scoring after "
            "seeing a level's outcome",
            "dropping an orbit because it is slow, or because its result is "
            "inconvenient",
            "reporting a partially completed level as a completed panel",
            "quoting the carried rows as though this run had recomputed them",
        ],
        "reporting_commitment": (
            "the completed panel size is reported whichever way the sign "
            "agreement comes out, and the Section 7.2 asymmetry paragraph is "
            "rewritten to the panel this run actually completes -- not to the "
            "one it was aiming at."),
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r42 prereg] {OUT.name} "
          f"sha256={payload['preregistration_sha256'][:16]}; carrying "
          f"{len(rec['rows'])} rows, sealed digest {sealed[:16]} verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
