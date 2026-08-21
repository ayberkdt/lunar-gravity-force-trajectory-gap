#!/usr/bin/env python3
"""Write the hashed record of R52's departure from its own registration.

R52's registered class P carries two clauses: block B must return the same
per-level verdicts as block A, and its score must order in radial span at every
budget that decides. The first clause holds at every one of the sixteen cells.
The second does not, at two budgets.

No registered class covers that combination. Q is about magnitudes the blocks
do not share, and they do share them; R is about a level changing hands or a
turnover moving, and neither happened; S is about too few comparisons
resolving, and enough resolve. So the four classes do not partition the outcome
space this campaign landed in, which is a defect of the registration rather
than an ambiguity in the data.

This file records that, the arithmetic behind it, and the mechanism the
ordering failure is attributable to, and pins the registration, the verdict and
block A's verdict by digest. It asserts nothing new: every count is read back
out of r52_verdict.json and r51_verdict.json.

Usage:  python rev52_departure_record.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
sys.path.insert(0, str(Path(__file__).resolve().parent))

import rev10_sobol_confirmatory as base                       # noqa: E402

PREREG = METRICS / "r52_preregistration.json"
VERDICT = METRICS / "r52_verdict.json"
SIBLING = METRICS / "r51_verdict.json"
OUT = METRICS / "r52_registration_departure.json"


def main() -> int:
    for p in (PREREG, VERDICT, SIBLING):
        if not p.exists():
            raise SystemExit(f"{p.name} is missing; nothing to record")
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    v = json.loads(VERDICT.read_text(encoding="utf-8"))
    a = json.loads(SIBLING.read_text(encoding="utf-8"))

    if v["preregistration_sha256"] != prereg["preregistration_sha256"]:
        raise SystemExit("the verdict was scored against a different "
                         "registration than the one on disk")
    if v["p_ordering_clause_holds"]:
        raise SystemExit("the ordering clause holds; there is no departure "
                         "to record")

    cells = [(b, lv, c) for b, blk in v["by_budget"].items()
             for lv, c in blk["levels"].items()]
    agree = sum(1 for *_, c in cells if c.get("agrees_with_block_a"))

    # the narrowest level is where the ordering statistic breaks, and it breaks
    # against a floor: read out both blocks' scores there so the mechanism is
    # in the record rather than asserted in prose
    narrow = {}
    for b in v["budgets_read"]:
        lv = sorted(v["by_budget"][b]["levels"], key=float)[0]
        bs = v["by_budget"][b]["levels"][lv]["uncapped"]
        as_ = a["by_budget"][b]["levels"][lv]["uncapped"]
        narrow[b] = {
            "level_km": float(lv),
            "block_a_score": as_["score"], "block_a_tally":
                [as_["radial"], as_["fixed"]],
            "block_b_score": bs["score"], "block_b_tally":
                [bs["radial"], bs["fixed"]],
        }

    record = {
        "schema": "r52_registration_departure_v1",
        "created_utc": base.utc_now(),
        "departs_from": PREREG.name,
        "departs_from_declared_sha256": prereg["preregistration_sha256"],
        "departs_from_file_sha256": base.file_hash(PREREG),
        "scored_by": VERDICT.name,
        "scored_by_sha256": base.file_hash(VERDICT),
        "compared_against": SIBLING.name,
        "compared_against_sha256": base.file_hash(SIBLING),
        "outcome_returned": v["outcome"],
        "what_the_registration_says": prereg["outcomes"]["P_replicates"],
        "clause_that_holds": (
            f"the per-level verdicts match block A in {agree} of "
            f"{len(cells)} cells, the turnover falls between the same level "
            f"pair at every budget, and no wide cell differs between the "
            f"blocks by the registered factor of "
            f"{v['block_magnitude_factor']}"),
        "clause_that_fails": (
            "the span-ordering clause: Kendall tau of block B's ceiling-free "
            "level scores is "
            + ", ".join(f"{t:.3f}" for t in v["kendall_tau_uncapped"])
            + " at budgets " + ", ".join(v["budgets_read"])
            + f"; it is not unity at beta "
            + ", ".join(v["budgets_failing_p_ordering_clause"])),
        "why_no_class_fits": (
            "Q is about magnitudes the blocks do not share and they share "
            "them; R is about a level changing hands or a turnover moving and "
            "neither happened; S is about too few comparisons resolving and "
            "enough resolve. The four classes therefore do not partition the "
            "outcome space, which is a defect of this registration."),
        "mechanism": (
            "the failure is saturation at the narrowest level, not a "
            "disagreement about span. Block A's three narrow levels sit at a "
            "score of exactly -1.000, and Kendall tau ignores ties, so it "
            "returns unity. In block B a single orbit at the narrowest level "
            "favours the radial policy, lifting that level above the two "
            "above it and making the pair discordant. The parent campaign's "
            "registration names this artefact: where the constant degree wins "
            "essentially everything the score saturates against its own floor "
            "and a single orbit changing hands inverts the pair."),
        "narrowest_level_by_budget": narrow,
        "consequence_applied": (
            "the result is reported as replicating on the verdicts, the "
            "turnover and the magnitudes, and as failing the ordering clause "
            "at two budgets with the reason given. No class is selected after "
            "the fact to make the registration appear to cover it, and the "
            "ordering statistic is not quoted as evidence of replication."),
        "registered_classes_cover_result": v["registered_classes_cover_result"],
    }
    OUT.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    print(f"  outcome {record['outcome_returned']}, verdict clause holds "
          f"({agree}/{len(cells)}), ordering clause fails at beta "
          + ", ".join(v["budgets_failing_p_ordering_clause"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
