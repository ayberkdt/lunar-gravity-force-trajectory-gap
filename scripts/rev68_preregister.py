"""Freeze R68 (O59/O60) before any timed propagation.

What is fixed here and not revisited after the result is seen: the two arms,
the population, the budget, the scoring rule, the timing band, the first-guess
rule, the treatment of misses and ceiling-censored cells, and the reading of
each declared outcome.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r68_preregistration.json"


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r68-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for n in ("rev68_timing_full.py", "rev68_preregister.py",
              "rev65_timing_family.py", "rev48_interior_timing.py",
              "rev13_timing_match.py", "rev14_budget_trajectory.py"):
        inputs[f"python_codes/{n}"] = sha(HERE / n)
    for n in ("r12_kernel_cost_curve.json",
              "r14_trajectory_A_beta_1.00.json",
              "r14_trajectory_B_beta_1.00.json",
              "r18_span_sweep_A_beta_1.00.json",
              "r18_span_sweep_B_beta_1.00.json",
              "r48_interior_timing.json",
              "r65_timing_family.json"):
        p = METRICS / n
        if p.exists():
            inputs[f"metrics/{n}"] = sha(p)

    payload = {
        "schema": "r68_preregistration_v1",
        "campaign": "R68 (O59/O60): the full-design comparison at equal "
                    "measured kernel time",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "The principal population comparison is matched on a per-call "
            "quadratic proxy calibrated on the reference arc. Measured "
            "machine time has been matched only on declared fourteen-orbit "
            "panels -- (O13) for the radial endpoint, (O48)/(O57)/(O58) for "
            "the interior member -- and those panels are stated not to be "
            "population tallies. The open objection is therefore that the "
            "paper is framed around compute allocation while its principal "
            "comparison is not matched on compute as the machine spends it. "
            "This campaign runs the measured-time construction on the whole "
            "of both coverage designs.",
        "arms": {
            "endpoint": {
                "label": "(O59)",
                "member": "the budget-calibrated radial endpoint, k = 1.00, "
                          "as archived by R14 at beta = 1",
                "question": "does the constant degree keep its advantage when "
                            "the comparison is matched on measured kernel "
                            "time rather than on nominal per-call work, on "
                            "the full designs rather than on a panel?"},
            "interior": {
                "label": "(O60)",
                "member": "the sampled interior member, k = 0.50, as archived "
                          "by R18 at beta = 1",
                "question": "does the interior member's operation-count "
                            "advantage survive a measured-time match at "
                            "population scale, where (O48)/(O57) could only "
                            "report a fourteen-orbit panel with no resolved "
                            "majority?"},
        },
        "locked_choices": {
            "population": "every orbit of coverage designs A and B that "
                          "carries a scored member at beta = 1, read from "
                          "r18_span_sweep_{A,B}_beta_1.00.json: 64 + 64. No "
                          "orbit is selected, and no orbit is dropped after "
                          "its result is seen.",
            "budget": "beta = 1; the member policies are the archived ones, "
                      "unchanged and not recalibrated",
            "scoring_rule": "the reference-inclusive envelope rule used "
                            "everywhere else, errors read at the tighter "
                            "level, resolved iff the absolute error "
                            "difference exceeds the summed envelope",
            "matching_level": "measured gravity-kernel time at the tighter "
                              "level, the level the errors are read at; this "
                              "is R44's level-consistency requirement applied "
                              "to time instead of to work",
            "timing_band": "0.95 <= T_constant/T_member <= 1.05, tighter than "
                           "the 0.90-1.10 band of the panel campaigns because "
                           "a population tally is more exposed to a "
                           "systematic band edge than fourteen cells are",
            "max_refine": 5,
            "first_guess": "N_1 = c^{-1}(T_member_measured / n_RHS(constant, "
                           "archived at the tighter level)). The quantity "
                           "matched is total kernel time and a constant "
                           "degree makes fewer calls than a switching member, "
                           "so inverting the member's mean per-call cost -- "
                           "the (O48)/(O58) rule -- starts systematically low. "
                           "This changes how many propagations a cell costs, "
                           "not what is being matched.",
            "refinement_update": "the first pass steps on the quadratic "
                                 "assumption the panel campaigns used; from "
                                 "the second pass on, the step exponent is "
                                 "the slope measured between that cell's own "
                                 "last two passes, clamped to [1.2, 3.5]. "
                                 "Like the first-guess rule this changes how "
                                 "many propagations a cell costs, not what is "
                                 "matched: the accepted degree is still the "
                                 "one whose measured time lands nearest unity.",
            "reporting": "both arms are reported whatever they return, and "
                         "each design's tally is reported separately; pooling "
                         "the two designs is not a licensed summary here.",
            "timing_match_miss": "a cell that still misses the band after "
                                 "refinement keeps its nearest integer match "
                                 "and is flagged; misses are counted in the "
                                 "record and in the manuscript sentence, not "
                                 "absorbed",
            "ceiling_censored": "a cell whose matched degree reaches the "
                                "orbit's adopted reference degree cannot be "
                                "scored against that reference; it is "
                                "reported as ceiling-censored and excluded "
                                "from the tallies rather than counted as a "
                                "win for either policy",
            "idle_machine": "the pipeline refuses to start while other python "
                            "processes are alive. That check runs once, "
                            "before the first cell of an invocation, and does "
                            "not re-arm; rev68_quiet_audit.py is the "
                            "after-the-fact test for a machine joined "
                            "mid-run.",
        },
        "declared_outcomes": {
            "_order": "evaluated per arm, in the order listed; the first that "
                      "holds is the one reported",
            "A_constant_majority": "the constant degree takes a resolved "
                                   "majority on both designs. Reading for the "
                                   "endpoint arm: the ordering the manuscript "
                                   "reports at equal nominal per-call work "
                                   "holds at equal measured time on the full "
                                   "population, and the compute-allocation "
                                   "framing is matched on compute as the "
                                   "machine spends it. Reading for the "
                                   "interior arm: the interior advantage is "
                                   "an operation-count result that does not "
                                   "survive hardware time, which is the "
                                   "bound (O57)/(O58) could only place on a "
                                   "panel.",
            "B_member_majority": "the member takes a resolved majority on "
                                 "both designs. Reading for the endpoint arm: "
                                 "the reversal reported at nominal per-call "
                                 "work does not survive a measured-time "
                                 "match, and the manuscript's principal claim "
                                 "must be restated as proxy-specific. Reading "
                                 "for the interior arm: the interior member "
                                 "is a deployable measured-time improvement "
                                 "and the operation-count caveat can be "
                                 "lifted.",
            "C_split": "the two designs disagree, or the resolved cells split "
                       "without a majority. Reading: the measured-time "
                       "statement is population dependent and is reported per "
                       "design, with no pooled claim.",
            "D_unresolved_dominated": "more than half the scored cells are "
                                      "undecided on both designs. Reading: "
                                      "measured time does not discriminate "
                                      "the pair at this population and the "
                                      "existing panel-level statement remains "
                                      "the ceiling of what can be said.",
        },
        "what_this_cannot_settle": [
            "the comparison stays single-machine and single-implementation: a "
            "measured-time match is a statement about this kernel on this "
            "CPU, not a hardware-independent runtime claim",
            "the member policies are the archived budget-calibrated ones, so "
            "the calibration still knows the radial dwell of the reference "
            "arc in advance; this campaign changes what is held equal, not "
            "how the schedule was calibrated",
            "the horizon stays seven days and the geometry stays inside the "
            "sampled factor box",
        ],
        "driver": "python_codes/rev68_timing_full.py",
        "records": ["metrics/r68_timing_full_endpoint.json",
                    "metrics/r68_timing_full_interior.json"],
        "numbering_note": "r66 is held by the allocation-anatomy figure and "
                          "r67 by the review-measurement record, so this "
                          "campaign takes r68; the observation labels are "
                          "(O59) for the endpoint arm and (O60) for the "
                          "interior arm.",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r68-prereg] wrote {OUT.name}: 2 arms, {len(inputs)} inputs "
          f"hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
