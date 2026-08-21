"""Freeze R69 (O61) before any propagation.

The schedules themselves are frozen first, by `rev69_reference_free.py plan`,
and this registration hashes that state file. That order matters here more
than usual: the claim under test is that the schedule can be built without
reading a reference arc, so the artifact that proves it is the schedule table,
built and sealed before anything was flown.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r69_preregistration.json"


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r69-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    state = METRICS / "r69_reference_free_state.json"
    if not state.exists():
        print("[r69-prereg] run `rev69_reference_free.py plan` first: the "
              "frozen schedules are what this registration seals")
        return 2

    inputs = {}
    for n in ("rev69_reference_free.py", "rev69_preregister.py",
              "rev68_timing_full.py", "rev14_budget_pareto.py",
              "rev18_span_sweep.py", "rev12_atallah.py"):
        inputs[f"python_codes/{n}"] = sha(HERE / n)
    for n in ("r69_reference_free_state.json", "r68_timing_full_interior.json",
              "r68_preregistration.json", "r12_kernel_cost_curve.json",
              "r14_trajectory_A_beta_1.00.json",
              "r14_trajectory_B_beta_1.00.json",
              "r18_span_sweep_A_beta_1.00.json",
              "r18_span_sweep_B_beta_1.00.json"):
        p = METRICS / n
        if p.exists():
            inputs[f"metrics/{n}"] = sha(p)

    payload = {
        "schema": "r69_preregistration_v1",
        "campaign": "R69 (O61): the interior candidate calibrated without a "
                    "reference trajectory",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O60) established that the exploratory k = 0.5 candidate keeps "
            "its advantage at equal measured kernel time on both full "
            "designs, but its budget is calibrated on the altitude history of "
            "an already-propagated reference arc. The manuscript states that "
            "limitation: the calibration knows the radial dwell in advance. "
            "If a high-fidelity reference propagation is a precondition, the "
            "schedule is an offline comparison rather than something a "
            "propagator can be handed. This campaign asks whether the same "
            "schedule can be built from the initial osculating elements "
            "alone.",
        "locked_choices": {
            "population": "every orbit of coverage designs A and B, 64 + 64, "
                          "the same population as (O60); no selection, no "
                          "post-hoc removal",
            "member": "k = 0.50, beta = 1, the (O60) candidate",
            "what_changes": "only the altitude history the two calibrations "
                            "read. Both the radial tolerance eps_A and the "
                            "family scale s_k are bisected on the two-body "
                            "history r(E) = a0 (1 - e0 cos E) of the initial "
                            "osculating elements, sampled on the same uniform "
                            "120 s grid through Kepler's equation, instead of "
                            "on the propagated reference arc.",
            "what_does_not_change": "the family construction, the budget "
                                    "definition beta = <N^2>/N_crit^2, "
                                    "N_crit from the orbit's own perilune "
                                    "altitude, the constant comparator "
                                    "degree, the 10 km binning, the floor and "
                                    "the reference-degree cap",
            "comparator": "a constant degree refined over integers until its "
                          "measured gravity-kernel time falls within "
                          "0.95-1.05 of the member's, the (O60) protocol "
                          "unchanged, serial on an idle machine",
            "scoring": "errors at the tighter level against the same archived "
                       "reference as every other campaign, under the "
                       "reference-inclusive envelope rule. The reference is "
                       "used to score, never to build: the claim is about "
                       "constructing the schedule without foreknowledge.",
            "member_envelope": "the member has no second tolerance level in "
                               "this campaign, so its envelope is reused from "
                               "the archived k = 0.5 member at the same "
                               "budget and is recorded as reused rather than "
                               "measured. This is the one quantity in the "
                               "resolution threshold that is not measured "
                               "here.",
            "budget_adherence": "reported three ways, all against the same "
                                "target beta N_crit^2: the predicted work the "
                                "bisection met on the Kepler history, the "
                                "work the frozen table spends when sampled on "
                                "the altitude history actually flown, and the "
                                "call-weighted mean N^2 of the propagation. "
                                "The second is the primary one, because it is "
                                "the paper's own sampled convention.",
            "budget_gate": 0.05,
            "timing_match_miss": "kept, flagged, and included in the primary "
                                 "tallies, as in (O59)/(O60)",
            "ceiling_censored": "excluded from the tallies and reported",
        },
        "declared_outcomes": {
            "_order": "evaluated in the order listed; the first that holds is "
                      "the one reported",
            "A_ordering_and_budget_hold": "the member keeps a resolved "
                                          "majority on both designs and the "
                                          "median absolute sampled budget "
                                          "miss is at most 5 per cent on "
                                          "both. Reading: the schedule can be "
                                          "set from the orbital elements "
                                          "alone; a precomputed high-fidelity "
                                          "reference trajectory is not "
                                          "required, and the limitation the "
                                          "manuscript states is lifted for "
                                          "this member.",
            "B_ordering_holds_budget_degrades": "the member keeps a resolved "
                                                "majority on both designs but "
                                                "the median budget miss "
                                                "exceeds 5 per cent. Reading: "
                                                "the ordering is robust to "
                                                "reference-free calibration "
                                                "while the budget guarantee "
                                                "is not, so the comparison "
                                                "stays honest only if the "
                                                "achieved budget is reported "
                                                "with it.",
            "C_one_design_loses": "one design loses the majority. Reading: "
                                  "reference-free calibration is not "
                                  "population-independent and the result is "
                                  "reported per design with no pooled claim.",
            "D_both_lose": "both designs lose the majority. Reading: the "
                           "(O60) advantage depends on the foreknowledge in "
                           "the calibration, which is a finding about the "
                           "policy's deployability and is reported as such.",
        },
        "what_this_cannot_settle": [
            "the calibration is still made once, before the arc; it is not a "
            "causal recalibration that updates as the orbit evolves",
            "N_crit still comes from the design perilune rather than from an "
            "onboard estimate",
            "scoring still uses an archived reference trajectory",
            "the horizon stays seven days and the comparison stays "
            "single-machine",
        ],
        "driver": "python_codes/rev69_reference_free.py",
        "records": ["metrics/r69_reference_free.json",
                    "metrics/r69_reference_free_state.json"],
        "numbering_note": "r68 is held by the full-design measured-time "
                          "campaign; this one takes r69 and the observation "
                          "label (O61).",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r69-prereg] wrote {OUT.name}: {len(inputs)} inputs hashed, "
          f"schedules sealed before propagation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
