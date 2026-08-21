"""Freeze R65 (O58) before any timed propagation.

Five things are fixed here and are not revisited after the result is seen: the
three k values, the panel, the scoring rule, the timing band, and the rule that
all three members are reported.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r65_preregistration.json"

K_VALUES = ["0.25", "0.50", "0.75"]


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r65-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for n in ("rev65_timing_family.py", "rev65_preregister.py",
              "rev64_interior_timing_tighter.py", "rev48_interior_timing.py",
              "rev13_timing_match.py"):
        inputs[f"python_codes/{n}"] = sha(HERE / n)
    for n in ("r48_interior_timing_selection.json",
              "r64_interior_timing_tighter.json",
              "r12_kernel_cost_curve.json",
              "r18_span_sweep_A_beta_1.00.json",
              "r18_span_sweep_B_beta_1.00.json"):
        p = METRICS / n
        if p.exists():
            inputs[f"metrics/{n}"] = sha(p)

    payload = {
        "schema": "r65_preregistration_v1",
        "campaign": "R65 (O58): the sampled interior family under "
                    "level-consistent measured time",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O48) and (O57) compare one member, k = 0.5, chosen because it "
            "was the most frequent sampled minimum of the seven-day nominal "
            "sweep rather than because the family has an optimum there. "
            "(O57)'s negative result is therefore a statement about k = 0.5 "
            "and not about interior allocation: it does not rule out a "
            "different concentration holding a measured-time advantage. This "
            "campaign asks the family-level question.",
        "locked_choices": {
            "k_values": K_VALUES,
            "panel": "(O57)'s fourteen orbits, themselves (O48)'s, read from "
                     "r48_interior_timing_selection.json; no orbit is "
                     "selected here",
            "budget": "beta = 1, the nominal family construction unchanged",
            "scoring_rule": "the reference-inclusive envelope rule used "
                            "everywhere else, errors at the tighter level, "
                            "resolved iff M_res > 1",
            "timing_band": "0.90 <= T_fix/T_k <= 1.10",
            "reporting": "all three members are reported. No per-orbit argmin "
                         "over k is taken and no single deployable k is "
                         "claimed from this campaign; selecting the "
                         "favourable member after the fact would make this an "
                         "oracle rather than a family sweep.",
        },
        "construction": {
            "members": "k = 0.25 and k = 0.75 are propagated here; the "
                       "k = 0.50 column is reused from (O57) unchanged and "
                       "marked as reused in the record",
            "comparator": "each member has its own constant degree matched on "
                          "measured kernel time at the tighter tolerance, "
                          "refined over integer degrees until the band is met "
                          "or the integer step cannot improve it",
            "timing_match_miss": "a cell that still misses the band after "
                                 "refinement keeps its nearest integer match "
                                 "and is flagged; (O57) had two such cells "
                                 "from a single-pass refinement and this "
                                 "campaign reports them rather than absorbing "
                                 "them",
            "idle_machine": "every timed stage refuses to start while other "
                            "python processes are alive",
        },
        "declared_outcomes": {
            "_order": "evaluated in the order listed; the first that holds is "
                      "the one reported",
            "A_other_member_wins": "k = 0.25 or k = 0.75 takes a resolved "
                                   "majority against its own time-matched "
                                   "comparator. Reading: (O57)'s negative "
                                   "result is specific to k = 0.5 and the "
                                   "measured-time optimum sits elsewhere in "
                                   "the family.",
            "B_family_shows_no_advantage": "no member takes a resolved "
                                           "majority, the resolved cells "
                                           "favouring the constant degree or "
                                           "splitting. Reading: the sampled "
                                           "interior family shows no "
                                           "measured-time advantage on this "
                                           "panel, which is a stronger "
                                           "statement than (O57)'s and is the "
                                           "one the manuscript then makes.",
            "C_scattered_optimum": "different orbits favour different members "
                                   "with no member ahead overall. Reading: a "
                                   "single global radius-only concentration "
                                   "parameter is not sufficient and the "
                                   "optimum is geometry dependent; this is "
                                   "reported as an observation about the "
                                   "family, not as a deployable rule.",
            "D_unresolved_dominated": "more than half of every member's cells "
                                      "undecided. Reading: measured time does "
                                      "not discriminate the family on "
                                      "fourteen orbits and enlarging the "
                                      "timing panel is not worth its cost.",
        },
        "driver": "python_codes/rev65_timing_family.py",
        "numbering_note": "records carry the prefix r65 because r58 is held "
                          "by the equal-work endpoint control; the "
                          "observation label is (O58).",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r65-prereg] wrote {OUT.name}: k = {', '.join(K_VALUES)}, "
          f"{len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
