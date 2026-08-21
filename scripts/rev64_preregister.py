"""Freeze R64 (O57) before any timed propagation.

The panel is not chosen here: it is (O48)'s fourteen orbits, read from that
campaign's selection record. Budget and interior index are (O48)'s. The single
change under test is the level the timing is made at.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r64_preregistration.json"


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r64-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for n in ("rev64_interior_timing_tighter.py", "rev64_preregister.py",
              "rev48_interior_timing.py", "rev13_timing_match.py"):
        inputs[f"python_codes/{n}"] = sha(HERE / n)
    for n in ("r48_interior_timing_selection.json", "r48_interior_timing.json",
              "r12_kernel_cost_curve.json",
              "r18_span_sweep_A_beta_1.00.json",
              "r18_span_sweep_B_beta_1.00.json"):
        p = METRICS / n
        if p.exists():
            inputs[f"metrics/{n}"] = sha(p)

    payload = {
        "schema": "r64_preregistration_v1",
        "campaign": "R64 (O57): the measured-kernel-time comparator matched "
                    "at the tolerance the errors are scored at",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O48) matches the comparator on kernel time measured at the "
            "tight level while every error in the comparison is read at the "
            "tighter level. That is the level inconsistency (O42) removed for "
            "realized work, and it is why (O48) bounds cost without "
            "superseding the operation-count comparison. This campaign moves "
            "the timing to the scoring level and changes nothing else.",
        "locked_choices": {
            "panel": "(O48)'s fourteen orbits, read from "
                     "r48_interior_timing_selection.json; no orbit is "
                     "selected here",
            "member": "k = 0.5 at beta = 1, as in (O48)",
            "not_revisited": "the panel is not enlarged to the full 128 "
                             "orbits in this campaign. Whether that is worth "
                             "doing is a question this campaign answers, not "
                             "one it prejudges.",
        },
        "construction": {
            "first_pass": "degree inverted on the member's TIGHTER call "
                          "histogram against the measured cost curve, where "
                          "(O48) used the tight histogram",
            "timed_stages": "the member's contention-free re-run, the "
                            "first-pass comparator and the refined comparator "
                            "are all timed at the tighter level",
            "refinement": "second-pass degree from the measured total kernel "
                          "time ratio with c(N) ~ N^2 locally, capped at the "
                          "orbit's adopted reference degree",
            "envelope": "the refined comparator is also propagated at the "
                        "tight level, untimed, because the envelope needs the "
                        "pair",
            "idle_machine": "every timed stage refuses to start while other "
                            "python processes are alive",
        },
        "declared_outcomes": {
            "A_interior_holds": "a resolved majority for the member under "
                                "level-consistent measured time. Reading: the "
                                "operation-count caveat loosens and the "
                                "interior result can be stated on measured "
                                "time for this panel.",
            "B_constant_wins": "a resolved majority for the constant degree. "
                               "Reading: the interior advantage is an "
                               "operation-count result and the manuscript "
                               "says so on measured time as well as on the "
                               "proxy.",
            "C_unresolved_dominated": "more than half the panel undecided. "
                                      "Reading: measured time does not "
                                      "decide this comparison on fourteen "
                                      "orbits, and enlarging the timing panel "
                                      "to 128 orbits is not worth its cost "
                                      "unless the comparator degree moves "
                                      "materially.",
            "exclusivity_note": "these three are evaluated in order and the "
                                "first that holds is the one reported. "
                                "(O56)'s classes overlapped because they were "
                                "written without this rule.",
        },
        "driver": "python_codes/rev64_interior_timing_tighter.py",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r64-prereg] wrote {OUT.name}: {len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
