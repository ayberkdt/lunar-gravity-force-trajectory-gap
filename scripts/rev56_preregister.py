"""Freeze the R56 (O56) panel, budget, policies and outcome classes before any
propagation.

Three choices are locked here and are not revisited after the result is seen:
the panel is the eight Design-A orbits that already carry a sixty-day
reference, the budget is beta = 1, and the interior policy is k = 0.5. If the
result suggests a different k would have done better, that is a different
campaign, not an amendment to this one.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r56_preregistration.json"

PANEL = [0, 1, 2, 3, 5, 6, 7, 8]


def sha(p: Path) -> dict:
    b = p.read_bytes()
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b)}


def main() -> int:
    if OUT.exists():
        print(f"[r56-prereg] {OUT.name} already exists; refusing to rewrite")
        return 1
    inputs = {}
    for name in ("rev56_longarc_interior.py", "rev56_preregister.py",
                 "rev20_span_longarc.py", "rev18_span_sweep.py",
                 "rev17_longarc60.py"):
        inputs[f"python_codes/{name}"] = sha(HERE / name)
    for name in ("r17_longarc60.json", "r20_span_longarc.json",
                 "r18_span_sweep_A_beta_1.00.json"):
        inputs[f"metrics/{name}"] = sha(METRICS / name)

    payload = {
        "schema": "r56_preregistration_v1",
        "campaign": "R56 (O56): the interior member at sixty days, "
                    "recalibrated on the sixty-day arc and matched at the "
                    "scoring tolerance",
        "registered_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "registered_before_any_propagation": True,
        "motivation":
            "(O30) reused the frozen seven-day degree table over a sixty-day "
            "arc, where the member spends a median 1.22 times the constant "
            "endpoint's per-call budget, and compared the two on the nominal "
            "per-call accounting. Its negative result therefore cannot "
            "distinguish 'the interior advantage does not survive the "
            "horizon' from 'a seven-day allocation is the wrong allocation "
            "for a sixty-day problem'. This campaign removes both confounds "
            "at once: the member is recalibrated on the sixty-day reference "
            "epochs, and the comparator is matched on realized total "
            "quadratic work read at the tolerance the errors are scored at, "
            "the convention (O42), (O53), (O54) and (O55) established.",
        "locked_choices": {
            "panel": "the eight Design-A orbits carrying a sixty-day "
                     f"reference, sobol indices {PANEL}; no population is "
                     "selected in this campaign",
            "budget": "beta = 1 only",
            "policies": "k = 0 (constant comparator) and k = 0.5 (the "
                        "interior member the manuscript already carries)",
            "not_revisited": "if the result suggests another k or another "
                             "budget would have done better, that is a "
                             "separate campaign and not an amendment here",
        },
        "construction": {
            "recalibration": "the k = 0.5 table is rebuilt with its scale "
                             "bisected so that <N_k^2> equals beta * "
                             "N_crit^2 over the sixty-day reference epochs, "
                             "tight level, the same convention (O25) applies "
                             "on seven-day epochs; only the arc changes",
            "comparator": "N* = round(N_0 sqrt(W_k/W_0)) with both works at "
                          "the tighter level, W_k from the recalibrated "
                          "member's own telemetry and W_0 from the archived "
                          "(O30) constant endpoint; the achieved ratio is "
                          "measured from the propagated runs, not assumed",
            "propagation": "both policies at both tolerance levels over 60 "
                           "days against the archived (O27) references, the "
                           "R17 contract unchanged",
            "censoring": "a comparator at or above the adopted reference "
                         "degree is censored and reported, never clamped",
            "scoring": "the reference-inclusive envelope rule used "
                       "everywhere else, errors read at the tighter level",
        },
        "declared_outcomes": {
            "A_interior_holds": "a resolved majority for k = 0.5 with median "
                                "E_fix/E_0.5 > 1. Reading: the (O30) result "
                                "was substantially an artifact of seven-day "
                                "calibration, and the manuscript's horizon "
                                "caveat becomes a statement that allocation "
                                "requires horizon-consistent calibration.",
            "B_constant_wins": "a resolved majority for the constant degree "
                               "with median ratio < 1. Reading: the interior "
                               "optimum itself moves with the horizon, which "
                               "the state-transition account of "
                               "Section VII predicts, and the constructive "
                               "interior claim is narrowed to the seven-day "
                               "horizon deliberately rather than by default.",
            "C_mixed": "the resolved subset is directionless. Reading: no "
                       "horizon statement is supported either way and the "
                       "manuscript says so.",
            "D_unresolved_dominated": "more than four of the eight orbits "
                                      "undecided. Reading: the sixty-day "
                                      "numerical envelope cannot resolve this "
                                      "comparison, which is itself a useful "
                                      "stopping result: further sixty-day "
                                      "populations would not decide it "
                                      "either.",
            "note": "all four are publishable and all four change the same "
                    "places: the sixty-day paragraph of Section VIII, the "
                    "horizon bound in the limitations, and (O30)'s entry in "
                    "the experiment contract.",
        },
        "plan_validation":
            "the plan subcommand was run on all eight orbits before this "
            "registration was written: 8 eligible, 0 missing a sixty-day "
            "reference, and the recalibrated tables meet the budget on the "
            "sixty-day epochs to within 0.07 percent.",
        "driver": "python_codes/rev56_longarc_interior.py",
        "inputs_as_read": inputs,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r56-prereg] wrote {OUT.name}: panel {PANEL}, "
          f"{len(inputs)} inputs hashed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
