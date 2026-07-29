"""R15 amendment to the frozen R14 pre-registration.

Written after R14 completed and before any R15 aggregate result was inspected.
It does not alter the R14 protocol, which stands as issued and as hashed; it adds
the rules the R15 items need, chiefly the failure-accounting convention that must
exist before the budget grid is ever pushed below the range R14 propagated.

Usage:  python rev15_preregister_amendment.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PARENT = METRICS / "r14_preregistration.json"
OUTPUT = METRICS / "r15_preregistration_amendment.json"

PAYLOAD = {
    "schema": "r15_preregistration_amendment_v1",
    "amends": "r14_preregistration.json",
    "status": ("additive amendment; the R14 protocol is unchanged and its hash "
               "remains valid. Written after R14 completed and before any R15 "
               "aggregate result was inspected."),

    "why": ("A second design audit raised four items the R14 protocol did not "
            "cover: whether the fixed comparator is the strongest available one, "
            "whether the budget calibration may use the reference trajectory, "
            "whether the force-defect statistic is resolved by the 120-s output "
            "grid, and how a policy failure would be counted at budgets below "
            "those R14 propagated."),

    "R15-H failure accounting": {
        "frozen_before": "any propagation at beta < 0.50",
        "status_in_R14": ("empirically moot: all 819 R14 trajectory sidecars "
                          "report status complete, with no surface impact and no "
                          "short arc, including every orbit at beta = 0.50"),
        "rules": [
            "a surface impact is a policy failure, recorded as such and never dropped",
            "survival counts are reported per budget beside every error statistic",
            "no full-arc RMS is synthesized for an arc that did not complete",
            ("where survival differs between the two policies, a "
             "common-survival-time error is reported in addition to, never "
             "instead of, the full-arc comparison"),
            "population medians state the size of the surviving subset explicitly",
        ],
    },

    "R15-A comparator hierarchy": {
        "F_near": "argmin_N |N^2 - B|, the R14 comparator; may overspend the budget",
        "F_sat": "max{N : N^2 <= B}, the budget-saturating degree; never overspends",
        "F_oracle": ("argmin over {N : N^2 <= B} of the measured seven-day error, "
                     "on a fixed ladder of offsets below F_sat declared in advance "
                     "as [0,1,2,3,4,6,8,12,16,24]"),
        "reporting_rule": ("F_oracle is a post-hoc lower envelope over the fixed "
                           "family and is never presented as a selectable policy. "
                           "Claims of the form 'better than any constant degree at "
                           "this budget' require F_oracle; claims of the form "
                           "'better than the constant degree at this budget' "
                           "require only F_sat."),
        "panel": "16 orbits, 8 per design, spread over perilune, at beta = 1",
    },

    "R15-B calibration hierarchy": {
        "oracle": ("bisection on the archived truth altitude history; the R14 "
                   "method, retained as the upper bound on budget adherence"),
        "pilot": ("bisection on the altitude history of an N = 40, "
                  "loose-tolerance pilot arc"),
        "kepler": ("bisection on the two-body altitude history implied by the "
                   "initial state alone; no propagation"),
        "reporting_rule": ("the realized budget is always measured on the true "
                           "arc, and the calibration error is reported separately "
                           "from the integrator's departure from the nominal "
                           "per-call match"),
    },

    "R15-D cadence convergence": {
        "cadences": [480.0, 240.0, 120.0, 30.0, 10.0],
        "coarse_by": "decimation of the archived 120-s grid",
        "fine_by": "re-propagation of the truth with a finer output grid",
        "acceptance": ("the 120-s grid is declared resolved if the median defect "
                       "ratio moves by less than 10% and no orbit changes which "
                       "policy has the smaller defect. Otherwise the finest "
                       "converged cadence is adopted and the force-level sweep is "
                       "recomputed on it."),
    },

    "unchanged_from_R14": [
        "the budget grid and its definition as a nominal per-call quadratic budget",
        "the truth-inclusive pairwise resolution rule",
        "the censoring convention at the degree cap",
        "the numerical contract (vector tight and tighter, 60-s max step, 120-s grid)",
        "the prohibition on calling an unresolved comparison a tie or a win",
        "the prohibition on the words optimal, near-optimal and best allocation",
    ],

    "prohibited": [
        "revising any R14 result to a different answer; R15 items may only "
        "strengthen a comparator, remove an information advantage, or test a "
        "sampling assumption",
        "reporting an R15 item in place of the R14 result it tests; both are kept",
        "extending the budget grid below 0.50 before the failure rules above are "
        "in force",
    ],
}


def main() -> int:
    payload = dict(PAYLOAD)
    payload["created_utc"] = base.utc_now()
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    payload["parent_protocol_sha256"] = parent["protocol_sha256"]
    payload["source"] = base.provenance()
    payload["amendment_sha256"] = base.object_hash(
        {k: v for k, v in payload.items()
         if k not in ("created_utc", "source", "amendment_sha256")})
    base.atomic_json(OUTPUT, payload)
    print(f"[amendment] {OUTPUT.name}")
    print(f"  parent    = {payload['parent_protocol_sha256'][:16]}")
    print(f"  amendment = {payload['amendment_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
