"""Amendment to the R14 registration: one Phase-A calibration point at beta = 0.62.

The R25 amendment (r25_preregistration_amendment.json, 736a71ab078152fe) already
fixed the wording for a bisection step at beta = 0.62 and declared its three
outcomes before any number existed. It could not be executed, and the reason is
the subject of this file.

That amendment states the protocol as "rev14_budget_trajectory, then
rev18_span_sweep, then rev19_equal_total_work, each called with --beta 0.62 and
otherwise unchanged". That is incomplete. rev14_budget_trajectory.build_specs
does not compute a budget: it reads the calibrated Atallah tolerance and the
constant comparator degree for that budget out of the frozen Phase-A record
r14_budget_pareto.json, whose grid is {0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00}.
A budget outside that grid raises KeyError before a single orbit is propagated.
Executing the bisection therefore requires extending a frozen calibration grid,
which the parent R14 registration treats as a registration act and not a
scheduling one, and which the manuscript already says out loud in (O33).

What is being extended, and what is not. The Phase-A calibration is a
deterministic, integration-noise-free construction: for each orbit it bisects
eps_A in log space until the 10-km-binned degree history consumes
<N_A^2> = beta * N_crit^2 on the archived truth epochs, and it picks the integer
constant degree nearest the same work. It has no free parameter, no stopping
choice and no outcome in it. Running it at a new beta adds a point to a curve;
it does not add a decision. Every rule it applies -- bisection target 1%, floor
2, cap at the adopted truth degree, censoring above cap, no interpolated degree
-- is inherited verbatim from the parent registration and none is relaxed here.

Where this departs from the parent registration, stated plainly. R14's
adaptive_extension_rule permits added trajectory budgets only at grid values
bracketing a crossing, and its prohibited list includes "quoting a crossover
budget more precisely than the sampled grid supports". Both bracketing grid
values, 0.50 and 0.75, have been propagated on both designs, so that rule is
satisfied and spent. This run goes beyond it. It is therefore not a
pre-registered result and is not reported as one: the manuscript's headline
bracket stays the pre-registered (0.50, 0.75], and the 0.62 point is reported as
a declared post-hoc localization, labelled as such wherever it appears.

Why the archive cannot be edited. r14_budget_pareto.json is pinned by sha256 in
three sealed manifests (R14, R18, R21) at bbce69671cea4ec0. Merging a new budget
into it would break all three integrity gates at once. The 0.62 calibration is
written to its own record instead, and the trajectory driver is pointed at that
record in the parent process only, where build_specs runs; no archived script is
modified, since those are pinned too.

The self-check that makes the separate record admissible. The extension script
recomputes an archived budget (beta = 0.75, both designs, all 128 orbits) with
the same worker and requires it to reproduce the archived calibration exactly --
tolerance, comparator degree, achieved work and censoring flag, every orbit. If
the reproduction fails anywhere, no 0.62 record is written. The point of the
check is that a separate file computed by the same code on the same inputs is
the same object, and that claim is tested rather than asserted.

Outcomes are not restated here. They are the E/F/G of the R25 amendment, fixed
before any 0.62 number existed, and this file does not touch them.

Usage:  python rev28_calibration_amendment.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "r28_calibration_amendment.json"
R14_PREREG = METRICS / "r14_preregistration.json"
R25_AMEND = METRICS / "r25_preregistration_amendment.json"
PARETO = METRICS / "r14_budget_pareto.json"

BETA = 0.62


def main() -> int:
    r14 = json.loads(R14_PREREG.read_text(encoding="utf-8"))
    r25 = json.loads(R25_AMEND.read_text(encoding="utf-8"))

    payload = {
        "schema": "r28_calibration_amendment_v1",
        "created_utc": base.utc_now(),
        "campaign": "R28 -- Phase-A calibration extension enabling the R25 (O33) "
                    "bisection step at beta = 0.62",
        "beta": BETA,

        "amends": {
            "file": R14_PREREG.name,
            "protocol_sha256": r14["protocol_sha256"],
            "frozen_budget_grid": r14["budget_grid"],
            "executes": R25_AMEND.name,
            "executes_sha256": r25["amendment_sha256"],
        },

        "blocking_fact": (
            "rev14_budget_trajectory.build_specs reads the calibrated Atallah "
            "tolerance and the constant comparator degree out of "
            "r14_budget_pareto.json rather than computing them, so a beta "
            "outside the frozen grid raises KeyError before any orbit is "
            "propagated. The R25 amendment's protocol clause is incomplete for "
            "this reason, and this file supplies the missing step."),

        "written_after_seeing": (
            "the complete beta = 0.50 / 0.75 / 1.00 realized-work results on "
            "both designs, which bracket the sign change at (0.50, 0.75]. This "
            "amendment is post hoc and is not presented as blind."),

        "what_is_added": (
            "one deterministic, integration-noise-free Phase-A calibration "
            "point. For each of the 128 orbits it bisects eps_A in log space "
            "until the 10-km-binned degree history consumes <N_A^2> = beta * "
            "N_crit^2 on the archived truth epochs, and selects the integer "
            "constant degree nearest the same work. The construction carries no "
            "free parameter, no stopping choice and no outcome."),

        "rules_inherited_verbatim": {
            "calibration_rule": r14["calibration_rule"],
            "fixed_comparator_rule": r14["fixed_comparator_rule"],
            "censoring_rule": r14["censoring_rule"],
            "resolution_rule": r14["resolution_rule"],
            "numerical_contract": r14["numerical_contract"],
            "note": "none of these is relaxed, reparameterized or reinterpreted",
        },

        "departure_from_parent_registration": (
            "R14's adaptive_extension_rule permits added trajectory budgets only "
            "at grid values bracketing a crossing, and its prohibited list "
            "includes 'quoting a crossover budget more precisely than the "
            "sampled grid supports'. Both bracketing grid values, 0.50 and 0.75, "
            "are already propagated on both designs, so that rule is satisfied "
            "and spent. This run goes beyond it and is therefore not a "
            "pre-registered result."),

        "reporting_commitment": (
            "the manuscript's headline bracket remains the pre-registered "
            "(0.50, 0.75] on both designs. The beta = 0.62 point is reported as "
            "a declared post-hoc localization, labelled as such wherever it "
            "appears, alongside the pre-registered bracket rather than in place "
            "of it. Both halves of the bracket are equally reportable and the "
            "wording for each is already fixed in the R25 amendment's outcomes "
            "E, F and G. Nothing recorded at 0.50, 0.75, 1.00, 1.25 or 1.50 is "
            "touched, recomputed or restated by this run."),

        "single_draw_commitment": (
            "beta = 0.62 is the only calibration point this amendment adds. If "
            "its result is inconvenient, no further bisection is run and no "
            "third value is drawn: the bracket is reported at whatever width "
            "this step leaves it. A truncated chain is reported with its "
            "completion fraction rather than as a located crossing."),

        "archive_integrity": {
            "archived_calibration_record": PARETO.name,
            "archived_sha256": base.file_hash(PARETO),
            "pinned_in_manifests": ["r14_final_experiment_manifest.json",
                                    "r18_final_experiment_manifest.json",
                                    "r21_final_experiment_manifest.json"],
            "rule": (
                "the archived record is not edited and no archived script is "
                "modified, because both are sha256-pinned in sealed manifests. "
                "The new budget is written to its own record, "
                "r28_budget_pareto_beta_0.62.json, and rev14_budget_trajectory "
                "is pointed at it in the parent process only, where build_specs "
                "runs. Worker processes receive their specs through the task "
                "payload and never read either record."),
        },

        "admissibility_self_check": (
            "before writing any beta = 0.62 record, the extension script "
            "recomputes an archived budget (beta = 0.75, both designs, all 128 "
            "orbits) with the same worker and requires exact reproduction of the "
            "archived calibration -- tolerance, comparator degree, achieved work "
            "and censoring flag, every orbit. If reproduction fails anywhere, no "
            "0.62 record is written and the campaign does not start."),

        "outcomes": {
            "unchanged_from": R25_AMEND.name,
            "note": ("outcomes E (crossing in the upper half), F (crossing in "
                     "the lower half) and G (unresolved or split) were fixed "
                     "before any 0.62 number existed and are not restated, "
                     "reweighted or reordered here"),
        },
    }
    payload["amendment_sha256"] = base.object_hash(
        {k: v for k, v in payload.items() if k != "created_utc"})
    base.atomic_json(OUT, payload)
    print(f"[r28 amend] {OUT.name}")
    print(f"  amendment_sha256 = {payload['amendment_sha256']}")
    print(f"  amends R14 protocol {payload['amends']['protocol_sha256'][:16]}")
    print(f"  executes R25 amendment {payload['amends']['executes_sha256'][:16]}")
    print(f"  frozen grid = {payload['amends']['frozen_budget_grid']}")
    print(f"  adding beta = {BETA} as declared post-hoc")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
