"""R50 amendment: two further budgets, declared before any ladder has run.

The registration fixed four budgets, 1.00, 0.75, 0.62 and 0.50, chosen so the
ladder is read at the budget the geometry sentence is written at and below it.
That grid is one-sided for the question this population asks. The widest level
is expected to carry its crossing below the grid, as the SH sub-box already
does; the narrowest may not cross inside it at all. With a grid that stops at
1.00 the narrow end could only be reported as "does not cross by the declared
budget", which is a statement about the grid rather than about the span.

Budgets 1.25 and 1.50 make the reading two-sided: the crossing is then located
from below the low end to above the high one, on measured cells rather than on
the edge of a grid. Both are members of the frozen budget grid the whole paper
uses, and the standard calibration already computes them, so nothing about the
protocol changes.

This is written while the base of block A is still propagating. No ladder has
run, no calibration exists, and no comparison from this population has been
computed, so the extension is declared blind in the sense that matters: it
cannot have been chosen to suit a number. The state of the campaign at the
moment of writing is recorded below from disk rather than asserted.

The two budgets are conditional on the clock, in the registered order, and after
the four declared budgets on both blocks. A budget the clock does not reach is
reported as declared and not run.

Usage:  python rev50_budget_extension_amendment.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PREREG = METRICS / "r50_preregistration.json"
OUT = METRICS / "r50_budget_extension_amendment.json"

ADDED = [1.25, 1.50]


def campaign_state() -> dict:
    """What exists on disk right now, read rather than claimed."""
    state = {}
    for name, key in (("span_ladder_a", "RS1"), ("span_ladder_b", "RS2")):
        conv = METRICS / f"r50_{name}_convergence.json"
        state[name] = {
            "prepass_rows": (METRICS / f"r50_{name}_rows.json").exists(),
            "operating_point": (METRICS
                                / f"r50_{name}_operating_point.json").exists(),
            "base_indexed_complete": (
                conv.exists()
                and bool(json.loads(conv.read_text(encoding="utf-8"))
                         .get("complete"))),
            "calibration": (METRICS
                            / f"r50_budget_pareto_{name}.json").exists(),
            "ladder_records": sorted(
                p.name for p in
                METRICS.glob(f"r19_equal_total_work_{key}_beta_*.json")),
            "trajectory_records": sorted(
                p.name for p in
                METRICS.glob(f"r14_trajectory_{key}_beta_*.json")),
        }
    return state


def main() -> int:
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    state = campaign_state()

    laddered = [r for v in state.values() for r in v["ladder_records"]]
    calibrated = [k for k, v in state.items() if v["calibration"]]
    if laddered or calibrated:
        raise SystemExit(
            "a ladder record or a calibration already exists; this amendment "
            "declares itself blind and that claim would be false: "
            f"{laddered or calibrated}")

    payload = {
        "schema": "r50_budget_extension_amendment_v1",
        "created_utc": base.utc_now(),
        "amends": PREREG.name,
        "amends_sha256": prereg["preregistration_sha256"],
        "amends_file_sha256": base.file_hash(PREREG),
        "added_budgets": ADDED,
        "registered_budgets": prereg["budgets"],
        "reason": (
            "the registered grid is one-sided for this population. The widest "
            "level is expected to carry its crossing below 0.50, as the SH "
            "sub-box already does, and the narrowest may not cross inside the "
            "grid at all; with a grid stopping at 1.00 the narrow end could "
            "only be reported as not crossing by the declared budget, which is "
            "a statement about the grid rather than about the radial span. "
            "1.25 and 1.50 make the crossing locatable from both ends on "
            "measured cells."),
        "why_this_is_not_post_hoc": (
            "no ladder record, no trajectory record and no calibration exists "
            "for either block at the moment this file is written; the state "
            "below is read from disk. The extension cannot have been chosen to "
            "suit a result because no result of this population exists."),
        "protocol_unchanged": (
            "1.25 and 1.50 are members of the frozen budget grid in "
            "r14_budget_pareto.json, which the whole paper uses, and the "
            "standard calibration computes them for every population. Nothing "
            "about the ladder, the resolution rule or the work match changes."),
        "conditionality": (
            "run after the four registered budgets on both blocks, in the "
            "order 1.25 then 1.50, and only if the clock allows. A budget the "
            "clock does not reach is reported as declared and not run."),
        "reporting": (
            "these two budgets are marked as an amendment wherever they are "
            "quoted, in the same way the post-hoc 0.62 budget is marked "
            "elsewhere in the paper, and the distinction between the four "
            "registered budgets and these two is kept in the table rather than "
            "flattened into one grid."),
        "campaign_state_at_amendment": state,
        "source": base.provenance(),
    }
    payload["amendment_sha256"] = base.object_hash(payload)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r50] {OUT.name} sha256={payload['amendment_sha256'][:16]}")
    print(f"      added {ADDED} to the registered {prereg['budgets']}")
    for name, v in state.items():
        print(f"      {name}: rows={v['prepass_rows']} op={v['operating_point']} "
              f"base={v['base_indexed_complete']} cal={v['calibration']} "
              f"ladders={len(v['ladder_records'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
