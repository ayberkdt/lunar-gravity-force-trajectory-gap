"""The declared departure of R53 from its own registered window.

R53 registered a cutoff and then ran past it. That is the whole of the
departure, and it is written down here for the same reason R51 and R52 wrote
theirs: a reader cannot tell an outcome-independent extension from a
result-conditional one by looking at the artefacts, so the campaign has to say
which it was and hand over the timestamps that let the reader check.

The extension has a history and the history matters, because what it produced
changed on 15 August. Three relaunches on 14 August produced nothing: the
equatorial cell finished its trajectory stage and then broke a pool worker on
every attempt at its span stage, so no ladder record was completed and the
column stood at the five cells that were complete inside the registered window.
A fourth relaunch, on 15 August, completed both remaining cells. That is the
part a reader has to be told plainly: **two of the seven cells ran outside the
registered window, after a five-cell verdict had been sealed, and the verdict
was then recomputed over all seven.**

What makes the extension reportable rather than disqualifying is that it could
not have selected an outcome. The cell list, its order and the scoring rule were
fixed in `r53_preregistration.json` before the first propagation; no population,
orbit, parameter, budget or reference degree was added; the extension resumed a
registered list rather than choosing among alternatives; and the two cells that
ran were the two the registration had already named as sixth and seventh. What
changed was the clock, and only the clock. None of that is visible in the
artefacts, which is what makes the declaration necessary rather than optional.

The registered outcome class Z, a cell the clock does not reach, no longer
fires: every declared cell ran. The defect it exposed is recorded anyway,
because it was real while it applied and because the registration's classes
still do not partition the ways this campaign could stop.

Usage:  python rev53_write_departure.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def sha(name: str) -> str:
    return hashlib.sha256((METRICS / name).read_bytes()).hexdigest()


def main() -> int:
    prereg = json.loads(
        (METRICS / "r53_preregistration.json").read_text(encoding="utf-8"))
    verdict = json.loads(
        (METRICS / "r53_verdict.json").read_text(encoding="utf-8"))

    outside = ["SE", "SF"]

    payload = {
        "schema": "r53_registration_departure_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "departs_from": "r53_preregistration.json",
        "departs_from_declared_sha256": prereg["preregistration_sha256"],
        "departs_from_file_sha256": sha("r53_preregistration.json"),
        "scored_by": "r53_verdict.json",
        "scored_by_sha256": sha("r53_verdict.json"),
        "outcome_returned": verdict["outcome"],

        "departure": (
            "the registered window was extended four times. The first three "
            "extensions, on 14 August, completed no cell. The fourth, on 15 "
            "August, completed the two remaining cells, and the verdict was "
            "recomputed over all seven. Two of the seven cells therefore ran "
            "outside the registered window and after a five-cell verdict had "
            "been sealed."),
        "what_the_registration_says":
            prereg["scheduled_window"],
        "what_happened": [
            {"local": "2026-08-14T04:51", "event":
             "supervisor started inside the registered window; the two "
             "operational-elliptical cells were already on disk"},
            {"local": "2026-08-14T06:55", "event": "design C complete"},
            {"local": "2026-08-14T07:31", "event": "high-apolune complete"},
            {"local": "2026-08-14T09:10", "event":
             "polar complete; the equatorial cell was not started, at 390 MB "
             "free on the metrics volume against the supervisor's 500 MB "
             "floor; the campaign stopped with five cells"},
            {"local": "2026-08-14T09:10", "event":
             "verdict sealed at five cells, W_monotone, with the equatorial "
             "and frozen-like cells recorded as declared and not run"},
            {"local": "2026-08-14T12:11", "event":
             "supervisor relaunched with cutoff 16:00, past the registered "
             "10:45, after disk was freed by moving a raw tree to another "
             "volume"},
            {"local": "2026-08-14T13:09", "event":
             "the equatorial cell's trajectory stage completed, 64 rows, no "
             "failures"},
            {"local": "2026-08-14T15:58", "event":
             "its span stage killed a pool worker abruptly after 226.4 min; "
             "rc=1, no record written"},
            {"local": "2026-08-14T15:59", "event":
             "second attempt at 11 workers, same failure in 0.5 min"},
            {"local": "2026-08-14T16:02", "event":
             "third attempt at 6 workers, same failure in 0.4 min; the "
             "campaign was stopped with five cells and the column as sealed"},
            {"local": "2026-08-15T00:47", "event":
             "fourth relaunch, cutoff 09:30, at 4 workers. The abrupt worker "
             "deaths were diagnosed as contention on this machine rather than "
             "as a property of the cell: the pool held at 4 where it had "
             "broken at 11, 8 and 6. The equatorial cell was resumed from its "
             "completed trajectory record rather than restarted, so no "
             "propagated orbit was recomputed and no archived record was "
             "overwritten"},
            {"local": "2026-08-15T01:54", "event":
             "equatorial cell complete, span and work-matched stages written"},
            {"local": "2026-08-15T03:53", "event":
             "frozen-like cell complete in 119.0 min; the campaign finished "
             "with all seven cells, 5.5 h inside its own cutoff"},
            {"local": "2026-08-15T03:53", "event":
             "verdict recomputed over seven cells, W_monotone unchanged; the "
             "five-cell verdict is kept beside it as "
             "r53_verdict.sealed_20260814.json"},
        ],
        "what_the_extension_produced": (
            "two ladder records, the near-equatorial and frozen-like cells, "
            "and a recomputed verdict. Both cells are quoted in the main text "
            "and drawn in the regime map. The column reported is therefore no "
            "longer the one that stood inside the registered window, and every "
            "number taken from those two cells is a number produced outside "
            "it."),
        "why_it_is_declared_anyway": (
            "the cell list, its order and the scoring rule were fixed before "
            "the first propagation, and the extension added no population, "
            "orbit, parameter, budget or reference degree. The two cells that "
            "ran were the two the registration had already named sixth and "
            "seventh, resumed rather than re-selected, so the extension could "
            "not have chosen which cells or which outcome to report. What "
            "changed was the clock. A reader cannot verify that from the "
            "artefacts, which is what makes the declaration necessary rather "
            "than optional."),

        "second_departure": (
            "outcome Z of the registration is defined by the clock, 'a cell "
            "the clock does not reach'. While it applied, it did not fit: the "
            "equatorial cell was stopped inside the registered window by a "
            "resource floor, 390 MB free against a 500 MB floor, not by the "
            "clock. Z no longer fires, since every declared cell ran, but the "
            "registration's classes still do not partition the ways this "
            "campaign could stop, and a later campaign reusing these classes "
            "should say so before it runs rather than after."),

        "cells_reported": verdict["cells_run"],
        "cells_run_outside_the_registered_window": outside,
        "cells_declared_and_not_run": verdict["cells_declared_and_not_run"],
        "verdict_recomputed_after_the_extension": True,
        "verdict_before_the_extension": "r53_verdict.sealed_20260814.json",
        "verdict_before_the_extension_sha256":
            sha("r53_verdict.sealed_20260814.json"),
        "outcome_changed_by_the_extension": False,
        "registered_classes_cover_result": False,
    }

    payload["departure_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    out = METRICS / "r53_registration_departure.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] {out.name}")
    print(f"  outcome {payload['outcome_returned']}, "
          f"{len(verdict['cells_run'])} cells, "
          f"{len(outside)} of them outside the registered window")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
