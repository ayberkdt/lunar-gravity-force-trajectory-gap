"""Supervisor for the JGCD pre-submission campaigns.

Runs the three registered campaigns to a wall clock, in the order their value
to the submission was declared, and then keeps the machine busy with a queue of
extensions rather than idling. Three rules, all learned the hard way on this
project:

  * a stage that cannot finish inside the remaining window is not started,
    because a half-propagated population supports no statement at all;
  * a stage counts as done only when its artifact parses, not when its file
    exists -- a run that died mid-write once left a zero-byte JSON that an
    existence test happily read as success;
  * the queue is ordered before the run, so what gets done when time runs out
    is decided by declared priority and not by whatever happened to be next.

Usage:
    python revJ_supervisor.py --stop-at "2026-08-09T09:30" --workers 11
    python revJ_supervisor.py --status
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"
LOG = HERE / "rJ_supervisor.log"
PROGRESS = METRICS / "rJ_supervisor_progress.json"

# Each stage: (name, argv, artifact proving completion, minutes needed, env)
STAGES = [
    ("J1.base",   ["revJ1_crossfield.py", "base"],   "rJ1_cases",       90, {}),
    ("J1.radial", ["revJ1_crossfield.py", "radial"], "rJ1_radial.json", 75, {}),
    ("J1.score",  ["revJ1_crossfield.py", "score"],  "rJ1_score.json",   8, {}),
    # Cheap and load-bearing: it is what lets J1/J2/J3 be read against the
    # manuscript's own instrument instead of a lookalike written this week.
    ("fidelity",  ["revJ_fidelity.py"],     "rJ_fidelity_check.json",   35, {}),
    # Re-scores the archived populations under the campaign's own definitions
    # so the two gravity solutions are compared under one convention, not two.
    ("compare",   ["revJ_compare.py"],      "rJ_field_comparison.json",  45, {}),
    ("J2.run",    ["revJ2_fullforce.py",  "run"],    "rJ2_cases",      170, {}),
    ("J2.score",  ["revJ2_fullforce.py",  "score"],  "rJ2_score.json",  20, {}),
    ("J3.run",    ["revJ3_tolerance.py",  "run"],    "rJ3_cases",      110, {}),
    ("J3.score",  ["revJ3_tolerance.py",  "score"],  "rJ3_score.json",   8, {}),
    # --- extensions, in declared priority order --------------------------
    # A second budget is the cheapest way to show that neither replication is
    # a property of the one budget the main text anchors on.
    # E1 ran at the original population size and was invalidated when the
    # population doubled, so it is queued again -- but low, below the two
    # controls and the primary-field ladder. Its absence costs the budget
    # ladder a row; their absence would cost the text a retraction.
    ("E2.J2beta050", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_beta_0.50.json", 170, {"JCAMP_BETA": "0.50"}),
    ("E2.J2beta050.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_beta_0.50.json", 20, {"JCAMP_BETA": "0.50"}),
    # The force-metric budget grid on the second solution. No propagation: it
    # runs on the J1 reference trajectories that already exist.
    ("E3.J1pareto", ["revJ1_crossfield.py", "pareto"],
     "rJ1_budget_pareto.json", 40, {}),
    # The budget dependence of the force ranking on the primary solution, under
    # the campaign's statistic, so the two fields can be compared like with
    # like on where the reversal turns on.
    ("E3b.primarygrid", ["revJ_budget_grid.py"], "rJ_budget_grid_primary.json",
     60, {}),
    # The full-dynamics population is the weakest count of the three campaigns,
    # so it is enlarged before anything else is broadened. The rule is nested:
    # the first 48 by Sobol index contains the first 24 already on disk.
    ("E4.J2n48", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_n48.json", 150, {"JCAMP_J2_ORBITS": "48"}),
    ("E4.J2n48.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_n48.json", 25, {"JCAMP_J2_ORBITS": "48"}),
    # Then double the cross-solution population by continuing its own sequence.
    ("E5.J1extend", ["revJ1_crossfield.py", "extend"],
     "rJ1_extension.json", 260, {}),
    # And carry the tolerance control onto the second archived population.
    ("E6.J3designB", ["revJ3_tolerance.py", "run"],
     "rJ3_run_complete_designB.json", 120, {"JCAMP_J3_DESIGN": "B"}),
    ("E6.J3designB.score", ["revJ3_tolerance.py", "score"],
     "rJ3_score_designB.json", 10, {"JCAMP_J3_DESIGN": "B"}),
    # At the budget where the force metric is a coin flip (beta ~ 0.75, from
    # E3), what does the trajectory metric say? If it is still decisive there,
    # the two metrics disagree hardest exactly where one of them stops
    # discriminating -- which is the sharpest form of the paper's claim.
    ("E7.J1beta075.base", ["revJ1_crossfield.py", "base"],
     "rJ1_base_complete_beta_0.75.json", 45, {"JCAMP_BETA": "0.75"}),
    ("E7.J1beta075", ["revJ1_crossfield.py", "radial"],
     "rJ1_radial_beta_0.75.json", 60, {"JCAMP_BETA": "0.75"}),
    ("E7.J1beta075.score", ["revJ1_crossfield.py", "score"],
     "rJ1_score_beta_0.75.json", 10, {"JCAMP_BETA": "0.75"}),
    ("E7b.J2beta075", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_beta_0.75.json", 60, {"JCAMP_BETA": "0.75"}),
    ("E7b.J2beta075.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_beta_0.75.json", 20, {"JCAMP_BETA": "0.75"}),
    # The whole confirmatory design under full dynamics, not a subset of it.
    ("E8.J2n64", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_n64.json", 90, {"JCAMP_J2_ORBITS": "64"}),
    ("E8.J2n64.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_n64.json", 30, {"JCAMP_J2_ORBITS": "64"}),
    # And a second population under full dynamics, so the full-dynamics result
    # rests on two independent designs the way the gravity-only one does.
    ("E9.J2designB", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_n64_designB.json", 240,
     {"JCAMP_J2_ORBITS": "64", "JCAMP_J2_DESIGN": "B"}),
    ("E9.J2designB.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_n64_designB.json", 30,
     {"JCAMP_J2_ORBITS": "64", "JCAMP_J2_DESIGN": "B"}),
    # Ordered by what each run removes rather than by what it adds. The two
    # controls come first because each deletes a sentence from the limitations
    # text; the primary-field ladder next, because it stops the sharpest claim
    # in the section from resting on one gravity solution; the extra budget
    # rows last, because they only lengthen a table that already makes its
    # point.
    #
    # Closes the caveat this campaign opened: the archived policy does not
    # spend the declared budget under the new dynamics. This is the other arm,
    # calibrated on the full-dynamics reference so that it does. Only the
    # radial trajectories are new; reference and comparator are shared.
    ("E13.J2recal", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_n64_recal.json", 75,
     {"JCAMP_J2_ORBITS": "64", "JCAMP_J2_RECAL": "1"}),
    ("E13.J2recal.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_n64_recal.json", 25,
     {"JCAMP_J2_ORBITS": "64", "JCAMP_J2_RECAL": "1"}),
    # Closes "one area-to-mass ratio": the same comparison with SRP eight times
    # larger. Every trajectory is new here, including the reference, because
    # the added force is felt by all of them.
    ("E14.J2highAM", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_am0.08.json", 90, {"JCAMP_J2_AOVERM": "0.08"}),
    ("E14.J2highAM.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_am0.08.json", 25, {"JCAMP_J2_AOVERM": "0.08"}),
    # The same budget ladder on the primary solution, from archived
    # trajectories. No integration: it re-scores runs that already exist.
    ("E12.primaryladder", ["revJ_ladder_primary.py"], "rJ_ladder_primary.json",
     40, {}),
    # Refills the half-budget row at the enlarged population, so every row of
    # the ladder is drawn from the same 64 orbits.
    ("E1.J1beta050.base", ["revJ1_crossfield.py", "base"],
     "rJ1_base_complete_beta_0.50.json", 25, {"JCAMP_BETA": "0.50"}),
    ("E1.J1beta050", ["revJ1_crossfield.py", "radial"],
     "rJ1_radial_beta_0.50.json", 45, {"JCAMP_BETA": "0.50"}),
    ("E1.J1beta050.score", ["revJ1_crossfield.py", "score"],
     "rJ1_score_beta_0.50.json", 10, {"JCAMP_BETA": "0.50"}),
    # A budget above the crossing, so the ladder on the second solution runs
    # from "no reversal" through the crossing to "reversal, growing".
    ("E10.J1beta125.base", ["revJ1_crossfield.py", "base"],
     "rJ1_base_complete_beta_1.25.json", 45, {"JCAMP_BETA": "1.25"}),
    ("E10.J1beta125", ["revJ1_crossfield.py", "radial"],
     "rJ1_radial_beta_1.25.json", 60, {"JCAMP_BETA": "1.25"}),
    ("E10.J1beta125.score", ["revJ1_crossfield.py", "score"],
     "rJ1_score_beta_1.25.json", 10, {"JCAMP_BETA": "1.25"}),
    ("E11.J2beta125", ["revJ2_fullforce.py", "run"],
     "rJ2_run_complete_beta_1.25.json", 60, {"JCAMP_BETA": "1.25"}),
    ("E11.J2beta125.score", ["revJ2_fullforce.py", "score"],
     "rJ2_score_beta_1.25.json", 20, {"JCAMP_BETA": "1.25"}),
    # Last, and deliberately re-entrant: the comparison re-scores only when a
    # campaign record is newer than it, and the tables are rebuilt every pass.
    # That way the final artifacts always reflect what actually finished, even
    # though the queue is re-entered many times before the window closes.
    ("final.compare", ["revJ_compare.py", "--stale-check"],
     "rJ_never_satisfied.tex", 45, {"JCAMP_COMPARE_TAG": "_final"}),
    ("final.tables", ["revJ_tables.py"], "rJ_never_satisfied.tex", 5, {}),
]

# A directory artifact is judged by how many case records it holds, because a
# directory exists from the moment the first worker starts writing into it.
# "name+k" means k more records than the stage before it left behind.
EXPECTED_CASES = {"rJ1_cases": 32 * 2 * 2, "rJ2_cases": 24 * 3 * 2,
                  "rJ3_cases": 16 * 3 * 2, "rJ2_cases+96": 24 * 3 * 2 + 96}


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def remaining_min(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def artifact_done(artifact: str) -> bool:
    """Readable and complete, not merely present."""
    if artifact in EXPECTED_CASES:
        path = METRICS / artifact.split("+")[0]
        if not path.is_dir():
            return False
        n = sum(1 for _ in path.glob("*/*.json"))
        return n >= EXPECTED_CASES[artifact]
    path = METRICS / artifact
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    if "complete" in payload:
        return bool(payload["complete"])
    return bool(payload.get("rows") or payload.get("counts"))


def record(name: str, rc: int, minutes: float) -> None:
    p = (json.loads(PROGRESS.read_text(encoding="utf-8"))
         if PROGRESS.exists() and PROGRESS.stat().st_size
         else {"schema": "rJ_supervisor_progress_v1", "stages": []})
    p["stages"].append({"stage": name, "rc": rc, "minutes": round(minutes, 1),
                        "finished_local": datetime.now().isoformat(
                            timespec="seconds")})
    p["updated_local"] = datetime.now().isoformat(timespec="seconds")
    PROGRESS.write_text(json.dumps(p, indent=2), encoding="utf-8")


def run_stage(name: str, argv: list[str], workers: int, stop_at: datetime,
              extra_env: dict) -> int:
    import os
    cmd = [sys.executable, *argv, "--workers", str(workers)]
    if len(argv) > 1 and argv[1] in ("base", "run", "radial"):
        cmd += ["--deadline", stop_at.isoformat(timespec="seconds")]
    env = {**os.environ, **extra_env}
    log(f"START {name}: {' '.join(argv)}"
        + (f"  env {extra_env}" if extra_env else ""))
    t0 = time.time()
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT, env=env)
    dt = (time.time() - t0) / 60.0
    log(f"END   {name}: rc={rc} in {dt:.1f} min")
    record(name, rc, dt)
    return rc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--skip", default="",
                    help="comma-separated stage names to treat as done")
    a = ap.parse_args()

    if a.status:
        for name, _argv, artifact, _need, _env in STAGES:
            print(f"{name:22s} {'done' if artifact_done(artifact) else 'pending':8s}"
                  f"  {artifact}")
        return 0

    if not a.stop_at:
        raise SystemExit("--stop-at is required")
    stop_at = datetime.fromisoformat(a.stop_at)
    skip = {s for s in a.skip.split(",") if s}
    log(f"supervisor: window {remaining_min(stop_at):.0f} min, "
        f"{a.workers} workers")

    for name, argv, artifact, need, extra_env in STAGES:
        if name in skip:
            log(f"SKIP  {name}: named on the command line")
            continue
        if artifact_done(artifact):
            log(f"SKIP  {name}: {artifact} already complete")
            continue
        left = remaining_min(stop_at)
        if left < need:
            log(f"STOP  {name}: needs {need} min, {left:.0f} left; "
                "refusing to start what cannot finish")
            break
        rc = run_stage(name, argv, a.workers, stop_at, extra_env)
        if rc != 0:
            log(f"WARN  {name} returned {rc}; continuing to the next stage so "
                "one failure does not idle the machine")
    log(f"supervisor: done, {remaining_min(stop_at):.0f} min left in window")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
