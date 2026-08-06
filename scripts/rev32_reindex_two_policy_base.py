"""Assemble the index of a two-policy base, which the pinned driver cannot.

What went wrong, stated plainly, because the fix only makes sense next to it.
The stratum bases were run with four of the six convergence-study policies
switched off, since the budget ladder reads only the truth and the
critical-degree comparator and dropping the other four is a third of the
propagation cost. The propagation honoured that. The *index assembly* did not:
rev11_full_convergence.orbit_summary computes the convergence study's
schedule-versus-comparator table from a hard-coded list of the four policies
that were switched off, so it raised on every orbit, every base was written with
complete=False, and every supervisor then did the correct thing with an
incomplete base -- it refused to calibrate or propagate a budget on it.

So the trajectories are all on disk and are exactly the ones the ladder needs.
What is missing is the index that says so.

This rebuilds that index from the sidecars already written. It propagates
nothing. orbit_summary is reproduced with the comparison block removed rather
than the block being made conditional in the pinned driver, because that driver
is sha256-pinned in the R11 manifest; everything else in the summary, including
the truth-inclusive envelope the resolution rule uses, is the archived
computation on the archived arrays.

An index built this way carries an explicit marker: the comparisons dictionary
is empty and the record says why. A stratum result quoted from a schedule
policy would then fail loudly rather than silently read a zero.

Usage:
    python rev32_reindex_two_policy_base.py --stratum polar
    python rev32_reindex_two_policy_base.py --registry r31 \
        --stratum operational_elliptical
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def _from_argv(flag: str, default: str | None = None) -> str:
    for i, arg in enumerate(sys.argv):
        if arg == flag and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if arg.startswith(flag + "="):
            return arg.split("=", 1)[1]
    if default is None:
        raise SystemExit(f"{flag} is required")
    return default


import population_registry as registry                      # noqa: E402

REGISTRY = _from_argv("--registry", "r30")
STRATUM = _from_argv("--stratum")
SPEC = registry.spec(REGISTRY, STRATUM)
ROWS = METRICS / f"{REGISTRY}_{STRATUM}_rows.json"
OUTPUT = METRICS / f"{REGISTRY}_{STRATUM}_convergence.json"
SMOKE_OUTPUT = METRICS / f"{REGISTRY}_{STRATUM}_convergence_smoke.json"
TREE = f"stratum_{STRATUM}_convergence"

os.environ["R11_TREE"] = TREE
os.environ["R11_CORRECTED"] = str(ROWS)
os.environ["R11_OUTPUT"] = str(OUTPUT)
os.environ["R11_SMOKE_OUTPUT"] = str(SMOKE_OUTPUT)

import rev10_sobol_confirmatory as base                     # noqa: E402
import rev11_full_convergence as fc                         # noqa: E402

fc.CORRECTED = ROWS
fc.OUTPUT = OUTPUT
fc.SMOKE_OUTPUT = SMOKE_OUTPUT
fc.CASE_ROOT = METRICS / "r11_cases" / TREE
fc.RAW_ROOT = METRICS / "r11_raw" / TREE
fc.POLICIES = ("truth", "fixed_critical")
fc.COMPARED = ("fixed_critical",)

NOTE = ("two-policy base: the four schedule policies of the convergence study "
        "were not propagated for this population, so the schedule-versus-"
        "comparator comparisons do not exist and this dictionary is empty by "
        "construction rather than by failure. No result is quoted from them. "
        "Declared in the registration before propagation.")


def orbit_summary(row: dict) -> dict:
    """fc.orbit_summary with the comparison block removed, nothing else."""
    index = row["sobol_index"]
    data = {}
    for policy in fc.POLICIES:
        for level in fc.LEVELS:
            data[(policy, level)] = fc.read_meta(index, policy, level, False)
    truth_self = base.common_error(
        data[("truth", "tight")][1], data[("truth", "tight")][2],
        data[("truth", "tighter")][1], data[("truth", "tighter")][2],
    )["pos_rms_m"]
    policies = {}
    for policy in fc.COMPARED:
        errors = {}
        for level in fc.LEVELS:
            errors[level] = base.common_error(
                data[(policy, level)][1], data[(policy, level)][2],
                data[("truth", level)][1], data[("truth", level)][2])
        self_diff = base.common_error(
            data[(policy, "tight")][1], data[(policy, "tight")][2],
            data[(policy, "tighter")][1], data[(policy, "tighter")][2],
        )["pos_rms_m"]
        policies[policy] = {
            "errors_against_same_tolerance_truth": errors,
            "self_difference_rms_m": self_diff,
            "truth_inclusive_envelope_m": self_diff + truth_self,
            "status": data[(policy, "tight")][0]["status"],
        }
    return {
        "sobol_index": index, "name": row["name"],
        "adopted_truth_degree": row["adopted_truth_degree"],
        "design_point": {k: row["design_point"][k] for k in
                         ("hp_km", "ha_km", "incl_deg", "eccentricity")},
        "n_work": row["n_work"], "n_critical": row["n_critical"],
        "truth_self_difference_rms_m": truth_self,
        "truth_status": data[("truth", "tight")][0]["status"],
        "policies": policies,
        "comparisons": {},
        "comparisons_absent_because": NOTE,
        "trajectory_records": [
            {"policy": p, "level": l,
             "status": data[(p, l)][0]["status"],
             "config_sha256": data[(p, l)][0]["config_sha256"],
             "raw_sha256": data[(p, l)][0]["raw_sha256"],
             "telemetry": data[(p, l)][0]["telemetry"]}
            for p in fc.POLICIES for l in fc.LEVELS],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--registry", default="r30")
    a = ap.parse_args()

    if not ROWS.exists():
        print(f"[abort] {ROWS.name} missing")
        return 2
    rows = json.loads(ROWS.read_text(encoding="utf-8"))["rows"]
    wall0 = time.perf_counter_ns()
    started = base.utc_now()

    summaries, failures = [], []
    for row in sorted(rows, key=lambda r: r["sobol_index"]):
        try:
            summaries.append(orbit_summary(row))
        except Exception as exc:                            # noqa: BLE001
            failures.append({"index": row["sobol_index"], "policy": "summary",
                             "level": "-", "status": "summary_incomplete",
                             "message": f"{type(exc).__name__}: {exc}"})
    complete = len(summaries) == len(rows) and not failures
    fc.write_index(summaries, complete, not complete, started, wall0, None,
                   False, failures, 0)

    d = json.loads(OUTPUT.read_text(encoding="utf-8"))
    d["base_scope"] = NOTE
    d["reindexed_by"] = "rev32_reindex_two_policy_base.py"
    base.atomic_json(OUTPUT, d)

    print(f"[{STRATUM}] {len(summaries)}/{len(rows)} orbits indexed, "
          f"{len(failures)} failures, complete={complete}")
    for f in failures[:5]:
        print(f"  !! {f['index']:03d} {f['message']}")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
