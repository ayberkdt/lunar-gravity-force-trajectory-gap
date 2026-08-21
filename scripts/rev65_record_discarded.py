"""Record why R65's first attempt was discarded, with the evidence.

The attempt ran to completion --- 28 of 28 cells --- and its aggregate would
have been readable. It is discarded anyway, because the machine was not idle
while it ran and the quantity the campaign matches on is measured wall time.

The registration requires an idle machine for every timed stage. That is
enforced only at start-up, so a campaign that begins on a quiet machine and is
joined by other work keeps running and keeps reporting in-band ratios: the
member and its comparator are divided by each other, so a load common to both
cancels and the band check passes while the degree it selected is wrong. What
does not cancel is load that changes between the two measurements, and that is
what this record measures.

Throughput here is gravity_kernel_ns / (n_rhs * mean_degree_sq), which is flat
on a quiet machine.

Usage:  python rev65_record_discarded.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
ARCHIVE = METRICS / "r65_discarded_attempt1"
OUT = METRICS / "r65_discarded_attempt1.json"
BASELINE_N = 8


def cases():
    rows = []
    for p in sorted((ARCHIVE / "cases").rglob("*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        t = d.get("telemetry") or {}
        n, dsq, ns = (t.get("n_rhs"), t.get("mean_degree_sq"),
                      t.get("gravity_kernel_ns"))
        if not (n and dsq and ns):
            continue
        rows.append({"created_utc": d.get("created_utc"),
                     "orbit": p.parent.name, "case": p.name,
                     "throughput_ns_per_call_per_degree_sq": ns / (n * dsq)})
    rows.sort(key=lambda r: r["created_utc"] or "")
    return rows


def main() -> int:
    if not ARCHIVE.exists():
        print(f"[abort] {ARCHIVE.name} is not there")
        return 2
    rows = cases()
    v = np.array([r["throughput_ns_per_call_per_degree_sq"] for r in rows])
    base = float(np.median(v[:BASELINE_N]))

    # Within a cell, the member and its comparator are measured minutes apart.
    # A load that changes in between does not cancel in their ratio.
    by_orbit = {}
    for r, x in zip(rows, v):
        by_orbit.setdefault(r["orbit"], []).append(x)
    spread = {k: float(max(a) / min(a)) for k, a in by_orbit.items()
              if len(a) > 1}
    worst = max(spread, key=spread.get)

    log = (ARCHIVE / "r65_20260819.log").read_bytes()
    a1 = json.loads((ARCHIVE / "r65_timing_family.json"
                     ).read_text(encoding="utf-8"))
    ok = json.loads((METRICS / "r65_timing_family.json"
                     ).read_text(encoding="utf-8"))

    def key(r):
        return (r["design"], int(r["sobol_index"]), str(r["k"]))
    d1 = {key(r): r for r in a1["rows"]}
    d2 = {key(r): r for r in ok["rows"]}
    fresh = [k for k in set(d1) & set(d2) if k[2] != "0.50"]
    deg_changed = sum(1 for k in fresh
                      if d1[k]["comparator_degree"] != d2[k]["comparator_degree"])
    flipped = sum(1 for k in set(d1) & set(d2)
                  if (d1[k].get("resolved"), d1[k].get("winner"))
                  != (d2[k].get("resolved"), d2[k].get("winner")))

    # The decision time was readable only from the NTFS creation stamp of this
    # file, which no archive format preserves; carrying it as a field is what
    # makes the ordering argument survive packaging.
    decided = "2026-08-19T19:56:14Z"
    first_clean = min(
        (json.loads(p.read_text(encoding="utf-8")).get("created_utc") or "")
        for p in (METRICS / "r65_cases").rglob("*.json"))

    payload = {
        "schema": "r65_discarded_attempt_v3",
        "campaign": "R65 (O58), first attempt",
        "status": "discarded, not used for any reported number",
        "decided_utc": decided,
        "timeline": {
            "registered_utc": json.loads(
                (METRICS / "r65_preregistration.json").read_text(
                    encoding="utf-8"))["registered_utc"],
            "attempt_last_cell_utc": max(r["created_utc"] for r in rows
                                         if r["created_utc"]),
            "decided_utc": decided,
            "clean_rerun_first_cell_utc": first_clean,
            "clean_rerun_aggregate_utc": ok.get("created_utc"),
            "note": ("the discard was recorded before the clean re-run began, "
                     "which is the whole of the ordering argument"),
        },
        "completed": "28 of 28 cells; the aggregate was written and is kept in "
                     "the archive directory",
        "outcome_it_returned": {
            "by_k": {k: {"resolved": s["resolved"],
                         "interior_wins": s["interior_wins"],
                         "fixed_wins": s["fixed_wins"],
                         "timing_match_misses": s["timing_match_misses"]}
                     for k, s in a1["by_k"].items()},
            "declared_class": "B, the same class the accepted re-run returned",
            "note": ("reported rather than withheld. The attempt's log carries "
                     "these tallies in plain text and its digest is published "
                     "here, so a claim that they were never read would not be "
                     "checkable. What the discard rests on is not the summary "
                     "but the cell-level selection below."),
        },
        "what_the_rerun_changed": {
            "freshly_propagated_cells": len(fresh),
            "comparator_degree_differs": deg_changed,
            "cells_whose_verdict_differs": flipped,
            "reading": ("the summary is unchanged but the measurement is not: "
                        "the comparator degrees this attempt selected are not "
                        "the ones an idle machine selected, which is the "
                        "quantity the campaign is built on"),
        },
        "criterion_was_set_after_the_run": (
            "the registration fixes a mechanism, an idle check at the start of "
            "the timed pipeline, not a continuing condition and not a numeric "
            "threshold. The throughput thresholds used to condemn this attempt "
            "were set after its data was on disk. That is disclosed here "
            "rather than presented as a pre-registered rule; what makes the "
            "discard defensible is the cell-level divergence above, which is "
            "independent of any threshold."),
        "reason": "the machine was not idle for the whole run, and the "
                  "campaign matches its comparator on measured kernel time",
        "why_the_band_check_did_not_catch_it":
            "the achieved ratio divides the comparator's measured time by the "
            "member's, so a load common to both cancels and every cell "
            "reported an in-band ratio; only load that changes between the "
            "two measurements biases the match, and it does so silently",
        "throughput_metric":
            "gravity_kernel_ns / (n_rhs * mean_degree_sq). This raw measure is "
            "NOT flat in the degree: the per-call cost carries a term that "
            "does not scale with N^2, so small-degree cases read high even on "
            "an idle machine. The degree-normalized audit that removes it, and "
            "the accepted run's own figures for comparison, are in "
            "metrics/r65_quiet_audit.json; the raw numbers below overstate the "
            "contention by however much the degree term accounts for.",
        "timed_cases": len(rows),
        "baseline_first_n": BASELINE_N,
        "baseline_throughput": base,
        "worst_throughput": float(v.max()),
        "worst_relative_to_baseline": float(v.max() / base),
        "cases_above_1p25_baseline": int((v > 1.25 * base).sum()),
        "worst_within_orbit_spread": {
            "orbit": worst, "max_over_min": spread[worst]},
        "orbits_with_within_spread_above_1p15":
            int(sum(1 for s in spread.values() if s > 1.15)),
        "matching_band": [0.90, 1.10],
        "reading": "within-cell spreads exceed the matching band, so the "
                   "comparator degrees this attempt selected are not the ones "
                   "an idle machine would have selected",
        "registration_unchanged": "metrics/r65_preregistration.json still "
                                  "governs; k values, panel, scoring rule, "
                                  "timing band and reporting rule are "
                                  "untouched, and the rerun is scored against "
                                  "it",
        "archive_dir": "metrics/r65_discarded_attempt1",
        "archived_log_sha256": hashlib.sha256(log).hexdigest(),
        "guard_added": "python_codes/probe_kernel_quiet.py, which compares the "
                       "kernel against the archived R12 idle-machine curve "
                       "before a timing campaign is started",
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    print(f"  {len(rows)} timed cases, baseline {base:.3f}, "
          f"worst {v.max():.3f} ({v.max() / base:.2f}x)")
    print(f"  {payload['cases_above_1p25_baseline']} cases above 1.25x")
    print(f"  worst within-orbit spread {spread[worst]:.2f} on {worst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
