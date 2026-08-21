"""Was the machine quiet while R68 measured? Audited after the run.

The construction, the metric and its degree correction are R65's; this file
only points them at the R68 case trees and keeps the R65 panels alongside as
the reference spread an idle machine produced. Nothing is propagated.

A population campaign runs for hours rather than for one evening, so the
per-arm split matters here: a load that arrived halfway through would show as
a wide normalized spread in one arm and not the other.

Usage:  python rev68_quiet_audit.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev65_quiet_audit as q65

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "r68_quiet_audit.json"

RUNS = {
    "r68_endpoint": "r68_cases/endpoint",
    "r68_interior": "r68_cases/interior",
    "r65_reference": "r65_cases",
    "r48_reference": "r48_cases",
}


def main() -> int:
    payload = {
        "schema": "r68_quiet_audit_v1",
        "purpose": ("post-hoc contention audit of the R68 population timing "
                    "campaign, from telemetry already in the case records; "
                    "no propagation is repeated"),
        "metric": ("gravity_kernel_ns / (n_rhs * mean_degree_sq), and the "
                   "same quantity divided by a fitted a + b/N^2 to remove the "
                   "degree dependence the raw measure carries"),
        "reference_runs": ("r65 and r48 are the accepted panel campaigns; "
                           "their spreads are what this machine produced when "
                           "it was quiet"),
        "runs": {},
    }
    for name, sub in RUNS.items():
        rows = q65.cases(sub)
        if not rows:
            print(f"[skip] {name}: no timed cases under {sub}")
            continue
        payload["runs"][name] = q65.audit(rows)

    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    hdr = (f"{'run':<20}{'n':>5}{'raw max':>9}{'norm max':>10}"
           f"{'>1.25x':>8}{'corr':>8}")
    print(hdr)
    print("-" * len(hdr))
    for name, a in payload["runs"].items():
        print(f"{name:<20}{a['timed_cases']:>5}"
              f"{a['raw']['worst_within_orbit_spread']:>9.2f}"
              f"{a['normalized']['worst_within_orbit_spread']:>10.2f}"
              f"{a['normalized']['cases_above_1p25_baseline']:>8d}"
              f"{a['degree_model']['corr_raw_throughput_vs_degree']:>8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
