"""Run rev14_budget_trajectory at beta = 0.62 against the R28 calibration record.

One job, and it exists because of one line: build_specs reads the Phase-A
calibration out of a module constant, and beta = 0.62 lives in a separate record
(r28_budget_pareto_beta_0.62.json) because the archived one is sha256-pinned in
three sealed manifests and must not be edited. Redirecting the constant is the
whole of the difference; every other input, tolerance and output path is the
driver's own.

Why redirecting a constant is safe here rather than merely convenient. The
Windows spawn trap that governs this codebase is that workers re-import their
module and lose any parent-only patch, so anything a worker reads at import time
must be set through the environment instead. This patch is not subject to it:
build_specs runs once in the parent, its output is placed in each task payload,
and no worker ever opens either calibration record. The patch is applied inside
main() so that a spawned child importing this module inherits nothing.

Usage:
    python rev28_r14_beta062.py --design A --workers 11 --deadline 2026-07-31T23:00:00+03:00
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import rev14_budget_trajectory as r14

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CALIBRATION = METRICS / "r28_budget_pareto_beta_0.62.json"
AMENDMENT = METRICS / "r28_calibration_amendment.json"
BETA = 0.62


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B"), required=True)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    if not CALIBRATION.exists():
        print(f"[abort] {CALIBRATION.name} missing; run "
              f"rev28_budget_pareto_extension.py first")
        return 2
    cal = json.loads(CALIBRATION.read_text(encoding="utf-8"))
    if cal.get("admissibility_check", {}).get("result") != \
            "exact_on_all_orbits_both_designs":
        print(f"[abort] {CALIBRATION.name} did not pass its reproduction check")
        return 2
    if not AMENDMENT.exists():
        print(f"[abort] {AMENDMENT.name} missing")
        return 2

    # the parent-only redirect; build_specs is the only reader and it runs here
    r14.PARETO = CALIBRATION
    specs = r14.build_specs(a.design, BETA)
    print(f"[r28-r14] design {a.design} beta {BETA} from {CALIBRATION.name}")
    print(f"  calibration {cal['amendment_sha256'][:16]}, "
          f"{len(specs)} orbits, "
          f"{sum(s['censored'] for s in specs.values())} censored")
    print(f"  status: {cal['status']}")

    return r14.run(a.design, BETA, a.workers,
                   r14.parse_deadline(a.deadline), a.limit)


if __name__ == "__main__":
    raise SystemExit(main())
