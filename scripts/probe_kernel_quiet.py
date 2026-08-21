"""Is the machine quiet enough to run a timing campaign?

The timing campaigns match a comparator degree on measured kernel time inside a
0.90--1.10 band, so a few per cent of CPU contention is enough to move a cell
across the band edge. R65's first attempt was lost that way: its throughput
drifted from the early baseline to 1.86 times it while other work ran.

This probe measures the same kernel at three degrees and compares each with the
archived R12 cost curve, which was taken on an idle machine. A ratio near unity
means the machine is quiet; anything above about 1.10 means it is not, and the
campaign should wait.

Usage:  python probe_kernel_quiet.py
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from rev3_common import load_model, kernel_args, warmup
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba

ROOT = Path(__file__).resolve().parents[1]
CURVE = ROOT / "metrics" / "r12_kernel_cost_curve.json"
LOG = ROOT / "metrics" / "r65_quiet_probe_log.json"
DEGREES = (60, 120, 240)
REPEATS = 4000
R_M = 1.8e6
TOLERANCE = 1.10


def per_call_ns(args, degree: int) -> float:
    x, y, z = R_M, 0.0, 0.0
    for _ in range(200):
        sh_accel_fixed_numba(x, y, z, degree, *args)
    s = np.empty(REPEATS)
    for i in range(REPEATS):
        t0 = time.perf_counter_ns()
        sh_accel_fixed_numba(x, y, z, degree, *args)
        s[i] = time.perf_counter_ns() - t0
    return float(np.median(s))


def main() -> int:
    ref = {r["degree"]: r["per_call_ns_median"]
           for r in json.loads(CURVE.read_text(encoding="utf-8"))["rows"]}
    model = load_model()
    args = kernel_args(model)
    warmup(model, args)

    print(f"{'N':>5}{'now [ns]':>12}{'archive [ns]':>14}{'ratio':>8}")
    ratios, measured = [], {}
    for n in DEGREES:
        if n not in ref:
            continue
        now = per_call_ns(args, n)
        ratio = now / ref[n]
        ratios.append(ratio)
        measured[str(n)] = {"now_ns": now, "archive_ns": ref[n],
                            "ratio": ratio}
        print(f"{n:>5}{now:>12.0f}{ref[n]:>14.0f}{ratio:>8.3f}")

    worst = max(ratios)
    # A precondition that leaves no artifact cannot be audited afterwards, and
    # the campaign this guard was written for was lost to exactly that kind of
    # unrecorded assumption. Every run appends.
    rec = {"utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "degrees": list(DEGREES), "repeats": REPEATS,
           "reference_curve": "metrics/r12_kernel_cost_curve.json",
           "per_degree": measured, "worst_ratio": worst,
           "tolerance": TOLERANCE, "verdict": "quiet" if worst <= TOLERANCE
           else "busy"}
    log = json.loads(LOG.read_text(encoding="utf-8")) if LOG.exists() else []
    log.append(rec)
    LOG.write_text(json.dumps(log, indent=2), encoding="utf-8")

    print(f"\nworst ratio {worst:.3f} against the archived idle-machine curve")
    print(f"[written] {LOG.name} ({len(log)} probes recorded)")
    if worst > TOLERANCE:
        print("[busy] the machine is not quiet enough for a timing campaign")
        return 1
    print("[quiet] safe to start a timing campaign")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
