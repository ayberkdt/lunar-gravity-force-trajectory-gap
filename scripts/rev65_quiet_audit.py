"""Was the machine quiet while a timing campaign measured? Audited after the run.

The band check cannot answer this. It divides the comparator's measured time by
the member's, so a load common to both cancels and every cell reports an in-band
ratio while the degree it selected is wrong. What does not cancel is load that
changes between the two measurements, and that is what this audit measures, from
telemetry the case records already carry.

Raw throughput is gravity_kernel_ns / (n_rhs * mean_degree_sq). It is not flat
in the degree: the per-call cost carries a term that does not scale with N^2
(setup, the low-order recursions, cache behaviour), so small-degree cases read
high on this measure even on an idle machine. Fitting

    thr(N) = a + b / N^2

over the run's own cases absorbs that, and the residual ratio thr/thr_fit is what
a load shows up in. Both the raw and the normalized spread are reported, because
the discarded R65 attempt was condemned on the raw one and the correction changes
how much of its spread was load.

Usage:  python rev65_quiet_audit.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "r65_quiet_audit.json"

RUNS = {
    "r65_accepted": "r65_cases",
    "r65_discarded_attempt1": "r65_discarded_attempt1/cases",
    "r64_reference": "r64_cases",
    "r48_reference": "r48_cases",
}


def cases(sub: str):
    root = METRICS / sub
    out = []
    if not root.exists():
        return out
    for p in sorted(root.rglob("*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        t = d.get("telemetry") or {}
        n, dsq, ns = (t.get("n_rhs"), t.get("mean_degree_sq"),
                      t.get("gravity_kernel_ns"))
        if n and dsq and ns:
            out.append({"utc": d.get("created_utc") or "", "orbit": p.parent.name,
                        "thr": ns / (n * dsq), "deg_sq": float(dsq)})
    out.sort(key=lambda r: r["utc"])
    return out


def audit(rows) -> dict:
    thr = np.array([r["thr"] for r in rows])
    inv = 1.0 / np.array([r["deg_sq"] for r in rows])
    # least squares on thr = a + b * (1/N^2)
    A = np.vstack([np.ones_like(inv), inv]).T
    (a, b), *_ = np.linalg.lstsq(A, thr, rcond=None)
    fit = a + b * inv
    norm = thr / fit
    corr = float(np.corrcoef(thr, np.sqrt(1.0 / inv))[0, 1])

    def spread(v):
        by = {}
        for r, x in zip(rows, v):
            by.setdefault(r["orbit"], []).append(x)
        s = [max(z) / min(z) for z in by.values() if len(z) > 1]
        return (float(max(s)), float(np.median(s))) if s else (1.0, 1.0)

    raw_max, raw_med = spread(thr)
    nrm_max, nrm_med = spread(norm)
    base = float(np.median(norm[:8]))
    return {
        "timed_cases": len(rows),
        "degree_model": {"a": float(a), "b": float(b),
                         "corr_raw_throughput_vs_degree": corr,
                         "residual_rms_over_median": float(
                             np.std(thr - fit) / np.median(thr))},
        "raw": {"worst_within_orbit_spread": raw_max,
                "median_within_orbit_spread": raw_med,
                "worst_over_opening_baseline": float(
                    thr.max() / np.median(thr[:8]))},
        "normalized": {"worst_within_orbit_spread": nrm_max,
                       "median_within_orbit_spread": nrm_med,
                       "cases_above_1p25_baseline": int(
                           (norm > 1.25 * base).sum()),
                       "worst_over_opening_baseline": float(norm.max() / base)},
    }


def main() -> int:
    payload = {
        "schema": "r65_quiet_audit_v1",
        "purpose": ("post-hoc contention audit of the timing campaigns, from "
                    "telemetry already in the sealed case records; no "
                    "propagation is repeated"),
        "metric": ("gravity_kernel_ns / (n_rhs * mean_degree_sq), and the same "
                   "quantity divided by a fitted a + b/N^2 to remove the "
                   "degree dependence the raw measure carries"),
        "why_the_band_check_is_blind": (
            "the achieved ratio divides the comparator's measured time by the "
            "member's, so a load common to both cancels; only load that "
            "changes between the two measurements biases the match"),
        "runs": {},
    }
    for name, sub in RUNS.items():
        rows = cases(sub)
        if not rows:
            print(f"[skip] {name}: no timed cases under {sub}")
            continue
        payload["runs"][name] = audit(rows)

    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}")
    hdr = f"{'run':<26}{'n':>5}{'raw max':>9}{'norm max':>10}{'>1.25x':>8}{'corr':>8}"
    print(hdr)
    print("-" * len(hdr))
    for name, a in payload["runs"].items():
        print(f"{name:<26}{a['timed_cases']:>5}"
              f"{a['raw']['worst_within_orbit_spread']:>9.2f}"
              f"{a['normalized']['worst_within_orbit_spread']:>10.2f}"
              f"{a['normalized']['cases_above_1p25_baseline']:>8d}"
              f"{a['degree_model']['corr_raw_throughput_vs_degree']:>8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
