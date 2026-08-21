"""SHA-256 integrity manifest for R67: two review-round re-readings.

R67 is a derived product, not a campaign. It propagates nothing, re-scores
nothing and adds no comparison the archived records do not already carry. A
review round asked two questions the manuscript had answered qualitatively:

  * does the span-versus-switches contrast at the declared budget survive a
    rank statistic, or is it a property of the linear correlation the R14
    mechanism block happens to record? It does not survive: the rank
    association with span is -0.09 and -0.06 against -0.37 and -0.27 in
    Pearson, and it is smaller than the rank association with the switch
    count. The main text now claims a weak association and points at the
    paired apolune ladder for the population evidence.
  * are the 28 undecided orbits of the forced-variational panel spread across
    the perilune range, as the supplement said? They are not: the undecided
    rate is 14 of 43 in the panel's lowest perilune tertile against 7 of 42
    and 7 of 43 above it.

The generator refuses to write unless its Pearson values reproduce the R14
record exactly, because a disagreement there would mean it joined the rows
differently and the rank values would be measuring something else.

Every input is an existing sealed record and stays owned by its own manifest;
they are named here and re-indexed nowhere, following the convention R59, R60
and R66 use for reused inputs.

Usage:  python rev67_finalize_manifest.py
"""

from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SCRIPTS = ROOT / "python_codes"

REUSED = [
    "r14_budget_pareto.json",
    "r14_descriptives.json",
    "r42_panel_verdict.json",
]
SCRIPTS_OWNED = ["rev67_review_measurements.py"]
RECORD = "r67_review_measurements.json"


def entry(path: Path) -> dict:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"sha256": h, "bytes": path.stat().st_size}


def main() -> int:
    payload = {
        "schema": "r67_final_experiment_manifest_v1",
        "created_utc": datetime.datetime.now(datetime.timezone.utc)
                        .isoformat().replace("+00:00", "Z"),
        "scope": ("R67: the rank form of the R14 span/switch association and "
                  "the perilune distribution of the R42 panel's undecided "
                  "comparisons, both read from sealed records."),
        "why": ("Two sentences of the manuscript rested on a linear "
                "correlation and on a qualitative description. Both are now "
                "measured, and both measurements changed what the sentences "
                "say."),
        "propagation": ("none. No orbit is propagated and no comparison is "
                        "re-scored. The generator refuses to write unless its "
                        "Pearson values reproduce the R14 mechanism block "
                        "exactly."),
        "trajectory_records": {
            "owner_manifests": ["r14"],
            "note": ("r14_trajectory_{A,B}_beta_*.json supply rho_budget and "
                     "are read only. Trajectory records are a partition, so "
                     "they are named here and re-indexed nowhere."),
        },
        "reused_inputs": {n: entry(METRICS / n) for n in REUSED},
        "scripts": {s: entry(SCRIPTS / s) for s in SCRIPTS_OWNED},
        "records": {f"metrics/{RECORD}": entry(METRICS / RECORD)},
    }
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["manifest_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    out = METRICS / "r67_final_experiment_manifest.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] metrics/{out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
