"""SHA-256 integrity manifest for the R22 supplementary tables.

R22 was the one campaign whose manifest had no generator: the JSON was written
by hand, so it carried neither a creation timestamp nor the numerical-kernel
identity that the other thirteen record, and its two generated tables were
filed under ``result_json`` rather than ``generated_tables``. This script
reproduces it on the twelve-field schema the rest of the package uses.

R22 propagates nothing. It reads the frozen R18 span-sweep summaries at every
declared budget and the R20 sixty-day sidecar telemetry in place, and emits two
supplementary tables plus one cost summary. The trajectories it reads stay
indexed under R18 and R20, so nothing is duplicated here.

Usage:  python rev22_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r22_final_experiment_manifest.json"

SCRIPTS = ["rev22_supplement_tables.py", "rev22_finalize_manifest.py"]
RESULT_JSON = ["r22_longarc_cost_summary.json"]
TABLES = ["r18_by_k_table.tex", "r20_longarc_detail_table.tex"]

# Every archived record the tables are built from. The five span-sweep budgets
# are the ones rev22_supplement_tables.py iterates over; listing fewer would
# understate what the by-k table depends on.
REUSED = ["r18_span_sweep_A_beta_0.50.json",
          "r18_span_sweep_B_beta_0.50.json",
          "r18_span_sweep_A_beta_1.00.json",
          "r18_span_sweep_B_beta_1.00.json",
          "r18_span_sweep_A_beta_1.50.json",
          "r20_span_longarc.json"]


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path) -> dict:
    out = {}
    for n in names:
        p = base / n
        out[n] = ({"sha256": sha(p), "bytes": p.stat().st_size}
                  if p.exists() else {"missing": True})
    return out


def main() -> int:
    payload = {
        "schema": "r22_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R22: two supplementary tables derived from archived R18 and "
                  "R20 records. Propagates nothing; reads the frozen "
                  "span-sweep summaries and the sixty-day sidecar telemetry "
                  "and emits the by-k budget aggregate and the sixty-day "
                  "per-orbit error/envelope/cost table."),
        "relationship_to_r18_r20": (
            "derived only. The by-k aggregate is read from the five R18 "
            "span-sweep summaries and the sixty-day table from the R20 record "
            "and its k_* sidecars, all in place. No trajectory is propagated "
            "or re-indexed here, so R18 and R20 remain the sole index of the "
            "arcs behind these tables."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "no kernel call: derived tables only"},
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": {
            "n_sidecars": 0, "n_raw_arrays": 0,
            "note": "derived tables only; reads the R20 sidecars in place"},
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "reused_inputs")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
