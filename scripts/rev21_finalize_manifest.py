"""Manifest for R21, the gradient-degree sensitivity check of the variational solve.

R21 propagates nothing. It re-reads archived reference trajectories and
evaluates the gravity gradient twice at each sampled epoch, so its manifest
records the driver script, the archived inputs it reads, and the two outputs it
writes. The convention is that of the R16 manifest, which is likewise a
field-level check with no trajectory tree of its own.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r21_final_experiment_manifest.json"

SCRIPTS = ["rev21_gradient_sensitivity.py", "rev21_finalize_manifest.py"]
RESULTS = ["r21_gradient_sensitivity.json", "r21_gradient_sensitivity_table.tex"]
REUSED = ["r14_variational_budget.json", "r14_budget_pareto.json"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def index(names, base: Path) -> dict:
    """Record every expected file, present or not. An earlier version dropped
    absent names from the dict, so a file that was never produced left no trace
    at all rather than a record the integrity check could fail on."""
    out = {}
    for n in names:
        p = base / n
        out[n] = ({"sha256": sha256(p), "bytes": p.stat().st_size}
                  if p.exists() else {"missing": True})
    return out


def main() -> int:
    payload = {
        "schema": "r21_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": "R21: relative truncation of the reference gravity gradient at "
                 "the variational solve's degree-120 evaluation, and the forcing "
                 "it misrepresents, on the eight-orbit panel of the fixed-budget "
                 "variational experiment. No trajectory is propagated: the "
                 "archived reference trajectories of the earlier campaigns are "
                 "re-read and the gradient is evaluated twice per sampled epoch.",
        "relationship_to_r14": "reads the R14 variational panel for the orbit "
                               "list, adopted truth degrees and predicted "
                               "displacements, and the R14 budget-pareto record "
                               "for the per-policy defect RMS; reuses the "
                               "archived truth trajectories indexed under the "
                               "earlier manifests and duplicates none of them.",
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23"},
        "reused_inputs": index(REUSED, METRICS),
        "scripts": index(SCRIPTS, CODE),
        "result_json": index(RESULTS, METRICS),
        "trajectory_tree": {"n_sidecars": 0, "n_raw_arrays": 0,
                            "note": "field-level check; propagates nothing"},
    }
    # Same canonical form as the other twelve sealed manifests. This one used
    # indent=1 with sort_keys, so its digest could not be reproduced by the
    # package's own convention and read as a broken seal.
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r21-manifest] wrote {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  scripts {len(payload['scripts'])}, results {len(payload['result_json'])}, "
          f"reused {len(payload['reused_inputs'])}")
    missing = [k for sec in ("scripts", "result_json", "reused_inputs")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
