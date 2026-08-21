"""SHA-256 integrity manifest for R57: the gradient-degree audit, completed.

R39 registered a sixteen-orbit audit and solved thirteen before its compute
ran out; the three it left are the most expensive in the subset, all at
reference degree 900. R57 solves exactly those three and carries R39's
thirteen forward byte for byte, so the completion record is a superset of the
panel record and nothing is re-solved or re-scored.

Ownership follows that split. The completion record and the solver belong
here. The audit table stays owned by the R39 manifest, which has always
indexed it, even though it is now drawn from this record: rev39_table.py
prefers the completion file when it exists, and R39's manifest is re-sealed
for the new digest as a generated-table refresh.

R39's own manifest is left as sealed. Its `completion` block reports thirteen
of sixteen because that is what R39 finished; the state of the audit as a
whole is this manifest's `supersedes` block.

Usage:  python rev57_finalize_manifest.py
"""

from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SCRIPTS = ROOT / "python_codes"

RESULT = "r57_gradient_degree_completion.json"
SCRIPT = "rev57_gradient_complete.py"
REUSED = ["r39_gradient_degree_panel.json", "r39_preregistration.json"]


def entry(path: Path) -> dict:
    return {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size}


def main() -> int:
    rec = json.loads((METRICS / RESULT).read_text(encoding="utf-8"))
    s = rec["summary"]
    payload = {
        "schema": "r57_final_experiment_manifest_v1",
        "created_utc": datetime.datetime.now(datetime.timezone.utc)
                        .isoformat().replace("+00:00", "Z"),
        "scope": ("R57: the three unfinished orbits of the (O40) gradient-degree "
                  "audit, solved at reference degree 900, with R39's thirteen "
                  "carried forward unchanged."),
        "propagation": ("three augmented variational solves, one per orbit. No "
                        "propagated trajectory is produced and no comparison is "
                        "re-scored; the carried rows are copied byte for byte."),
        "supersedes": {
            "record": "r39_gradient_degree_panel.json",
            "why": ("R39 reports 13 of 16 registered orbits and band L as a "
                    "partial band. With these three solved the audit covers "
                    "all 16 and every declared band is complete; the R39 "
                    "manifest is left sealed at its own final state."),
            "orbits": s["orbits"],
            "resolved": s["resolved"],
            "side_changes_resolved": s["side_changes_resolved"],
            "side_changes_unresolved": s["side_changes_unresolved"],
            "bands_complete": sorted(s["bands_complete"]),
        },
        "reused_inputs": {n: entry(METRICS / n) for n in REUSED},
        "scripts": {SCRIPT: entry(SCRIPTS / SCRIPT)},
        "result_json": {RESULT: entry(METRICS / RESULT)},
        "generated_tables_note": ("r39_gradient_degree_table.tex is drawn from "
                                  "this record but stays indexed by the R39 "
                                  "manifest, which has always owned it."),
    }
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["manifest_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    out = METRICS / "r57_final_experiment_manifest.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] metrics/{out.name}: {s['orbits']} orbits, "
          f"{s['side_changes_resolved']} resolved side changes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
