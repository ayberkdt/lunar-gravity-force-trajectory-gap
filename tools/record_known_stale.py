"""Rebuild known_stale_digests.json from the current manifests and scripts.

    python tools/record_known_stale.py

Use this only after deliberately resolving or accepting a driver difference.
It records every mismatch that exists right now, so running it blindly would
silently absorb new drift, which is exactly what the audit exists to catch.
``audit_manifest_digests.py`` fails on anything not already recorded here.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"
SCRIPTS = ROOT / "scripts"
OUT = ROOT / "known_stale_digests.json"

WHAT_THIS_IS = (
    "Driver scripts whose current bytes differ from the digest their campaign "
    "manifest recorded, i.e. they were edited after that manifest was frozen. "
    "The differences are substantive (hundreds to thousands of bytes), so the "
    "digests are NOT refreshed: doing so would assert that the current script "
    "produced the archived results, which has not been demonstrated. The "
    "resolution is to re-run the affected campaign with the current driver, or "
    "to restore the driver to the version that ran, and then remove its entry "
    "here."
)

EFFECT_ON_CI = (
    "audit_manifest_digests.py fails on any mismatch that is NOT listed here, "
    "so new drift breaks the build. Listed entries are reported as known and "
    "do not fail. Run with --strict to fail on these as well."
)


def main() -> int:
    entries = []
    for manifest_path in sorted(METRICS.glob("r*_final_experiment_manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for name, meta in (manifest.get("scripts") or {}).items():
            if not isinstance(meta, dict) or not meta.get("sha256"):
                continue
            script = SCRIPTS / Path(name).name
            if not script.exists():
                continue
            current = script.read_bytes()
            observed = hashlib.sha256(current).hexdigest()
            if observed == meta["sha256"]:
                continue
            entries.append({
                "manifest": manifest_path.name,
                "script": Path(name).name,
                "recorded_sha256": meta["sha256"],
                "recorded_bytes": meta.get("bytes"),
                "observed_sha256": observed,
                "observed_bytes": len(current),
            })

    previous = []
    if OUT.exists():
        previous = json.loads(OUT.read_text(encoding="utf-8")).get("entries", [])
    known = {(e["manifest"], e["script"]) for e in previous}
    added = [e for e in entries if (e["manifest"], e["script"]) not in known]

    payload = {
        "schema": "known_stale_digests_v1",
        "recorded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "what_this_is": WHAT_THIS_IS,
        "effect_on_ci": EFFECT_ON_CI,
        "entries": entries,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {OUT.name}: {len(entries)} entries "
          f"({len(previous)} before, {len(added)} newly absorbed)")
    for entry in added:
        print(f"  [new] {entry['manifest']}: {entry['script']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
