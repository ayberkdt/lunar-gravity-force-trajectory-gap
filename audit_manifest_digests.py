"""Check every driver script against the digest its campaign manifest records.

    python audit_manifest_digests.py

Each campaign manifest under ``metrics/`` records a SHA-256 for the driver
scripts that produced it. This walks all fourteen manifests and reports any
script whose current bytes differ, which is the signal that a driver was edited
after its manifest was finalized and that the manifest needs refreshing with
the campaign's own ``revNN_finalize_manifest.py``.

Exit status is 0 when every recorded digest matches.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
METRICS = ROOT / "metrics"
SCRIPTS = ROOT / "scripts"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    matched = 0
    missing: list[tuple[str, str]] = []
    mismatched: list[tuple[str, str, str, str]] = []

    manifests = sorted(METRICS.glob("r*_final_experiment_manifest.json"))
    for manifest_path in manifests:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for name, meta in (manifest.get("scripts") or {}).items():
            recorded = meta.get("sha256") if isinstance(meta, dict) else None
            if not recorded:
                continue
            script = SCRIPTS / Path(name).name
            if not script.exists():
                missing.append((manifest_path.name, name))
                continue
            actual = sha256(script)
            if actual == recorded:
                matched += 1
            else:
                mismatched.append(
                    (manifest_path.name, name, recorded, actual))

    print(f"manifests scanned : {len(manifests)}")
    print(f"digests matching  : {matched}")
    print(f"scripts missing   : {len(missing)}")
    print(f"digests differing : {len(mismatched)}")

    for manifest_name, script_name in missing:
        print(f"  [missing]  {manifest_name}: {script_name}")
    for manifest_name, script_name, recorded, actual in mismatched:
        print(f"  [differs]  {manifest_name}: {script_name}")
        print(f"             recorded {recorded}")
        print(f"             actual   {actual}")

    return 0 if not missing and not mismatched else 1


if __name__ == "__main__":
    raise SystemExit(main())
