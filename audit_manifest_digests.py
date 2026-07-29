"""Check every driver script against the digest its campaign manifest records.

    python audit_manifest_digests.py            # fail on undocumented drift
    python audit_manifest_digests.py --strict   # fail on any mismatch at all

Each campaign manifest under ``metrics/`` records a SHA-256 and a byte count
for the driver scripts that produced it. This walks all fourteen manifests and
compares them against the scripts as they stand.

A mismatch means a driver was edited after its manifest was finalized. Those
already known are listed in ``known_stale_digests.json`` with the recorded and
observed values; they are reported but do not fail the default run, so that any
*new* drift is what breaks the build. ``--strict`` ignores that list and fails
on every mismatch, which is the state the archive should reach before release:
each entry is resolved by re-running the affected campaign with the current
driver, or restoring the driver to the version that ran, and then deleting its
entry.

Exit status is 0 when nothing unexpected was found.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
METRICS = ROOT / "metrics"
SCRIPTS = ROOT / "scripts"
KNOWN_STALE = ROOT / "known_stale_digests.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_known() -> dict[tuple[str, str], dict]:
    if not KNOWN_STALE.exists():
        return {}
    payload = json.loads(KNOWN_STALE.read_text(encoding="utf-8"))
    return {(e["manifest"], e["script"]): e for e in payload.get("entries", [])}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="fail on documented mismatches too")
    args = parser.parse_args()

    known = {} if args.strict else load_known()

    matched = 0
    missing: list[tuple[str, str]] = []
    documented: list[tuple[str, str]] = []
    unexpected: list[tuple[str, str, str, str]] = []

    manifests = sorted(METRICS.glob("r*_final_experiment_manifest.json"))
    for manifest_path in manifests:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for name, meta in (manifest.get("scripts") or {}).items():
            recorded = meta.get("sha256") if isinstance(meta, dict) else None
            if not recorded:
                continue
            script_name = Path(name).name
            script = SCRIPTS / script_name
            if not script.exists():
                missing.append((manifest_path.name, script_name))
                continue
            actual = sha256(script)
            if actual == recorded:
                matched += 1
            elif (manifest_path.name, script_name) in known:
                documented.append((manifest_path.name, script_name))
            else:
                unexpected.append(
                    (manifest_path.name, script_name, recorded, actual))

    print(f"manifests scanned      : {len(manifests)}")
    print(f"digests matching       : {matched}")
    print(f"scripts missing        : {len(missing)}")
    print(f"documented differences : {len(documented)}")
    print(f"undocumented           : {len(unexpected)}")

    for manifest_name, script_name in missing:
        print(f"  [missing]     {manifest_name}: {script_name}")
    for manifest_name, script_name in documented:
        print(f"  [known-stale] {manifest_name}: {script_name}")
    for manifest_name, script_name, recorded, actual in unexpected:
        print(f"  [DRIFT]       {manifest_name}: {script_name}")
        print(f"                recorded {recorded}")
        print(f"                actual   {actual}")

    if documented and not args.strict:
        print(f"\n{len(documented)} known difference(s); see "
              f"known_stale_digests.json. Run with --strict to treat them as "
              f"failures.")

    return 1 if (missing or unexpected) else 0


if __name__ == "__main__":
    raise SystemExit(main())
