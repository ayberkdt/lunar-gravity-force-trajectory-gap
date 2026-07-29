"""Check the archive against the digests its campaign manifests record.

    python audit_manifest_digests.py            # fail on undocumented drift
    python audit_manifest_digests.py --strict   # fail on any mismatch at all

Each campaign manifest under ``metrics/`` records a SHA-256 for the driver
scripts that produced it and for the tables, figures and result records they
wrote. This walks all fourteen manifests and compares every one of those
against the files as they stand.

A mismatch means the file changed after its manifest was finalized. Those
already known are listed in ``known_stale_digests.json`` with their recorded
and observed values, and with any evidence gathered about whether the current
driver still yields the archived numbers; they are reported but do not fail
the default run, so that any *new* drift is what breaks the build. ``--strict``
ignores that list and fails on every mismatch, which is the state the archive
should reach before an archival release.

Exit status is 0 when nothing unexpected was found.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
METRICS = ROOT / "metrics"
KNOWN_STALE = ROOT / "known_stale_digests.json"

# manifest section -> directories the artifact may live in
SECTIONS = {
    "scripts": ("scripts",),
    "generated_tables": ("metrics",),
    "generated_figures": ("figures", "metrics"),
    "result_json": ("metrics",),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_known() -> set[tuple[str, str, str]]:
    if not KNOWN_STALE.exists():
        return set()
    payload = json.loads(KNOWN_STALE.read_text(encoding="utf-8"))
    return {(e["manifest"], e.get("category", "scripts"),
             e.get("file", e.get("script", "")))
            for e in payload.get("entries", [])}


def locate(name: str, roots: tuple[str, ...]) -> Path | None:
    for directory in roots:
        candidate = ROOT / directory / Path(name).name
        if candidate.exists():
            return candidate
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="fail on documented mismatches too")
    args = parser.parse_args()

    known = set() if args.strict else load_known()

    matched: dict[str, int] = {}
    missing: list[tuple[str, str, str]] = []
    documented: list[tuple[str, str, str]] = []
    unexpected: list[tuple[str, str, str, str, str]] = []

    manifests = sorted(METRICS.glob("r*_final_experiment_manifest.json"))
    for manifest_path in manifests:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for section, roots in SECTIONS.items():
            for name, meta in (manifest.get(section) or {}).items():
                recorded = meta.get("sha256") if isinstance(meta, dict) else None
                if not recorded:
                    continue
                short = Path(name).name
                path = locate(name, roots)
                if path is None:
                    missing.append((manifest_path.name, section, short))
                    continue
                actual = sha256(path)
                if actual == recorded:
                    matched[section] = matched.get(section, 0) + 1
                elif (manifest_path.name, section, short) in known:
                    documented.append((manifest_path.name, section, short))
                else:
                    unexpected.append(
                        (manifest_path.name, section, short, recorded, actual))

    total_matched = sum(matched.values())
    print(f"manifests scanned      : {len(manifests)}")
    print(f"digests matching       : {total_matched}")
    for section in sorted(matched):
        print(f"    {section:20s} {matched[section]}")
    print(f"files missing          : {len(missing)}")
    print(f"documented differences : {len(documented)}")
    print(f"undocumented           : {len(unexpected)}")

    for manifest_name, section, short in missing:
        print(f"  [missing]     {manifest_name}: {section}/{short}")
    for manifest_name, section, short in documented:
        print(f"  [known-stale] {manifest_name}: {section}/{short}")
    for manifest_name, section, short, recorded, actual in unexpected:
        print(f"  [DRIFT]       {manifest_name}: {section}/{short}")
        print(f"                recorded {recorded}")
        print(f"                actual   {actual}")

    if documented and not args.strict:
        print(f"\n{len(documented)} known difference(s); see "
              f"known_stale_digests.json and docs/DIGEST_STATUS.md. "
              f"Run with --strict to treat them as failures.")

    return 1 if (missing or unexpected) else 0


if __name__ == "__main__":
    raise SystemExit(main())
