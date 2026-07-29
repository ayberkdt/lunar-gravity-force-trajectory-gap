"""Rebuild known_stale_digests.json from the current manifests and files.

    python tools/record_known_stale.py

Covers every digest a campaign manifest records: driver scripts, generated
tables and figures, and result records. Any note already written against an
entry is carried over, so evidence about a difference is not lost when the
file is regenerated.

Use this only after deliberately resolving or accepting a difference. It
records every mismatch that exists right now, so running it blindly would
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
OUT = ROOT / "known_stale_digests.json"

# manifest section -> directories to look the artifact up in
SECTIONS = {
    "scripts": ("scripts",),
    "generated_tables": ("metrics",),
    "generated_figures": ("figures", "metrics"),
    "result_json": ("metrics",),
}

WHAT_THIS_IS = (
    "Files whose current bytes differ from the digest their campaign manifest "
    "recorded, meaning they changed after that manifest was frozen. The "
    "digests are NOT refreshed: re-hashing would assert that the current file "
    "is the one the manifest describes, which is a separate claim from whether "
    "it still yields the same numbers. Where re-running the driver has been "
    "shown to regenerate the archived output byte-for-byte, that evidence is "
    "recorded in the entry's 'reproduction' field and in "
    "docs/DIGEST_STATUS.md."
)

EFFECT_ON_CI = (
    "audit_manifest_digests.py fails on any mismatch that is NOT listed here, "
    "so new drift breaks the build. Listed entries are reported as known and "
    "do not fail. Run with --strict to fail on these as well."
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def locate(name: str, roots: tuple[str, ...]) -> Path | None:
    for directory in roots:
        candidate = ROOT / directory / Path(name).name
        if candidate.exists():
            return candidate
    return None


def main() -> int:
    previous: dict[tuple[str, str, str], dict] = {}
    if OUT.exists():
        for entry in json.loads(OUT.read_text(encoding="utf-8")).get("entries", []):
            key = (entry["manifest"], entry.get("category", "scripts"),
                   entry["file"] if "file" in entry else entry.get("script", ""))
            previous[key] = entry

    entries = []
    for manifest_path in sorted(METRICS.glob("r*_final_experiment_manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for section, roots in SECTIONS.items():
            for name, meta in (manifest.get(section) or {}).items():
                recorded = meta.get("sha256") if isinstance(meta, dict) else None
                if not recorded:
                    continue
                path = locate(name, roots)
                if path is None:
                    continue
                observed = sha256(path)
                if observed == recorded:
                    continue
                key = (manifest_path.name, section, Path(name).name)
                entry = {
                    "manifest": manifest_path.name,
                    "category": section,
                    "file": Path(name).name,
                    "recorded_sha256": recorded,
                    "recorded_bytes": meta.get("bytes"),
                    "observed_sha256": observed,
                    "observed_bytes": path.stat().st_size,
                }
                carried = previous.get(key, {}).get("reproduction")
                if carried:
                    entry["reproduction"] = carried
                entries.append(entry)

    added = [e for e in entries
             if (e["manifest"], e["category"], e["file"]) not in previous]

    payload = {
        "schema": "known_stale_digests_v2",
        "recorded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "what_this_is": WHAT_THIS_IS,
        "effect_on_ci": EFFECT_ON_CI,
        "entries": entries,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    by_category: dict[str, int] = {}
    for entry in entries:
        by_category[entry["category"]] = by_category.get(entry["category"], 0) + 1
    print(f"wrote {OUT.name}: {len(entries)} entries "
          f"({len(previous)} before, {len(added)} newly absorbed)")
    for category, count in sorted(by_category.items()):
        print(f"  {category:20s} {count}")
    for entry in added:
        print(f"  [new] {entry['manifest']}: {entry['category']}/{entry['file']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
