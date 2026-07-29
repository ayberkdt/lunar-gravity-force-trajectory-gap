"""Report every environment that appears in the archived provenance records.

    python environments/show_environments.py

The drivers wrote interpreter and package versions into the result records as
they ran. This groups those records by the environment they recorded, which is
what the table in this directory's README is built from. Run it after adding a
campaign to check the README is still complete.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"

FIELDS = ("python", "numpy", "scipy", "numba", "tudatpy")


def scan(text: str, field: str) -> str | None:
    for pattern in (rf'"{field}_version"\s*:\s*"([^"]+)"',
                    rf'"{field}"\s*:\s*"([^"]+)"',
                    rf'"[a-z_]*_{field}"\s*:\s*"([^"]+)"'):
        found = re.search(pattern, text)
        if found:
            return found.group(1).split()[0]
    return None


def main() -> int:
    groups: dict[tuple, list[str]] = defaultdict(list)
    for path in sorted(METRICS.glob("*.json")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        found = {f: scan(text, f) for f in FIELDS}
        if not any(found[f] for f in ("python", "numpy", "scipy")):
            continue
        groups[tuple(found[f] for f in FIELDS)].append(path.name)

    kernels = set()
    for path in sorted(METRICS.glob("r*_final_experiment_manifest.json")):
        manifest = json.loads(path.read_text(encoding="utf-8"))
        kernel = manifest.get("numerical_kernel") or {}
        if kernel.get("lunaris_commit"):
            kernels.add((kernel.get("lunaris_release_tag"),
                         kernel["lunaris_commit"]))

    for signature, files in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        labelled = ", ".join(f"{name} {value}"
                             for name, value in zip(FIELDS, signature) if value)
        print(f"\n{labelled}")
        print(f"  {len(files)} record(s)")
        campaigns = sorted({name.split("_")[0] for name in files})
        print(f"  campaigns: {' '.join(campaigns)}")

    print("\nPinned source snapshot(s):")
    for tag, commit in sorted(kernels):
        print(f"  {tag or '(untagged)'}  {commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
