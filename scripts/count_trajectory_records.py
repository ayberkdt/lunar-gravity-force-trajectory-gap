"""Count the distinct trajectory records the campaign manifests index.

The supplement quotes this number beside its non-duplication claim. The full
integrity gate recomputes it, but only after hashing every recorded file,
which is minutes of I/O over ten gigabytes. The number itself is a property of
what the manifests *say*, not of the bytes, so this reads the manifests alone
and reuses the gate's own parser and its trajectory pattern rather than
restating either.

Usage:  python count_trajectory_records.py
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import check_manifest_integrity as gate  # noqa: E402


def main() -> int:
    paths = (sorted(glob.glob(str(gate.METRICS
                                  / "r[0-9]*_final_experiment_manifest.json")))
             + sorted(glob.glob(str(gate.METRICS
                                    / "rJ_final_experiment_manifest.json"))))
    owners: dict[str, list[str]] = {}
    for p in paths:
        payload = json.loads(Path(p).read_text(encoding="utf-8"))
        out: list = []
        gate.collect(payload, "root", out)
        seen = set()
        for name, _digest, _missing, _section in out:
            norm = name.replace("\\", "/")
            if not norm.startswith("metrics/"):
                norm = f"metrics/{norm}"
            if gate.TRAJECTORY.match(norm) and norm not in seen:
                seen.add(norm)
                owners.setdefault(norm, []).append(Path(p).name)
    shared = {k: v for k, v in owners.items() if len(v) > 1}
    print(f"manifests: {len(paths)}")
    print(f"distinct trajectory records: {len(owners)}")
    print(f"records under more than one manifest: {len(shared)}")
    for k in sorted(shared)[:5]:
        print(f"    {k} <- {', '.join(shared[k])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
