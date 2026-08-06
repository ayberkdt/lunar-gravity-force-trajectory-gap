"""One place that answers "what populations does this registration declare?"

Two registrations carry propagated populations that are neither design A, B nor
C. R30 declares five geometry strata under a ``strata`` map; R31 declares one
out-of-box elliptical population under a single ``design``. The stratum drivers
work on either, so they need one accessor rather than a branch in each of them.

The two registrations are deliberately not merged into one shape. An R30
stratum is a sub-box of a range that was already sampled and can only locate an
existing claim; R31 widens the range and is a scope extension. Giving them one
registration would make that distinction a field rather than a document, and it
is the kind of distinction that gets lost when it is a field.

Every entry carries the same three things the drivers need: the frozen design
file, its hash, and the design key that appears in trajectory filenames.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def registration(registry: str) -> dict:
    return json.loads(
        (METRICS / f"{registry}_preregistration.json").read_text(
            encoding="utf-8"))


def populations(registry: str) -> dict:
    """name -> {file, design_sha256, design_key}, whichever shape it is in."""
    reg = registration(registry)
    if "strata" in reg:
        return dict(reg["strata"])
    design = reg["design"]
    return {reg["population"]: {"file": design["file"],
                                "design_sha256": design["design_sha256"],
                                "design_key": design["design_key"]}}


def spec(registry: str, name: str) -> dict:
    pops = populations(registry)
    if name not in pops:
        raise SystemExit(f"{name} is not declared in "
                         f"{registry}_preregistration.json: {sorted(pops)}")
    return pops[name]
