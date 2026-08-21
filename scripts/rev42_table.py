"""LaTeX table for the completed forced-variational panel (R42).

Same two columns of evidence the R37 table carries -- which panel orbits the
resolution rule leaves undecided, and which the calibration channel places
outside the band the original eight-orbit panel occupied -- read from
metrics/r42_panel_verdict.json. The table builder is R37's, imported and
pointed at the new verdict, so the presentation cannot drift from the one the
manuscript already shows.

Usage:  python rev42_table.py
"""

from __future__ import annotations

from pathlib import Path

import rev37_table as builder

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

builder.VERDICT = METRICS / "r42_panel_verdict.json"
builder.OUT = METRICS / "r42_panel_extension_table.tex"

if __name__ == "__main__":
    raise SystemExit(builder.main())
