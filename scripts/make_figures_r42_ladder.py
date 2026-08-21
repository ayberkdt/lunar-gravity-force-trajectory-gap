"""Instrument-ladder figure with the completed panel on its third rung.

make_figures_r34.py forms the two R37 filenames inside ``load_variational``
itself, so the variational column of ``fig_instrument_ladder.pdf`` is drawn from
the 80-orbit level even after the chain was completed to 128. This script does
for the figure what rev42_instrument_ladder.py did for the table and
make_figures_r42.py did for the parity figure: it builds a metrics view in a
scratch directory --- hard links to the real files, so nothing is duplicated on
disk and nothing can be written back --- puts the R42 record and verdict in it
under the two R37 names, and calls the sealed maker unchanged.

The substitution is only which scored record the third column draws. The panel
membership rule, the column set, the jitter seed, the median marks and the
figure itself are make_figures_r34.py's. The other three columns read R14
records that the substitution does not touch, so they are bit-identical to the
sealed figure.

Usage:  python make_figures_r42_ladder.py
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import make_figures_r34 as maker

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SUBSTITUTE = {
    "r37_variational_extension.json": "r42_variational_completion.json",
    "r37_panel_verdict.json": "r42_panel_verdict.json",
}


def main() -> int:
    view = Path(tempfile.mkdtemp(prefix="r42_ladder_view_"))
    try:
        for src in METRICS.glob("*.json"):
            if src.name in SUBSTITUTE:
                continue
            try:
                os.link(src, view / src.name)
            except OSError:
                shutil.copy2(src, view / src.name)
        for name, real in SUBSTITUTE.items():
            shutil.copy2(METRICS / real, view / name)
        maker.METRICS = view
        maker.figure_ladder()
    finally:
        maker.METRICS = METRICS
        shutil.rmtree(view, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
