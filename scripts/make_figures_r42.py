"""Panel parity figure with the completed panel, without editing the sealed maker.

make_figures_r19.py reads the panel from the two R37 filenames, which it forms
inside the plotting function itself; there is no constant to rebind. Rather than
copy the plotting code, this script builds a metrics view in a scratch directory
--- hard links to the real files, so nothing is duplicated on disk and nothing
can be written back --- and puts the R42 record and verdict in it under the two
R37 names. The maker is then pointed at that view and called unchanged.

The substitution is only which scored record the figure draws. The panel
membership rule, the open-marker convention for undecided orbits and the figure
itself are make_figures_r19.py's.

Usage:  python make_figures_r42.py
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import make_figures_r19 as maker

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SUBSTITUTE = {
    "r37_variational_extension.json": "r42_variational_completion.json",
    "r37_panel_verdict.json": "r42_panel_verdict.json",
}


def main() -> int:
    view = Path(tempfile.mkdtemp(prefix="r42_metrics_view_"))
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
        # maker.main(), not variational_parity() directly: the style
        # application lives in main(), and bypassing it drew this figure in
        # matplotlib defaults next to figures in the paper's style.
        maker.main()
    finally:
        maker.METRICS = METRICS
        shutil.rmtree(view, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
