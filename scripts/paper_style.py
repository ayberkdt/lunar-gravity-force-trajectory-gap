"""Shared publication figure style: Computer Modern (usetex), Okabe-Ito
colors, inward ticks, subtle grid. Import and call apply() before plotting."""

from __future__ import annotations

import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette
C1 = "#0072B2"  # blue
C2 = "#E69F00"  # orange
C3 = "#009E73"  # bluish green
C4 = "#CC79A7"  # reddish purple
C5 = "#D55E00"  # vermillion
C6 = "#56B4E9"  # sky blue
GRAY = "0.45"


def apply() -> None:
    use_tex = os.environ.get("PAPER_USE_TEX", "1") != "0"
    plt.rcParams.update({
        "text.usetex": use_tex,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage[T1]{fontenc}",
        "font.family": "serif",
        "font.size": 9.0,
        "axes.labelsize": 9.0,
        "legend.fontsize": 8.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.3,
        "lines.markersize": 4.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.45,
        "ytick.minor.width": 0.45,
        "axes.grid": True,
        "grid.color": "0.5",
        "grid.alpha": 0.22,
        "grid.linewidth": 0.45,
        "legend.frameon": False,
        "legend.handlelength": 1.9,
        "legend.borderaxespad": 0.4,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })
