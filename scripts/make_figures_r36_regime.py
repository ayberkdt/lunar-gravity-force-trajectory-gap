"""R36: the budget-geometry regime map.

One figure answering the operational question the paper's separate results only
answer jointly: for a given orbit population and a given declared per-call
budget, which allocation wins?

  panel (a)  budget-calibrated radial endpoint against its equal-budget
             constant degree            (R14 propagated ladder)
  panel (b)  the interior member k = 0.5 against a constant degree matched on
             realized total quadratic work   (R19)

Colour is log10 of the median error ratio E_fixed / E_policy over the
population, so positive (blue) favours the varying-degree policy and negative
(red) favours the constant degree, matching the convention already used by
fig_doe_regime.pdf. Cell text is the resolved win-loss count.

Three properties of the evidence are drawn rather than left to the caption:

  * rows are populations, never pooled. R30's registration forbids pooling the
    strata with each other or with the coverage designs, because the sub-boxes
    overlap. A, B and C are three independent draws of one factor box; SH is a
    sub-box of that same box; OE lies outside it. The row grouping and the rule
    between groups say so, and no marginal is ever drawn.
  * cells whose comparison is mostly undecided by the resolution rule are
    hatched, and cells with almost nothing resolved are drawn blank: a median
    ratio computed where one comparison of 57 resolves is not a result.
  * beta = 0.62 is a declared post-hoc localization under the R28 amendment,
    which commits to labelling it as such wherever it appears. Its column is
    marked and separated from the pre-registered grid.

Reads only frozen campaign records; runs no propagation.

Usage:  python make_figures_r36_regime.py
"""

from __future__ import annotations

import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import paper_style
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

paper_style.apply()

ROOT = pathlib.Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"


# Populations, ordered by the radial span of the range they sample. The order
# is a reading aid, not a partition: see the module docstring.
# OEU is not a sixth population. It is OE repropagated with its reference
# degree raised so the calibrated schedule is never clipped at it, and it sits
# directly under its parent because the pair is the reading, not either row
# alone. Where the ceiling bound, the capped row credits the radial policy with
# a defect that vanished by construction.
ROWS = [("A", "A"), ("B", "B"), ("C", "C"), ("SH", "SH"), ("OE", "OE"),
        ("OEU", "OE$^{\\ast}$")]
GROUPS = [(0, 3, "coverage box", r"$31$--$150\times181$--$599$"),
          (3, 4, "sub-box of it", r"$80$--$149\times521$--$600$"),
          (4, 6, "outside it", r"$80$--$119\times701$--$2479$")]

# beta = 3.00 is dropped from the plotted grid: it carries a single cell, on
# design A in panel (a), and the resolution rule leaves that one undecided, so
# drawing a column for it adds an empty band and no verdict. The caption says so.
BETAS = ["0.50", "0.62", "0.75", "1.00", "1.25", "1.50"]
POSTHOC = {"0.62"}

MIN_RESOLVED = 6
HATCH_ABOVE = 0.5


def panel_a():
    """Radial endpoint against its equal-budget constant degree."""
    out = {}
    for key, _ in ROWS:
        for p in sorted(METRICS.glob(f"r14_trajectory_{key}_beta_*.json")):
            beta = re.search(r"beta_(\d\.\d+)", p.name).group(1)
            s = json.loads(p.read_text(encoding="utf-8"))["summary"]
            out[(key, beta)] = {
                "rho": s["rho_budget"]["median"],
                "policy_wins": s["resolved_atallah_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "resolved": s["resolved"],
                "orbits": s["orbits"],
            }
    return out


def panel_b():
    """Interior member k = 0.5 against a work-matched constant degree."""
    out = {}
    for key, _ in ROWS:
        for p in sorted(METRICS.glob(f"r19_equal_total_work_{key}*.json")):
            m = re.search(r"beta_(\d\.\d+)", p.name)
            beta = m.group(1) if m else "1.00"
            s = json.loads(p.read_text(encoding="utf-8"))["summary"]
            out[(key, beta)] = {
                "rho": s["median_rho"],
                "policy_wins": s["resolved_interior_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "resolved": s["resolved"],
                "orbits": s["orbits"],
            }
    return out


def draw(ax, data, title, cmap, norm, show_xlabel=True):
    n_r, n_b = len(ROWS), len(BETAS)

    # the post-hoc budget gets a quiet background band instead of rules
    for j, beta in enumerate(BETAS):
        if beta in POSTHOC:
            ax.add_patch(plt.Rectangle((j, 0), 1, n_r, facecolor="0.955",
                                       edgecolor="none", zorder=0))

    for i, (key, _) in enumerate(ROWS):
        for j, beta in enumerate(BETAS):
            cell = data.get((key, beta))
            if cell is None:
                # never propagated: left blank, so the eye sees only evidence
                ax.add_patch(plt.Rectangle((j + .06, i + .06), .88, .88,
                                           facecolor="none", edgecolor="0.86",
                                           lw=0.45, ls=(0, (1.6, 1.6)),
                                           zorder=2))
                continue
            undecided = 1.0 - cell["resolved"] / max(cell["orbits"], 1)
            blank = cell["resolved"] < MIN_RESOLVED
            val = np.log10(max(cell["rho"], 1e-12))
            face = "white" if blank else cmap(norm(val))
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=face,
                                       edgecolor="white", lw=1.1, zorder=2))
            if blank:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none",
                                           edgecolor="0.7", lw=0.5,
                                           hatch="/////", zorder=3))
                ax.text(j + .5, i + .5, "n.r.", ha="center", va="center",
                        fontsize=6.0, color="0.45", zorder=4)
                continue
            if undecided > HATCH_ABOVE:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none",
                                           edgecolor="0.42", lw=0.0,
                                           hatch="/////", alpha=0.40, zorder=3))
            dark = abs(norm(val) - 0.5) > 0.32
            ax.text(j + .5, i + .5,
                    f"{cell['policy_wins']}--{cell['fixed_wins']}",
                    ha="center", va="center", fontsize=6.4, zorder=4,
                    color="white" if dark else "0.12")

    ax.set_xlim(0, n_b)
    ax.set_ylim(n_r, 0)
    ax.set_xticks(np.arange(n_b) + 0.5)
    ax.set_xticklabels([rf"${b}^{{\dagger}}$" if b in POSTHOC else f"${b}$"
                        for b in BETAS], fontsize=8)
    ax.set_yticks(np.arange(n_r) + 0.5)
    ax.set_yticklabels([lab for _, lab in ROWS], fontsize=8)
    if show_xlabel:
        ax.set_xlabel(r"declared per-call budget $\beta$", labelpad=2)
    ax.set_title(title, fontsize=8.2, pad=4)
    ax.grid(False)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    # the rows above and below are different populations, not one ordered family
    for start, stop, _, _ in GROUPS[:-1]:
        ax.plot([0, n_b], [stop, stop], color="0.3", lw=0.9,
                solid_capstyle="butt", zorder=5, clip_on=False)


def main() -> int:
    a, b = panel_a(), panel_b()
    norm = TwoSlopeNorm(vmin=-2.2, vcenter=0.0, vmax=2.2)
    cmap = plt.get_cmap("RdBu")

    fig = plt.figure(figsize=(6.46, 2.85), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[0.235, 1.0], hspace=0.10,
                          wspace=0.02)
    lab_axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    axes = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

    draw(axes[0], a,
         r"(a) budget-calibrated radial endpoint vs.\ equal-budget "
         r"constant degree", cmap, norm, show_xlabel=False)
    draw(axes[1], b,
         r"(b) interior member $k{=}0.5$ vs.\ work-matched constant degree "
         r"\emph{(exploratory)}", cmap, norm, show_xlabel=True)

    for ax_lab in lab_axes:
        ax_lab.set_xlim(0, 1)
        ax_lab.set_ylim(len(ROWS), 0)
        ax_lab.axis("off")
        for start, stop, name, geom in GROUPS:
            mid = (start + stop) / 2
            ax_lab.text(0.98, mid - 0.14, name, ha="right", va="center",
                        fontsize=6.4, color="0.2")
            ax_lab.text(0.98, mid + 0.28, geom, ha="right", va="center",
                        fontsize=5.8, color="0.45")
            ax_lab.plot([0.995, 0.995], [start + .12, stop - .12],
                        color="0.75", lw=0.8, clip_on=False)

    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                      ax=axes, location="right", fraction=0.028, pad=0.010,
                      extend="both", ticks=[-2, -1, 0, 1, 2])
    cb.set_label(r"constant $\leftarrow$ $\log_{10}(E_{\mathrm{fixed}}/E_{\mathrm{policy}})$ $\rightarrow$ varying",
                 fontsize=7.0, labelpad=3)
    cb.ax.tick_params(labelsize=7)
    # direction words sit outside the bar's own axes, as its title and its
    # x-label, so they cannot land on the extend triangles or the tick labels
    # direction words sit outside the bar's own axes, as its title
    # and its x-label, so they cannot land on the extend triangles
    cb.outline.set_linewidth(0.4)

    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "fig_regime_map.pdf")
    plt.close(fig)
    print("[written] figures/fig_regime_map.pdf")

    rec = {"schema": "r36_regime_map_v1",
           "note": ("median rho is E_fixed/E_policy over the whole population; "
                    "win-loss counts are resolved comparisons only. Rows are "
                    "separate populations and are never pooled."),
           "plotted_budgets": BETAS,
           "omitted_from_plot": {
               "3.00": ("one cell only (design A, panel a) and the resolution "
                        "rule leaves it undecided; it is in the record and in "
                        "r14_trajectory_A_beta_3.00.json")},
           "min_resolved_for_colour": MIN_RESOLVED,
           "hatch_above_undecided_fraction": HATCH_ABOVE,
           "post_hoc_budgets": sorted(POSTHOC),
           "panels": {"a_radial_vs_fixed": {f"{k}|{v}": val
                                            for (k, v), val in a.items()},
                      "b_interior_vs_workmatched": {f"{k}|{v}": val
                                                    for (k, v), val in b.items()}}}
    (METRICS / "r36_regime_map.json").write_text(
        json.dumps(rec, indent=2), encoding="utf-8")
    print("[written] metrics/r36_regime_map.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
