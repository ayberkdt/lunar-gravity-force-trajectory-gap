"""R36: the budget-geometry regime map.

One figure answering the operational question the paper's separate results only
answer jointly: for a given orbit population and a given declared per-call
budget, which allocation wins?

  panel (a)  budget-calibrated radial endpoint against its equal-budget
             constant degree            (R14 propagated ladder)
  panel (b)  the interior member k = 0.5 against a constant degree matched on
             realized total quadratic work   (R19, tight-level match; the
             scoring-tolerance rematch R44 covers designs A and B only and
             is not drawn, so each panel carries one accounting)

Colour is log10 of the median error ratio E_fixed / E_policy over the
population, so positive (blue) favours the varying-degree policy and negative
(red) favours the constant degree, matching the convention already used by
fig_doe_regime.pdf. Cell text is the resolved win-loss count.

Three properties of the evidence are drawn rather than left to the caption:

  * rows are populations, never pooled. R30's registration forbids pooling the
    strata with each other or with the coverage designs, because the sub-boxes
    overlap. A, B and C are three independent draws of one factor box; SH, SL,
    SP, SE and SF are sub-boxes of that same box, the first two restricting a
    radial coordinate and the other three an angle; OE lies outside it. SH and
    SL carry a group each rather than a shared one, because they restrict
    opposite ends of the radial range and a single geometry label could not
    describe both. The row grouping and the rules between groups say so, and no
    marginal is ever drawn.

    The span-ladder rows are the one exception and are marked as one. They are
    apolune levels of a single paired population rather than separate
    populations, and each pools the two identity blocks of that population,
    which its own registration licenses because the blocks are the same design
    at the same levels. What is still never pooled is one level with another:
    the levels are what the ladder varies, and averaging them would erase the
    only quantity it measures.
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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

paper_style.apply()

ROOT = pathlib.Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"

# Two renderings of one map. The annotated form carries both panels, all
# fourteen rows and the resolved win--loss count in every cell, and is the
# supplement's. The plain form is the main text's: panel (a) alone, without
# the three angle sub-boxes and without the counts, because at main-text width
# a second exploratory panel and sixty-five small counts cost more reading
# than they repay. Colour scale, cell states and marks are identical.
ANNOTATE = True


# Populations, ordered by the radial span of the range they sample. The order
# is a reading aid, not a partition: see the module docstring.
# OEU is not a sixth population. It is OE repropagated with its reference
# degree raised so the calibrated schedule is never clipped at it, and it sits
# directly under its parent because the pair is the reading, not either row
# alone. Where the ceiling bound, the capped row credits the radial policy with
# a defect that vanished by construction.
# Row order carries the reading. The three angle sub-boxes sit together
# because they share a verdict; the span sub-box then sits directly above the
# out-of-box population because those two are the rows that move the crossing,
# and putting them adjacent lets the eye follow the radial-span argument down
# the figure instead of hopping over the angle block.
ROWS = [("A", "A"), ("B", "B"), ("C", "C"),
        ("SP", "polar"), ("SE", "equatorial"), ("SF", "frozen-like"),
        ("SH", "high span"), ("SL", "low perilune"),
        ("OE", "wide elliptic"), ("OEU", "wide ell., ceiling-free"),
        ("RS300", "$300$"), ("RS600", "$600$"),
        ("RS1200", "$1200$"), ("RS2400", "$2400$")]
# The ladder rows carried a level-wide ceiling mark while only one identity
# block had been re-propagated without the ceiling. Both blocks now have been,
# at the four registered budgets, so those cells are drawn ceiling-free and the
# mark moves to where the ceiling still stands: the two amended budgets, which
# have no ceiling-free counterpart on either block. It is a per-cell mark now
# rather than a per-row one, because within a row the two readings differ by
# column.

# The paired ladder: one population read at four prescribed apolunes, with
# every other element of an identity held fixed across them. Its rows are read
# down the column, which is the axis the ladder varies; the populations above
# are read across.
R50_BLOCKS = {"RS1": "r50_span_ladder_a_design_frozen.json",
              "RS2": "r50_span_ladder_b_design_frozen.json"}
# The ceiling-free re-propagation of each block (O50 for A, O51 for B). It
# covers the registered grid only, so the amended budgets fall back to the
# capped parent and are marked where they are drawn. The uncapped designs carry
# the parent's sobol indices unchanged, so the level index is the parent's.
UNCAPPED = {"RS1": "RS1U", "RS2": "RS2U"}
CEILING_FREE_BETAS = {"1.00", "0.75", "0.62", "0.50"}


def ladder_source(parent: str, beta: str) -> str:
    return UNCAPPED[parent] if beta in CEILING_FREE_BETAS else parent


def is_capped_cell(row_key: str, beta: str) -> bool:
    return row_key.startswith("RS") and beta not in CEILING_FREE_BETAS
# The five strata are sub-boxes of one box but not of one kind, and the
# distinction is the reading: SH restricts the apolune end of the radial range
# and SL the perilune end, while SP, SE and SF restrict an angle and leave both
# radial ranges at the box. They therefore carry a group label each rather than
# one averaged range.
GROUPS = [(0, 3, "coverage box", r"$31$--$150\times181$--$599$"),
          (3, 6, "sub-boxes in angle", r"$31$--$150\times183$--$598$"),
          (6, 7, "sub-box in span", r"$80$--$149\times521$--$600$"),
          (7, 8, "sub-box in perilune", r"$30$--$50\times182$--$594$"),
          (8, 10, "outside it", r"$80$--$119\times701$--$2479$"),
          (10, 14, "paired ladder, apolune", r"$82$--$120\times$ level, km")]

# Groups are separated by a gap rather than by a rule drawn across the cells.
# The gap is what carries the figure's structure: the reader is meant to see
# four population classes, not ten rows. At 0.42 the blocks read as one field
# with faint seams, which is the main reason this panel takes longer to decode
# than the paper's other figures.
GROUP_GAP = 0.80


def _group_of(i):
    for g, (start, stop, _, _) in enumerate(GROUPS):
        if start <= i < stop:
            return g
    return len(GROUPS) - 1


def _y(i):
    """Top edge of row i once the inter-group gaps are inserted."""
    return i + _group_of(i) * GROUP_GAP


TOTAL_H = len(ROWS) + (len(GROUPS) - 1) * GROUP_GAP

# The main-text figure omits the three angle sub-boxes. This is a presentation
# choice with a reason, not a convenience: the paper's own reading of those
# rows is that restricting inclination or the argument of perilune leaves both
# verdicts in place, so they carry a null result that one clause of prose
# states completely, while the figure in the main text exists to show where the
# crossing moves. Dropping them removes a whole group and three of fourteen
# rows from a panel a reader has to decode twice.
#
# Only the drawing is filtered. panel_a() and panel_b() still read every
# population, the record still contains every cell, and the annotated
# supplement figure still draws all fourteen rows, so nothing here decides what
# evidence exists -- only what the main text's copy of the map shows. The
# caption has to say so, and it does.
MAIN_DROP = {"SP", "SE", "SF"}


def _select_plot_rows() -> None:
    """Narrow ROWS/GROUPS to what the main-text figure draws."""
    global ROWS, GROUPS, TOTAL_H
    if ANNOTATE:
        return
    kept, groups, i = [], [], 0
    for start, stop, name, geom in GROUPS:
        block = [r for r in ROWS[start:stop] if r[0] not in MAIN_DROP]
        if not block:
            continue
        kept.extend(block)
        groups.append((i, i + len(block), name, geom))
        i += len(block)
    ROWS, GROUPS = kept, groups
    TOTAL_H = len(ROWS) + (len(GROUPS) - 1) * GROUP_GAP

# budgets above 1.50 are dropped from the plotted grid: design A was never
# propagated at 2.00, and the three cells that exist there (B at 2.00, both
# designs at 3.00) sit below or barely at the colour floor except B at 2.00,
# whose row the main text's trajectory table carries. The caption says so.
BETAS = ["0.50", "0.62", "0.75", "1.00", "1.25", "1.50"]
POSTHOC = {"0.62"}

MIN_RESOLVED = 6
HATCH_ABOVE = 0.5


def _ladder_levels():
    """sobol index -> apolune level, per block, from the frozen designs."""
    idx = {}
    for key, fname in R50_BLOCKS.items():
        d = json.loads((METRICS / fname).read_text(encoding="utf-8"))
        idx[key] = {o["sobol_index"]: o["apolune_level_km"]
                    for o in d["orbits"]}
    return idx


def ladder_cells(kind):
    """Per-level cells of the paired ladder, pooled over its two blocks.

    The level records do not exist as files: a block's record is one population
    of sixty-four orbits carrying all four levels, so the cells are accumulated
    from the per-orbit rows rather than read from a summary block.
    """
    idx = _ladder_levels()
    acc = {}
    for parent in R50_BLOCKS:
        for beta in BETAS:
            src = ladder_source(parent, beta)
            p = METRICS / ((f"r14_trajectory_{src}_beta_{beta}.json")
                           if kind == "a"
                           else f"r19_equal_total_work_{src}_beta_{beta}.json")
            if not p.exists():
                continue
            for r in json.loads(p.read_text(encoding="utf-8"))["rows"]:
                level = idx[parent][r["sobol_index"]]
                if kind == "a":
                    c = r["comparison"]
                    rho = c.get("rho_budget")
                    won = (c.get("resolved_winner")
                           or c.get("raw_winner")) == "atallah"
                    resolved = bool(c["resolved"])
                else:
                    rho = r.get("rho_workmatched")
                    won = r["winner"] == "interior"
                    resolved = bool(r["resolved"])
                cell = acc.setdefault((f"RS{level:.0f}", beta),
                                      {"rho": [], "policy_wins": 0,
                                       "fixed_wins": 0, "resolved": 0,
                                       "orbits": 0})
                cell["orbits"] += 1
                if rho is not None:
                    cell["rho"].append(float(rho))
                if resolved:
                    cell["resolved"] += 1
                    cell["policy_wins" if won else "fixed_wins"] += 1
    return {k: {"rho": float(np.median(c["rho"])) if c["rho"] else 1.0,
                "policy_wins": c["policy_wins"], "fixed_wins": c["fixed_wins"],
                "resolved": c["resolved"], "orbits": c["orbits"]}
            for k, c in acc.items()}


def posthoc_cells_scored():
    """(population, beta) pairs the post-hoc budget campaign actually scored.

    The 0.62 column of the seven R53 populations is drawn from that campaign's
    verdict rather than from whatever files happen to be on disk. A cell whose
    chain did not complete leaves a trajectory record behind, and globbing
    would draw it: the campaign's own supervisor refuses to quote such a cell,
    and so must the figure. Designs A and B carry 0.62 from (O34) and are not
    governed by this record.
    """
    p = METRICS / "r53_verdict.json"
    if not p.exists():
        return None
    v = json.loads(p.read_text(encoding="utf-8"))
    governed = set(v["report"])
    scored = set(v["cells_run"])
    return governed, scored


def _posthoc_blocked(key, beta, gate):
    if gate is None or beta != "0.62":
        return False
    governed, scored = gate
    return key in governed and key not in scored


def panel_a():
    """Radial endpoint against its equal-budget constant degree."""
    out = {}
    gate = posthoc_cells_scored()
    for key, _ in ROWS:
        for p in sorted(METRICS.glob(f"r14_trajectory_{key}_beta_*.json")):
            beta = re.search(r"beta_(\d\.\d+)", p.name).group(1)
            if _posthoc_blocked(key, beta, gate):
                continue
            s = json.loads(p.read_text(encoding="utf-8"))["summary"]
            out[(key, beta)] = {
                "rho": s["rho_budget"]["median"],
                "policy_wins": s["resolved_atallah_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "resolved": s["resolved"],
                "orbits": s["orbits"],
            }
    out.update(ladder_cells("a"))
    return out


def panel_b():
    """Interior member k = 0.5 against a work-matched constant degree."""
    out = {}
    gate = posthoc_cells_scored()
    for key, _ in ROWS:
        # anchor the population key: an unanchored glob let OE also match the
        # OEU files, whose later sort position overwrote every OE cell of
        # panel (b) with the uncapped population's numbers
        for p in sorted(METRICS.glob(f"r19_equal_total_work_{key}.json")) + \
                sorted(METRICS.glob(f"r19_equal_total_work_{key}_beta_*.json")):
            m = re.search(r"beta_(\d\.\d+)", p.name)
            beta = m.group(1) if m else "1.00"
            if _posthoc_blocked(key, beta, gate):
                continue
            s = json.loads(p.read_text(encoding="utf-8"))["summary"]
            out[(key, beta)] = {
                "rho": s["median_rho"],
                "policy_wins": s["resolved_interior_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "resolved": s["resolved"],
                "orbits": s["orbits"],
            }
    out.update(ladder_cells("b"))
    return out


def _text_colour(rgba):
    """Black or white by the cell's own luminance, not by a norm threshold.

    The RdBu ramp is dark at both ends and pale in the middle, so a symmetric
    distance-from-centre test mislabels the light blues. Relative luminance
    reads the colour actually drawn.
    """
    r, g, b = rgba[:3]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if lum < 0.45 else "0.10"


def draw(ax, data, title, cmap, norm, show_xlabel=True):
    """title=None draws the panel bare: the main-text figure lets the caption
    name the comparison instead of printing it twice."""
    n_b = len(BETAS)

    # The post-hoc budget is marked by the dagger on its tick label alone. An
    # earlier version tinted the whole column, which collided with the tint
    # used for cells that were never propagated: at a glance the two greys
    # were the same grey, so a declared-post-hoc column read as missing data.

    for i, (key, _) in enumerate(ROWS):
        y = _y(i)
        for j, beta in enumerate(BETAS):
            cell = data.get((key, beta))
            if cell is None:
                # never propagated: a flat tint reads as absence at a glance,
                # where the previous dotted outline read as noise
                ax.add_patch(plt.Rectangle((j, y), 1, 1, facecolor="0.945",
                                           edgecolor="white", lw=1.1,
                                           zorder=2))
                continue
            undecided = 1.0 - cell["resolved"] / max(cell["orbits"], 1)
            blank = cell["resolved"] < MIN_RESOLVED
            val = np.log10(max(cell["rho"], 1e-12))
            face = (1, 1, 1, 1) if blank else cmap(norm(val))
            ax.add_patch(plt.Rectangle((j, y), 1, 1, facecolor=face,
                                       edgecolor="white", lw=1.1, zorder=2))
            if blank:
                ax.text(j + .5, y + .5, "n.r.", ha="center", va="center",
                        fontsize=7.0 if ANNOTATE else 7.6, color="0.45", zorder=4)
                continue
            if undecided > HATCH_ABOVE:
                ax.add_patch(plt.Rectangle((j, y), 1, 1, facecolor="none",
                                           edgecolor=(1, 1, 1, 0.85), lw=0.0,
                                           hatch="////", zorder=3))
            if ANNOTATE:
                label = f"{cell['policy_wins']}--{cell['fixed_wins']}"
                if is_capped_cell(key, beta):
                    label += r"$^{\mathrm{c}}$"
            elif is_capped_cell(key, beta):
                label = "c"
            else:
                continue
            ax.text(j + .5, y + .5, label,
                    ha="center", va="center",
                    fontsize=7.2 if ANNOTATE else 8.0, zorder=4,
                    color=_text_colour(face))

    ax.set_xlim(0, n_b)
    ax.set_ylim(TOTAL_H, -0.35)
    ax.set_xticks(np.arange(n_b) + 0.5)
    ax.set_xticklabels([rf"${b}^{{\dagger}}$" if b in POSTHOC else f"${b}$"
                        for b in BETAS], fontsize=8 if ANNOTATE else 9)
    ax.set_yticks([_y(i) + 0.5 for i in range(len(ROWS))])
    ax.set_yticklabels([lab for _, lab in ROWS],
                       fontsize=8 if ANNOTATE else 9)
    if show_xlabel:
        ax.set_xlabel(r"declared per-call budget $\beta$", labelpad=3)
    if title is not None:
        ax.set_title(title, fontsize=8.2, pad=5)
    ax.grid(False)
    ax.minorticks_off()
    ax.tick_params(which="both", length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)


def main() -> int:
    # Data first, then the row filter: the panels and the record are built
    # from every population regardless of which ones this pass draws.
    a, b = panel_a(), panel_b()
    _select_plot_rows()
    norm = TwoSlopeNorm(vmin=-2.2, vcenter=0.0, vmax=2.2)
    cmap = plt.get_cmap("RdBu")

    # The main-text rendering draws panel (a) alone. Panel (b) is exploratory
    # and its full result already lives in the equal-work table and the prose
    # that reads it, so at main-text width the second panel repeated a result
    # the text carries while doubling what the reader has to decode. The
    # annotated supplement rendering keeps both panels: no evidence leaves the
    # paper, only the main text's copy of the map narrows to the claim it
    # supports.
    n_panels = 2 if ANNOTATE else 1

    # Height follows the row count rather than being fixed, so that adding a
    # population does not silently squash the cells of every other one. The
    # two-panel per-row height and fixed overhead are the ones the nine-row
    # version was drawn at. The single-panel rendering spends more height per
    # row on purpose: at 0.117 in the cells read as stripes, and a taller cell
    # is what lets the map read as a matrix at a glance.
    if ANNOTATE:
        height = 1.19 + 2 * TOTAL_H * 0.117
    else:
        height = 0.94 + TOTAL_H * 0.155
    fig = plt.figure(figsize=(6.46, height), constrained_layout=True)
    gs = fig.add_gridspec(n_panels, 2, width_ratios=[0.235, 1.0], hspace=0.10,
                          wspace=0.02)
    lab_axes = [fig.add_subplot(gs[k, 0]) for k in range(n_panels)]
    axes = [fig.add_subplot(gs[k, 1]) for k in range(n_panels)]

    if ANNOTATE:
        draw(axes[0], a,
             r"(a) radial endpoint vs.\ constant degree, equal per-call budget",
             cmap, norm, show_xlabel=False)
        draw(axes[1], b,
             r"(b) interior member $k{=}0.5$ vs.\ work-matched constant "
             r"\emph{(exploratory)}", cmap, norm, show_xlabel=True)
    else:
        # No in-figure title: the caption names the comparison.
        draw(axes[0], a, None, cmap, norm, show_xlabel=True)

    # Both panels carry the same rows, so the population names and their
    # geometry ranges are written once, beside panel (a). Repeating the whole
    # block beside panel (b) doubled the text in the margin and invited the
    # reader to check whether the two lists agreed, which is not a question
    # the figure should raise.
    # The single-panel rendering drops the geometry ranges from the margin:
    # the prose that reads the figure defines every row, the ranges live in
    # the populations table, and at main-text size the extra grey line per
    # group was most of what made the margin read as clutter. The annotated
    # supplement rendering keeps them.
    PLAIN_NAMES = {"paired ladder, apolune": "paired ladder,\napolune (km)"}
    for which, ax_lab in enumerate(lab_axes):
        ax_lab.set_xlim(0, 1)
        ax_lab.set_ylim(TOTAL_H, -0.35)
        ax_lab.axis("off")
        for start, stop, name, geom in GROUPS:
            top, bot = _y(start), _y(stop - 1) + 1.0
            mid = 0.5 * (top + bot)
            # a bracket spanning only its own rows, so the grouping is read
            # from the margin and no rule ever crosses a cell
            ax_lab.plot([0.985, 0.985], [top + .10, bot - .10],
                        color="0.62", lw=0.9, solid_capstyle="butt")
            if which:
                continue
            if ANNOTATE:
                ax_lab.text(0.940, mid - 0.32, name, ha="right", va="center",
                            fontsize=7.4 if ANNOTATE else 8.0, color="0.10",
                            fontweight="semibold")
                ax_lab.text(0.940, mid + 0.32, geom, ha="right", va="center",
                            fontsize=6.0 if ANNOTATE else 6.5, color="0.52")
            else:
                ax_lab.text(0.940, mid, PLAIN_NAMES.get(name, name),
                            ha="right", va="center",
                            fontsize=7.4 if ANNOTATE else 8.0,
                            color="0.10", fontweight="semibold",
                            linespacing=1.25)

    # One panel gets its colorbar underneath, horizontal and shortened: the
    # tall right-hand bar was sized for two stacked panels and would dominate
    # a single one. The two-panel rendering keeps the shared right-hand bar.
    if ANNOTATE:
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                          ax=axes, location="right", fraction=0.030, pad=0.012,
                          extend="both", ticks=[-2, -1, 0, 1, 2])
    else:
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                          ax=axes, location="bottom", fraction=0.050,
                          pad=0.022, shrink=0.55, aspect=34,
                          extend="both", ticks=[-2, -1, 0, 1, 2])
    # "median" belongs in the label, not the caption: the colour is a
    # population median, and the figure is the artefact that has to say so.
    # The plain rendering puts the bar under the panel, so its label is as wide
    # as the figure and anything longer is cut by the media box: at full length
    # the last word printed as "median favors vary". The arrow sides there drop
    # the repeated "median", which the leading term already supplies.
    arrows = ("median favors constant $\\leftarrow$ $\\rightarrow$ median favors varying"
              if ANNOTATE else
              "favors constant $\\leftarrow$ $\\rightarrow$ favors varying")
    cb.set_label(r"$\log_{10}$ median $E_{\mathrm{fix}}/E_{\mathrm{policy}}$: " + arrows,
                 fontsize=6.6 if ANNOTATE else 8.8,
                 labelpad=4 if ANNOTATE else 2)
    cb.ax.minorticks_off()
    cb.ax.tick_params(which="both", labelsize=7 if ANNOTATE else 7.6,
                      length=2.5, width=0.4,
                      pad=1.5)
    cb.outline.set_linewidth(0.4)

    # three cell states, named in the figure so the caption need not carry them
    legend_handles = [
        Patch(facecolor="0.90", edgecolor="white", hatch="////",
              label="over half undecided"),
        # "no scored cell" rather than "not propagated": one of these cells,
        # the near-equatorial sub-box at the post-hoc budget, has 64 propagated
        # trajectories on disk and no scored comparison, because the two stages
        # that score it never finished. Calling it unpropagated would print a
        # claim the campaign's own departure record contradicts.
        Patch(facecolor="0.945", edgecolor="white", label="no scored cell"),
        Patch(facecolor="white", edgecolor="0.75", lw=0.5,
              label="too few resolved (n.r.)"),
    ]
    # The two superscripts carried by the row and tick labels were defined in
    # the body text only, which left the panel unreadable on its own and let
    # the dagger collide with the one in Table 6. They are named here instead.
    # The plain rendering keeps the same five definitions but in the shortest
    # wording that still says what each mark asserts; the annotated rendering
    # keeps the full wording, whose qualifications the supplement's reader is
    # the one who needs.
    if ANNOTATE:
        mark_handles = [
            Line2D([], [], linestyle="none",
                   label=r"$^{\mathrm{c}}$ no ceiling-free run; capped values "
                         r"shown"),
            Line2D([], [], linestyle="none",
                   label=r"$^{\dagger}$ post hoc where added; registered for "
                         r"the ladder"),
        ]
        cells_legend = axes[-1].legend(
            handles=legend_handles, loc="upper center",
            bbox_to_anchor=(0.5, -0.20), ncol=3, frameon=False,
            fontsize=6.3, handlelength=1.5, handleheight=0.9,
            columnspacing=1.6, handletextpad=0.5)
        axes[-1].legend(handles=mark_handles, loc="upper center",
                        bbox_to_anchor=(0.5, -0.30), ncol=2, frameon=False,
                        fontsize=6.3, handlelength=0.0, handletextpad=0.0,
                        columnspacing=2.0)
        axes[-1].add_artist(cells_legend)
    else:
        mark_handles = [
            Line2D([], [], linestyle="none",
                   label=r"c: capped"),
            Line2D([], [], linestyle="none",
                   label=r"$\dagger$: post hoc"),
        ]
        # One panel frees the full figure width, so the cell states and the
        # two superscripts share a single strip under the colorbar instead of
        # stacking two rows. constrained_layout reserves the strip's height.
        fig.legend(handles=legend_handles + mark_handles,
                   loc="outside lower center", ncol=5, frameon=False,
                   fontsize=8.1 if ANNOTATE else 8.7, handlelength=1.1, handleheight=0.9,
                   columnspacing=0.9, handletextpad=0.35)

    FIGURES.mkdir(exist_ok=True)
    name = "fig_regime_map.pdf" if ANNOTATE else "fig_regime_map_main.pdf"
    # Without this the writer stamps the current time into the PDF, so an
    # unchanged figure hashes differently on every run and the manifest digest
    # reads stale for a reason that is not the content.
    fig.savefig(FIGURES / name, metadata={"CreationDate": None})
    plt.close(fig)
    print(f"[written] figures/{name}")

    rec = {"schema": "r36_regime_map_v1",
           "note": ("median rho is E_fixed/E_policy over the whole population; "
                    "win-loss counts are resolved comparisons only. Rows are "
                    "separate populations and are never pooled, with one "
                    "declared exception: the RS rows are apolune levels of the "
                    "single paired ladder population, each pooling that "
                    "population's two identity blocks, which its registration "
                    "licenses because the blocks are the same design at the "
                    "same levels. No level is ever pooled with another."),
           "ladder_rows": {
               "levels_km": [300, 600, 1200, 2400],
               "blocks_pooled": sorted(R50_BLOCKS),
               "amended_budgets": ["1.25", "1.50"],
               "amendment_note": ("for the RS rows only, 1.25 and 1.50 are an "
                                  "amendment declared before any ladder of "
                                  "that population had run; for every other "
                                  "row they are part of the frozen grid"),
               "source": "per-orbit rows, aggregated by level; no summary "
                         "block exists per level",
               "ceiling_free_budgets": sorted(CEILING_FREE_BETAS),
               "ceiling_free_sources": sorted(UNCAPPED.values()),
               "capped_cells_marked_c": ("the amended budgets, which have no "
                                         "ceiling-free re-propagation on "
                                         "either block and are drawn from the "
                                         "capped parents")},
           "plotted_budgets": BETAS,
           "omitted_from_plot": {
               "2.00": ("design B only in panel (a), 5-10 of 15 resolved; "
                        "design A was never propagated there, so no column "
                        "is drawn; the row is in "
                        "r14_trajectory_B_beta_2.00.json and the main "
                        "text's trajectory table"),
               "3.00": ("two cells in panel (a): design A resolves one "
                        "comparison (0-1), design B two (1-1); records in "
                        "r14_trajectory_A_beta_3.00.json and "
                        "r14_trajectory_B_beta_3.00.json")},
           "min_resolved_for_colour": MIN_RESOLVED,
           "hatch_above_undecided_fraction": HATCH_ABOVE,
           "post_hoc_budgets": sorted(POSTHOC),
           "panels": {"a_radial_vs_fixed": {f"{k}|{v}": val
                                            for (k, v), val in a.items()},
                      "b_interior_vs_workmatched": {f"{k}|{v}": val
                                                    for (k, v), val in b.items()}}}
    if ANNOTATE:
        (METRICS / "r36_regime_map.json").write_text(
            json.dumps(rec, indent=2), encoding="utf-8")
        print("[written] metrics/r36_regime_map.json")
    return 0


if __name__ == "__main__":
    import sys
    # --plain draws the main-text rendering: panel (a) alone, without the
    # per-cell counts. The two renderings come from two invocations, so they can be run
    # against different states of the disk and once were: on 14 August the
    # main-text figure was rebuilt while a cell the supplement's figure had
    # correctly left blank was half-written, and for four hours the printed
    # paper showed a cell its own text said did not exist. The record is
    # written by the annotated run only, which means it agrees with the
    # supplement's figure and says nothing about the main text's. Both are on
    # the integrity checker's derivative list for that reason; run them
    # together or the check is what catches the drift.
    if "--plain" in sys.argv[1:]:
        ANNOTATE = False
    raise SystemExit(main())
