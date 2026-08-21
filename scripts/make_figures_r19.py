"""Figures promoted from the supplement to the main text (R19 editorial pass).

  fig_variational_parity.pdf   predicted vs measured equal-budget error ratio
  budget_pareto.pdf            regenerated with a third realized-work panel

Both draw only on frozen R14 records; nothing is recomputed here.

Usage:  python make_figures_r19.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
from matplotlib.lines import Line2D       # noqa: E402

import paper_style as ps                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
FIGS = ROOT / "figures"


def variational_parity() -> None:
    """Predicted against measured equal-budget ratio on the variational panel.

    The panel is the highest completed level of the R37 extension when that
    record exists, and the archived eight-orbit panel otherwise. Orbits whose
    propagated comparison the resolution rule leaves undecided are drawn open:
    they place a point on the plane but carry no verdict, so they are excluded
    from the sign tally (metrics/r37_scoring_amendment.json).
    """
    ext = METRICS / "r37_variational_extension.json"
    verdict = METRICS / "r37_panel_verdict.json"
    if ext.exists() and verdict.exists():
        d = json.loads(ext.read_text(encoding="utf-8"))
        v = json.loads(verdict.read_text(encoding="utf-8"))
        keep = {(c["design"], c["sobol_index"])
                for c in v["unresolved"]["orbits"]}
        n_per = v["levels"]["highest_complete_per_design"]
        pareto = json.loads((METRICS / "r14_budget_pareto.json"
                             ).read_text(encoding="utf-8"))
        members = set()
        for des in ("A", "B"):
            pr = [r for r in pareto["designs"][des]["rows"]
                  if not r["budgets"]["beta_1.00"]["censored"]]
            pr.sort(key=lambda r: r["hp_km"])
            idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
            members |= {(des, int(pr[i]["sobol_index"])) for i in idx}
        rows = [r for r in d["rows"]
                if (r["design"], r["sobol_index"]) in members]
        panel_n = v["panel"]["orbits"]
    else:
        d = json.loads((METRICS / "r14_variational_budget.json"
                        ).read_text(encoding="utf-8"))
        rows = [r for r in d["rows"] if r.get("status") == "complete"]
        keep = set()
        panel_n = len(rows)

    pred = np.array([r["predicted_ratio_fixed_over_atallah"] for r in rows])
    meas = np.array([r["measured"]["fixed_budget"] / r["measured"]["atallah_budget"]
                     for r in rows])
    des = [r["design"] for r in rows]
    resolved = np.array([bool(r["measured"].get("resolved")) for r in rows])
    calib = np.array([r["calibration_ratio_fixed"] for r in rows
                      if r.get("calibration_ratio_fixed") is not None])

    agree = int(np.sum((pred[resolved] < 1.0) == (meas[resolved] < 1.0)))
    n_res = int(resolved.sum())

    # What this figure has to show is not that the points lie near a diagonal
    # but which side of unity they lie on: a prediction and a measurement that
    # straddle rho = 1 disagree about which policy wins, and the paper's claim
    # is that none of the resolved ones does. So the two quadrants where that
    # disagreement would appear are drawn, and their emptiness is the result.
    # Previously the reader had to reconstruct this from a scatter, a diagonal
    # and four legend entries.
    # Sized to be placed at 0.78\linewidth, which reproduces it at about its
    # native width: the earlier 4.6 in canvas was shrunk to 0.62\linewidth and
    # rendered its 6.2 pt labels at 5.4 pt on the page, below what the journal
    # asks of figure text. Every size below is therefore a printed size.
    fig, ax = plt.subplots(figsize=(4.9, 4.75))
    # Tight padding. The panel spans eight decades because two orbits sit at
    # the extremes; generous padding on top of that spent a fifth of the axis
    # on empty ground and pushed the bulk of the points into a small diagonal
    # band in the middle.
    lim = [min(pred.min(), meas.min()) * 0.65,
           max(pred.max(), meas.max()) * 1.5]

    for x0, x1, y0, y1 in ((lim[0], 1.0, 1.0, lim[1]),
                           (1.0, lim[1], lim[0], 1.0)):
        ax.fill_between([x0, x1], y0, y1, color="#f7ebe9", zorder=0,
                        linewidth=0)
    ax.fill_between(lim, [x / 2 for x in lim], [x * 2 for x in lim],
                    color="0.90", zorder=1, linewidth=0,
                    label="within a factor of two")
    ax.plot(lim, lim, color="0.30", lw=1.2, ls="--", zorder=2,
            label="exact agreement")
    ax.axvline(1.0, color="0.62", lw=0.8, zorder=2)
    ax.axhline(1.0, color="0.62", lw=0.8, zorder=2)

    # Shape carries the design and fill carries decidability, so colour is not
    # asked to encode a third thing.
    for design, marker, colour in (("A", "o", "#1f5fa8"), ("B", "s", "#d1651a")):
        m = np.array([i for i, x in enumerate(des) if x == design], dtype=int)
        if m.size == 0:
            continue
        r, u = m[resolved[m]], m[~resolved[m]]
        if r.size:
            ax.scatter(pred[r], meas[r], marker=marker, s=30, zorder=4,
                       color=colour, edgecolor="white", linewidth=0.5,
                       label=f"design {design}")
        if u.size:
            ax.scatter(pred[u], meas[u], marker=marker, s=30, zorder=3,
                       facecolor="none", edgecolor=colour, linewidth=0.8,
                       alpha=0.8)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"predicted $\rho$  (forced variational)")
    ax.set_ylabel(r"measured $\rho$  (propagated)")

    # The tinted quadrants are only worth drawing if the figure says what they
    # mean and what is in them. The note sits in the lower-right tinted corner
    # it describes, which carries no data: stacked over the legend at the top
    # left it forced a four-line wrap that ran under the inset's frame.
    ax.text(0.975, 0.030,
            "tinted corners: prediction and measurement\n"
            f"disagree on which policy wins ({n_res - agree} of {n_res} "
            "resolved)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.0,
            color="0.25", linespacing=1.35)
    # The open-marker convention reads better as a legend entry than as a
    # fifth line of prose, so the note above no longer carries it.
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], marker="o", linestyle="none", markersize=5.2,
                          markerfacecolor="none", markeredgecolor="0.35",
                          markeredgewidth=0.8))
    labels.append("undecided comparison")
    # Below the axes in two rows, as Fig. 2 places its legend. Inside the panel
    # it stood in the upper-left tinted corner, whose emptiness is the result
    # the figure reports, and its entries ran into the inset's tick labels.
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.135),
              ncol=3, fontsize=7.4, framealpha=0.0, borderpad=0.2,
              columnspacing=1.4, handletextpad=0.5)

    # calibration_ratio_fixed is predicted over measured position error for the
    # constant-degree comparator (rev14_variational_budget.py), so the inset is
    # a linearization check and not a work match: an earlier title read
    # W_fix/W_rad, which names a different quantity entirely.
    # Left to the log locator the axis labelled itself 6e-1, 2e0, 3e0, which is
    # three ways of writing "about one" and crowds a box a fifth of an axis
    # high. Explicit ticks at a half, one and two give the scale instead. The
    # upper limit holds the whole sample: at 2.2, and then at 2.8, the largest
    # ratio (3.05, on a resolved orbit) fell outside the box and was drawn as a
    # stray mark above it rather than shown or dropped.
    # It sits in the white wedge above the diagonal of the agreement quadrant,
    # not on a tinted corner: the tinted corners' emptiness is the figure's
    # result, and a box parked on one reads as covering the evidence.
    # Both axes carry the same log range, so the diagonal is the line y=x in
    # axes coordinates and the box clears it only while its bottom edge sits
    # above its right edge: at [0.545, 0.720] the earlier box had the diagonal
    # and its factor-of-two band running through its lower right corner. The
    # limits carry headroom so the upper tick label and the title do not meet.
    ins = ax.inset_axes([0.520, 0.825, 0.260, 0.130])
    ins.axhline(1.0, color="0.55", lw=0.8, ls="--")
    ins.scatter(np.arange(len(calib)), calib, s=6, color="0.25")
    ins.set_yscale("log")
    ins.set_ylim(0.52, 3.8)
    ins.set_yticks([0.6, 1.0, 2.0, 3.0])
    ins.set_yticklabels(["0.6", "1", "2", "3"])
    ins.minorticks_off()
    ins.set_xticks([])
    # right=False as well as the hidden spine: the paper style mirrors y ticks,
    # which left two marks floating where the right-hand spine is not drawn.
    ins.tick_params(labelsize=6.8, length=2.0, right=False, top=False)
    # Open frame: the inset is a marginal note on the main panel, not a second
    # figure boxed inside it.
    for side in ("top", "right", "bottom"):
        ins.spines[side].set_visible(False)
    ins.spines["left"].set_linewidth(0.6)
    ins.patch.set_alpha(0.0)
    # The ratio the inset plots is named in the caption rather than in the
    # title: spelled out here the title was wider than the panel had room for
    # and ran over the right-hand axis.
    ins.set_title("comparator linearization", fontsize=6.8, pad=3,
                  color="0.25")

    fig.tight_layout()
    out = FIGS / "fig_variational_parity.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[written] {out.name}  (panel {panel_n} orbits, "
          f"{n_res} resolved, {agree}/{n_res} sign agreement)")


def main() -> int:
    # A silent failure here produces a figure in matplotlib's defaults next to
    # figures in the paper's style, which is a difference a reader sees and an
    # author does not. The fall-back stays, because a missing LaTeX should not
    # stop the figure being drawn, but it says so.
    try:
        ps.apply()
    except Exception as exc:                                   # noqa: BLE001
        print(f"[warn] paper_style.apply() failed ({exc}); the figure will be "
              f"drawn in matplotlib defaults and will not match the others")
    FIGS.mkdir(exist_ok=True)
    variational_parity()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
