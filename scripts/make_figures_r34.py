"""R34 figures: the instrument ladder and the paired per-orbit record.

Two figures, both built entirely from frozen R14 records (no propagation):

  fig_instrument_ladder.pdf
      One column per instrument, one point per orbit, showing how the
      fixed-over-radial ratio moves as each instrument carries more of the
      dynamics.  Replaces the prose ladder of counts.

  fig_beta1_per_orbit.pdf
      The paired seven-day comparison at beta = 1, one bar per orbit, sorted,
      with each orbit's own summed resolution envelope drawn as the band that
      decides it.  Replaces the walls of tallies in the main text.
"""

from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import paper_style
import matplotlib.pyplot as plt

paper_style.apply()

ROOT = pathlib.Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"

BETA_KEY = "beta_1.00"
C_FIXED = paper_style.C1
C_RADIAL = paper_style.C5
GRAY = paper_style.GRAY


def load_force(design):
    """Per-orbit fixed-over-radial ratios for the two force-level instruments."""
    pareto = json.loads((METRICS / "r14_budget_pareto.json").read_text())
    defect, proxy = [], []
    for row in pareto["designs"][design]["rows"]:
        entry = row.get("budgets", {}).get(BETA_KEY)
        if entry is None or entry.get("censored"):
            continue
        rad, fix = entry["atallah"]["defect"], entry["fixed"]["defect"]
        defect.append(fix["defect_rms_m_s2"] / rad["defect_rms_m_s2"])
        proxy.append(
            abs(fix["in_track_displacement_proxy_m"])
            / abs(rad["in_track_displacement_proxy_m"])
        )
    return np.array(defect), np.array(proxy)


def load_variational():
    """Predicted ratios on the forced-variational panel.

    The panel is the highest completed level of the R37 extension when that
    record exists, and the archived eight-orbit panel otherwise, so this column
    and the corresponding row of the instrument-ladder table cannot disagree.
    """
    ext = METRICS / "r37_variational_extension.json"
    verdict = METRICS / "r37_panel_verdict.json"
    if ext.exists() and verdict.exists():
        v = json.loads(verdict.read_text())
        n_per = v["levels"]["highest_complete_per_design"]
        pareto = json.loads((METRICS / "r14_budget_pareto.json").read_text())
        members = set()
        for des in ("A", "B"):
            pr = [r for r in pareto["designs"][des]["rows"]
                  if not r["budgets"]["beta_1.00"]["censored"]]
            pr.sort(key=lambda r: r["hp_km"])
            idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
            members |= {(des, int(pr[i]["sobol_index"])) for i in idx}
        rec = json.loads(ext.read_text())
        return np.array([r["predicted_ratio_fixed_over_atallah"]
                         for r in rec["rows"]
                         if (r["design"], r["sobol_index"]) in members])
    rec = json.loads((METRICS / "r14_variational_budget.json").read_text())
    return np.array(
        [r["predicted_ratio_fixed_over_atallah"] for r in rec["rows"]]
    )


def load_propagated(design):
    """Per-orbit errors and the summed envelope that decides each comparison.

    Scored at the ``tight`` level with the tight-to-tighter self-difference as
    the envelope, which is the convention of the R14 budget campaign: this
    reproduces its reported verdict counts and median ratio exactly.
    """
    traj = json.loads(
        (METRICS / f"r14_trajectory_{design}_beta_1.00.json").read_text()
    )
    e_fix, e_rad, env = [], [], []
    for row in traj["rows"]:
        rad = row["policies"]["atallah_budget"]
        fix = row["policies"]["fixed_budget"]
        e_rad.append(rad["error_tight"]["pos_rms_m"])
        e_fix.append(fix["error_tight"]["pos_rms_m"])
        env.append(
            rad["truth_inclusive_envelope_m"] + fix["truth_inclusive_envelope_m"]
        )
    return np.array(e_fix), np.array(e_rad), np.array(env)


def _strip(ax, x, values, color, rng):
    """A jittered one-dimensional scatter with its median marked."""
    jitter = rng.uniform(-0.13, 0.13, size=values.size)
    ax.scatter(
        x + jitter, values, s=7, color=color, alpha=0.55,
        edgecolors="none", zorder=3,
    )
    med = float(np.median(values))
    ax.plot(
        [x - 0.28, x + 0.28], [med, med],
        color="black", lw=1.4, zorder=5, solid_capstyle="butt",
    )
    return med


def figure_ladder():
    rng = np.random.default_rng(20260803)
    fig, ax = plt.subplots(figsize=(5.4, 3.2))

    defect_a, proxy_a = load_force("A")
    defect_b, proxy_b = load_force("B")
    var = load_variational()
    ratios = []
    for design in ("A", "B"):
        e_fix, e_rad, _ = load_propagated(design)
        ratios.append(e_fix / e_rad)
    prop = np.concatenate(ratios)

    columns = [
        (np.concatenate([defect_a, defect_b]), "Defect RMS", C_RADIAL),
        (np.concatenate([proxy_a, proxy_b]), r"Proxy $D_I$", C_RADIAL),
        (var, "Variational", C_FIXED),
        (prop, "Propagated", C_FIXED),
    ]

    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1e8)
    ax.axhspan(1e-5, 1.0, color=C_FIXED, alpha=0.05, zorder=0)
    ax.axhline(1.0, color=GRAY, lw=0.8, ls="--", zorder=2)

    labels = []
    for i, (values, label, color) in enumerate(columns):
        med = _strip(ax, i, values, color, rng)
        n_radial = int((values > 1.0).sum())
        labels.append(label)
        ax.annotate(
            rf"{n_radial}/{values.size}",
            xy=(i, 1.02), xycoords=("data", "axes fraction"),
            ha="center", va="bottom", fontsize=8, color=GRAY,
            annotation_clip=False,
        )
        ax.annotate(
            rf"${med:.3g}$",
            xy=(i + 0.30, med), ha="left", va="center",
            fontsize=7.5, color="black",
        )

    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.55, len(columns) - 0.30)
    ax.set_ylabel(r"fixed-over-radial ratio")
    for y, txt, va in ((3e6, "radial better", "top"),
                       (3e-5, "constant degree better", "bottom")):
        ax.text(
            len(columns) - 0.42, y, txt, fontsize=7.5, color=GRAY,
            ha="right", va=va, style="italic",
        )
    ax.set_title(
        r"orbits favouring radial, of those each instrument scores",
        fontsize=8, color=GRAY, pad=13,
    )
    fig.savefig(FIGURES / "fig_instrument_ladder.pdf")
    plt.close(fig)


def figure_per_orbit():
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)

    for ax, design in zip(axes, ("A", "B")):
        e_fix, e_rad, env = load_propagated(design)
        ratio = e_fix / e_rad
        # The envelope, expressed on the same ratio axis: the factor by which
        # the two errors must differ before the comparison is decidable.
        gap = np.abs(e_fix - e_rad)
        resolved = gap > env
        order = np.argsort(ratio)

        x = np.arange(ratio.size)
        r = ratio[order]
        res = resolved[order]

        ax.axhline(1.0, color=GRAY, lw=0.8, ls="--", zorder=2)
        ax.vlines(
            x[~res], 1.0, r[~res],
            color=GRAY, lw=1.6, alpha=0.45, zorder=3,
        )
        favours_fixed = r < 1.0
        for mask, color in (
            (res & favours_fixed, C_FIXED),
            (res & ~favours_fixed, C_RADIAL),
        ):
            ax.vlines(x[mask], 1.0, r[mask], color=color, lw=1.6, zorder=4)

        ax.set_yscale("log")
        ax.set_xlabel(f"design~{design} orbits, sorted")
        ax.set_xlim(-1.5, ratio.size + 0.5)
        n_fix = int((res & favours_fixed).sum())
        n_rad = int((res & ~favours_fixed).sum())
        ax.set_title(
            rf"{n_fix} constant \textbar\ {n_rad} radial \textbar\ "
            rf"{int((~res).sum())} undecided",
            fontsize=8, color=GRAY,
        )

    axes[0].set_ylabel(r"$\rho_{\mathrm{budget}}=E_{\mathrm{fixed}}/E_{\mathrm{radial}}$")
    handles = [
        plt.Line2D([], [], color=C_FIXED, lw=2, label="resolved, constant better"),
        plt.Line2D([], [], color=C_RADIAL, lw=2, label="resolved, radial better"),
        plt.Line2D([], [], color=GRAY, lw=2, alpha=0.45, label="undecided"),
    ]
    axes[1].legend(handles=handles, loc="upper left", fontsize=7.2)
    fig.savefig(FIGURES / "fig_beta1_per_orbit.pdf")
    plt.close(fig)


if __name__ == "__main__":
    figure_ladder()
    figure_per_orbit()
    print("wrote fig_instrument_ladder.pdf and fig_beta1_per_orbit.pdf")
