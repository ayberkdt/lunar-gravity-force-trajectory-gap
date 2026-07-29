"""The O25 budget figure: force-defect curves and propagated operating points.

Two rows, one column per independent population:
  top    -- integration-noise-free truncation force defect against budget,
  bottom -- propagated seven-day trajectory RMS error against budget.

Population medians with 10-90 percentile bands. Budgets at which any orbit is
censored (the comparator would need degrees above the adopted truth degree) are
drawn with open markers and joined by a dashed segment, because their medians
are taken over a smaller subset. The accuracy-target operating point of the
archived benchmark is drawn at its own measured budget, and the interval in
which the two force curves cross is shaded. The trajectory panels are
operating points at the budgets that were propagated, not a continuous frontier.

Usage:
    python make_figures_r14.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import paper_style as ps

ps.apply()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, NullLocator

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"
GRID = [0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00]

AT_C, FX_C, OR_C = ps.C1, ps.C5, ps.C3


def key_of(b):
    return f"beta_{b:.2f}"


def band(ax, x, med, lo, hi, color, label, censored):
    x = np.asarray(x, float)
    med = np.asarray(med, float)
    ax.fill_between(x, lo, hi, color=color, alpha=0.13, lw=0)
    cen = np.asarray(censored, bool)
    # solid where every orbit contributes, dashed once any is censored
    for i in range(len(x) - 1):
        style = "--" if (cen[i] or cen[i + 1]) else "-"
        ax.plot(x[i:i + 2], med[i:i + 2], style, color=color, lw=1.3, zorder=3)
    ax.plot(x[~cen], med[~cen], "o", color=color, ms=4.0, zorder=4, label=label)
    if cen.any():
        ax.plot(x[cen], med[cen], "o", mfc="white", mec=color, mew=1.1,
                ms=4.0, zorder=4)


def crossing_interval(betas, ratios):
    """Bracketing budgets where the median defect ratio crosses unity."""
    for i in range(len(betas) - 1):
        if (ratios[i] - 1.0) * (ratios[i + 1] - 1.0) < 0.0:
            return betas[i], betas[i + 1]
    return None


def main() -> int:
    pareto = json.loads((METRICS / "r14_budget_pareto.json").read_text())
    oracle = None
    if (METRICS / "r14_oracle.json").exists():
        oracle = json.loads((METRICS / "r14_oracle.json").read_text())
    traj = {}
    for p in sorted(METRICS.glob("r14_trajectory_*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        traj.setdefault(d["design"], {})[float(d["beta"])] = d

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), sharex=True)
    for col, des in enumerate(("A", "B")):
        s = pareto["designs"][des]["summary"]
        keys = [key_of(b) for b in GRID]
        live = [k for k in keys if s[k].get("orbits")]
        bs = [s[k]["beta_requested"] for k in live]
        cen = [s[k]["censored"] > 0 for k in live]

        # ---------------------------------------------------------- force panel
        ax = axes[0, col]
        for who, color, lab in (("atallah", AT_C, "budget-calibrated radial"),
                                ("fixed", FX_C, "constant degree")):
            st = [s[k][f"{who}_defect_rms_m_s2"] for k in live]
            band(ax, bs, [q["median"] for q in st], [q["p10"] for q in st],
                 [q["p90"] for q in st], color, lab, cen)
        if oracle and des in oracle["designs"]:
            os_ = oracle["designs"][des]["summary"]
            ok = [k for k in keys if os_.get(k, {}).get("orbits")]
            if ok:
                ax.plot([os_[k]["beta_requested"] for k in ok],
                        [os_[k]["oracle_defect_rms_m_s2"]["median"] for k in ok],
                        ":", color=OR_C, lw=1.3, zorder=3)
                ax.plot([os_[k]["beta_requested"] for k in ok],
                        [os_[k]["oracle_defect_rms_m_s2"]["median"] for k in ok],
                        "s", color=OR_C, ms=3.4, zorder=4,
                        label="allocation benchmark (O26)")
        o = s["original"]
        bo = pareto["designs"][des]["beta_original"]["median"]
        ax.plot([bo], [o["atallah_defect_rms_m_s2"]["median"]], "*",
                color=AT_C, ms=10.0, mec="0.2", mew=0.5, zorder=6,
                label="accuracy-target point")
        ratios = [s[k]["R_a_defect_rms"]["median"] for k in live]
        cx = crossing_interval(bs, ratios)
        if cx:
            ax.axvspan(cx[0], cx[1], color=ps.GRAY, alpha=0.16, lw=0, zorder=0)
        ax.axvline(1.0, color="0.35", lw=0.7, ls=(0, (4, 3)), zorder=1)
        ax.set_yscale("log")
        # the 10th-percentile band collapses many decades at the top budgets;
        # clip to the range the medians actually occupy plus two decades
        meds = [s[k][f"{w}_defect_rms_m_s2"]["median"]
                for k in live for w in ("atallah", "fixed")]
        ax.set_ylim(min(meds) / 30.0, max(meds) * 12.0)
        ax.set_title(rf"Design {des}", fontsize=9.0)
        if col == 0:
            ax.set_ylabel(r"truncation defect $|\Delta a|_{\mathrm{RMS}}$ [m\,s$^{-2}$]")

        # ----------------------------------------------------- trajectory panel
        ax = axes[1, col]
        if des in traj and traj[des]:
            tb = sorted(traj[des])
            for who, color, lab in (("atallah_error_m", AT_C, "budget-calibrated radial"),
                                    ("fixed_error_m", FX_C, "constant degree")):
                st = [traj[des][b]["summary"][who] for b in tb]
                tc = [len(traj[des][b]["censored"]) > 0 for b in tb]
                band(ax, tb, [q["median"] for q in st], [q["p10"] for q in st],
                     [q["p90"] for q in st], color, lab, tc)
            # keep the axis on the range where both policies are visible, and
            # annotate rather than silently clip whatever falls outside it
            vis = [traj[des][b]["summary"][w]["median"]
                   for b in tb if b >= 0.7
                   for w in ("atallah_error_m", "fixed_error_m")]
            if vis:
                lo, hi = min(vis) / 25.0, max(vis) * 25.0
                ax.set_ylim(lo, hi)
                for b in tb:
                    for w, cc in (("atallah_error_m", AT_C),
                                  ("fixed_error_m", FX_C)):
                        v = traj[des][b]["summary"][w]["median"]
                        if v > hi:
                            ax.annotate(f"{v:.0f}\\,m", xy=(b, hi),
                                        xytext=(0, -3), textcoords="offset points",
                                        ha="center", va="top", fontsize=6.0,
                                        color=cc,
                                        arrowprops=dict(arrowstyle="-|>", color=cc,
                                                        lw=0.8, shrinkA=0, shrinkB=0,
                                                        mutation_scale=6),
                                        annotation_clip=False)
            for b in tb:
                sm = traj[des][b]["summary"]
                ax.annotate(f"{sm['unresolved']}/{sm['orbits']}",
                            xy=(b, sm["fixed_error_m"]["median"]),
                            xytext=(0, -10), textcoords="offset points",
                            ha="center", fontsize=6.2, color="0.35")
            if cx:
                ax.axvspan(cx[0], cx[1], color=ps.GRAY, alpha=0.16, lw=0, zorder=0)
        ax.axvline(1.0, color="0.35", lw=0.7, ls=(0, (4, 3)), zorder=1)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xticks(GRID)
        ax.set_xticklabels([f"{b:g}" for b in GRID])
        # a log x-axis otherwise draws its own decade labels through ours
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        axes[0, col].xaxis.set_minor_locator(NullLocator())
        ax.set_xlim(0.44, 3.45)
        axes[0, col].set_xlim(0.44, 3.45)
        ax.set_xlabel(r"gravity-work budget $\beta = W/N_{\mathrm{crit}}^{2}$")
        if col == 0:
            ax.set_ylabel(r"seven-day position RMS [m]")

    handles = [Line2D([], [], color=AT_C, marker="o", ms=4, lw=1.3,
                      label="budget-calibrated radial rule"),
               Line2D([], [], color=FX_C, marker="o", ms=4, lw=1.3,
                      label="constant degree, same nominal per-call budget"),
               Line2D([], [], color=AT_C, marker="*", ms=9, lw=0, mec="0.2",
                      label="archived accuracy-target operating point"),
               Line2D([], [], color="0.5", marker="o", mfc="white", ms=4, ls="--",
                      label="budget with censored orbits")]
    if oracle:
        handles.insert(2, Line2D([], [], color=OR_C, marker="s", ms=3.4, ls=":",
                                 lw=1.3, label="allocation benchmark (O26)"))
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.legend(handles=handles, loc="lower center",
               ncol=2 if len(handles) > 4 else 3,
               bbox_to_anchor=(0.5, 0.0), fontsize=7.6)
    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "budget_pareto.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[written] {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
