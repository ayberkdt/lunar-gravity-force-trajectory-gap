"""Main-text figure for the published-rule benchmark (R12/R13).

Three panels that carry the section's argument without repeating its numbers:

  (a) ordered seven-day error ratios of the Atallah rule against the two fixed
      comparators, both populations, with resolved comparisons marked. The rule
      beats the critical fixed degree almost everywhere; against its own
      work-matched degree the ratios cluster around unity and few resolve.
  (b) why so few resolve: the resolution margin M_res against the error gap,
      with the resolution threshold as a diagonal. Almost every matched-work
      pair sits an order of magnitude inside its own numerical envelope, and the
      third-tolerance retest moves the contested orbits along the arrow.
  (c) the noise-free measurement: truncation acceleration defect of the rule
      against that of its work-matched fixed degree, one point per orbit. Every
      orbit of both populations lies below the diagonal.

Sources: metrics/r12_atallah_campaign{,_designB}.json,
         metrics/r13_resolution_diagnosis.json, metrics/r13_ultratight.json,
         metrics/r13_force_defect.json
Product: figures/fig_atallah_benchmark.pdf
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import paper_style as ps

ROOT = Path(__file__).resolve().parents[1]
ps.apply()

camp = {d: json.loads((ROOT / f"metrics/r12_atallah_campaign{s}.json").read_text())
        for d, s in (("A", ""), ("B", "_designB"))}
diag = json.loads((ROOT / "metrics/r13_resolution_diagnosis.json").read_text())
ultra = json.loads((ROOT / "metrics/r13_ultratight.json").read_text())
defect = json.loads((ROOT / "metrics/r13_force_defect.json").read_text())

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))

# ---------------------------------------------------------------- (a) ratios
ax = axes[0]
for d, marker in (("A", "o"), ("B", "^")):
    rows = camp[d]["rows"]
    for key, color, label in (
            ("atallah_vs_fixed_critical", ps.C1, "vs critical fixed"),
            ("atallah_vs_fixed_work_atallah", ps.C2, "vs work-matched fixed")):
        pairs = sorted(((x["comparisons"][key]["rho"],
                         x["comparisons"][key]["resolved"]) for x in rows
                        if x["comparisons"][key]["rho"]), key=lambda z: z[0])
        r = [p_[0] for p_ in pairs]
        res = [p_[1] for p_ in pairs]
        x = np.linspace(0, 100, len(r))
        ax.plot([xx for xx, rr in zip(x, res) if not rr],
                [yy for yy, rr in zip(r, res) if not rr], marker,
                color=color, markersize=2.0, markerfacecolor="none",
                markeredgewidth=0.4, linestyle="none",
                label=None)
        ax.plot([xx for xx, rr in zip(x, res) if rr],
                [yy for yy, rr in zip(r, res) if rr], marker,
                color=color, markersize=2.4, markeredgewidth=0,
                linestyle="none",
                label=(label if d == "A" else None))
ax.axhline(1.0, color=ps.GRAY, lw=0.8, ls="--")
ax.set_yscale("log")
ax.set_xlabel(r"orbits, ordered [\%]")
ax.set_ylabel(r"$\rho$ (comparator / rule)")
ax.set_title("(a) seven-day error ratio", fontsize=8.5)
ax.legend(loc="lower right", frameon=False, handletextpad=0.3,
          labelspacing=0.25, borderpad=0.2)

# ------------------------------------------------------- (b) resolution margin
ax = axes[1]
for d, marker in (("A", "o"), ("B", "^")):
    rows = diag["designs"][d]["rows"]
    ax.plot([r["gap_m"] for r in rows], [r["m_res"] for r in rows], marker,
            color=ps.C2, markersize=2.6, alpha=0.6, markeredgewidth=0,
            label=f"design {d}")
before = [r["previous"]["m_res"] for r in ultra["rows"]]
after = [r["comparison"]["m_res"] for r in ultra["rows"]]
gaps = [r["comparison"]["absolute_error_difference_m"] for r in ultra["rows"]]
for g, b, a in zip(gaps, before, after):
    ax.annotate("", xy=(g, a), xytext=(g, b),
                arrowprops=dict(arrowstyle="-|>", lw=0.45, color=ps.C4,
                                alpha=0.75, mutation_scale=5,
                                shrinkA=0, shrinkB=0))
ax.axhline(1.0, color=ps.GRAY, lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"error gap $|E_{\mathrm{At}}-E_{\mathrm{fix}}|$ [m]")
ax.set_ylabel(r"$M_{\mathrm{res}}$")
ax.set_title("(b) matched-work resolution", fontsize=8.5)
ax.legend(loc="lower right", frameon=False, handletextpad=0.4)

# ------------------------------------------------------------- (c) force defect
ax = axes[2]
lo, hi = 1e-12, 1e-6
for d, marker in (("A", "o"), ("B", "^")):
    rows = defect["designs"][d]["rows"]
    ax.plot([r["fixed_work"]["defect_rms_m_s2"] for r in rows],
            [r["atallah"]["defect_rms_m_s2"] for r in rows], marker,
            color=ps.C3, markersize=2.6, alpha=0.6, markeredgewidth=0,
            label=f"design {d}")
ax.plot([lo, hi], [lo, hi], color=ps.GRAY, lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel(r"fixed degree $|\Delta a|$ [m\,s$^{-2}$]")
ax.set_ylabel(r"rule $|\Delta a|$ [m\,s$^{-2}$]")
ax.set_title("(c) truncation force defect", fontsize=8.5)
ax.legend(loc="upper left", frameon=False, handletextpad=0.4)

fig.tight_layout(pad=0.4)
out = ROOT / "figures" / "fig_atallah_benchmark.pdf"
fig.savefig(out)
print(f"[written] {out.relative_to(ROOT)}")
