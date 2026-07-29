"""Round-2 figures: potential-blend energy drift and GGGRX transfer."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import paper_style
from paper_style import C1, C2, C3, C4

BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
FIGS = BASE / "figures"
FIGS.mkdir(exist_ok=True)

paper_style.apply()


def _load(name):
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


# ---- Fig: energy drift of the four blending policies
d = _load("r2_potential_blend.json")
rev = np.array(d["t_rev_series"])
res = d["results"]
style = {"fixed": (C1, "-", r"fixed $N=120$"),
         "switch": (C4, "-.", "discrete switch"),
         "blend_accel": (C2, "-", "acceleration blend"),
         "blend_pot": (C3, "--", "potential blend (corrected)")}
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(5.2, 4.6), sharex=True,
                              gridspec_kw={"height_ratios": [1.15, 1.0]})
for name, (color, ls, lbl) in style.items():
    e = np.array(res[name]["energy_rel_series"])
    ax.plot(rev, e, ls, color=color, lw=1.35, label=lbl)
    ax2.plot(rev, np.abs(e) + 1e-18, ls, color=color, lw=1.25)
ax.axhline(0.0, color="0.25", lw=0.6)
ax.set_ylabel(r"Signed $\Delta E/E_0$")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax2.set_yscale("log")
ax2.set_ylim(1e-12, 5e-5)
ax2.set_xlabel("Revolutions")
ax2.set_ylabel(r"$|\Delta E/E_0|$")
fig.savefig(FIGS / "fig_blend_energy.pdf")
plt.close(fig)

# ---- Fig: GGGRX vs JGGRX spectrum overlay
g = _load("r2_gggrx_transfer.json")
j = _load("r1_spectrum_pfit.json")
ng = np.array(g["spectrum_arrays"]["n"]); sg = np.array(g["spectrum_arrays"]["sigma_coeff_rms"])
nj = np.array(j["spectrum_arrays"]["n"]); sj = np.array(j["spectrum_arrays"]["sigma_coeff_rms"])
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(5.0, 4.4), sharex=True,
                              gridspec_kw={"height_ratios": [1.6, 1.0]})
ax.loglog(nj, sj, "-", color=C1, lw=1.1, label=r"JGGRX\_1800F (JPL)")
ax.loglog(ng, sg, "--", color=C2, lw=0.9,
          dashes=(3.5, 2.5), label=r"GGGRX\_1200L (GSFC)")
ax.set_ylabel(r"Per-coefficient RMS $\sigma_n$")
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.text(0.97, 0.95, r"both: $\hat p_{\mathrm{spec}}=2.13$",
        transform=ax.transAxes, ha="right", va="top", fontsize=8)
common_n, ij, ig = np.intersect1d(nj.astype(int), ng.astype(int),
                                  assume_unique=True, return_indices=True)
valid = common_n >= 2
ratio = sg[ig[valid]] / sj[ij[valid]] - 1.0
ax2.semilogx(common_n[valid], 100.0 * ratio, color=C3, lw=0.9)
ax2.axhline(0.0, color="0.25", lw=0.6)
ax2.set_xlabel(r"Degree $n$")
ax2.set_ylabel(r"$100(\sigma_n^{\rm GG}/\sigma_n^{\rm JG}-1)$ [\%]")
fig.savefig(FIGS / "fig_gggrx_spectrum.pdf")
plt.close(fig)

print("round-2 figures written")
