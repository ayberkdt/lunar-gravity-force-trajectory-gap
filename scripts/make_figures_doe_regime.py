"""Regime map for the stratified 24-orbit designed geometry matrix.

Three panels over the (perilune, apolune) design plane: the equal-work ratio
rho_work, the critical-altitude ratio rho_crit (both best-of-three-schedules,
values > 1 favor scheduling), and the gravity-kernel time saving of the best
schedule.  Marker size encodes inclination.  Source: Stage-2 7-day matrix.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

import paper_style as ps

ROOT = Path(__file__).resolve().parents[1]
ps.apply()

data = json.loads((ROOT / "metrics/r7_doe_matrix_stage2.json").read_text())
rows = data["rows"]
scheds = ("sched_emp", "sched_down", "sched_up")

hp = np.array([r["hp_km"] for r in rows])
ha = np.array([r["ha_km"] for r in rows])
inc = np.array([r["incl_deg"] for r in rows])
best_work = np.array([max(r["ratios"][s]["rho_work"] for s in scheds) for r in rows])
best_crit = np.array([max(r["ratios"][s]["rho_crit"] for s in scheds) for r in rows])
saving = np.array([
    100.0 * r["ratios"][max(scheds, key=lambda s: r["ratios"][s]["rho_work"])]
    ["grav_time_saving_vs_crit"] for r in rows])

# Inclination -> marker size (28..120 pt^2)
size = 28.0 + (inc - inc.min()) / max(inc.max() - inc.min(), 1.0) * 92.0

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), constrained_layout=True)

for ax, values, title, is_ratio in (
    (axes[0], best_work, r"$\rho_{\mathrm{work}}$", True),
    (axes[1], best_crit, r"$\rho_{\mathrm{crit}}$", True),
    (axes[2], saving, r"grav.\ time saving [\%]", False),
):
    if is_ratio:
        lv = np.log10(values)
        m = float(np.max(np.abs(lv))) * 1.02
        norm = TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m)
        sc = ax.scatter(hp, ha, c=lv, s=size, cmap="RdBu", norm=norm,
                        edgecolors="0.25", linewidths=0.4)
        cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label(r"$\log_{10}$ " + title + r" ($>0$: sched.\ wins)")
    else:
        sc = ax.scatter(hp, ha, c=values, s=size, cmap="viridis",
                        edgecolors="0.25", linewidths=0.4)
        cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label(title)
    ax.set_xlabel(r"Perilune altitude [km]")
    ax.set_ylabel(r"Apolune altitude [km]")

# Inclination size legend on the first panel
for iv in (10, 50, 90):
    s = 28.0 + (iv - inc.min()) / max(inc.max() - inc.min(), 1.0) * 92.0
    axes[0].scatter([], [], s=s, c="0.6", edgecolors="0.25", linewidths=0.4,
                    label=rf"$i={iv}^\circ$")
axes[0].legend(title="inclination", loc="upper right", labelspacing=0.8,
               borderpad=0.5, fontsize=6.5, title_fontsize=7)

fig.savefig(ROOT / "figures/fig_doe_regime.pdf")
plt.close(fig)
print("[written] figures/fig_doe_regime.pdf")
