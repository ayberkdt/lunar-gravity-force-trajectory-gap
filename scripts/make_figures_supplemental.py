from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import paper_style as ps


ROOT = Path(__file__).resolve().parents[1]
ps.apply()

data = json.loads((ROOT / "metrics/supplemental_pstar_uncertainty.json").read_text())
p = np.asarray(data["p_grid"])
sse = np.asarray(data["sse"])
lo, hi = data["bootstrap"]["p_quantiles_2p5_50_97p5"][0], data["bootstrap"]["p_quantiles_2p5_50_97p5"][2]

fig, ax = plt.subplots(figsize=(4.4, 2.7))
ax.plot(p, sse / len(data["altitudes_km"]), color=ps.C1)
ax.axvspan(lo, hi, color=ps.C2, alpha=0.20, label="central 95% altitude-resampling range")
ax.axvline(data["best_p"], color=ps.C5, linestyle="--", label=rf"minimum $p^*={data['best_p']:.3f}$")
ax.set_xlabel(r"Tail-budget exponent $p$")
ax.set_ylabel(r"Mean squared degree mismatch")
ax.legend()
fig.tight_layout()
fig.savefig(ROOT / "figures/fig_pstar_objective.pdf")
plt.close(fig)
