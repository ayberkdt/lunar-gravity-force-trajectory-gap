"""Generate all paper figures (PDF) from the JSON outputs in metrics/."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
FIGS = BASE / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 9.5,
    "font.family": "serif",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
})

C1, C2, C3, C4 = "#1b4965", "#bb7b18", "#5c8a3c", "#8a3c5c"


def _load(name: str) -> dict:
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


# ---- Fig: truncation criteria N_min(h)
d = _load("e2_truncation_criteria.json")
rows = d["rows"]
h = [r["altitude_km"] for r in rows]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.plot(h, [r["empirical_nmin_tail1e2"] for r in rows], "o-", color="black",
        lw=1.8, ms=4.5, label="Empirical (JGGRX spectrum)", zorder=5)
ax.plot(h, [r["kaula_p1_7_tail1e2"] for r in rows], "s--", color=C1,
        lw=1.4, ms=4, label=r"Spectrum-weighted, $p=1.7$")
ax.plot(h, [r["kaula_p2_0_tail1e2"] for r in rows], "^--", color=C3,
        lw=1.4, ms=4, label=r"Spectrum-weighted, $p=2.0$ (Kaula)")
ax.plot(h, [r["attenuation_only_floor1e3"] for r in rows], "d--", color=C2,
        lw=1.4, ms=4, label=r"Attenuation-only, $10^{-3}$ floor")
ax.set_xlabel("Altitude $h$ [km]")
ax.set_ylabel(r"Recommended truncation degree $N_{\min}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([20, 30, 50, 80, 100, 150, 200, 300])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_yticks([25, 50, 100, 200, 400, 600])
ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
# Explicit major ticks are not enough on a log axis: matplotlib still labels
# the minor decade subdivisions, and over a narrow range those labels collide
# with the major ones and print as one run of digits.
for _axis in (ax.get_xaxis(), ax.get_yaxis()):
    _axis.set_minor_formatter(mticker.NullFormatter())
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_truncation.pdf")
plt.close(fig)

# ---- Fig: kernel timing
d = _load("e3_kernel_timing.json")
sweep = d["degree_sweep"]
N = np.array([r["degree"] for r in sweep], dtype=float)
t = np.array([r["best_us"] for r in sweep])
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.loglog(N, t, "o-", color=C1, lw=1.6, ms=4, label="Serial Kahan kernel (measured)")
ref = t[-1] * (N / N[-1]) ** 2
ax.loglog(N, ref, ":", color="gray", lw=1.2, label=r"$\mathcal{O}(N^2)$ reference slope")
ax.set_xlabel(r"Truncation degree $N$")
ax.set_ylabel(r"Evaluation time [$\mu$s]")
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_timing.pdf")
plt.close(fig)

# ---- Fig: band shares
d = _load("e1_band_shares_60_100_nmax300.json")
res = d["results"]
h = [r["altitude_km"] for r in res]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.semilogy(h, [r["share_of_pert_total"]["band_2_60"] for r in res], "o-",
            color=C1, lw=1.6, ms=4, label=r"$2 \leq n \leq 60$")
ax.semilogy(h, [r["share_of_pert_total"]["band_61_100"] for r in res], "s-",
            color=C2, lw=1.6, ms=4, label=r"$61 \leq n \leq 100$")
ax.semilogy(h, [r["share_of_pert_total"]["tail_gt_100"] for r in res], "^-",
            color=C3, lw=1.6, ms=4, label=r"$n > 100$")
ax.axhline(1e-2, color="gray", ls=":", lw=1)
ax.text(155, 1.25e-2, "1% threshold", color="gray", fontsize=8)
ax.set_xlabel("Altitude $h$ [km]")
ax.set_ylabel("Band RMS share of perturbation")
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_band_share.pdf")
plt.close(fig)

# ---- Fig: switching jump vs quantization step
d = _load("e5_switch_jump.json")
rows = d["rows"]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
for h_km, color, marker in [(50.0, C1, "o"), (100.0, C2, "s")]:
    sel = [r for r in rows if r["altitude_km"] == h_km]
    q = [r["step"] for r in sel]
    j = [100.0 * r["jump_over_pert"] for r in sel]
    ax.loglog(q, j, marker + "-", color=color, lw=1.6, ms=4.5,
              label=f"$h = {h_km:.0f}$ km")
ax.set_xlabel(r"Quantization step $q$ (degrees)")
ax.set_ylabel("RMS jump / perturbation RMS [%]")
ax.set_xticks([5, 10, 25, 50])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_switch_jump.pdf")
plt.close(fig)

# ---- Fig: orbit-level error vs degree (E6)
d = _load("e6_orbit_mapping.json")
rows = d["rows"]
N = [r["degree"] for r in rows]
e = [r["rms_m"] for r in rows]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.semilogy(N, e, "o-", color=C1, lw=1.7, ms=4.5, zorder=5,
            label="24 h RMS position error")
marks = [(56, C3, r"$p=2.0$ pick"), (70, "black", r"$p=1.7$ pick"),
         (124, C2, "attenuation pick")]
import matplotlib.transforms as mtransforms
blend = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
for x, color, lbl in marks:
    ax.axvline(x, color=color, ls=":", lw=1.1)
    ax.text(x + 2.0, 0.965, lbl, rotation=90, rotation_mode="anchor",
            ha="right", va="bottom", fontsize=7.5, color=color,
            transform=blend)
ax.set_xlabel(r"Truncation degree $N$")
ax.set_ylabel("Position error vs. $N=300$ truth [m]")
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.savefig(FIGS / "fig_orbit_error.pdf")
plt.close(fig)

# ---- Fig: eccentric-orbit cost/accuracy trade (E7)
d = _load("e7_adaptive_orbit.json")
pts = []
for key in ("rows", "rows_tight_budget", "rows_floor60"):
    for r in d.get(key, []):
        pts.append((r["run"], r["cost_rel_fixed138"], r["rms_m"]))
label_map = {
    "fixed_138": (r"fixed $N=138$", C1, "o"),
    "fixed_106": (r"fixed $N=106$", C1, "s"),
    "adaptive_p17_q10": (r"schedule $\varepsilon=10^{-2}$, floor 37", C2, "^"),
    "adaptive_p17_q10_eps0.001": (r"schedule $\varepsilon=10^{-3}$, floor 37", C2, "v"),
    "adaptive_p17_q10_eps0.01_floor60": (r"schedule $\varepsilon=10^{-2}$, floor 60", C3, "D"),
    "adaptive_p17_q10_eps0.001_floor60": (r"schedule $\varepsilon=10^{-3}$, floor 60", C3, "P"),
}
fig, ax = plt.subplots(figsize=(4.9, 3.4))
for run, cost, rms in pts:
    lbl, color, marker = label_map[run]
    ax.semilogy(cost, rms, marker, color=color, ms=7, label=lbl)
ax.set_xlabel(r"Relative cost (mean $N^2$; fixed $N=138$ = 1)")
ax.set_ylabel("24 h RMS position error [m]")
ax.set_xlim(0.25, 1.1)
ax.legend(frameon=False, fontsize=7.5, loc="upper right")
fig.savefig(FIGS / "fig_orbit_pareto.pdf")
plt.close(fig)

print("figures written:", sorted(p.name for p in FIGS.glob("*.pdf")))
