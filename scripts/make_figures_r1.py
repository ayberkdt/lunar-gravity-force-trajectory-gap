"""Revision-1 figure generation from metrics/r1_*.json."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import paper_style
from paper_style import C1, C2, C3, C4, C5

BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
FIGS = BASE / "figures"
FIGS.mkdir(exist_ok=True)

paper_style.apply()


def _load(name: str) -> dict:
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


# ---------------- Fig: spectrum with fits and regularization boundary
d = _load("r1_spectrum_pfit.json")
n = np.array(d["spectrum_arrays"]["n"])
sig = np.array(d["spectrum_arrays"]["sigma_coeff_rms"])
f600 = d["fit_10_600"]
ffull = d["fit_2_1800"]
fig, (ax, axr) = plt.subplots(2, 1, figsize=(5.4, 4.9), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1.2],
                                           "hspace": 0.12})
ax.loglog(n, sig, "-", color="0.55", lw=0.7, label=r"JGGRX\_1800F $\sigma_n$")
xs = np.array([10.0, 600.0])
ax.loglog(xs, 10 ** (f600["logK"] - f600["p"] * np.log10(xs)), "-", color=C1,
          lw=1.8, label=rf"OLS fit $n\in[10,600]$: $p={f600['p']:.2f}$")
xs2 = np.array([2.0, 1800.0])
ax.loglog(xs2, 10 ** (ffull["logK"] - ffull["p"] * np.log10(xs2)), "--",
          color=C3, lw=1.4, label=rf"OLS fit $n\in[2,1800]$: $p={ffull['p']:.2f}$")
ax.axvspan(600, 1800, color="0.85", alpha=0.5)
ax.text(950, sig[5], "constrained\n($n>600$)", fontsize=8, color="0.35",
        ha="center")
ax.set_ylabel(r"Per-coefficient RMS $\sigma_n$")
ax.legend(frameon=False, fontsize=8, loc="lower left")
# residual panel (10..600 fit)
mask = (n >= 10)
resid = np.log10(sig[mask]) - (f600["logK"] - f600["p"] * np.log10(n[mask]))
axr.semilogx(n[mask], resid, "-", color=C1, lw=0.7)
axr.axhline(0, color="k", lw=0.6)
axr.axvspan(600, 1800, color="0.85", alpha=0.5)
axr.set_xlabel(r"Degree $n$")
axr.set_ylabel("Fit residual [dex]")
fig.savefig(FIGS / "fig_spectrum.pdf")
plt.close(fig)

# ---------------- Fig: degree-RMS verification
d = _load("r1_degree_rms_verification.json")
rows = [r for r in d["rows"] if r["altitude_km"] == 50.0]
degs = [r["degree"] for r in rows]
fig, ax = plt.subplots(figsize=(4.9, 3.0))
ax.plot(degs, [r["vector_ratio"] for r in rows], "o-", color=C1,
        label="total vector")
ax.plot(degs, [r["radial_ratio"] for r in rows], "s--", color=C2,
        label="radial only")
ax.axhline(1.0, color="k", lw=0.7)
ax.set_xlabel(r"Degree $n$")
ax.set_ylabel("Sampled RMS / analytic RMS")
ax.set_ylim(0.9, 1.1)
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_degree_rms_check.pdf")
plt.close(fig)

# ---------------- Fig: truncation criteria (updated, vector forms)
d = _load("r1_spectrum_pfit.json")
rows = d["criteria_rows"]
h = [r["altitude_km"] for r in rows]
pdense = _load("supplemental_pstar_uncertainty.json")
pstar_by_h = dict(zip(pdense["altitudes_km"],
                      pdense["predicted_nmin_at_best_p"]))
pstar_by_h.update({20.0: 328, 30.0: 219})
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.plot(h, [r["emp_vector_eps1e2"] for r in rows], "o-", color="black",
        lw=1.8, ms=4.5, label="Empirical (JGGRX spectrum)", zorder=5)
ax.plot(h, [r["proxy_vec_p1_7"] for r in rows], "s--", color=C1, lw=1.4,
        ms=4, label=r"Proxy, $p=1.7$")
ax.plot(h, [pstar_by_h[float(x)] for x in h], "P--", color=C4, lw=1.2,
        ms=4, label=r"Proxy, $p_{\mathrm{fit}}=1.76$")
ax.plot(h, [r["proxy_vec_p2_0"] for r in rows], "^--", color=C3, lw=1.4,
        ms=4, label=r"Proxy, $p=2.0$ (Kaula)")
ax.plot(h, [r["atten_1e3"] for r in rows], "d--", color=C2, lw=1.4, ms=4,
        label=r"Conventional attenuation, $f=10^{-3}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([20, 30, 50, 80, 100, 150, 200, 300])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_yticks([25, 50, 100, 200, 400, 600])
ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_xlabel("Altitude $h$ [km]")
ax.set_ylabel(r"Recommended truncation degree $N_{\min}$")
ax.legend(frameon=False, fontsize=7.5)
fig.savefig(FIGS / "fig_truncation.pdf")
plt.close(fig)

# ---------------- Fig: band shares with CIs
d = _load("r1_band_shares.json")
rows = d["rows"]
h = [r["altitude_km"] for r in rows]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
for key, color, marker, lbl in (
        ("b2_60", C1, "o", r"$2 \leq n \leq 60$"),
        ("b61_100", C2, "s", r"$61 \leq n \leq 100$"),
        ("tail", C3, "^", r"$n > 100$")):
    y = np.array([r[key]["share"] for r in rows])
    lo = np.array([r[key]["ci95"][0] for r in rows])
    hi = np.array([r[key]["ci95"][1] for r in rows])
    ax.errorbar(h, y, yerr=[y - lo, hi - y], fmt=marker + "-", color=color,
                lw=1.5, ms=4, capsize=2, label=lbl)
ax.set_yscale("log")
ax.axhline(1e-2, color="gray", ls=":", lw=1)
ax.text(150, 1.25e-2, r"1\% threshold", color="gray", fontsize=8)
ax.set_xlabel("Altitude $h$ [km]")
ax.set_ylabel("Normalized band-difference RMS")
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_band_share.pdf")
plt.close(fig)

# ---------------- Fig: timing (median + IQR)
d = _load("r1_kernel_timing.json")
sweep = d["degree_sweep"]
N = np.array([r["degree"] for r in sweep], dtype=float)
med = np.array([r["median_us"] for r in sweep])
q1 = np.array([r["q1_us"] for r in sweep])
q3 = np.array([r["q3_us"] for r in sweep])
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.loglog(N, med, "o-", color=C1, lw=1.6, ms=4,
          label="Serial Kahan kernel (median)")
ax.fill_between(N, q1, q3, color=C1, alpha=0.25, lw=0, label="IQR over blocks")
ref = med[-1] * (N / N[-1]) ** 2
ax.loglog(N, ref, ":", color="gray", lw=1.2, label=r"$N^2$ reference slope")
ax.set_xlabel(r"Truncation degree $N$")
ax.set_ylabel(r"Evaluation time [$\mu$s]")
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_timing.pdf")
plt.close(fig)

# ---------------- Fig: switching jump with CIs
d = _load("r1_switch_jump.json")
rows = d["rows"]
fig, ax = plt.subplots(figsize=(4.9, 3.2))
for h_km, color, marker in [(50.0, C1, "o"), (100.0, C2, "s")]:
    sel = [r for r in rows if r["altitude_km"] == h_km]
    q = np.array([r["step"] for r in sel], dtype=float)
    y = np.array([100.0 * r["jump_over_pert"] for r in sel])
    lo = np.array([100.0 * r["ci95"][0] for r in sel])
    hi = np.array([100.0 * r["ci95"][1] for r in sel])
    ax.errorbar(q, y, yerr=[y - lo, hi - y], fmt=marker + "-", color=color,
                lw=1.5, ms=4.5, capsize=2, label=f"$h = {h_km:.0f}$ km")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([5, 10, 25, 50])
ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
ax.set_xlabel(r"Quantization step $q$ (degrees)")
ax.set_ylabel(r"RMS jump / perturbation RMS [\%]")
ax.legend(frameon=False, fontsize=8)
fig.savefig(FIGS / "fig_switch_jump.pdf")
plt.close(fig)

# ---------------- Fig: 24 h truncation mapping, both tolerance sets
d = _load("r1_orbit_mapping.json")
rows = d["rows"]
conv = {r["degree"]: r for r in d["convergence_control"]["rows"]}
N6 = [r["degree"] for r in rows]
eB = [r["pos_rms_m"] for r in rows]
fig, ax = plt.subplots(figsize=(4.9, 3.4))
ax.semilogy(N6, eB, "o-", color=C1, lw=1.7, ms=4.5,
            label=r"baseline tol. (rtol $10^{-11}$)")
NT = sorted(conv)
ax.semilogy(NT, [conv[k]["rms_tight_tol"] for k in NT], "s--", color=C2,
            lw=1.5, ms=4.5, label=r"tight tol. (rtol $10^{-12}$)")
floor = d["convergence_control"]["integrator_floor_rms_m"]
ax.axhline(floor, color="gray", ls=":", lw=1)
ax.text(22, floor * 0.80, "baseline-tolerance integration floor",
        fontsize=7.5, color="gray", ha="left", va="top")
blend = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
for x, color, lbl in [(56, C3, r"$p=2.0$ pick"), (68, "black", r"$p^{*}=1.76$ pick"),
                      (124, C2, "attenuation pick")]:
    ax.axvline(x, color=color, ls=":", lw=1.0)
    ax.text(x, 0.985, lbl, rotation=90, ha="right", va="top",
            fontsize=7.5, color=color, transform=blend,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
ax.set_xlabel(r"Truncation degree $N$")
ax.set_ylabel("24 h RMS position error [m]")
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.savefig(FIGS / "fig_orbit_error.pdf")
plt.close(fig)

print("static r1 figures written")

# ---------------- Fig: 7-day long-arc error growth (tight tolerance)
d = _load("r1_longarc_tight.json")
runs = [("fixed_138", C1, "-", r"fixed $N=138$"),
        ("fixed_106", C1, "--", r"fixed $N=106$"),
        ("sched_dwell_alt", C2, "-", r"schedule $\varepsilon=10^{-3}$, floor 60, down"),
        ("sched_dwell_up", C3, "-", r"schedule, up-quantized"),
        ("sched_mindwell600", C4, "-", "schedule, 600 s min dwell")]
fig, ax = plt.subplots(figsize=(5.2, 3.5))
for run, color, ls, lbl in runs:
    s = d["error_series"][run]
    t = np.array(s["t_s"]) / 86400.0
    e = np.sqrt(np.array(s["radial_m"])**2 + np.array(s["in_track_m"])**2 +
                np.array(s["cross_track_m"])**2)
    # decimate to a smooth envelope: max over 0.1-day windows
    nbin = int(t[-1] / 0.1)
    tb, eb = [], []
    for k in range(nbin):
        m = (t >= k * 0.1) & (t < (k + 1) * 0.1)
        if m.any():
            tb.append(t[m].mean())
            eb.append(e[m].max())
    ax.semilogy(tb, eb, ls, color=color, lw=1.5, label=lbl)
ax.set_xlabel("Time [days]")
ax.set_ylabel("Position error envelope vs. $N=300$ truth [m]")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
fig.savefig(FIGS / "fig_longarc_growth.pdf")
plt.close(fig)

# ---------------- Fig: schedule profile + dwell (multi-panel)
d2 = _load("r1_longarc.json")
prof = d2["profile"]
tp = np.array(prof["t_s"]) / 3600.0
fig, axs = plt.subplots(3, 1, figsize=(5.4, 5.6), sharex=True,
                        gridspec_kw={"height_ratios": [1.4, 1.0, 1.0]})
axs[0].plot(tp, prof["altitude_km"], color=C1, lw=1.4)
axs[0].set_ylabel("Altitude [km]")
axs[1].step(tp, prof["degree"], where="post", color=C2, lw=1.4)
axs[1].set_ylabel("Scheduled $N$")
sw = _load("r1_switch_instrumentation.json")
case = sw["cases"]["perilune_start/scheduled"]["series"]
ts = np.array(case["t_s"]) / 3600.0
axs[2].plot(ts, case["h_s"], ".-", color=C3, lw=0.8, ms=2.5,
            label="accepted step size")
att = np.array(case["attempts"])
rej = att > 1
axs[2].plot(ts[rej], np.array(case["h_s"])[rej], "x", color=C5, ms=5,
            label="step with rejected attempt(s)")
degs = np.array(case["degree"])
swi = np.flatnonzero(degs[1:] != degs[:-1]) + 1
for x in ts[swi]:
    axs[2].axvline(x, color="0.8", lw=0.6, zorder=0)
axs[2].set_ylabel("Step size [s]")
axs[2].set_xlabel("Time [h]")
axs[2].set_xlim(0, float(ts[-1]))
axs[2].set_ylim(0, 445)
axs[2].legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2,
              columnspacing=1.2, handletextpad=0.5)
fig.savefig(FIGS / "fig_schedule_panels.pdf")
plt.close(fig)

# ---------------- Fig: measured cost vs accuracy (7-day, tight)
rows = {r["run"]: r for r in d["rows"]}
label_map = {
    "fixed_138": (r"fixed $N=138$", C1, "o"),
    "fixed_106": (r"fixed $N=106$", C1, "s"),
    "sched_dwell_alt": (r"sched.\ $\varepsilon{=}10^{-3}$ fl.60 down", C2, "^"),
    "sched_naive_alt": (r"sched.\ $\varepsilon{=}10^{-2}$ fl.37 down", C2, "v"),
    "sched_dwell_up": (r"sched.\ up-quantized", C3, "D"),
    "sched_coarse3": (r"sched.\ coarse 3-level", C4, "P"),
    "sched_mindwell600": (r"sched.\ 600 s min dwell", C4, "X"),
}
fig, ax = plt.subplots(figsize=(5.2, 3.5))
for run, (lbl, color, marker) in label_map.items():
    r = rows[run]
    ax.semilogy(r["grav_s"], r["pos_rms_m"], marker, color=color, ms=7,
                label=lbl)
ax.set_xlabel("Measured gravity-kernel wall time [s] (7-day arc)")
ax.set_ylabel("7-day RMS position error [m]")
ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=1)
fig.savefig(FIGS / "fig_pareto_walltime.pdf")
plt.close(fig)

# ---------------- Fig: switch zoom (step behavior near one switch)
case = sw["cases"]["apolune_start/scheduled"]["series"]
ts = np.array(case["t_s"]); hs = np.array(case["h_s"])
degs = np.array(case["degree"]); att = np.array(case["attempts"])
swi = np.flatnonzero(degs[1:] != degs[:-1]) + 1
mid = swi[len(swi) // 3]
t0 = ts[mid]
m = np.abs(ts - t0) <= 1500.0
fig, (a1, a2) = plt.subplots(2, 1, figsize=(5.0, 3.8), sharex=True)
a1.step((ts[m] - t0) / 60.0, degs[m], where="post", color=C2, lw=1.5)
a1.set_ylabel("Active $N$")
a2.plot((ts[m] - t0) / 60.0, hs[m], "o-", color=C3, lw=1.2, ms=3.5)
rej = m & (att > 1)
a2.plot((ts[rej] - t0) / 60.0, hs[rej], "x", color=C5, ms=6,
        label="rejected attempt(s)")
a2.axvline(0, color="0.7", lw=0.8)
a2.set_xlabel("Time relative to switch [min]")
a2.set_ylabel("Step size [s]")
a2.legend(frameon=False, fontsize=7.5)
fig.savefig(FIGS / "fig_switch_zoom.pdf")
plt.close(fig)

print("long-arc figures written")
