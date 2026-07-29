"""Publication figures for the publication-readiness (R3) controls."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import paper_style
from paper_style import C1, C2, C3, C4, C5, C6, GRAY

BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
FIGS = BASE / "figures"
FIGS.mkdir(exist_ok=True)
paper_style.apply()


def load(name):
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


# Dense fixed-degree sweep: fidelity, in-track endpoint, timing, coherence.
d = load("r3_degree_sweep.json")
styles = {
    "phaseA": (C1, "o", "polar, perilune start"),
    "phaseB": (C2, "s", "polar, apolune start"),
    "inc60": (C3, "^", r"$i=60^\circ$, perilune start"),
}
fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.3), sharex=True)
for key, (color, marker, label) in styles.items():
    rows = d["results"][key]["rows"]
    n = np.array([r["degree"] for r in rows])
    rms = np.array([r["pos_rms_m"] for r in rows])
    itr = np.array([r["ric_final_m"]["in_track"] for r in rows])
    grav = np.array([r["grav_s"] for r in rows])
    coh = np.array([d["omitted_band_track_projection"][key][str(int(x))]
                    ["coherence_abs_mean_over_rms"] for x in n])
    axs[0, 0].plot(n, rms, color=color, marker=marker, label=label)
    axs[0, 1].plot(n, itr, color=color, marker=marker)
    axs[1, 0].plot(n, grav, color=color, marker=marker)
    axs[1, 1].plot(n, coh, color=color, marker=marker)
axs[0, 0].set_yscale("log")
axs[0, 0].set_ylabel("7-day RMS position error [m]")
axs[0, 0].legend(loc="upper right", fontsize=7.4)
axs[0, 1].axhline(0.0, color="0.25", lw=0.6)
axs[0, 1].set_ylabel("Final in-track error [m]")
axs[1, 0].set_ylabel("Measured gravity time [s]")
axs[1, 1].set_ylabel(r"$|\overline{a_I}|/{\rm RMS}(a_I)$")
for ax in axs[1, :]:
    ax.set_xlabel(r"Fixed truncation degree $N$")
for tag, ax in zip(("(a)", "(b)", "(c)", "(d)"), axs.flat):
    ax.text(0.02, 0.95, tag, transform=ax.transAxes, va="top")
fig.tight_layout()
fig.savefig(FIGS / "fig_degree_sweep.pdf")
plt.close(fig)


# Multi-geometry schedule matrix.
m = load("r3_longarc_matrix.json")
case_titles = {
    "M_phaseB": r"$50\times300$ km polar, apolune start",
    "M_inc60": r"$50\times300$ km, $i=60^\circ$",
    "M_100x300": r"$100\times300$ km polar",
    "M_moonpa": r"$50\times300$ km polar, DE440/MOON\_PA",
}
run_order = ["fixed_138", "fixed_106", "sched_down", "sched_up",
             "sched_mindwell600", "sched_emp"]
run_label = {"fixed_138": "fixed 138", "fixed_106": "fixed 106",
             "sched_down": "down", "sched_up": "up",
             "sched_mindwell600": "dwell", "sched_emp": "empirical"}
colors = [C1, C6, C5, C3, C4, C2]
# sharey is deliberately off: the four geometries span very different error
# ranges, so a shared axis clipped the smallest bars of the lower-range panels.
# Each panel is scaled to its own data, with headroom for the rotated labels.
fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.1), sharey=False)
for ax, (case, title) in zip(axs.flat, case_titles.items()):
    by = {r["run"]: r for r in m["cases"][case]["rows"]}
    y = [by[k]["pos_rms_m"] for k in run_order]
    x = np.arange(len(run_order))
    ax.bar(x, y, color=colors, edgecolor="0.25", linewidth=0.35)
    ax.set_yscale("log")
    ax.set_title(title, fontsize=8.5)
    ax.set_xticks(x, [run_label[k] for k in run_order], rotation=32,
                  ha="right")
    ax.set_ylabel("7-day RMS position error [m]")
    for xx, yy in zip(x, y):
        ax.text(xx, yy * 1.08, f"{yy:.1f}", ha="center", va="bottom",
                fontsize=6.6, rotation=90)
    ax.set_ylim(bottom=min(y) / 2.0, top=max(y) * 4.0)
fig.tight_layout()
fig.savefig(FIGS / "fig_longarc_matrix.pdf")
plt.close(fig)


# Event-aligned direct switching statistics.
s = load("r3_switch_direct.json")
fig, axs = plt.subplots(2, 2, figsize=(7.1, 4.8), sharex=True)
for col, phase in enumerate(("perilune_start", "apolune_start")):
    aligned = s["cases"][f"{phase}/scheduled"]["event_aligned"]
    for direction, color, label in (("down", C5, "downward switch"),
                                    ("up", C3, "upward switch")):
        z = aligned[direction]
        x = np.array(z["bin_center_s"]) / 60.0
        med = np.array([np.nan if v is None else v for v in z["median_step_s"]])
        q1 = np.array([np.nan if v is None else v for v in z["q1_step_s"]])
        q3 = np.array([np.nan if v is None else v for v in z["q3_step_s"]])
        rp = np.array([np.nan if v is None else v for v in z["rejection_probability"]])
        axs[0, col].plot(x, med, color=color, label=label)
        axs[0, col].fill_between(x, q1, q3, color=color, alpha=0.18, linewidth=0)
        axs[1, col].plot(x, rp, color=color)
    axs[0, col].axvline(0.0, color="0.2", lw=0.7)
    axs[1, col].axvline(0.0, color="0.2", lw=0.7)
    axs[0, col].set_title(phase.replace("_", " "), fontsize=8.5)
    axs[1, col].set_xlabel("Time from switch [min]")
axs[0, 0].set_ylabel("Accepted step [s]\nmedian and IQR")
axs[1, 0].set_ylabel("Rejection probability")
axs[0, 0].legend(loc="upper right", fontsize=7.3)
fig.tight_layout()
fig.savefig(FIGS / "fig_switch_aggregate.pdf")
plt.close(fig)


# Rotation interpolation parity: matrix and robust principal-angle metrics.
r = load("r3_rotation_parity.json")
fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.8))
markers = {"nodes": "o", "offgrid": "s", "stagelike": "^"}
set_colors = {"nodes": C1, "offgrid": C2, "stagelike": C5}
for probe in ("nodes", "offgrid", "stagelike"):
    rows = [z for z in r["rows"] if z["probe_set"] == probe]
    x = np.arange(len(rows))
    labels = [f'{int(z["duration_days"])}d/{int(z["table_dt_s"])}s' for z in rows]
    axs[0].plot(x, [z["frobenius_max"] for z in rows], marker=markers[probe],
                color=set_colors[probe], label=probe)
    axs[1].plot(x, [z["principal_angle_max_rad"] for z in rows],
                marker=markers[probe], color=set_colors[probe], label=probe)
for ax, ylabel in zip(axs, (r"Maximum $\|R_{\rm tab}-R_{\rm SPICE}\|_F$",
                            "Maximum principal angle [rad]")):
    ax.set_yscale("log")
    ax.set_xticks(np.arange(4), labels, rotation=25)
    ax.set_xlabel("Window / table cadence")
    ax.set_ylabel(ylabel)
axs[0].legend(fontsize=7.5)
fig.tight_layout()
fig.savefig(FIGS / "fig_rotation_parity.pdf")
plt.close(fig)

print("round-3 figures written")
