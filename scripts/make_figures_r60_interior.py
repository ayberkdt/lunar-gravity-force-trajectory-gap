"""R60: the allocation-interior existence figure (main text, Section VIII).

Draws the one result Section VIII exists to establish, from the archived
(O28) family records and nothing else:

  panel (a)  per-orbit normalized error E_k / E_0 across the five sampled
             members, medians with p10--p90 intervals, designs A and B
  panel (b)  where each orbit's lowest raw sampled error falls, counted
             per k from the record's own archived argmin

Sources (frozen; this script propagates nothing and rescored nothing):

  metrics/r18_span_sweep_A_beta_1.00.json
  metrics/r18_span_sweep_B_beta_1.00.json

Two conventions are inherited from the record rather than chosen here:

  * the error field is `error_m`, the level the record's own `best_k`
    argmin was taken at. This script recomputes the argmin from `error_m`
    and REFUSES to draw if any orbit disagrees with the archived `best_k`,
    so the figure cannot silently adopt a different scoring convention.
  * the per-k counts of panel (b) are the archived `best_k_counts`; the
    recount is a check, not the source.

The raw argmin is drawn as the record defines it: no resolution rule is
applied and no cell is called a tie, because (O28) defined `best_k` on raw
errors. Resolved endpoint comparisons stay in the prose and Table 3.

Usage:  python make_figures_r60_interior.py
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

K_ORDER = ["0.00", "0.25", "0.50", "0.75", "1.00"]
DESIGNS = (("A", "o", "#1f5fa8"), ("B", "s", "#d1651a"))
INTERIOR = {"0.25", "0.50", "0.75"}


def load(design: str) -> dict:
    p = METRICS / f"r18_span_sweep_{design}_beta_1.00.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    if d["beta"] != 1.0:
        raise SystemExit(f"[refuse] {p.name}: beta is {d['beta']}, not 1.0")
    return d


def ratios(d: dict) -> np.ndarray:
    """Per-orbit E_k / E_0 on the record's own error field, orbits x k."""
    out = []
    for r in d["rows"]:
        e0 = r["entries"]["0.00"]["error_m"]
        out.append([r["entries"][k]["error_m"] / e0 for k in K_ORDER])
    return np.asarray(out)


def checked_counts(d: dict) -> dict[str, int]:
    """Archived best_k_counts, refused unless the raw argmin reproduces it."""
    recount: dict[str, int] = {k: 0 for k in K_ORDER}
    for r in d["rows"]:
        errs = {k: r["entries"][k]["error_m"] for k in K_ORDER}
        argmin = min(errs, key=errs.get)
        if argmin != r["best_k"]:
            raise SystemExit(f"[refuse] {d['design']}/{r['name']}: archived "
                             f"best_k {r['best_k']} != argmin(error_m) "
                             f"{argmin}; conventions have diverged")
        recount[argmin] += 1
    archived = {k: int(v) for k, v in d["summary"]["best_k_counts"].items()}
    for k in K_ORDER:
        if recount[k] != archived.get(k, 0):
            raise SystemExit(f"[refuse] {d['design']}: recount {recount} != "
                             f"archived best_k_counts {archived}")
    return recount


def main() -> int:
    data = {des: load(des) for des, _, _ in DESIGNS}

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(6.0, 2.35), constrained_layout=True)

    # (a) the sampled response, normalized per orbit so the population's
    # spread in absolute difficulty cannot hide the interior shape
    x = [float(k) for k in K_ORDER]
    for des, marker, colour in DESIGNS:
        rr = ratios(data[des])
        med = np.median(rr, axis=0)
        p10 = np.percentile(rr, 10, axis=0)
        p90 = np.percentile(rr, 90, axis=0)
        ax_a.fill_between(x, p10, p90, color=colour, alpha=0.14, lw=0)
        ax_a.plot(x, med, marker=marker, ms=4.5, lw=1.2, color=colour,
                  label=f"design {des}", markeredgecolor="white",
                  markeredgewidth=0.5)
    ax_a.axhline(1.0, color="0.55", lw=0.8, ls="--", zorder=0)
    ax_a.set_yscale("log")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax_a.set_xlabel(r"allocation concentration $k$", labelpad=2)
    ax_a.set_ylabel(r"$E_k/E_{k=0}$, per orbit")
    ax_a.legend(loc="upper left", fontsize=6.4, frameon=False,
                borderpad=0.2, handletextpad=0.5, labelspacing=0.3)
    ax_a.set_title(r"(a) median and p10--p90 of the normalized error",
                   fontsize=7.6, pad=4)
    ax_a.tick_params(labelsize=7.5)

    # (b) where the raw sampled minimum falls, per orbit
    counts = {des: checked_counts(data[des]) for des, _, _ in DESIGNS}
    width = 0.34
    xs = np.arange(len(K_ORDER))
    ax_b.axvspan(0.5, 3.5, color="0.93", zorder=0)
    for i, (des, _, colour) in enumerate(DESIGNS):
        vals = [counts[des][k] for k in K_ORDER]
        ax_b.bar(xs + (i - 0.5) * width, vals, width=width, color=colour,
                 alpha=0.88, label=f"design {des}", zorder=2)
        for xx, v in zip(xs + (i - 0.5) * width, vals):
            ax_b.text(xx, v + 0.8, str(v), ha="center", va="bottom",
                      fontsize=6.2, color="0.25")
    n_int = {des: sum(counts[des][k] for k in K_ORDER if k in INTERIOR)
             for des, _, _ in DESIGNS}
    n_orb = {des: sum(counts[des].values()) for des, _, _ in DESIGNS}
    # upper-left, clear of the k=0.5 bars and their count labels
    ax_b.text(0.03, 0.965,
              "interior minima:\n"
              + ", ".join(f"{n_int[d]}/{n_orb[d]} ({d})"
                          for d, _, _ in DESIGNS),
              transform=ax_b.transAxes, ha="left", va="top",
              fontsize=6.6, color="0.20", linespacing=1.35)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax_b.set_xlabel(r"allocation concentration $k$", labelpad=2)
    ax_b.set_ylabel("orbits with lowest sampled error")
    ax_b.set_ylim(0, 33)
    ax_b.set_title(r"(b) location of the raw sampled minimum",
                   fontsize=7.6, pad=4)
    ax_b.tick_params(labelsize=7.5)
    ax_b.legend(loc="upper right", fontsize=6.4, frameon=False,
                borderpad=0.2, handletextpad=0.5, labelspacing=0.3)

    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "fig_allocation_interior.pdf"
    fig.savefig(out)
    plt.close(fig)
    for des, _, _ in DESIGNS:
        print(f"[checked] design {des}: interior minima "
              f"{n_int[des]}/{n_orb[des]}, per-k {counts[des]}")
    print(f"[written] figures/{out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
