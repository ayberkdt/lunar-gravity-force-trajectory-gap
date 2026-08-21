"""Allocation anatomy: what the two equal-budget policies actually do along one orbit.

Section VI reports the endpoint comparison as a population verdict, and Section
VII explains it dynamically. Between the two the reader has no picture of the
allocation itself: the degree history the budget buys, and the trajectory
response that history produces. This figure supplies exactly that, and nothing
that is not already archived.

  panel (a)  altitude over the first two revolutions
  panel (b)  the two degree histories over the same window, from the archived
             binned schedule and the archived critical degree
  panel (c)  signed in-track displacement against the shared reference over the
             full seven-day arc, for both policies

The figure is explanatory, not evidential. The orbit is chosen by GEOMETRY
ALONE, by the deterministic rule

    argmin over design A of  |log(hp / median hp)| + |log(ha / median ha)|

evaluated on the design's own perilune and apolune coordinates, so no outcome
and no judgement enters the selection. It returns sobolA_032 with a score of
0.029 against 0.130 for the runner-up. The population evidence stays in the
budget-Pareto figure and Table 1.

Sources (frozen; this script propagates nothing and rescores nothing):

  metrics/r14_trajectory_A_beta_1.00.json                 design point, verdict
  metrics/r14_cases/A_beta_1.00/<orbit>/atallah_budget_tight.json
                                                          binned degree table
  <r11 raw>/convergence/<orbit>/truth_tight.npz           shared reference
  <r11 raw>/convergence/<orbit>/fixed_critical_tight.npz  constant degree
  <r14 raw>/A_beta_1.00/<orbit>/atallah_budget_tight.npz  radial history

The comparator at beta = 1 is the R11 critical-degree run, which is what the
(O14) record reuses; the R10 baseline runs of the same name are the earlier
scalar-tolerance campaign and are NOT interchangeable with it. Three refusals
guard against drawing anything the record does not say:

  * the output grids of the three runs must coincide;
  * the position RMS recomputed here must reproduce the archived
    atallah_error_m and fixed_error_m;
  * the mean squared degree of the reconstructed schedule must meet the
    campaign's one-percent work-match target against N_crit^2.

Usage:  python make_figures_allocation_anatomy.py
"""

from __future__ import annotations

import json
import os
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

# The run store was consolidated onto D:; the in-tree paths are kept as the
# fallback so the script still runs from a fresh checkout of the package.
RAW_ROOTS = [pathlib.Path(r"D:/Makale_Kosular/metrics"), METRICS]
R14_RAW_ROOTS = [pathlib.Path(r"D:/Makale_Kosular/raw_offload/r14_raw"),
                 METRICS / "r14_raw"]

RECORD = METRICS / "r14_trajectory_A_beta_1.00.json"
# GL1800F's own reference radius, as read from the coefficient-file header by
# every campaign and archived with the band shares.
RADIUS_SOURCE = METRICS / "e1_band_shares_60_100_nmax300.json"

# Fig. 3 (make_figures_r14.py) sets the convention for this policy pair:
# AT_C, FX_C = ps.C1, ps.C5. Two adjacent figures drawing the same two
# policies must not swap it, so the constants are taken from there rather
# than chosen here.
C_RAD = paper_style.C1   # budget-calibrated radial history
C_FIX = paper_style.C5   # constant degree
C_ALT = "0.35"

WORK_MATCH_TOL = 0.01     # the campaign's declared work-match target
ERROR_TOL = 1e-9          # relative, against the archived scalars


def first(paths, *parts):
    for base in paths:
        p = base.joinpath(*parts)
        if p.exists():
            return p
    raise SystemExit("[refuse] not found in any run store: %s" % (parts,))


def reference_radius_m() -> float:
    d = json.loads(RADIUS_SOURCE.read_text(encoding="utf-8"))
    r = float(d["gravity_model"]["reference_radius_m"])
    if not 1.7e6 < r < 1.8e6:
        raise SystemExit("[refuse] reference radius %r is not the Moon's" % r)
    return r


def select_orbit(rows) -> int:
    """Geometry-only choice: closest to the design medians in log(hp), log(ha)."""
    hp = np.array([r["design_point"]["hp_km"] for r in rows])
    ha = np.array([r["design_point"]["ha_km"] for r in rows])
    score = np.abs(np.log(hp / np.median(hp))) + np.abs(np.log(ha / np.median(ha)))
    return int(np.argmin(score))


def ric(sol, ref):
    """Radial / in-track / cross-track components of (sol - ref), reference frame."""
    d = (sol[:3] - ref[:3]).T
    r = ref[:3].T
    v = ref[3:].T
    er = r / np.linalg.norm(r, axis=1, keepdims=True)
    ec = np.cross(r, v)
    ec = ec / np.linalg.norm(ec, axis=1, keepdims=True)
    ei = np.cross(ec, er)
    return (np.sum(d * er, 1), np.sum(d * ei, 1), np.sum(d * ec, 1))


def perilune_indices(h: np.ndarray) -> np.ndarray:
    """Output epochs that are strict local minima of altitude, plus the start."""
    interior = np.where((h[1:-1] < h[:-2]) & (h[1:-1] <= h[2:]))[0] + 1
    return np.concatenate(([0], interior)) if h[0] < h[1] else interior


def main() -> int:
    rec = json.loads(RECORD.read_text(encoding="utf-8"))
    if rec["beta"] != 1.0 or rec["design"] != "A":
        raise SystemExit("[refuse] %s is not the design-A beta = 1 record" % RECORD.name)
    rows = rec["rows"]
    j = select_orbit(rows)
    row = rows[j]
    name = row["name"]
    dp = row["design_point"]

    cfg = json.loads((METRICS / "r14_cases" / "A_beta_1.00" / name /
                      "atallah_budget_tight.json").read_text(encoding="utf-8"))["config"]
    n_crit = int(cfg["n_critical"])
    if n_crit != int(row["fixed_degree"]) or not row["reuse_fixed_critical"]:
        raise SystemExit("[refuse] %s: the beta = 1 comparator is not the "
                         "critical degree in this record" % name)

    ref_p = first(RAW_ROOTS, "r11_raw", "convergence", name, "truth_tight.npz")
    fix_p = first(RAW_ROOTS, "r11_raw", "convergence", name, "fixed_critical_tight.npz")
    rad_p = first(R14_RAW_ROOTS, "A_beta_1.00", name, "atallah_budget_tight.npz")

    ref, fix, rad = (np.load(p) for p in (ref_p, fix_p, rad_p))
    t = np.asarray(ref["t_s"], float)
    if not (np.allclose(t, fix["t_s"]) and np.allclose(t, rad["t_s"])):
        raise SystemExit("[refuse] %s: the three runs do not share an output grid" % name)
    Y, F, A = (np.asarray(d["state_si"], float) for d in (ref, fix, rad))

    # --- refusal 1: the recomputed errors must be the record's own ---
    e_fix = float(np.sqrt(np.mean(np.linalg.norm(F[:3] - Y[:3], axis=0) ** 2)))
    e_rad = float(np.sqrt(np.mean(np.linalg.norm(A[:3] - Y[:3], axis=0) ** 2)))
    for got, want, label in ((e_fix, row["comparison"]["fixed_error_m"], "fixed"),
                             (e_rad, row["comparison"]["atallah_error_m"], "radial")):
        if abs(got - want) > ERROR_TOL * max(abs(want), 1e-30):
            raise SystemExit("[refuse] %s: recomputed %s error %.6f != archived "
                             "%.6f" % (name, label, got, want))

    R = reference_radius_m()
    h_ref = (np.linalg.norm(Y[:3], axis=0) - R) / 1e3
    h_pol = (np.linalg.norm(A[:3], axis=0) - R) / 1e3

    # --- the radial history, read off the archived binned schedule ---
    table = {float(k): int(v) for k, v in cfg["atallah_degree_table"].items()}
    edges = np.array(sorted(table))
    degrees = np.array([table[e] for e in edges])
    idx = np.clip(np.searchsorted(edges, h_pol, side="right") - 1, 0, len(edges) - 1)
    n_rad = degrees[idx]

    # --- refusal 2: the reconstruction must still hold the declared budget ---
    beta_hat = float(np.mean(n_rad.astype(float) ** 2)) / float(n_crit ** 2)
    if abs(beta_hat - 1.0) > WORK_MATCH_TOL:
        raise SystemExit("[refuse] %s: schedule reconstructed from the archived "
                         "table gives beta = %.4f, outside the campaign's "
                         "%.0f%% work-match target" % (name, beta_hat,
                                                       100 * WORK_MATCH_TOL))

    # --- window: the first two revolutions, located from the altitude series ---
    peri = perilune_indices(h_ref)
    if len(peri) < 3:
        raise SystemExit("[refuse] %s: fewer than three perilune passages found" % name)
    end = int(peri[2])
    w = slice(0, end + 1)
    t_h = t[w] / 3600.0

    _, ei_fix, _ = ric(F, Y)
    _, ei_rad, _ = ric(A, Y)
    t_d = t / 86400.0

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        3, 1, figsize=(5.5, 3.5), constrained_layout=True,
        gridspec_kw={"height_ratios": [0.8, 0.85, 1.3]})

    # (a) where the orbit is
    ax_a.plot(t_h, h_ref[w], color=C_ALT, lw=1.1)
    ax_a.set_ylabel("altitude [km]")
    ax_a.set_xlim(t_h[0], t_h[-1])
    ax_a.set_title("(a) altitude over two revolutions", fontsize=7.6, pad=3)
    ax_a.tick_params(labelsize=7.5, labelbottom=False)

    # (b) what each policy spends there
    ax_b.step(t_h, n_rad[w], where="post", color=C_RAD, lw=1.1,
              label="radial history")
    ax_b.axhline(n_crit, color=C_FIX, lw=1.2,
                 label=r"constant degree $N_{\mathrm{crit}}=%d$" % n_crit)
    ax_b.set_ylabel("degree $N$ (log)")
    ax_b.set_xlabel("time [h]", labelpad=2)
    ax_b.set_xlim(t_h[0], t_h[-1])
    ax_b.set_yscale("log")
    ax_b.set_yticks([30, 60, 120, 240])
    ax_b.set_yticklabels(["30", "60", "120", "240"])
    ax_b.set_yticks([], minor=True)
    ax_b.set_ylim(24, 560)
    ax_b.set_title("(b) equal nominal per-call budget, two ways to spend it",
                   fontsize=7.6, pad=3)
    ax_b.tick_params(labelsize=7.5)
    ax_b.legend(loc="upper center", fontsize=6.4, frameon=False,
                borderpad=0.2, handletextpad=0.6, labelspacing=0.25, ncol=2)

    # (c) what the trajectory does with it
    ax_c.axhline(0.0, color="0.6", lw=0.7, zorder=0)
    ax_c.plot(t_d, ei_rad, color=C_RAD, lw=1.0, label="radial history")
    ax_c.plot(t_d, ei_fix, color=C_FIX, lw=1.0, label="constant degree")
    ax_c.set_xlabel("time [d]", labelpad=2)
    ax_c.set_ylabel("in-track displacement [m]")
    ax_c.set_xlim(0, t_d[-1])
    ax_c.set_title("(c) signed in-track displacement from the shared reference",
                   fontsize=7.6, pad=3)
    ax_c.tick_params(labelsize=7.5)
    ax_c.legend(loc="upper left", fontsize=6.4, frameon=False,
                borderpad=0.2, handletextpad=0.6, labelspacing=0.25)
    ax_c.text(0.02, 0.60,
              r"arc RMS %.1f\,m against %.1f\,m" % (e_rad, e_fix),
              transform=ax_c.transAxes, ha="left", va="top",
              fontsize=6.6, color="0.20")

    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "fig_allocation_anatomy.pdf"
    # Without this the writer stamps the current time into the PDF, so an
    # unchanged figure hashes differently on every run and the manifest
    # digest reads stale for a reason that is not the content.
    fig.savefig(out, metadata={"CreationDate": None})
    plt.close(fig)

    print("[selected by geometry] %s: perilune %.1f km, apolune %.1f km, "
          "i = %.1f deg" % (name, dp["hp_km"], dp["ha_km"], dp["incl_deg"]))
    print("[checked] errors reproduce the record: radial %.4f m, constant "
          "%.4f m, rho = %.4f" % (e_rad, e_fix, e_fix / e_rad))
    print("[checked] reconstructed schedule holds beta = %.5f, degrees "
          "%d-%d against N_crit = %d" % (beta_hat, n_rad.min(), n_rad.max(), n_crit))
    print("[written] figures/%s" % out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
