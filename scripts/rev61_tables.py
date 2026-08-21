"""Tables for R61 (O42-ext): the scoring-tolerance rematch on the third
coverage design and the five geometry strata.

Two tables, because two different questions are being answered.

  * ``r61_equal_work_table.tex`` is the evidence: every registered cell, the
    tight-level tally beside the level-consistent one, so a reader can see
    which cells the convention moves and by how much.
  * ``r61_bracket_shift_table.tex`` is the claim: one row per population, the
    crossing bracket under each convention. This is the table Section IX.B's
    sentence should cite, because that sentence is about brackets and not
    about cells.

Three conventions, each carried over from rev44_tables.py for a reason:

Paired baseline. The tight-level counterpart of an R61 cell is the R19 cell,
not the R30 ladder. R19 holds the same quantity equal at the other level;
the ladder is a different comparator, and pairing against it would report a
comparator change as a convention change.

Like-for-like ratios. ``median_of_ratios`` is the statistic comparable with
median rho; ``ratio_of_medians`` is emitted alongside under its own name so
the two cannot be silently swapped.

Boundary cells. A tally that flips the leading side when the resolution
multiple moves over 0.5, 1 and 2 is not a verdict, and two of this
campaign's flips are narrow (17-18 and 18-21). Every cell is re-tallied at
the three cuts, nothing repropagated, and a cell whose leading side is not
stable across them is marked as a boundary cell and is reported as one
wherever its bracket is quoted.

Usage:  python rev61_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

POPS = [("C", "Design C"), ("SL", "Low perilune"), ("SP", "Near-polar"),
        ("SE", "Near-equatorial"), ("SF", "Frozen-like"),
        ("SH", "High apolune")]
BETAS = (0.50, 0.75, 1.00)
CUTS = (0.5, 1.0, 2.0)


def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def load(prefix: str, key: str, beta: float):
    p = METRICS / f"{prefix}_{key}_{beta_tag(beta)}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def lead(interior: int, fixed: int) -> str:
    if interior > fixed:
        return "interior"
    if fixed > interior:
        return "fixed"
    return "tie"


def retally(rows: list[dict]) -> dict:
    """Resolution-cut re-count. Nothing is repropagated: only the multiple of
    the recorded summed envelope required to call a comparison resolved
    changes."""
    out = {}
    for cut in CUTS:
        w = l = u = 0
        for r in rows:
            if not r.get("work_matched_error_m"):
                continue
            diff = r["work_matched_error_m"] - r["interior_error_m"]
            thr = cut * r["resolution_threshold_m"]
            if diff > thr:
                w += 1
            elif -diff > thr:
                l += 1
            else:
                u += 1
        out[f"{cut:g}"] = {"interior": w, "fixed": l, "unresolved": u,
                           "lead": lead(w, l)}
    return out


def is_boundary(rt: dict) -> bool:
    return len({v["lead"] for v in rt.values()}) > 1


def bracket(leads: dict) -> str:
    """The crossing bracket over the sampled grid, stated in the paper's own
    convention: the interval between the largest budget the constant degree
    leads at and the smallest budget above it the interior member leads at.

    A grid of three budgets can only ever locate the crossing to one of these
    answers, and 'below'/'above' are reported as such rather than as an
    interval the grid never sampled.
    """
    ordered = sorted(leads)
    if leads[ordered[0]] == "interior":
        return f"below {ordered[0]:.2f}"
    for lo, hi in zip(ordered, ordered[1:]):
        if leads[lo] == "fixed" and leads[hi] == "interior":
            return f"({lo:.2f}, {hi:.2f}]"
    if leads[ordered[-1]] == "fixed":
        return f"above {ordered[-1]:.2f}"
    return "not located"


def collect() -> dict:
    out: dict = {}
    for key, name in POPS:
        cells, leads_new, leads_old = {}, {}, {}
        for beta in BETAS:
            new = load("r61_equal_work_tighter", key, beta)
            old = load("r19_equal_total_work", key, beta)
            if not new or not old:
                continue
            sn, so = new["summary"], old["summary"]
            rt = retally(new["rows"])
            ln = lead(sn["resolved_interior_wins"], sn["resolved_fixed_wins"])
            lo = lead(so["resolved_interior_wins"], so["resolved_fixed_wins"])
            leads_new[beta], leads_old[beta] = ln, lo
            rows = new["rows"]
            ratios = [r["rho_workmatched"] for r in rows
                      if r.get("rho_workmatched")]
            e_int = [r["interior_error_m"] for r in rows
                     if r.get("interior_error_m")]
            e_fix = [r["work_matched_error_m"] for r in rows
                     if r.get("work_matched_error_m")]
            shift = [r["work_matched_degree"] - r["constant_endpoint_degree"]
                     for r in rows if r.get("work_matched_degree")]
            tight = [r["achieved_work_ratio_tight"] for r in rows
                     if r.get("achieved_work_ratio_tight")]
            cells[f"{beta:.2f}"] = {
                "beta": beta,
                "level_consistent": {
                    "orbits": sn["orbits"], "resolved": sn["resolved"],
                    "interior": sn["resolved_interior_wins"],
                    "fixed": sn["resolved_fixed_wins"],
                    "unresolved": sn["unresolved"],
                    "median_rho": sn["median_rho"], "lead": ln,
                    "achieved_work_ratio_tighter":
                        sn["achieved_work_ratio_tighter"],
                    "achieved_work_ratio_tight_median": float(
                        np.median(tight)),
                    "median_of_ratios": float(np.median(ratios)),
                    "ratio_of_medians": float(
                        np.median(e_fix) / np.median(e_int)),
                    "comparator_degree_shift_median": float(np.median(shift)),
                    "censored": sum(1 for r in rows if r.get("censored")),
                    "orbits_where_interior_loses": sorted(
                        r["sobol_index"] for r in rows
                        if r.get("winner") == "fixed"),
                },
                "tight_level": {
                    "orbits": so["orbits"], "resolved": so["resolved"],
                    "interior": so["resolved_interior_wins"],
                    "fixed": so["resolved_fixed_wins"],
                    "unresolved": so["unresolved"],
                    "median_rho": so["median_rho"], "lead": lo,
                },
                "resolution_cut_retally": rt,
                "boundary_cell": is_boundary(rt),
                "lead_changed": ln != lo,
            }
        if not cells:
            continue
        b_new, b_old = bracket(leads_new), bracket(leads_old)
        out[key] = {
            "population": name, "design_key": key, "cells": cells,
            "bracket_tight_level": b_old,
            "bracket_level_consistent": b_new,
            "bracket_moved": b_new != b_old,
            "boundary_cells": [b for b, c in cells.items()
                               if c["boundary_cell"]],
        }
    return out


def evidence_table(data: dict) -> str:
    lines = ["\\begin{tabular}{@{}l c r r r r r@{}}", "\\toprule",
             "population & $\\beta$ & match & interior & fixed & unres. & "
             "median $\\rho$ \\\\", "\\midrule"]
    for i, (key, d) in enumerate(data.items()):
        if i:
            lines.append("\\midrule")
        first = True
        for tag, c in d["cells"].items():
            for which, label in (("tight_level", "tight level"),
                                 ("level_consistent", "scoring tol.")):
                v = c[which]
                head = d["population"] if first else ""
                beta = f"{c['beta']:.2f}" if label == "tight level" else ""
                mark = ("\\,$^{\\dagger}$" if which == "level_consistent"
                        and c["boundary_cell"] else "")
                lines.append(
                    f"{head} & {beta} & {label}{mark} & {v['interior']} & "
                    f"{v['fixed']} & {v['unresolved']} & "
                    f"{v['median_rho']:.2f} \\\\")
                first = False
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def tex_bracket(b: str) -> str:
    """Math mode around the interval only. 'below 0.50' set whole in math
    mode typesets the word as five italic variables."""
    for word in ("below", "above"):
        if b.startswith(word):
            return f"{word} ${b[len(word):].strip()}$"
    if b.startswith("("):
        return f"${b}$"
    return b


def bracket_table(data: dict) -> str:
    lines = ["\\begin{tabular}{@{}l l l c@{}}", "\\toprule",
             "population & tight-level match & scoring-tolerance match & "
             "moves \\\\", "\\midrule"]
    for key, d in data.items():
        mark = "yes" if d["bracket_moved"] else "no"
        dag = "\\,$^{\\dagger}$" if d["boundary_cells"] else ""
        lines.append(f"{d['population']} & "
                     f"{tex_bracket(d['bracket_tight_level'])} & "
                     f"{tex_bracket(d['bracket_level_consistent'])}{dag} & "
                     f"{mark} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    data = collect()
    if not data:
        print("[r61-tables] no cells found")
        return 1
    (METRICS / "r61_equal_work_table.tex").write_text(
        evidence_table(data), encoding="utf-8")
    (METRICS / "r61_bracket_shift_table.tex").write_text(
        bracket_table(data), encoding="utf-8")
    (METRICS / "r61_manuscript_descriptives.json").write_text(
        json.dumps(data, indent=2), encoding="utf-8")
    print("[written] r61_equal_work_table.tex, r61_bracket_shift_table.tex, "
          "r61_manuscript_descriptives.json")
    moved = [d["population"] for d in data.values() if d["bracket_moved"]]
    for key, d in data.items():
        flags = (" BOUNDARY " + ",".join(d["boundary_cells"])
                 if d["boundary_cells"] else "")
        print(f"  {d['population']:<17}{d['bracket_tight_level']:>14} -> "
              f"{d['bracket_level_consistent']:<14}"
              f"{'MOVES' if d['bracket_moved'] else 'same':<7}{flags}")
    print(f"  brackets moved: {len(moved)} of {len(data)} "
          f"({', '.join(moved) if moved else 'none'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
