"""Tables for the R23 controls on the constructive claim.

Emits, from the records and never by hand:

  r23_oracle_table.tex          the interior member against the constant family
                                represented by its best member, next to the
                                single nominated degree it is currently
                                compared with
  r23_manuscript_descriptives.json
                                every number the manuscript quotes from R23,
                                including the R19 beta = 0.5 result when it
                                exists

Usage:  python rev23_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def load(name: str):
    p = METRICS / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def oracle_table(d: dict) -> str:
    s = d["summary"]
    o, sat = s["vs_oracle"], s["vs_saturating"]
    lines = ["\\begin{tabular}{@{}l r r r r@{}}", "\\toprule",
             "Constant comparator & interior & comparator & unres. & "
             "median $\\rho$ \\\\",
             "\\midrule",
             f"nominated, budget-saturating & {sat['interior_wins']} & "
             f"{sat['sat_wins']} & {sat['unresolved']} & "
             f"{sat['median_rho']:.2f} \\\\",
             f"best under budget (post-hoc oracle) & {o['interior_wins']} & "
             f"{o['oracle_wins']} & {o['unresolved']} & "
             f"{o['median_rho']:.2f} \\\\",
             "\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def per_orbit_table(d: dict) -> str:
    lines = ["\\begin{tabular}{@{}l r r r r r r@{}}", "\\toprule",
             "orbit & $h_p$ (km) & $N_{\\mathrm{sat}}$ & "
             "$N_{\\mathrm{oracle}}$ & $E_{\\mathrm{int}}$ (m) & "
             "$E_{\\mathrm{oracle}}$ (m) & verdict \\\\",
             "\\midrule"]
    for r in d["rows"]:
        v = r["winner_vs_oracle"]
        mark = {"interior": "I", "oracle": "O", None: "U"}[v]
        lines.append(
            f"{r['design']}{r['sobol_index']:03d} & {r['hp_km']:.0f} & "
            f"{r['n_sat']} & {r['n_oracle']} & "
            f"{r['interior_error_m']:.3g} & {r['oracle_error_m']:.3g} & "
            f"{mark} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def ultra_table(d: dict) -> str:
    """Does the beta = 1 separation survive a decade of tolerance?"""
    lines = ["\\begin{tabular}{@{}l r r r r r@{}}", "\\toprule",
             "Group & orbits & resolved & interior & fixed & median "
             "$M_{\\mathrm{res}}$ \\\\",
             "\\midrule"]
    for group in ("resolved", "borderline", "all"):
        s = d["summary"].get(group)
        if not s:
            continue
        label = {"resolved": "previously resolved",
                 "borderline": "previously borderline",
                 "all": "whole panel"}[group]
        lines.append(
            f"{label} & {s['orbits']} & {s['resolved_after']} & "
            f"{s['interior_wins']} & {s['fixed_wins']} & "
            f"{s['m_res_median_before']:.2f} $\\rightarrow$ "
            f"{s['m_res_median_after']:.2f} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def cost_cells(d: dict) -> list:
    """One entry per (model degree, degree), aggregated over every session.

    `rows` holds a single session; `all_sessions` holds all of them. Reading
    `rows` when `all_sessions` exists is how a three-session panel gets quoted
    as one session, so the aggregation is done here once and everything
    downstream reads these cells.
    """
    src = d.get("all_sessions") or d["rows"]
    cells: dict = {}
    for r in src:
        cells.setdefault((r["model_degree"], r["degree"]), []).append(r)
    out = []
    for (model, degree), rs in sorted(cells.items()):
        pc = [x["per_call_protocol"]["per_call_ns_median"] for x in rs]
        bl = [x["block_protocol"]["per_call_ns_median"] for x in rs]
        out.append({
            "model_degree": model, "degree": degree, "sessions": len(rs),
            "per_call_ns": float(np.median(pc)),
            "block_ns": float(np.median(bl)),
            "block_over_percall": float(np.median(
                [x["block_over_percall"] for x in rs])),
            # session-to-session spread, so no timing is quoted without it
            "per_call_max_over_min": max(pc) / min(pc),
            "block_max_over_min": max(bl) / min(bl),
        })
    return out


def cost_table(d: dict) -> str:
    """The two archived protocols measured side by side, over all sessions."""
    lines = ["\\begin{tabular}{@{}r r r r r r@{}}", "\\toprule",
             "$N$ & model & per-call ($\\mu$s) & block ($\\mu$s) & "
             "block/per-call & session spread \\\\",
             "\\midrule"]
    for c in cost_cells(d):
        if c["degree"] not in (60, 120, 300, 900):
            continue
        spread = max(c["per_call_max_over_min"], c["block_max_over_min"])
        lines.append(
            f"{c['degree']} & {c['model_degree']} & "
            f"{c['per_call_ns']/1e3:.1f} & {c['block_ns']/1e3:.1f} & "
            f"{c['block_over_percall']:.2f} & "
            f"$\\times${spread:.2f} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def _report_ultra(d: dict, out: dict) -> None:
    (METRICS / "r23_ultra_table.tex").write_text(ultra_table(d),
                                                 encoding="utf-8")
    pc = d["panel_completeness"]
    # decided both times and changed sides -- not "was undecided, now decided"
    flips = [f"{r['design']}{r['sobol_index']:03d}" for r in d["rows"]
             if r["winner_before"] is not None
             and r["winner_after"] is not None
             and r["winner_after"] != r["winner_before"]]
    dissolved = [f"{r['design']}{r['sobol_index']:03d}" for r in d["rows"]
                 if r["resolved_before"] and not r["resolved_after"]]
    gained = [f"{r['design']}{r['sobol_index']:03d}" for r in d["rows"]
              if not r["resolved_before"] and r["resolved_after"]]
    out["ultratight_span"] = {
        "summary": d["summary"], "panel_completeness": pc,
        "verdict_flips": flips,
        "separations_dissolved": dissolved,
        "separations_gained": gained,
    }
    print("[written] r23_ultra_table.tex")
    print(f"  panel: {pc['aggregated']} aggregated, {pc['missing']} missing")
    for group in ("resolved", "borderline", "all"):
        s = d["summary"].get(group)
        if not s:
            continue
        print(f"  {group:11s}: {s['orbits']:3d} orbits, resolved "
              f"{s['resolved_after']:3d} (interior {s['interior_wins']}, "
              f"fixed {s['fixed_wins']}), M_res median "
              f"{s['m_res_median_before']:.2f} -> {s['m_res_median_after']:.2f}"
              f", sign stable {s['gain_sign_stable'] if 'gain_sign_stable' in s else s['gap_sign_stable']}/{s['orbits']}")
    print(f"  verdict flips: {flips or 'none'}")
    print(f"  separations that dissolved: {dissolved or 'none'}")
    print(f"  separations that appeared:  {gained or 'none'}")


def _report_cost(d: dict, out: dict) -> None:
    (METRICS / "r23_cost_table.tex").write_text(cost_table(d), encoding="utf-8")
    me = d.get("model_degree_effect") or []
    # every measured cell, not one session's worth of them
    prot = [r["block_over_percall"]
            for r in (d.get("all_sessions") or d["rows"])]
    spread = d.get("session_to_session_spread") or []
    out["cost_curve_unified"] = {
        "n_sessions": d.get("n_sessions", 1),
        "n_cells": len(prot),
        "ratio_900_over_300_remeasured": d.get("ratio_900_over_300_remeasured"),
        "ratio_900_over_300_per_session": d.get(
            "ratio_900_over_300_per_session"),
        "at_quoted_degrees": d["archived_curve_comparison"]["at_quoted_degrees"],
        "block_over_percall": {
            "min": min(prot), "max": max(prot),
            "median": float(np.median(prot))},
        "session_spread": {
            k: {"median": float(np.median([s[k]["rel_spread"]
                                           for s in spread])),
                "max": max(s[k]["rel_spread"] for s in spread)}
            for k in ("per_call_protocol", "block_protocol")} if spread else {},
        "model_degree_effect": me,
        "cells": cost_cells(d),
    }
    print("[written] r23_cost_table.tex")
    print(f"  block/per-call ratio across the ladder: "
          f"{min(prot):.2f} to {max(prot):.2f} (median {np.median(prot):.2f})")
    if me:
        v = [x["model900_over_model300_percall"] for x in me]
        print(f"  model-900 over model-300, per-call: "
              f"{min(v):.2f} to {max(v):.2f}")
    for cell in d["archived_curve_comparison"]["at_quoted_degrees"]:
        parts = ", ".join(f"{k.replace('_ns', '')} {v/1e3:.1f}"
                          for k, v in cell.items() if k != "degree")
        print(f"  N={cell['degree']} us: {parts}")


def main() -> int:
    out: dict = {}

    d = load("r23_ultratight_span.json")
    if d:
        _report_ultra(d, out)
    else:
        print("[pending] r23_ultratight_span.json not written yet")

    # the three-session panel supersedes the single-session pilot; the pilot is
    # kept on disk rather than overwritten, but it is not what gets quoted
    d = (load("r23_cost_curve_reproducibility.json")
         or load("r23_cost_curve_unified.json"))
    if d:
        _report_cost(d, out)
    else:
        print("[pending] no cost curve record yet")

    d = load("r23_oracle_vs_interior.json")
    if d:
        chk = d["endpoint_reproduction_check"]
        (METRICS / "r23_oracle_table.tex").write_text(
            oracle_table(d), encoding="utf-8")
        (METRICS / "r23_oracle_per_orbit_table.tex").write_text(
            per_orbit_table(d), encoding="utf-8")
        gains = [r["oracle_gain_over_sat"] for r in d["rows"]]
        out["oracle_vs_interior"] = {
            "summary": d["summary"],
            "endpoint_reproduction_check": chk,
            "median_oracle_gain_over_saturating": float(np.median(gains)),
            "orbits_interior_loses_to_oracle": sorted(
                f"{r['design']}{r['sobol_index']:03d}" for r in d["rows"]
                if r["winner_vs_oracle"] == "oracle"),
            "orbits_interior_beats_oracle": sorted(
                f"{r['design']}{r['sobol_index']:03d}" for r in d["rows"]
                if r["winner_vs_oracle"] == "interior"),
        }
        print("[written] r23_oracle_table.tex, r23_oracle_per_orbit_table.tex")
        print(f"  reproduction of the archived k=0 endpoint: worst "
              f"{chk['worst_relative_difference']:.2%} -> "
              f"{'PASS' if chk['passes'] else 'FAIL'}")
        for label, key in (("nominated degree", "vs_saturating"),
                           ("post-hoc oracle", "vs_oracle")):
            b = d["summary"][key]
            comp = b.get("sat_wins", b.get("oracle_wins"))
            print(f"  interior vs {label}: {b['interior_wins']} / {comp}, "
                  f"unresolved {b['unresolved']}, "
                  f"median rho {b['median_rho']:.2f}")
    else:
        print("[pending] r23_oracle_vs_interior.json not written yet")

    r19 = load("r19_manuscript_descriptives.json")
    if r19:
        out["r19_by_budget"] = r19
        for key, v in r19.items():
            if not key.endswith("beta_0.50"):
                continue
            print(f"  R19 {key}: realized {v['interior_wins']}/"
                  f"{v['fixed_wins']}/{v['unresolved']}, median rho "
                  f"{v['median_rho_realized']:.2f}")

    if out:
        (METRICS / "r23_manuscript_descriptives.json").write_text(
            json.dumps(out, indent=2), encoding="utf-8")
        print("[written] r23_manuscript_descriptives.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
