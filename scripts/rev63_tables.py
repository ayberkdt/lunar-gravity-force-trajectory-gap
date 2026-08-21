"""Tables and outcome adjudication for R63 (O55).

Four apolune levels instead of (O54)'s two, so the statistic of interest is
the *shape* of the margin across levels, not only whether two levels differ.
Everything is per level and per block; nothing is pooled across either.

The class is computed from the rule frozen in r63_preregistration.json. The
registration declares that class to be a summary of the beta = 1.00 row that
must be read beside the full table, so this generator always prints the full
table with it.

Usage:  python rev63_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BLOCKS = [("RS1U", "Block A, ceiling free"),
          ("RS2U", "Block B, ceiling free")]
BETAS = (1.00, 0.75, 0.50)
LEVELS = (300.0, 600.0, 1200.0, 2400.0)
CUTS = (0.5, 1.0, 2.0)


def beta_tag(b: float) -> str:
    return f"beta_{b:.2f}"


def apolune_map(key: str, beta: float) -> dict[int, float]:
    span = json.loads((METRICS / f"r18_span_sweep_{key}_{beta_tag(beta)}.json"
                       ).read_text(encoding="utf-8"))
    return {int(r["sobol_index"]): float(r["ha_km"]) for r in span["rows"]}


def lead(i: int, f: int) -> str:
    return "interior" if i > f else ("fixed" if f > i else "tie")


def tally(rows, ha, level, cut=None):
    i = f = u = 0
    rhos = []
    for r in rows:
        if ha.get(int(r["sobol_index"])) != level:
            continue
        if r.get("rho_workmatched"):
            rhos.append(r["rho_workmatched"])
        if cut is None:
            if not r.get("resolved"):
                u += 1
            elif r.get("winner") == "interior":
                i += 1
            elif r.get("winner") == "fixed":
                f += 1
            continue
        if not r.get("work_matched_error_m"):
            continue
        diff = r["work_matched_error_m"] - r["interior_error_m"]
        thr = cut * r["resolution_threshold_m"]
        if diff > thr:
            i += 1
        elif -diff > thr:
            f += 1
        else:
            u += 1
    return {"interior": i, "fixed": f, "unresolved": u, "lead": lead(i, f),
            "margin": i - f,
            "median_rho": float(np.median(rhos)) if rhos else None}


def baseline(key, beta, ha, level):
    p = METRICS / f"r19_equal_total_work_{key}_{beta_tag(beta)}.json"
    if not p.exists():
        return None
    return tally(json.loads(p.read_text(encoding="utf-8"))["rows"], ha, level)


def collect() -> dict:
    out: dict = {}
    for key, label in BLOCKS:
        cells = {}
        for beta in BETAS:
            p = METRICS / f"r63_ladder_uncapped_{key}_{beta_tag(beta)}.json"
            if not p.exists():
                continue
            rows = json.loads(p.read_text(encoding="utf-8"))["rows"]
            ha = apolune_map(key, beta)
            per_level = {}
            for lv in LEVELS:
                cuts = {f"{c:g}": tally(rows, ha, lv, c) for c in CUTS}
                per_level[f"{lv:.0f}"] = {
                    "level_consistent": tally(rows, ha, lv),
                    "tight_level": baseline(key, beta, ha, lv),
                    "resolution_cut_retally": cuts,
                    "cut_sensitive":
                        len({v["lead"] for v in cuts.values()}) > 1,
                }
            margins = [per_level[f"{lv:.0f}"]["level_consistent"]["margin"]
                       for lv in LEVELS]
            cells[f"{beta:.2f}"] = {
                "beta": beta, "levels": per_level, "margins": margins,
                "monotone_non_decreasing": all(
                    b >= a for a, b in zip(margins, margins[1:])),
                "interior_leads_from_600": all(
                    per_level[f"{lv:.0f}"]["level_consistent"]["lead"]
                    == "interior" for lv in LEVELS[1:]),
            }
        if cells:
            out[key] = {"block": label, "cells": cells}
    return out


def adjudicate(data: dict) -> dict:
    rows = {k: v["cells"].get("1.00") for k, v in data.items()}
    if not all(rows.values()) or len(rows) < 2:
        return {"class": None,
                "reason": "the beta = 1.00 row is not complete on both blocks"}
    mono = all(c["monotone_non_decreasing"] for c in rows.values())
    from600 = all(c["interior_leads_from_600"] for c in rows.values())
    reversal = any(
        c["levels"][f"{lv:.0f}"]["level_consistent"]["lead"] == "fixed"
        for c in rows.values() for lv in (1200.0, 2400.0))
    if reversal:
        cls, why = "C_reversal", (
            "the constant degree leads again at 1200 or 2400 km in at least "
            "one block")
    elif mono and from600:
        cls, why = "A_widening", (
            "the margin is non-decreasing across all four levels in both "
            "blocks and the interior member leads at 600 km and above")
    else:
        cls, why = "B_saturating", (
            "the interior member leads at 600 km and above but the margin is "
            "not monotone above 600 km in at least one block")
    return {"class": cls, "reason": why,
            "margins_at_beta_1": {k: c["margins"] for k, c in rows.items()},
            "levels_km": list(LEVELS)}


def table(data: dict) -> str:
    lines = ["\\begin{tabular}{@{}l c r r r r r@{}}", "\\toprule",
             "block & $\\beta$ & $h_a$ (km) & match & interior & fixed & "
             "median $\\rho$ \\\\", "\\midrule"]
    for n, (key, d) in enumerate(data.items()):
        if n:
            lines.append("\\midrule")
        first = True
        for tag, c in d["cells"].items():
            for lv in LEVELS:
                v = c["levels"][f"{lv:.0f}"]
                for which, name in (("tight_level", "tight level"),
                                    ("level_consistent", "scoring tol.")):
                    t = v[which]
                    if t is None:
                        continue
                    head = d["block"] if first else ""
                    bcol = (f"{c['beta']:.2f}" if lv == LEVELS[0]
                            and name == "tight level" else "")
                    lcol = f"{lv:.0f}" if name == "tight level" else ""
                    mark = ("\\,$^{\\dagger}$" if name == "scoring tol."
                            and v["cut_sensitive"] else "")
                    rho = ("--" if t["median_rho"] is None
                           else f"{t['median_rho']:.2f}")
                    lines.append(
                        f"{head} & {bcol} & {lcol} & {name}{mark} & "
                        f"{t['interior']} & {t['fixed']} & {rho} \\\\")
                    first = False
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    data = collect()
    if not data:
        print("[r63-tables] no records yet")
        return 1
    verdict = adjudicate(data)
    (METRICS / "r63_ladder_uncapped_table.tex").write_text(
        table(data), encoding="utf-8")
    (METRICS / "r63_manuscript_descriptives.json").write_text(
        json.dumps({"blocks": data, "outcome": verdict}, indent=2),
        encoding="utf-8")
    print("[written] r63_ladder_uncapped_table.tex, "
          "r63_manuscript_descriptives.json")
    for key, d in data.items():
        for tag, c in d["cells"].items():
            cells = " | ".join(
                f"{lv:.0f}km {c['levels'][f'{lv:.0f}']['level_consistent']['interior']}"
                f"-{c['levels'][f'{lv:.0f}']['level_consistent']['fixed']}"
                f"{'*' if c['levels'][f'{lv:.0f}']['cut_sensitive'] else ''}"
                for lv in LEVELS)
            print(f"  {d['block']:<24} b={tag}  {cells}   margins "
                  f"{c['margins']}")
    print(f"  OUTCOME: {verdict['class']} -- {verdict['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
