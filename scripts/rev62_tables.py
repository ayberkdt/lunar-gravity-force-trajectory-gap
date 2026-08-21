"""Tables and outcome adjudication for R62 (O54).

The campaign's question is about a *difference between two apolune levels*,
so every statistic here is computed per level and never pooled across them.
Pooling would average away exactly the contrast the ladder was built to
isolate.

The outcome class is decided by the rule frozen in r62_preregistration.json
and is computed here rather than read off by eye, so the class this campaign
returns is a function of the records and not of who is reading them.

Usage:  python rev62_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BLOCKS = [("RS1", "span_ladder_a", "Block A"),
          ("RS2", "span_ladder_b", "Block B")]
BETAS = (0.50, 0.75, 1.00)
LEVELS = (300.0, 600.0)
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
    """Resolved interior--fixed at one apolune level. With cut=None the
    record's own resolution is used; with a cut the comparison is re-decided
    at that multiple of the recorded envelope, nothing repropagated."""
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
            "median_rho": float(np.median(rhos)) if rhos else None}


def baseline_tally(key: str, beta: float, ha, level):
    p = METRICS / f"r19_equal_total_work_{key}_{beta_tag(beta)}.json"
    if not p.exists():
        return None
    rows = json.loads(p.read_text(encoding="utf-8"))["rows"]
    return tally(rows, ha, level)


def collect() -> dict:
    out: dict = {}
    for key, pop, label in BLOCKS:
        cells = {}
        for beta in BETAS:
            p = METRICS / f"r62_ladder_interior_{key}_{beta_tag(beta)}.json"
            if not p.exists():
                continue
            rows = json.loads(p.read_text(encoding="utf-8"))["rows"]
            ha = apolune_map(key, beta)
            per_level = {}
            for lv in LEVELS:
                cuts = {f"{c:g}": tally(rows, ha, lv, c) for c in CUTS}
                base = baseline_tally(key, beta, ha, lv)
                now = tally(rows, ha, lv)
                per_level[f"{lv:.0f}"] = {
                    "level_consistent": now,
                    "tight_level": base,
                    "resolution_cut_retally": cuts,
                    "cut_sensitive": len({v["lead"] for v in cuts.values()}) > 1,
                }
            a, b = per_level[f"{LEVELS[0]:.0f}"], per_level[f"{LEVELS[1]:.0f}"]
            cells[f"{beta:.2f}"] = {
                "beta": beta,
                "levels": per_level,
                "leads_differ": (a["level_consistent"]["lead"]
                                 != b["level_consistent"]["lead"]),
                "margin_300": (a["level_consistent"]["interior"]
                               - a["level_consistent"]["fixed"]),
                "margin_600": (b["level_consistent"]["interior"]
                               - b["level_consistent"]["fixed"]),
            }
        if cells:
            out[key] = {"block": label, "population": pop, "cells": cells}
    return out


def adjudicate(data: dict) -> dict:
    """The frozen rule, applied to the beta = 0.50 row of both blocks."""
    rows = {k: v["cells"].get("0.50") for k, v in data.items()}
    if not all(rows.values()) or len(rows) < 2:
        return {"class": None,
                "reason": "the beta = 0.50 row is not complete on both blocks"}
    per_block = {}
    for k, c in rows.items():
        lo = c["levels"]["300"]
        hi = c["levels"]["600"]
        per_block[k] = {
            "lead_300": lo["level_consistent"]["lead"],
            "lead_600": hi["level_consistent"]["lead"],
            "cut_sensitive": lo["cut_sensitive"] or hi["cut_sensitive"],
            "margin_moves_toward_interior": c["margin_600"] > c["margin_300"],
        }
    textbook = all(v["lead_300"] == "fixed" and v["lead_600"] == "interior"
                   and not v["cut_sensitive"] for v in per_block.values())
    same_side = all(v["lead_300"] == v["lead_600"] for v in per_block.values())
    toward = all(v["margin_moves_toward_interior"]
                 for v in per_block.values())
    if textbook:
        cls, why = "A_shift_preserved", (
            "both blocks lead the constant degree at 300 km and the interior "
            "member at 600 km, and neither level is cut-sensitive")
    elif same_side:
        cls, why = "C_shift_absent", (
            "the leading side is the same at both levels in both blocks")
    else:
        cls, why = "B_shift_attenuated", (
            "the step moves the tally toward the interior member"
            + ("" if toward else " in only one block")
            + ", but the leading-side difference is missing or cut-sensitive "
              "in at least one block")
    return {"class": cls, "reason": why, "per_block": per_block}


def table(data: dict) -> str:
    lines = ["\\begin{tabular}{@{}l c r r r r r@{}}", "\\toprule",
             "block & $\\beta$ & $h_a$ (km) & match & interior & fixed & "
             "median $\\rho$ \\\\", "\\midrule"]
    for i, (key, d) in enumerate(data.items()):
        if i:
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
        print("[r62-tables] no records yet")
        return 1
    verdict = adjudicate(data)
    (METRICS / "r62_ladder_interior_table.tex").write_text(
        table(data), encoding="utf-8")
    (METRICS / "r62_manuscript_descriptives.json").write_text(
        json.dumps({"blocks": data, "outcome": verdict}, indent=2),
        encoding="utf-8")
    print("[written] r62_ladder_interior_table.tex, "
          "r62_manuscript_descriptives.json")
    for key, d in data.items():
        for tag, c in d["cells"].items():
            lo = c["levels"]["300"]["level_consistent"]
            hi = c["levels"]["600"]["level_consistent"]
            lo0 = c["levels"]["300"]["tight_level"]
            hi0 = c["levels"]["600"]["tight_level"]
            print(f"  {d['block']} beta={tag}  "
                  f"300km {lo0['interior']}-{lo0['fixed']} -> "
                  f"{lo['interior']}-{lo['fixed']} ({lo['lead']})   "
                  f"600km {hi0['interior']}-{hi0['fixed']} -> "
                  f"{hi['interior']}-{hi['fixed']} ({hi['lead']})")
    print(f"  OUTCOME: {verdict['class']} -- {verdict['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
