"""Tables for R68 (O59/O60): the full-design measured-time comparison.

Three tables, and each is deliberately narrow in what it licenses.

  r68_measured_time_summary_table  the registered primary result, one row per
                                   arm and design. This is the table that
                                   answers the runtime-practicality objection,
                                   so it carries the achieved time ratio next
                                   to the tally rather than in a footnote.
  r68_band_table                   every cell outside the prescribed band,
                                   with the side it missed on. Registered as
                                   reportable, not as excludable.
  r68_apolune_posthoc_table        the apolune stratification of the interior
                                   arm. Post-hoc: it was not registered, it is
                                   a diagnostic for why a perilune-stratified
                                   fourteen-orbit panel and a full design read
                                   differently, and the table caption has to
                                   say so wherever it is used.

Counts follow the registration: a cell outside the time-match band keeps its
nearest integer match and stays in the primary tally. The without-misses
tally is in r68_band_admissibility.json as a sensitivity and is not printed
as a result here.

Usage:  python rev68_tables.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rev68_timing_full as r68            # noqa: E402

ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"

SUMMARY = METRICS / "r68_measured_time_summary_table.tex"
BAND = METRICS / "r68_band_table.tex"
APOLUNE = METRICS / "r68_apolune_posthoc_table.tex"

APOLUNE_BANDS = ((0, 300), (300, 500), (500, 700))


def rows_of(arm: str) -> list:
    return json.loads(r68.out_path(arm).read_text(encoding="utf-8"))["rows"]


def fmt(x, digits=2):
    return "--" if x is None else f"{x:.{digits}f}"


def summary_table() -> None:
    lines = [r"\begin{tabular}{@{}l l r r r r r@{}}", r"\toprule",
             r"policy & design & orbits & resolved & variable--fixed & "
             r"median $\rho$ & median $T$ ratio \\", r"\midrule"]
    for arm, label in (("endpoint", r"radial endpoint, $k=1$"),
                       ("interior", r"interior member, $k=0.5$")):
        member = "radial" if arm == "endpoint" else "interior"
        rows = rows_of(arm)
        for j, design in enumerate(r68.DESIGNS):
            d = [r for r in rows if r["design"] == design]
            res = [r for r in d if r["resolved"]]
            rho = sorted(r["rho_constant_over_member"] for r in d
                         if r["rho_constant_over_member"])
            tr = [r["achieved_time_ratio"] for r in d]
            wins = sum(1 for r in res if r["winner"] == member)
            fix = sum(1 for r in res if r["winner"] == "constant")
            lines.append(
                f"{label if j == 0 else ''} & {design} & {len(d)} & "
                f"{len(res)} & {wins}--{fix} & "
                f"{fmt(float(np.median(rho)) if rho else None, 3)} & "
                f"{fmt(float(np.median(tr)), 3)} \\\\")
        if arm == "endpoint":
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {SUMMARY.name}")


def band_table() -> None:
    lines = [r"\begin{tabular}{@{}l r r r c l@{}}", r"\toprule",
             r"orbit & $h_p$ (km) & $N$ & $T$ ratio & side & verdict \\",
             r"\midrule"]
    n = 0
    for arm in ("endpoint", "interior"):
        member = "radial" if arm == "endpoint" else "interior"
        miss = [r for r in rows_of(arm) if r["timing_match_miss"]]
        if not miss:
            continue
        lines.append(r"\multicolumn{6}{@{}l}{" +
                     (r"radial endpoint, $k=1$" if arm == "endpoint"
                      else r"interior member, $k=0.5$") + r"} \\")
        for r in sorted(miss, key=lambda r: (r["design"], r["sobol_index"])):
            side = "high" if r["achieved_time_ratio"] > r68.BAND_HI else "low"
            verdict = (r["winner"] if r["resolved"] else "unresolved")
            verdict = {"radial": "radial", "interior": "interior",
                       "constant": "fixed"}.get(verdict, verdict)
            lines.append(
                f"{r['design']}{r['sobol_index']:03d} & {r['hp_km']:.1f} & "
                f"{r['comparator_degree']} & "
                f"{r['achieved_time_ratio']:.3f} & {side} & {verdict} \\\\")
            n += 1
    lines += [r"\bottomrule", r"\end{tabular}"]
    BAND.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {BAND.name} ({n} cells)")


def apolune_table() -> None:
    rows = rows_of("interior")
    lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
             r"$h_a$ (km) & orbits & resolved & interior & fixed & "
             r"median $\rho$ \\", r"\midrule"]
    for lo, hi in APOLUNE_BANDS:
        d = [r for r in rows if lo <= (r["ha_km"] or 0.0) < hi]
        res = [r for r in d if r["resolved"]]
        rho = sorted(r["rho_constant_over_member"] for r in d
                     if r["rho_constant_over_member"])
        lines.append(
            f"{lo}--{hi} & {len(d)} & {len(res)} & "
            f"{sum(1 for r in res if r['winner'] == 'interior')} & "
            f"{sum(1 for r in res if r['winner'] == 'constant')} & "
            f"{fmt(float(np.median(rho)) if rho else None, 2)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    APOLUNE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {APOLUNE.name} (post-hoc; caption must say so)")


def main() -> int:
    summary_table()
    band_table()
    apolune_table()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
