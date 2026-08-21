"""LaTeX table for the enlarged forced-variational panel (R37).

Two things the main text asserts and must be able to show: which panel orbits
the resolution rule leaves undecided, and which ones the calibration channel
places outside the band the original eight-orbit panel occupied. Both are read
from metrics/r37_panel_verdict.json, which is itself derived from the run's
own record.

Usage:  python rev37_table.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
VERDICT = METRICS / "r37_panel_verdict.json"
OUT = METRICS / "r37_panel_extension_table.tex"


def main() -> int:
    v = json.loads(VERDICT.read_text(encoding="utf-8"))
    panel, rec = v["panel"], v["record"]
    uncal = [c for c in v["outside_calibration_band"]
             if c["resolved"] and c.get("in_panel")]
    unres = [c for c in v["unresolved"]["orbits"] if c.get("in_panel")]

    lines = [r"\begin{tabular}{@{}l r r r r@{}}", r"\toprule",
             r"& orbits & resolved & sign agreement & calibration median \\",
             r"\midrule",
             rf"Panel (level {v['levels']['highest_complete_per_design']} "
             rf"per design) & {panel['orbits']} & {panel['resolved']} & "
             rf"{panel['sign_agreement']}/{panel['resolved']} & "
             rf"{panel['calibration_median']:.4f} \\",
             rf"Whole record & {rec['orbits']} & {rec['resolved']} & "
             rf"{rec['sign_agreement']}/{rec['resolved']} & "
             rf"{rec['calibration_median']:.4f} \\",
             r"\midrule",
             r"\multicolumn{5}{@{}l}{\emph{Resolved panel orbits outside the "
             r"$0.94$--$1.04$ calibration band}} \\"]
    if uncal:
        lines.append(r"\multicolumn{5}{@{}p{\linewidth}@{}}{" + ", ".join(
            rf"{c['design']}{c['sobol_index']:03d} "
            rf"($h_p={c['hp_km']:.0f}$~km, {c['calibration_ratio_fixed']:.3f})"
            for c in uncal) + r"} \\")
    else:
        lines.append(r"\multicolumn{5}{@{}p{\linewidth}@{}}{none} \\")
    lines += [r"\midrule",
              r"\multicolumn{5}{@{}l}{\emph{Panel orbits the resolution rule "
              r"leaves undecided}} \\"]
    if unres:
        lines.append(r"\multicolumn{5}{@{}p{\linewidth}@{}}{" + ", ".join(
            rf"{c['design']}{c['sobol_index']:03d} "
            rf"($h_p={c['hp_km']:.0f}$~km)" for c in unres) + r"} \\")
    else:
        lines.append(r"\multicolumn{5}{@{}p{\linewidth}@{}}{none} \\")
    lines += [r"\bottomrule", r"\end{tabular}"]

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {OUT.name}: panel {panel['orbits']} orbits, "
          f"{panel['sign_agreement']}/{panel['resolved']} resolved signs, "
          f"{len(uncal)} outside band, {len(unres)} undecided")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
