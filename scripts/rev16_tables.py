"""Tables for the R16 cross-solution / cross-body transfer test.

Emits, into metrics/:
  r16_transfer_table.tex        main-text summary, one row per field
  r16_transfer_detail_table.tex supplement, empirical vs proxy degrees

Usage:  python rev16_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SRC = METRICS / "r16_multibody_calibration.json"

# Rows are grouped by role so the reproducibility claim and the transfer claim
# stay visually separate; within a group the paper's own field comes first.
ORDER = ["JGGRX_1800F", "GRGM1200A", "GGGRX_1200L",
         "SHGJ180U", "JGMESS_160A", "JGMRO120D", "GOCO05c", "EGM96"]

PRETTY = {"JGGRX_1800F": "JGGRX\\_1800F", "GRGM1200A": "GRGM1200A",
          "GGGRX_1200L": "GGGRX\\_1200L", "SHGJ180U": "SHGJ180U",
          "JGMESS_160A": "JGMESS\\_160A", "JGMRO120D": "JGMRO120D",
          "GOCO05c": "GOCO05c", "EGM96": "EGM96"}


def load():
    d = json.loads(SRC.read_text(encoding="utf-8"))
    return d, {f["key"]: f for f in d["fields"]}


def summary_table() -> str:
    _, by_key = load()
    lines = ["\\begin{tabular}{@{}l l r r r r r@{}}", "\\toprule",
             "Field & Body & $N_{\\max}$ & $\\hat p_{\\mathrm{spec}}$ & "
             "$p_{\\mathrm{fit}}$ & $\\Delta_{\\mathrm{fit}}$ & "
             "$\\Delta_{p=2}$ \\\\",
             "\\midrule",
             "\\multicolumn{7}{@{}l}{\\emph{Same body, three independent "
             "solutions}}\\\\"]
    for key in ORDER:
        f = by_key[key]
        if key == "SHGJ180U":
            lines.append("\\midrule")
            lines.append("\\multicolumn{7}{@{}l}{\\emph{Other bodies}}\\\\")
        s = f["spectral_slope"]
        if f["degenerate"]:
            p_fit = "---"
            d_fit = "---"
        else:
            p_fit = f"{f['p_fit']:.3f}"
            d_fit = f"{f['rms_mismatch_p_fit']:.2f}"
        lines.append(
            f"{PRETTY[key]} & {f['body']} & {f['max_degree_in_file']} & "
            f"{s['p']:.3f} & {p_fit} & {d_fit} & "
            f"{f['rms_mismatch_p2']:.2f}\\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def detail_table() -> str:
    """Empirical and proxy degrees at four representative ratios."""
    d, by_key = load()
    ratios = d["band"]["ratios_R_over_r"]
    alts = d["band"]["moon_equivalent_altitudes_km"]
    picks = [alts.index(a) for a in (50.0, 100.0, 200.0, 300.0)]

    header = " & ".join(f"$R/r={ratios[i]:.4f}$" for i in picks)
    lines = ["\\begin{tabular}{@{}l l " + "r r " * len(picks) + "@{}}",
             "\\toprule",
             "Field & Body & " + " & ".join(
                 f"\\multicolumn{{2}}{{c}}{{{h}}}" for h in header.split(" & "))
             + "\\\\",
             "& & " + " & ".join(["$N^{\\mathrm{emp}}$", "proxy"] * len(picks))
             + "\\\\",
             "\\midrule"]
    for key in ORDER:
        f = by_key[key]
        cells = []
        for i in picks:
            row = f["criteria_rows"][i]
            proxy = "---" if f["degenerate"] else str(row["proxy_p_fit"])
            cells += [str(row["emp"]), proxy]
        lines.append(f"{PRETTY[key]} & {f['body']} & " + " & ".join(cells) + "\\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}; run rev16_multibody_calibration.py first")
    (METRICS / "r16_transfer_table.tex").write_text(summary_table(),
                                                    encoding="utf-8")
    (METRICS / "r16_transfer_detail_table.tex").write_text(detail_table(),
                                                           encoding="utf-8")
    print("[written] r16_transfer_table.tex, r16_transfer_detail_table.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
