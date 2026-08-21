"""LaTeX table for the gradient-degree audit of the enlarged panel (R39).

The manuscript's claim is about a side of ratio one, so the table is built
around that: for each audited orbit it gives the predicted ratio with the
gradient at degree 120, the same ratio with the gradient at the orbit's own
reference degree, the relative change between them, and whether the side
changed. Band membership travels with every row, because the registration
forbids reporting a partial band as a band and band L is one of four.

Read from metrics/r57_gradient_degree_completion.json when that record exists
and from metrics/r39_gradient_degree_panel.json otherwise. R57 is R39 carried
forward and finished: it copies R39's solved rows byte for byte and adds the
three the R39 campaign ran out of compute for, so the completion record is a
superset and the table drawn from it supersedes the partial one. Nothing is
re-solved or re-scored here.

Usage:  python rev39_table.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
COMPLETION = METRICS / "r57_gradient_degree_completion.json"
PANEL = METRICS / "r39_gradient_degree_panel.json"
OUT = METRICS / "r39_gradient_degree_table.tex"

BAND_LABEL = {"L": "L, $h_p<50$", "M": "M, $50$--$100$", "H": "H, $>100$",
              "X-fragile": "X, ratio near 1", "X-extreme": "X, ratio extreme"}
# R39 printed bands in this order and ascending perilune inside each. The
# completion record appends its three rows at the end rather than in place,
# so the order is imposed here instead of taken from the file.
BAND_ORDER = ["L", "M", "H", "X-fragile", "X-extreme"]


def fmt(x: float) -> str:
    """Ratios span 1e-4 to 1e5; keep three significant figures either way."""
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e4 or a < 1e-2:
        m, e = f"{x:.2e}".split("e")
        return f"${m}\\times 10^{{{int(e)}}}$"
    return f"${x:.3g}$"


def main() -> int:
    record = COMPLETION if COMPLETION.exists() else PANEL
    if not record.exists():
        print(f"[abort] {PANEL.name} missing")
        return 2
    rec = json.loads(record.read_text(encoding="utf-8"))
    comp = sorted(rec["comparison"],
                  key=lambda c: (BAND_ORDER.index(c["band"]), c["hp_km"]))
    st = rec["summary"]["band_status"]

    lines = [r"\begin{tabular}{@{}l l r r r r r c@{}}", r"\toprule",
             (r"band & orbit & $h_p$ & $N_{\mathrm{ref}}$ & res. & "
              r"$\rho$ at $G_{120}$ & $\rho$ at $G_{N_{\mathrm{ref}}}$ & "
              r"$\Delta$ \\"),
             r"\midrule"]
    last = None
    for c in comp:
        band = BAND_LABEL.get(c["band"], c["band"])
        show = band if c["band"] != last else ""
        last = c["band"]
        lines.append(
            f"{show} & {c['design']}{c['sobol_index']:03d} & "
            f"{c['hp_km']:.1f} & {c['adopted_truth_degree']} & "
            f"{'y' if c['measured_resolved'] else 'n'} & "
            f"{fmt(c['ratio_gradient_120'])} & "
            f"{fmt(c['ratio_gradient_reference'])} & "
            f"${c['relative_change']*100:+.1f}\\%$ \\\\")
    lines.append(r"\midrule")

    s = rec["summary"]
    if st:
        partial = [b for b, v in st.items() if not v["complete"]]
        done = [f"{b} ({v['solved']}/{v['declared']})"
                for b, v in st.items() if v["complete"]]
        note = ("Bands complete: " + ", ".join(done) + ". "
                + ("Incomplete: "
                   + ", ".join(f"{b} ({st[b]['solved']}/{st[b]['declared']})"
                               for b in partial)
                   + ", reported as a partial band and never as the band."
                   if partial else "All declared bands complete."))
    else:
        # the completion record carries a list of finished bands instead of
        # the partial-run bookkeeping, because none of them is partial now
        note = ("All declared bands complete: "
                + ", ".join(sorted(s["bands_complete"])) + ".")
    lines.append(r"\multicolumn{8}{@{}p{\linewidth}@{}}{\footnotesize " + note
                 + r"} \\")
    by_band = s.get("abs_relative_change_by_band")
    if by_band:
        lines.append(
            r"\multicolumn{8}{@{}p{\linewidth}@{}}{\footnotesize "
            r"Median $|\Delta|$ by band: "
            + ", ".join(f"{b} {by_band[b]*100:.2g}\\%"
                        for b in BAND_ORDER if b in by_band)
            + r".} \\")
    lines.append(
        r"\multicolumn{8}{@{}p{\linewidth}@{}}{\footnotesize Side changes: "
        f"{s['side_changes_resolved']} of {s['resolved']} resolved, "
        f"{s['side_changes_unresolved']} among the unresolved. "
        r"$|\Delta|$ median "
        f"{s['abs_relative_change']['median']*100:.2g}\\%, p90 "
        f"{s['abs_relative_change']['p90']*100:.1f}\\%, max "
        f"{s['abs_relative_change']['max']*100:.1f}\\%." + r"} \\")
    lines += [r"\bottomrule", r"\end{tabular}"]

    OUT.write_text("% auto-generated by rev39_table.py -- do not edit by hand\n"
                   + "\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {OUT.name}: {len(comp)} orbits, "
          f"{s['side_changes_resolved']} resolved side changes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
