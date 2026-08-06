"""R29 tables: the design-C ladder, and the three designs side by side.

The pinned generator rev19_tables.py is imported and used, not edited: it and
its output are sha256-pinned in the R19 and R25 manifests, and every quantity
that appears for designs A and B is still computed by that pinned code on the
pinned records. What is new here is design C and the assembly.

One accessor cannot be reused. rev19_tables.load leaves beta = 1 unsuffixed,
because that is where the R19 archive lives; design C writes beta = 1 with an
explicit suffix so that its records stay outside the subtree the R19 manifest
claims. The loader below knows that, and nothing else differs.

The post-hoc midpoint carries its dagger on design C exactly as it does on
designs A and B. The three-design table is restricted to the grid, because the
bracket statement is about the grid and a table that mixes the two invites the
reader to average them.

Usage:  python rev29_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev19_tables as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BETAS = (2.00, 1.50, 1.25, 1.00, 0.75, 0.62, 0.50)   # descending
GRID = (1.50, 1.25, 1.00, 0.75, 0.50)
POST_HOC = (0.62,)
MARK = "\\textsuperscript{$\\dagger$}"


def load(design: str, beta: float):
    if design == "C":
        p = METRICS / f"r19_equal_total_work_C_beta_{beta:.2f}.json"
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    return r19.load(design, beta)


def beta_cell(beta: float) -> str:
    return f"{beta:.2f}" + (MARK if beta in POST_HOC else "")


def rows_for(design: str, betas) -> list[str]:
    rows = []
    for beta in betas:
        d = load(design, beta)
        pc = r19.per_call(design, beta)
        if not d or not pc:
            continue
        s = d["summary"]
        a = s["achieved_work_ratio"]
        rows.append(f"{beta_cell(beta)} & nominal per call & {pc['interior']} & "
                    f"{pc['fixed']} & {pc['unresolved']} & "
                    f"{pc['median_of_ratios']:.2f} & target-matched by "
                    f"construction \\\\")
        rows.append(f" & realized total & {s['resolved_interior_wins']} & "
                    f"{s['resolved_fixed_wins']} & {s['unresolved']} & "
                    f"{s['median_rho']:.2f} & "
                    f"${a['median']:.3f}$ $[{a['min']:.3f},{a['max']:.3f}]$ \\\\")
    return rows


def designC_table() -> str:
    """Design C alone, in the columns designs A and B are reported in."""
    rows = rows_for("C", BETAS)
    lines = ["\\begin{tabular}{@{}l l r r r r r@{}}", "\\toprule",
             "$\\beta$ & budget held equal & interior & fixed & unres. & "
             "median $\\rho$ & work match \\\\", "\\midrule"]
    lines += rows
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def three_design_table() -> str:
    """The realized-work tally of all three designs on the frozen grid."""
    lines = ["\\begin{tabular}{@{}l r r r r r r r r r@{}}", "\\toprule",
             "& \\multicolumn{3}{c}{Design A} & \\multicolumn{3}{c}{Design B} & "
             "\\multicolumn{3}{c}{Design C} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
             "$\\beta$ & int. & fix. & $\\rho$ & int. & fix. & $\\rho$ & "
             "int. & fix. & $\\rho$ \\\\", "\\midrule"]
    for beta in GRID:
        cells = []
        any_row = False
        for design in ("A", "B", "C"):
            d = load(design, beta)
            if not d:
                cells += ["---", "---", "---"]
                continue
            any_row = True
            s = d["summary"]
            cells += [str(s["resolved_interior_wins"]),
                      str(s["resolved_fixed_wins"]),
                      f"{s['median_rho']:.2f}"]
        if any_row:
            lines.append(f"{beta:.2f} & " + " & ".join(cells) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    out: dict = {}
    for design in ("A", "B", "C"):
        for beta in BETAS:
            d = load(design, beta)
            pc = r19.per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            out[f"{design}_beta_{beta:.2f}"] = {
                "design": design, "beta": beta,
                "pre_registered": beta not in POST_HOC,
                "orbits": s["orbits"], "resolved": s["resolved"],
                "interior_wins": s["resolved_interior_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "unresolved": s["unresolved"],
                "median_rho_realized": s["median_rho"],
                "achieved_work_ratio": s["achieved_work_ratio"],
                "per_call": pc,
            }
    (METRICS / "r29_designC_table.tex").write_text(designC_table(),
                                                   encoding="utf-8")
    (METRICS / "r29_three_design_table.tex").write_text(three_design_table(),
                                                        encoding="utf-8")
    (METRICS / "r29_manuscript_descriptives.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("[written] r29_designC_table.tex, r29_three_design_table.tex, "
          "r29_manuscript_descriptives.json")
    for key, v in out.items():
        if v["design"] != "C":
            continue
        flag = "" if v["pre_registered"] else "  [POST HOC]"
        print(f"  {key}{flag}: realized {v['interior_wins']}/{v['fixed_wins']}/"
              f"{v['unresolved']} of {v['orbits']}, rho "
              f"{v['median_rho_realized']:.3f}, work "
              f"{v['achieved_work_ratio']['median']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
