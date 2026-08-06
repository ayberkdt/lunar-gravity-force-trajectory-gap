"""The equal-work table with the post-hoc beta = 0.62 midpoint folded in.

Why this is not an edit to rev19_tables.py. That generator and its output,
r19_equal_work_table.tex, are both sha256-pinned in the R19 and R25 manifests.
Editing either breaks two sealed integrity gates for a row that is not even
pre-registered. So the pinned generator is imported and used unchanged --- every
pre-registered row in this table is still computed by the pinned code, on the
pinned records --- and only the assembly is new.

What is new is one budget and one mark. The beta = 0.62 rows carry a dagger and
the caption says what the dagger means, because the R28 amendment commits to
labelling the midpoint wherever it appears rather than letting it read as one
more grid point. The pre-registered bracket is not restated here; it is the main
text's job and the main text keeps it.

The midpoint is placed in budget order rather than appended, because a ladder
read out of order is harder to check than a ladder with a marked row in it.

Usage:  python rev28_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev19_tables as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

# descending budget order, with the post-hoc midpoint in its place
BETAS = (1.50, 1.25, 1.00, 0.75, 0.62, 0.50)
POST_HOC = (0.62,)
MARK = "\\textsuperscript{$\\dagger$}"


def beta_cell(beta: float) -> str:
    return f"{beta:.2f}" + (MARK if beta in POST_HOC else "")


def table() -> str:
    """rev19_tables.table() with one extra budget and a mark on it."""
    blocks = []
    for design in ("A", "B"):
        rows = []
        for beta in BETAS:
            d = r19.load(design, beta)
            pc = r19.per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            a = s["achieved_work_ratio"]
            rows.append(
                f"{beta_cell(beta)} & nominal per call & {pc['interior']} & "
                f"{pc['fixed']} & {pc['unresolved']} & "
                f"{pc['median_of_ratios']:.2f} & target-matched by "
                f"construction \\\\")
            rows.append(
                f" & realized total & {s['resolved_interior_wins']} & "
                f"{s['resolved_fixed_wins']} & {s['unresolved']} & "
                f"{s['median_rho']:.2f} & "
                f"${a['median']:.3f}$ $[{a['min']:.3f},{a['max']:.3f}]$ \\\\")
        if rows:
            blocks.append((design, rows))

    lines = ["\\begin{tabular}{@{}c l l r r r r r@{}}", "\\toprule",
             "& $\\beta$ & budget held equal & interior & fixed & unres. & "
             "median $\\rho$ & work match \\\\",
             "\\midrule"]
    for i, (design, rows) in enumerate(blocks):
        if i:
            lines.append("\\midrule")
        label = (f"\\multirow{{{len(rows)}}}{{*}}"
                 f"{{\\rotatebox[origin=c]{{90}}{{\\emph{{Design {design}}}}}}}")
        lines.append(f"{label} & {rows[0]}")
        lines += [f" & {r}" for r in rows[1:]]
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def midpoint_table() -> str:
    """The post-hoc budget alone, for the supplement.

    The main text carries the pre-registered ladder and nothing else, so the
    midpoint needs a table of its own rather than two marked rows inside the
    registered one. Same columns, so the two read against each other directly.
    """
    lines = ["\\begin{tabular}{@{}c l r r r r r@{}}", "\\toprule",
             "Design & budget held equal & interior & fixed & unres. & "
             "median $\\rho$ & work match \\\\",
             "\\midrule"]
    for i, design in enumerate(("A", "B")):
        d = r19.load(design, 0.62)
        pc = r19.per_call(design, 0.62)
        if not d or not pc:
            continue
        if i:
            lines.append("\\midrule")
        s = d["summary"]
        a = s["achieved_work_ratio"]
        lines.append(
            f"\\multirow{{2}}{{*}}{{{design}}} & nominal per call & "
            f"{pc['interior']} & {pc['fixed']} & {pc['unresolved']} & "
            f"{pc['median_of_ratios']:.2f} & target-matched by construction \\\\")
        lines.append(
            f" & realized total & {s['resolved_interior_wins']} & "
            f"{s['resolved_fixed_wins']} & {s['unresolved']} & "
            f"{s['median_rho']:.2f} & "
            f"${a['median']:.3f}$ $[{a['min']:.3f},{a['max']:.3f}]$ \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    out: dict = {}
    for design in ("A", "B"):
        for beta in BETAS:
            d = r19.load(design, beta)
            pc = r19.per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            out[f"{design}_{r19.beta_tag(beta)}"] = {
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
    (METRICS / "r28_equal_work_table.tex").write_text(table(), encoding="utf-8")
    (METRICS / "r28_midpoint_table.tex").write_text(midpoint_table(),
                                                    encoding="utf-8")
    (METRICS / "r28_manuscript_descriptives.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("[written] r28_equal_work_table.tex, r28_midpoint_table.tex, "
          "r28_manuscript_descriptives.json")
    for key, v in out.items():
        flag = "" if v["pre_registered"] else "  [POST HOC]"
        print(f"  {key}{flag}: realized {v['interior_wins']}/{v['fixed_wins']}/"
              f"{v['unresolved']} of {v['orbits']}, rho "
              f"{v['median_rho_realized']:.3f}, work "
              f"{v['achieved_work_ratio']['median']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
