"""How much of the paper's tallies depends on where the resolution cut sits.

Every count the manuscript quotes for the budget comparisons follows the same
rule: a comparison is resolved when the error difference clears the summed
truth-inclusive envelope, that is when

    M_res = |E_a - E_b| / (env_a + env_b) > 1.

The cut at one is a convention, not a theorem, and a referee is entitled to ask
whether the counts sit on a cliff. This answers that from the archived records
alone. Nothing is propagated: every comparison already carries both errors and
its summed threshold, so the same rows are simply re-tallied at other cuts, and
the distribution of M_res itself is reported so the reader can see where the
comparisons actually lie rather than inferring it from three tallies.

A count that barely moves between cuts of 0.5 and 2 is not resting on the
convention. A count that swings is, and this table says which is which.

Usage:  python rev27_threshold_sensitivity.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CUTS = (0.5, 1.0, 2.0)


def records() -> list[tuple[str, float, Path]]:
    out = []
    for p in sorted(METRICS.glob("r19_equal_total_work_*.json")):
        m = re.fullmatch(r"r19_equal_total_work_([AB])(?:_beta_(\d+\.\d+))?\.json",
                         p.name)
        if not m:
            continue
        out.append((m.group(1), float(m.group(2) or 1.00), p))
    return sorted(out, key=lambda t: (t[0], -t[1]))


def m_res_values(path: Path) -> list[float]:
    d = json.loads(path.read_text(encoding="utf-8"))
    vals = []
    for r in d.get("rows", []):
        thr = r.get("resolution_threshold_m")
        ea, eb = r.get("interior_error_m"), r.get("work_matched_error_m")
        if not thr or ea is None or eb is None:
            continue
        vals.append((abs(eb - ea) / thr, "interior" if eb > ea else "fixed"))
    return vals


def tally(vals, cut: float) -> tuple[int, int, int]:
    interior = sum(1 for m, w in vals if m > cut and w == "interior")
    fixed = sum(1 for m, w in vals if m > cut and w == "fixed")
    return interior, fixed, len(vals) - interior - fixed


def main() -> int:
    blocks, summary = [], {}
    for design, beta, path in records():
        vals = m_res_values(path)
        if not vals:
            continue
        m = np.array([v for v, _ in vals])
        row = {"design": design, "beta": beta, "comparisons": len(vals),
               "m_res_median": float(np.median(m)),
               "m_res_p25": float(np.percentile(m, 25)),
               "m_res_p75": float(np.percentile(m, 75)),
               "near_cut_fraction": float(np.mean((m > 0.5) & (m < 2.0))),
               "tallies": {}}
        cells = []
        for cut in CUTS:
            i, f, u = tally(vals, cut)
            row["tallies"][f"cut_{cut}"] = {"interior": i, "fixed": f,
                                            "unresolved": u}
            cells.append(f"{i}--{f}")
        blocks.append((design, beta, row, cells))
        summary[f"{design}_beta_{beta:.2f}"] = row

    # Design label as a rotated cell spanning its own rows, so the group
    # boundary is visible rather than inferred from where the letter stops.
    lines = ["\\begin{tabular}{@{}c r r r r r@{}}", "\\toprule",
             "& $\\beta$ & $M_{\\mathrm{res}}>0.5$ & "
             "$M_{\\mathrm{res}}>1$ & $M_{\\mathrm{res}}>2$ & "
             "median $M_{\\mathrm{res}}$ \\\\",
             "\\midrule"]
    by_design: dict[str, list[str]] = {}
    for design, beta, row, cells in blocks:
        by_design.setdefault(design, []).append(
            f"{beta:.2f} & {cells[0]} & {cells[1]} & {cells[2]} & "
            f"{row['m_res_median']:.2f} \\\\")
    for i, (design, rows) in enumerate(by_design.items()):
        if i:
            lines.append("\\midrule")
        label = (f"\\multirow{{{len(rows)}}}{{*}}"
                 f"{{\\rotatebox[origin=c]{{90}}{{\\emph{{Design {design}}}}}}}")
        lines.append(f"{label} & {rows[0]}")
        lines += [f" & {r}" for r in rows[1:]]
    lines += ["\\bottomrule", "\\end{tabular}"]
    out = METRICS / "r27_threshold_sensitivity_table.tex"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    rec = METRICS / "r27_threshold_sensitivity.json"
    rec.write_text(json.dumps({
        "schema": "r27_threshold_sensitivity_v1",
        "what_this_is": (
            "the realized-work comparisons re-tallied at three resolution "
            "cuts, computed from the archived errors and envelopes; nothing "
            "is propagated and no error is recomputed"),
        "cuts": list(CUTS),
        "rows": summary}, indent=2), encoding="utf-8")

    print(f"[written] {out.name}, {rec.name}")
    print(f"{'design/beta':>14}  {'>0.5':>9} {'>1':>9} {'>2':>9}   "
          f"{'median':>7} {'near cut':>8}")
    for design, beta, row, cells in blocks:
        print(f"{design}@{beta:.2f}".rjust(14) +
              f"  {cells[0]:>9} {cells[1]:>9} {cells[2]:>9}   "
              f"{row['m_res_median']:>7.2f} {row['near_cut_fraction']:>7.0%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
