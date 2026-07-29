"""Budget-robustness table: the interpolation path at three declared budgets.

Reads every r18_span_sweep_{design}_beta_{beta}.json present and emits a single
table showing where the optimum sits, how much realized work it costs, and how
many comparisons the resolution rule can still decide as the budget grows.

Usage:  python rev18_budget_table.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
ORDER = [("A", "0.50"), ("B", "0.50"), ("A", "1.00"), ("B", "1.00"),
         ("A", "1.50"), ("B", "1.50")]


def load(design: str, beta: str):
    p = METRICS / f"r18_span_sweep_{design}_beta_{beta}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def main() -> int:
    lines = ["\\begin{tabular}{@{}l l r c r r r r@{}}", "\\toprule",
             "$\\beta$ & design & $n$ & best $k$ & span & error [m] & work & "
             "resolved \\\\",
             " & & & & median & median & median & vs.\\ both \\\\",
             "\\midrule"]
    desc = {}
    prev_beta = None
    for design, beta in ORDER:
        d = load(design, beta)
        if not d:
            continue
        s = d["summary"]
        by = s["by_k"]
        counts = {k: sum(1 for r in d["rows"] if r["best_k"] == k) for k in by}
        best = max(counts, key=counts.get)
        e = by[best]
        if beta != prev_beta:
            lines.append("\\addlinespace[2pt]")
            prev_beta = beta
        lines.append(
            f"{beta} & {design} & {s['orbits']} & {best} & "
            f"{e['median_span']:.2f} & {e['median_error_m']:.3f} & "
            f"$\\times{e['median_total_work_ratio_vs_constant']:.3f}$ & "
            f"{s['interior_best_resolved_against_both']} \\\\")
        desc[f"{design}_{beta}"] = {
            "orbits": s["orbits"], "best_k": best,
            "best_k_count": counts[best],
            "median_span": e["median_span"],
            "median_error_m": e["median_error_m"],
            "median_work_ratio": e["median_total_work_ratio_vs_constant"],
            "constant_error_m": by["0.00"]["median_error_m"],
            "radial_error_m": by["1.00"]["median_error_m"],
            "interior_best": s["orbits_with_interior_best"],
            "resolved_vs_both": s["interior_best_resolved_against_both"],
            "location_resolved": s["interior_best_location_resolved"]}
    lines += ["\\bottomrule", "\\end{tabular}"]
    (METRICS / "r18_budget_table.tex").write_text("\n".join(lines) + "\n",
                                                  encoding="utf-8")
    (METRICS / "r18_budget_descriptives.json").write_text(
        json.dumps(desc, indent=2), encoding="utf-8")
    print("[written] r18_budget_table.tex, r18_budget_descriptives.json")
    for k, v in desc.items():
        print(f"  {k}: best k={v['best_k']} on {v['best_k_count']}/"
              f"{v['orbits']}, err {v['median_error_m']:.3f} vs const "
              f"{v['constant_error_m']:.3f}, work "
              f"x{v['median_work_ratio']:.3f}, resolved-vs-both "
              f"{v['resolved_vs_both']}, location {v['location_resolved']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
