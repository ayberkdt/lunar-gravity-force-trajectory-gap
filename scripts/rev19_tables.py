"""Table for the R19 equal-realized-work comparison.

Emits metrics/r19_equal_work_table.tex and the descriptives the manuscript
quotes, so no number is transcribed by hand.

Two things this reports that the first version did not:

Budget scale. The comparison runs at beta = 1, 0.75 and 0.5, and the table
carries all three because the verdict is not the same at all three. Showing
only the anchor would show the claim everywhere except where it turns over.
The midpoint was added last, to locate a sign change that the two original
budgets could only bracket; a budget whose records are absent is skipped
rather than faked, so this generator is safe to run mid-campaign.

A like-for-like ratio. The manuscript contrasted the realized-work median rho,
which is a median of per-orbit ratios, against a per-call "sixfold", which is a
ratio of population medians. Those are different statistics and their quotient
overstates what changing the comparator did. Both accountings are therefore
summarized here by the median of per-orbit ratios, and the ratio-of-medians
figure is emitted alongside it under its own name so the two can never be
silently swapped again.

Usage:  python rev19_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BETAS = (1.50, 1.25, 1.00, 0.75, 0.50)
K = "0.50"


def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def _suffix(beta: float) -> str:
    return "" if abs(beta - 1.0) < 1e-12 else f"_{beta_tag(beta)}"


def load(design: str, beta: float):
    p = METRICS / f"r19_equal_total_work_{design}{_suffix(beta)}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def per_call(design: str, beta: float) -> dict | None:
    """The same comparison under the nominal per-call budget, for contrast."""
    p = METRICS / f"r18_span_sweep_{design}_{beta_tag(beta)}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    win = lose = unres = 0
    ratios, e_fix, e_int = [], [], []
    for r in d["rows"]:
        e0, ek = r["entries"]["0.00"], r["entries"].get(K, {})
        if not ek.get("error_m") or not e0.get("error_m"):
            continue
        thr = (e0.get("envelope_m") or 0.0) + (ek.get("envelope_m") or 0.0)
        diff = e0["error_m"] - ek["error_m"]
        if diff > thr:
            win += 1
        elif -diff > thr:
            lose += 1
        else:
            unres += 1
        ratios.append(e0["error_m"] / ek["error_m"])
        e_fix.append(e0["error_m"])
        e_int.append(ek["error_m"])
    if not ratios:
        return None
    return {
        "interior": win, "fixed": lose, "unresolved": unres,
        "median_of_ratios": float(np.median(ratios)),
        "ratio_of_medians": float(np.median(e_fix) / np.median(e_int)),
        "median_error_fixed_m": float(np.median(e_fix)),
        "median_error_interior_m": float(np.median(e_int)),
    }


def table() -> str:
    """Design label as a rotated spanning cell rather than a banner row.

    With four budgets and two accountings each, a design occupies eight rows.
    A full-width ``Design A`` row above them costs a line, breaks the rules of
    the table twice, and still leaves the reader to work out how far the group
    extends. A rotated label spanning the block says the same thing in the
    margin and lets the eye see the group boundary directly.
    """
    blocks = []
    for design in ("A", "B"):
        rows = []
        for beta in BETAS:
            d = load(design, beta)
            pc = per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            a = s["achieved_work_ratio"]
            rows.append(
                f"{beta:.2f} & nominal per call & {pc['interior']} & "
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


def main() -> int:
    out: dict = {}
    for design in ("A", "B"):
        for beta in BETAS:
            d = load(design, beta)
            pc = per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            rows = d["rows"]
            shift = [r["work_matched_degree"] - r["constant_endpoint_degree"]
                     for r in rows
                     if r.get("work_matched_degree")
                     and r.get("constant_endpoint_degree")]
            losers = sorted(r["sobol_index"] for r in rows
                            if r["winner"] == "fixed")
            out[f"{design}_{beta_tag(beta)}"] = {
                "design": design, "beta": beta,
                "orbits": s["orbits"], "resolved": s["resolved"],
                "interior_wins": s["resolved_interior_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "unresolved": s["unresolved"],
                "median_rho_realized": s["median_rho"],
                "achieved_work_ratio": s["achieved_work_ratio"],
                "comparator_degree_shift": {
                    "median": float(np.median(shift)),
                    "min": int(min(shift)), "max": int(max(shift))}
                if shift else None,
                "orbits_where_interior_loses": losers,
                "per_call": pc,
                "like_for_like_note": (
                    "median_of_ratios is the statistic comparable with "
                    "median_rho_realized; ratio_of_medians is the larger "
                    "figure the text previously quoted as 'sixfold'"),
            }
    (METRICS / "r19_equal_work_table.tex").write_text(table(), encoding="utf-8")
    (METRICS / "r19_manuscript_descriptives.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("[written] r19_equal_work_table.tex, "
          "r19_manuscript_descriptives.json")
    for key, v in out.items():
        pc = v["per_call"]
        print(f"  {key}: per-call {pc['interior']}/{pc['fixed']}/"
              f"{pc['unresolved']} (median-of-ratios "
              f"{pc['median_of_ratios']:.2f}, ratio-of-medians "
              f"{pc['ratio_of_medians']:.2f})")
        print(f"      realized {v['interior_wins']}/{v['fixed_wins']}/"
              f"{v['unresolved']}, median rho "
              f"{v['median_rho_realized']:.2f}, work ratio "
              f"{v['achieved_work_ratio']['median']:.3f}")
        print(f"      interior loses on: {v['orbits_where_interior_loses']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
