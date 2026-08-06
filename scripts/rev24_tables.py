"""LaTeX tables and manuscript numbers for the R24 campaigns.

Every number the manuscript quotes for R24 is emitted here from the campaign
records, so the prose and the tables cannot drift apart by hand-transcription.

Usage:  python rev24_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"


def _fmt(x, digits=2):
    if x is None:
        return "---"
    return f"{x:.{digits}f}"


def oracle_ultra_table() -> int:
    src = METRICS / "r24_oracle_ultra.json"
    if not src.exists():
        print("[r24 tables] r24_oracle_ultra.json missing; skipped")
        return 1
    d = json.loads(src.read_text(encoding="utf-8"))
    s, shrink = d["summary"], d["envelope_shrink_factor"]

    lines = [
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Comparator & \multicolumn{2}{c}{resolved} & interior & comparator"
        r" & unres. & flips \\",
        r"\cmidrule(lr){2-3}",
        r" & before & after & \multicolumn{2}{c}{after} & after & \\",
        r"\midrule",
    ]
    for tag, name in (("vs_oracle", "Best ladder member (oracle)"),
                      ("vs_sat", "Budget-saturating degree")):
        t = s[tag]
        unres = t["comparisons"] - t["resolved_after"]
        lines.append(
            f"{name} & {t['resolved_before']} & {t['resolved_after']} & "
            f"{t['interior_wins']} & {t['comparator_wins']} & {unres} & "
            f"{t['verdict_flips']} \\\\")
    lines += [
        r"\midrule",
        r"\multicolumn{7}{l}{\footnotesize Interior-member envelope, "
        r"tight-to-tighter over tighter-to-ultra: median "
        f"${_fmt(shrink['median'], 1)}$, range "
        f"${_fmt(shrink['min'], 2)}$--${_fmt(shrink['max'], 1)}$"
        r"} \\",
        r"\multicolumn{7}{l}{\footnotesize Panel: all "
        f"{s['vs_oracle']['comparisons']} comparisons, "
        f"{d['panel_completeness']['missing']} missing"
        r"} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    out = METRICS / "r24_oracle_ultra_table.tex"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[r24 tables] wrote {out.name}")
    return 0


def bin_control_table() -> int:
    src = METRICS / "r24_bin_control.json"
    if not src.exists():
        print("[r24 tables] r24_bin_control.json missing; skipped")
        return 1
    d = json.loads(src.read_text(encoding="utf-8"))
    s = d["summary"]
    inf = s["envelope_inflation_exact_over_binned"]
    n = s["comparisons"]
    lines = [
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"& binned & exact radius \\",
        r"\midrule",
        f"Resolved comparisons & {s['resolved_binned']} & "
        f"{s['resolved_exact']} \\\\",
        f"Won by the interior member & {s['interior_wins_binned']} & "
        f"{s['interior_wins_exact']} \\\\",
        f"Won by the constant comparator & "
        f"{s['resolved_binned'] - s['interior_wins_binned']} & "
        f"{s['resolved_exact'] - s['interior_wins_exact']} \\\\",
        r"\midrule",
        r"\multicolumn{3}{l}{\emph{What changed, and what did not}} \\",
        r"Verdicts reversed & \multicolumn{2}{c}{"
        f"{s['verdicts_changed']} of {n}" r"} \\",
        r"Resolution lost & \multicolumn{2}{c}{"
        f"{s['lost_resolution']} of {n}" r"} \\",
        r"Envelope, exact over binned & \multicolumn{2}{c}{"
        f"median ${_fmt(inf['median'], 1)}$ "
        f"(${_fmt(inf['min'], 2)}$--${_fmt(inf['max'], 1)}$)" r"} \\",
        r"Error, exact over binned & \multicolumn{2}{c}{"
        f"median ${_fmt(s['error_ratio_median'], 2)}$ "
        f"(${_fmt(s['error_ratio_min'], 2)}$--"
        f"${_fmt(s['error_ratio_max'], 2)}$)" r"} \\",
        r"Realized work, exact over binned & \multicolumn{2}{c}{"
        f"median ${_fmt(s['work_ratio_median'], 3)}$ "
        f"(${_fmt(s['work_ratio_min'], 3)}$--"
        f"${_fmt(s['work_ratio_max'], 3)}$)" r"} \\",
        r"\midrule",
        r"\multicolumn{3}{l}{\footnotesize Rebuilt binned thresholds reproduce "
        r"the archived ones to "
        f"${_fmt(s['worst_threshold_reproduction_rel_diff'], 1)}$ "
        r"relative difference, worst case.} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    out = METRICS / "r24_bin_control_table.tex"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[r24 tables] wrote {out.name}")
    return 0


def descriptives() -> int:
    """The prose numbers, so the text can be checked against the record."""
    out = {}
    src = METRICS / "r24_oracle_ultra.json"
    if src.exists():
        d = json.loads(src.read_text(encoding="utf-8"))
        rows = d["rows"]
        grew = [r for r in rows
                if r["interior_envelope_after_m"]
                > r["interior_envelope_before_m"]]
        new_or = [r for r in rows if r["resolved_vs_oracle_after"]
                  and not r["resolved_vs_oracle_before"]]
        undec = [r for r in rows if not r["resolved_vs_oracle_after"]]
        withdrawn = [r for r in rows if r["resolved_vs_sat_before"]
                     and not r["resolved_vs_sat_after"]]
        out["oracle_ultra"] = {
            "summary": d["summary"],
            "envelope_shrink": d["envelope_shrink_factor"],
            "newly_resolved_vs_oracle": [
                {"orbit": f"{r['design']}{r['sobol_index']:03d}",
                 "winner": r["winner_vs_oracle_after"],
                 "interior_error_m": r["interior_error_m"],
                 "oracle_error_m": r["oracle_error_m"]} for r in new_or],
            "envelope_grew": [
                {"orbit": f"{r['design']}{r['sobol_index']:03d}",
                 "before_m": r["interior_envelope_before_m"],
                 "after_m": r["interior_envelope_after_m"],
                 "interior_self_difference_m":
                     r["interior_self_difference_tighter_to_ultra_m"],
                 "truth_self_difference_m":
                     r["truth_self_difference_tighter_to_ultra_m"]}
                for r in grew],
            "still_undecided_vs_oracle": [
                {"orbit": f"{r['design']}{r['sobol_index']:03d}",
                 "interior_error_m": r["interior_error_m"],
                 "oracle_error_m": r["oracle_error_m"],
                 "interior_envelope_m": r["interior_envelope_after_m"]}
                for r in undec],
            "vs_sat_verdict_withdrawn": [
                {"orbit": f"{r['design']}{r['sobol_index']:03d}",
                 "was": r["winner_vs_sat_before"]} for r in withdrawn],
        }
    src = METRICS / "r24_bin_control.json"
    if src.exists():
        d = json.loads(src.read_text(encoding="utf-8"))
        out["bin_control"] = {
            "summary": d["summary"],
            "panel_completeness": d["panel_completeness"],
        }
    if not out:
        return 1
    p = METRICS / "r24_manuscript_descriptives.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[r24 tables] wrote {p.name}")
    return 0


def main() -> int:
    rc = 0
    rc |= oracle_ultra_table()
    rc |= bin_control_table()
    rc |= descriptives()
    return 0 if rc in (0, 1) else rc


if __name__ == "__main__":
    raise SystemExit(main())
