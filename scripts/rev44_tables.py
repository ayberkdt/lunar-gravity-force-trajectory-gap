"""Table for the R44 (O42) equal-realized-work comparison, matched at the
tighter level.

Emits metrics/r44_equal_work_table.tex and the descriptives the manuscript
quotes, so no number is transcribed by hand. This supersedes
rev19_tables.py as the source of the main text's tab:equal-work: the R19
match holds tight-level realized work equal while every error in the
comparison is read at the tighter level, and council review flagged that
level inconsistency. The R44 rows hold realized work equal at the level the
errors are scored at. The R19 record is untouched and remains re-tallied at
three resolution cuts in the supplement's threshold-sensitivity table.

Two conventions carried over from rev19_tables.py:

Budget scale. All propagated cells are carried because the verdict is not
the same at all of them; a budget whose records are absent is skipped rather
than faked, so this generator is safe to run mid-campaign.

Like-for-like ratios. median_of_ratios is the statistic comparable with
median rho; ratio_of_medians is emitted alongside under its own name so the
two can never be silently swapped. Under the R44 match the ratio of medians
is the *smaller* of the two in every cell, so no conservativity claim is
attached to either.

Usage:  python rev44_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BETAS_BY_DESIGN = {"A": (1.25, 1.00, 0.75, 0.50), "B": (1.50, 1.00, 0.75, 0.50)}
K = "0.50"
CUTS = (0.5, 1.0, 2.0)


def beta_tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def load(design: str, beta: float):
    p = METRICS / f"r44_equal_work_tighter_{design}_{beta_tag(beta)}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def per_call(design: str, beta: float) -> dict | None:
    """The same comparison under the nominal per-call budget, for contrast."""
    p = METRICS / f"r18_span_sweep_{design}_{beta_tag(beta)}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    win = lose = unres = 0
    ratios = []
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
    if not ratios:
        return None
    return {"interior": win, "fixed": lose, "unresolved": unres,
            "median_of_ratios": float(np.median(ratios))}


def member_overspend(design: str, beta: float) -> dict | None:
    """Member tighter-level realized work over the nominal constant degree's,
    read from the archived-telemetry pair each case config carries."""
    root = METRICS / f"r44_cases/{design}_workmatched_tighter_{beta_tag(beta)}"
    if not root.is_dir():
        return None
    ratios = []
    for p in sorted(root.glob("sobol*/fixed_tighter.json")):
        c = json.loads(p.read_text(encoding="utf-8"))["config"]
        ratios.append(c["target_total_quadratic_work_tighter"]
                      / c["constant_endpoint_total_work_tighter"])
    if not ratios:
        return None
    return {"median": float(np.median(ratios)),
            "min": float(min(ratios)), "max": float(max(ratios)),
            "n": len(ratios)}


def retally(rows: list[dict]) -> dict:
    """The resolution-cut re-count, nothing repropagated: only the multiple
    of the recorded summed envelope required to call a comparison resolved
    changes."""
    out = {}
    for cut in CUTS:
        w = l = u = 0
        for r in rows:
            if not r.get("work_matched_error_m"):
                continue
            diff = r["work_matched_error_m"] - r["interior_error_m"]
            thr = cut * r["resolution_threshold_m"]
            if diff > thr:
                w += 1
            elif -diff > thr:
                l += 1
            else:
                u += 1
        out[f"{cut:g}"] = {"interior": w, "fixed": l, "unresolved": u}
    return out


def table() -> str:
    blocks = []
    for design, betas in BETAS_BY_DESIGN.items():
        rows = []
        for beta in betas:
            d = load(design, beta)
            pc = per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            a = s["achieved_work_ratio_tighter"]
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
    for design, betas in BETAS_BY_DESIGN.items():
        for beta in betas:
            d = load(design, beta)
            pc = per_call(design, beta)
            if not d or not pc:
                continue
            s = d["summary"]
            rows = d["rows"]
            shift = [r["work_matched_degree"] - r["constant_endpoint_degree"]
                     for r in rows if r.get("work_matched_degree")]
            tight = [r["achieved_work_ratio_tight"] for r in rows
                     if r.get("achieved_work_ratio_tight")]
            e_int = [r["interior_error_m"] for r in rows
                     if r.get("interior_error_m")]
            e_fix = [r["work_matched_error_m"] for r in rows
                     if r.get("work_matched_error_m")]
            ratios = [r["rho_workmatched"] for r in rows
                      if r.get("rho_workmatched")]
            losers = sorted(r["sobol_index"] for r in rows
                            if r.get("winner") == "fixed")
            out[f"{design}_{beta_tag(beta)}"] = {
                "design": design, "beta": beta,
                "orbits": s["orbits"], "resolved": s["resolved"],
                "interior_wins": s["resolved_interior_wins"],
                "fixed_wins": s["resolved_fixed_wins"],
                "unresolved": s["unresolved"],
                "median_rho": s["median_rho"],
                "achieved_work_ratio_tighter": s["achieved_work_ratio_tighter"],
                "achieved_work_ratio_tight_median": float(np.median(tight)),
                "comparator_degree_shift": {
                    "median": float(np.median(shift)),
                    "min": int(min(shift)), "max": int(max(shift))},
                "censored": sum(1 for r in rows if r.get("censored")),
                "median_of_ratios": float(np.median(ratios)),
                "ratio_of_medians": float(np.median(e_fix) / np.median(e_int)),
                "resolution_cut_retally": retally(rows),
                "member_over_constant_work_tighter":
                    member_overspend(design, beta),
                "orbits_where_interior_loses": losers,
                "per_call": pc,
            }

    # at the declared budget, how the fixed-win sets compare with the nominal
    # match: the level change should not manufacture losses on new orbits
    for design in ("A", "B"):
        p18 = METRICS / f"r18_span_sweep_{design}_beta_1.00.json"
        key = f"{design}_beta_1.00"
        if not p18.exists() or key not in out:
            continue
        d18 = json.loads(p18.read_text(encoding="utf-8"))
        nominal_losers = set()
        for r in d18["rows"]:
            e0, ek = r["entries"]["0.00"], r["entries"].get(K, {})
            if not ek.get("error_m") or not e0.get("error_m"):
                continue
            thr = (e0.get("envelope_m") or 0.0) + (ek.get("envelope_m") or 0.0)
            if e0["error_m"] - ek["error_m"] < -thr:
                nominal_losers.add(r["sobol_index"])
        r44_losers = set(out[key]["orbits_where_interior_loses"])
        out[key]["nominal_fixed_wins"] = sorted(nominal_losers)
        out[key]["nominal_fixed_wins_recurring"] = sorted(
            nominal_losers & r44_losers)

    (METRICS / "r44_equal_work_table.tex").write_text(table(),
                                                      encoding="utf-8")
    (METRICS / "r44_manuscript_descriptives.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("[written] r44_equal_work_table.tex, "
          "r44_manuscript_descriptives.json")
    for key, v in out.items():
        print(f"  {key}: realized {v['interior_wins']}/{v['fixed_wins']}/"
              f"{v['unresolved']}, rho {v['median_rho']:.2f}, "
              f"work {v['achieved_work_ratio_tighter']['median']:.3f} "
              f"(tight-level {v['achieved_work_ratio_tight_median']:.3f}), "
              f"shift {v['comparator_degree_shift']['median']:.0f}, "
              f"censored {v['censored']}")
        ov = v["member_over_constant_work_tighter"]
        if ov:
            print(f"      member/N0 tighter work {ov['median']:.3f} "
                  f"[{ov['min']:.3f},{ov['max']:.3f}] over {ov['n']} orbits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
