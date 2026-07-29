"""Tables for the R17 sixty-day widened-geometry campaign.

Emits, into metrics/:
  r17_longarc60_table.tex        main text, one row per orbit
  r17_longarc60_growth_table.tex supplement, in-track growth by checkpoint

Both are built from metrics/r17_longarc60.json. Orbits that did not reach the
full arc are listed with their termination day and excluded from the numeric
columns, so the table reports the sample it actually has.

Usage:  python rev17_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
SRC = METRICS / "r17_longarc60.json"

CHECKPOINTS = (7, 14, 28, 42, 60)
SCHED_LABEL = {"sched_emp": "emp", "sched_down": "down", "sched_up": "up"}


def geometry_label(orb: dict) -> str:
    return (f"${orb['hp_km']:.0f}\\times{orb['ha_km']:.0f}$, "
            f"$i={orb['incl_deg']:.0f}^\\circ$")


def ordered_rows(payload: dict) -> list:
    """Perilune order, so the table can be read against the low-perilune
    mechanism the seven-day campaign identified."""
    return sorted((r for r in payload["rows"] if r.get("status") == "complete"),
                  key=lambda r: r["orbit"]["hp_km"])


def main_table(payload: dict) -> str:
    rows = ordered_rows(payload)
    lines = ["\\begin{tabular}{@{}l r r r r@{\\hspace{4pt}}l r c@{}}",
             "\\toprule",
             "Geometry & $N_{\\mathrm{crit}}$ & fixed crit.\\ & env. & "
             "\\multicolumn{2}{c}{best sched.} & $\\rho_{\\mathrm{crit}}$ & "
             "res. \\\\",
             " & & [m] & [m] & \\multicolumn{2}{c}{[m]} & & \\\\",
             "\\midrule"]
    for r in rows:
        orb = r["orbit"]
        if not r.get("reached_full_arc"):
            lines.append(
                f"{geometry_label(orb)} & \\multicolumn{{7}}{{c}}"
                f"{{reference orbit reached the surface on day "
                f"{r.get('arc_end_day', float('nan')):.1f}; excluded}}\\\\")
            continue
        pol = r["policies"]
        best = r["best_schedule_name"]
        e_crit = pol["fixed_crit"]["errors_against_same_tolerance_truth"][
            "tighter"]["pos_rms_m"]
        e_best = pol[best]["errors_against_same_tolerance_truth"][
            "tighter"]["pos_rms_m"]
        comp = r["comparisons"].get(f"{best}_vs_fixed_crit", {})
        env = comp.get("resolution_threshold_m")
        rho = comp.get("rho")
        res = ("yes" if comp.get("resolved") else "no")
        lines.append(
            f"{geometry_label(orb)} & {r['spec']['n_crit']} & "
            f"{e_crit:.1f} & {env:.1f} & {e_best:.1f} & "
            f"{SCHED_LABEL[best]} & {rho:.3f} & {res} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def growth_table(payload: dict) -> str:
    """Cumulative RMS is monotone by construction and is the growth measure;
    the instantaneous in-track value is signed and oscillates, so it is shown
    as its own row rather than mixed into the same series."""
    rows = [r for r in ordered_rows(payload) if r.get("reached_full_arc")]
    head = " & ".join(f"d{d}" for d in CHECKPOINTS)
    lines = ["\\begin{tabular}{@{}l l " + "r " * len(CHECKPOINTS) + "@{}}",
             "\\toprule",
             f"Geometry & quantity & {head} \\\\",
             "\\midrule"]
    for r in rows:
        orb = r["orbit"]
        best = r["best_schedule_name"]
        series = (
            (best, "pos_rms_m", f"sched.\\ {SCHED_LABEL[best]}, cum.\\ RMS"),
            (best, "in_track_at_checkpoint_m",
             f"sched.\\ {SCHED_LABEL[best]}, in-track"),
            ("fixed_crit", "pos_rms_m", "fixed crit., cum.\\ RMS"),
        )
        for policy, field, label in series:
            cps = r["policies"][policy]["checkpoints"]
            cells = []
            for d in CHECKPOINTS:
                key = f"d{d}"
                cells.append(f"{abs(cps[key][field]):.0f}"
                             if key in cps else "---")
            lines.append(f"{geometry_label(orb)} & {label} & "
                         + " & ".join(cells) + " \\\\")
        lines.append("\\addlinespace[2pt]")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def descriptives(payload: dict) -> dict:
    """Every number the manuscript quotes from this campaign, derived here
    rather than read off a table by hand."""
    rows = [r for r in payload["rows"]
            if r.get("status") == "complete" and r.get("reached_full_arc")]
    s = payload["summary"]

    per_orbit = []
    for r in rows:
        best = r["best_schedule_name"]
        c = r["comparisons"].get(f"{best}_vs_fixed_crit", {})
        per_orbit.append({
            "name": r["orbit"]["name"],
            "hp_km": r["orbit"]["hp_km"],
            "incl_deg": r["orbit"]["incl_deg"],
            "n_crit": r["spec"]["n_crit"],
            "best_schedule": best,
            "rho_crit_best": c.get("rho"),
            "resolved": c.get("resolved"),
            "winner": c.get("winner_if_resolved"),
        })
    per_orbit.sort(key=lambda x: x["hp_km"])

    # The two comparator families answer different questions and must not be
    # pooled: fixed_crit is the degree a user would pick for the orbit's
    # critical altitude, fixed_work is matched to the schedule's own mean
    # squared degree. Totalling them together once produced a misleading
    # "schedule wins" count that belonged entirely to the work-matched family.
    def totals(suffix: str) -> dict:
        cs = [c for key, c in s.get("comparisons", {}).items()
              if key.endswith(suffix)]
        return {"comparisons": sum(c["n_orbits"] for c in cs),
                "resolved": sum(c["resolved"] for c in cs),
                "resolved_schedule_wins": sum(
                    c["resolved_schedule_wins"] for c in cs),
                "resolved_fixed_wins": sum(
                    c["resolved_fixed_wins"] for c in cs),
                "unresolved": sum(c["unresolved"] for c in cs)}

    rhos = [p["rho_crit_best"] for p in per_orbit
            if p["rho_crit_best"] is not None]
    directions = sorted({p["best_schedule"] for p in per_orbit})

    return {
        "orbits_reaching_full_arc": s["orbits_reaching_60_days"],
        "orbits_attempted": s["orbits_attempted"],
        "orbits_terminated_early": s["orbits_terminated_early"],
        "early_termination": s["early_termination_days"],
        "orbits_with_work_matched_comparator": s.get(
            "orbits_with_work_matched_comparator"),
        "vs_critical_degree": totals("_vs_fixed_crit"),
        "vs_work_matched_degree": totals("_vs_fixed_work"),
        "best_schedule_rho_vs_critical": {
            "min": min(rhos) if rhos else None,
            "max": max(rhos) if rhos else None},
        "best_schedule_directions_observed": directions,
        "truth_self_difference_rms_m": s["truth_self_difference_rms_m"],
        "in_track_d60_over_d7": s["best_schedule_in_track_d60_over_d7"],
        "per_orbit_by_perilune": per_orbit,
    }


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}; run rev17_longarc60.py first")
    payload = json.loads(SRC.read_text(encoding="utf-8"))
    (METRICS / "r17_manuscript_descriptives.json").write_text(
        json.dumps(descriptives(payload), indent=2), encoding="utf-8")
    (METRICS / "r17_longarc60_table.tex").write_text(
        main_table(payload), encoding="utf-8")
    (METRICS / "r17_longarc60_growth_table.tex").write_text(
        growth_table(payload), encoding="utf-8")
    s = payload["summary"]
    print(f"[written] r17_longarc60_table.tex, r17_longarc60_growth_table.tex")
    print(f"  {s['orbits_reaching_60_days']} orbits reached 60 days of "
          f"{s['orbits_attempted']} attempted")
    for key, c in s.get("comparisons", {}).items():
        print(f"  {key:26s} n={c['n_orbits']:2d} resolved={c['resolved']:2d} "
              f"(sched {c['resolved_schedule_wins']}, "
              f"fixed {c['resolved_fixed_wins']}) "
              f"median rho={c['rho']['median'] if c['rho'] else None}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
