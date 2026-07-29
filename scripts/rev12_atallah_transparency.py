"""Transparency record for the perilune-tolerance-matched Atallah benchmark (R12).

The benchmark matches an error *target*, not a degree and not a delivered
accuracy: the acceleration tolerance handed to the published rule is the actual
worst-case acceleration truncation error that the critical fixed degree already
incurs at perilune, and the rule then chooses whatever degree drives its own
conservative analytical bound below that target. Because the bound is
conservative, the degree it selects at perilune is much higher than the critical
degree and the acceleration error it actually delivers there is much smaller
than the target. This script measures that gap explicitly, per orbit and for
both populations:

  N_crit                critical-altitude fixed degree
  N_At(r_p)             degree the rule uses in the perilune bin
  E_crit(r_p)           actual worst-case |a(N_crit) - a(N_truth)| at perilune
                        (this is the matched tolerance, by construction)
  E_At(r_p)             actual worst-case |a(N_At) - a(N_truth)| at perilune
  E_crit/E_At           how much extra real accuracy the conservative bound buys
  <N^2>_At/N_crit^2     the work it spends for it

It also records truth-degree cap contact. The rule is capped at the adopted
truth degree of each orbit; where it selects that cap in the perilune bin, its
model there coincides with the truth model, so its measured error is zero by
construction in those bins. The decision counts are therefore also recomputed on
the cap-free subgroup as a sensitivity check.

Usage:
    python rev12_atallah_transparency.py            # both populations
    python rev12_atallah_transparency.py --design A
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r12_atallah_transparency.json"
GRID = {"n_lat": 25, "n_lon": 48}  # same lat/lon grid as the campaign's matching

DESIGNS = {
    "A": {"case": METRICS / "r12_cases" / "atallah",
          "campaign": METRICS / "r12_atallah_campaign.json",
          "descriptives": METRICS / "r12_atallah_descriptives.json"},
    "B": {"case": METRICS / "r12_cases" / "atallah_designB",
          "campaign": METRICS / "r12_atallah_campaign_designB.json",
          "descriptives": METRICS / "r12_atallah_descriptives_designB.json"},
}

_MODELS: dict[int, tuple] = {}


def model_for(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def perilune_bin_degree(table: dict, hp_km: float, bin_km: float = 10.0) -> int:
    """Degree the flight-time selector returns at perilune (same clamping)."""
    tab = {float(k): int(v) for k, v in table.items()}
    hb = bin_km * math.floor(hp_km / bin_km)
    hb = min(max(tab), max(min(tab), hb))
    return tab[hb]


def per_orbit(design: str) -> list[dict]:
    cfgs = DESIGNS[design]
    campaign = json.loads(cfgs["campaign"].read_text())
    desc = json.loads(cfgs["descriptives"].read_text())
    cost = desc["per_orbit_cost"]
    rows = []
    for row in campaign["rows"]:
        index = int(row["sobol_index"])
        side = json.loads((cfgs["case"] / f"sobolA_{index:03d}" /
                           "atallah_tight.json").read_text())
        cfg, tel = side["config"], side["telemetry"]
        adopted = int(cfg["adopted_truth_degree"])
        n_crit = int(cfg["n_critical"])
        tol = float(cfg["atallah_tol_accel_m_s2"])
        hp_km = float(row["design_point"]["hp_km"])
        model, args = model_for(adopted)
        r_p = model.r_ref + hp_km * 1e3
        n_at = perilune_bin_degree(cfg["atallah_degree_table"], hp_km)
        at_cap = bool(n_at >= adopted)
        # actual worst-case acceleration truncation error of each degree at
        # perilune, against the same adopted truth degree used for the matching
        err_at = (0.0 if at_cap else
                  at.actual_truncation_error_max(model, args, r_p, n_at, adopted,
                                                 **GRID))
        counts = {int(k): int(v) for k, v in tel["degree_counts"].items()}
        n_rhs = sum(counts.values())
        c = cost.get(str(index), {})
        rows.append({
            "sobol_index": index, "hp_km": hp_km,
            "ha_km": float(row["design_point"]["ha_km"]),
            "adopted_truth_degree": adopted,
            "n_critical": n_crit,
            "n_atallah_perilune": n_at,
            "degree_ratio_atallah_over_critical": n_at / n_crit,
            "atallah_at_truth_cap_in_perilune_bin": at_cap,
            "rhs_fraction_at_truth_cap": sum(v for k, v in counts.items()
                                             if k >= adopted) / n_rhs,
            "tolerance_accel_m_s2": tol,
            "actual_error_critical_perilune_m_s2": tol,
            "actual_error_atallah_perilune_m_s2": err_at,
            "actual_error_ratio_critical_over_atallah": (
                None if err_at == 0.0 else tol / err_at),
            "atallah_bound_respected": bool(err_at <= tol),
            "work_ratio_atallah_over_critical": c.get("work_ratio_atallah_over_critical"),
            "atallah_mean_degree": c.get("atallah_mean_degree"),
            "n_work_atallah": c.get("n_work_atallah"),
        })
        print(f"  [{design}] idx={index:03d} N_crit={n_crit:3d} "
              f"N_At(r_p)={n_at:3d}{' (cap)' if at_cap else '     '} "
              f"tol={tol:.3e} E_At={err_at:.3e}", flush=True)
    return rows


def decision_split(design: str, rows: list[dict]) -> dict:
    """Decision counts restricted to orbits with and without truth-cap contact."""
    campaign = json.loads(DESIGNS[design]["campaign"].read_text())
    cap = {r["sobol_index"]: r["rhs_fraction_at_truth_cap"] > 0.0 for r in rows}
    out = {}
    for key in ("atallah_vs_fixed_critical", "atallah_vs_fixed_work_atallah"):
        for group, want in (("cap_contact", True), ("cap_free", False)):
            cs = [r["comparisons"][key] for r in campaign["rows"]
                  if cap.get(r["sobol_index"]) is want]
            res = [c for c in cs if c["resolved"]]
            rho = [c["rho"] for c in cs if c["rho"] is not None]
            out[f"{key}/{group}"] = {
                "pairs": len(cs), "resolved": len(res),
                "atallah_wins": sum(c["winner_if_resolved"] == "atallah" for c in res),
                "fixed_wins": sum(c["winner_if_resolved"] not in (None, "atallah")
                                  for c in res),
                "rho_median": float(np.median(rho)) if rho else None}
    return out


def stat(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def summarize(rows: list[dict]) -> dict:
    return {
        "orbits": len(rows),
        "n_critical": stat([r["n_critical"] for r in rows]),
        "n_atallah_perilune": stat([r["n_atallah_perilune"] for r in rows]),
        "degree_ratio": stat([r["degree_ratio_atallah_over_critical"] for r in rows]),
        "actual_error_ratio": stat(
            [r["actual_error_ratio_critical_over_atallah"] for r in rows]),
        "work_ratio": stat([r["work_ratio_atallah_over_critical"] for r in rows]),
        "orbits_at_truth_cap_in_perilune_bin": sum(
            r["atallah_at_truth_cap_in_perilune_bin"] for r in rows),
        "orbits_with_any_cap_contact": sum(
            r["rhs_fraction_at_truth_cap"] > 0.0 for r in rows),
        "rhs_fraction_at_cap_when_in_contact": stat(
            [r["rhs_fraction_at_truth_cap"] for r in rows
             if r["rhs_fraction_at_truth_cap"] > 0.0]),
        "bound_violations": sum(not r["atallah_bound_respected"] for r in rows),
    }


def _sci(value: float, digits: int = 2) -> str:
    """math-mode scientific notation without siunitx."""
    if value == 0.0:
        return "0"
    exp = int(math.floor(math.log10(abs(value))))
    mant = value / 10.0**exp
    return f"{mant:.{digits}f}\\times 10^{{{exp}}}"


def longtable(design: str, rows: list[dict]) -> str:
    body = []
    for r in rows:
        capped = r["atallah_at_truth_cap_in_perilune_bin"]
        err_at = ("0^{*}" if capped
                  else _sci(r["actual_error_atallah_perilune_m_s2"]))
        ratio = ("--" if r["actual_error_ratio_critical_over_atallah"] is None
                 else f"${_sci(r['actual_error_ratio_critical_over_atallah'], 1)}$")
        tol = _sci(r["tolerance_accel_m_s2"])
        body.append(
            f"    {r['sobol_index']:d} & {r['hp_km']:.0f} & {r['n_critical']:d} & "
            f"{r['n_atallah_perilune']:d} & {r['adopted_truth_degree']:d} & "
            f"{r['rhs_fraction_at_truth_cap'] * 100:.0f} & ${tol}$ & ${err_at}$ & "
            f"{ratio} & {r['work_ratio_atallah_over_critical']:.2f}\\\\")
    return f"""% auto-generated by rev12_atallah_transparency.py -- do not edit by hand
\\begin{{longtable}}{{r r r r r r l l r r}}
  \\caption{{Per-orbit matching record for the Atallah benchmark, design~{design}.
  $N_{{\\mathrm{{crit}}}}$ is the critical-altitude fixed degree,
  $N_{{\\mathrm{{At}}}}(r_p)$ the degree the published rule selects in the perilune
  bin, and $N_{{\\mathrm{{truth}}}}$ the adopted truth degree at which the rule is
  capped; $f_{{\\mathrm{{cap}}}}$ is the percentage of right-hand-side calls made at
  that cap. The tolerance $\\varepsilon_a$ handed to the rule is the actual
  worst-case acceleration truncation error of $N_{{\\mathrm{{crit}}}}$ at perilune,
  and $E_{{\\mathrm{{At}}}}$ is the actual worst-case error the selected degree
  delivers there (both against the same adopted truth degree, on the same
  $25\\times48$ lat/lon grid). Entries marked * are exactly zero because the rule
  sits at the truth-degree cap in that bin. The last two columns are the
  delivered accuracy gain $\\varepsilon_a/E_{{\\mathrm{{At}}}}$ and the quadratic
  work ratio $\\langle N^2\\rangle_{{\\mathrm{{At}}}}/N_{{\\mathrm{{crit}}}}^2$.}}
  \\label{{tab:atallah-matching-{design}}}\\\\
  \\toprule
  idx & $h_p$ [km] & $N_{{\\mathrm{{crit}}}}$ & $N_{{\\mathrm{{At}}}}(r_p)$ &
  $N_{{\\mathrm{{truth}}}}$ & $f_{{\\mathrm{{cap}}}}$ [\\%] & $\\varepsilon_a$ [m\\,s$^{{-2}}$] &
  $E_{{\\mathrm{{At}}}}$ [m\\,s$^{{-2}}$] & $\\varepsilon_a/E_{{\\mathrm{{At}}}}$ & $w$\\\\
  \\midrule
  \\endfirsthead
  \\toprule
  idx & $h_p$ [km] & $N_{{\\mathrm{{crit}}}}$ & $N_{{\\mathrm{{At}}}}(r_p)$ &
  $N_{{\\mathrm{{truth}}}}$ & $f_{{\\mathrm{{cap}}}}$ [\\%] & $\\varepsilon_a$ [m\\,s$^{{-2}}$] &
  $E_{{\\mathrm{{At}}}}$ [m\\,s$^{{-2}}$] & $\\varepsilon_a/E_{{\\mathrm{{At}}}}$ & $w$\\\\
  \\midrule
  \\endhead
  \\bottomrule
  \\endfoot
{chr(10).join(body)}
\\end{{longtable}}
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=("A", "B", "both"), default="both")
    ap.add_argument("--tables-only", action="store_true",
                    help="rebuild the LaTeX tables from the archived JSON")
    a = ap.parse_args()
    designs = ["A", "B"] if a.design == "both" else [a.design]
    payload = {"schema": "r12_atallah_transparency_v1",
               "created_utc": base.utc_now(),
               "grid": GRID, "designs": {}}
    if OUTPUT.exists():
        payload["designs"] = json.loads(OUTPUT.read_text()).get("designs", {})
    for d in designs:
        rows = (payload["designs"][d]["rows"] if a.tables_only
                else per_orbit(d))
        payload["designs"][d] = {"rows": rows, "summary": summarize(rows),
                                 "decision_split": decision_split(d, rows)}
        (METRICS / f"r12_atallah_matching_table_{d}.tex").write_text(
            longtable(d, rows), encoding="utf-8")
        print(f"[design {d}] {json.dumps(payload['designs'][d]['summary'], indent=2)}")
    base.atomic_json(OUTPUT, payload)
    print(f"[written] {OUTPUT.name} + per-design matching longtables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
