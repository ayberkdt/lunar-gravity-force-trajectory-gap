"""Why the matched-work Atallah comparisons are unresolved (R13 diagnosis).

The benchmark resolves a comparison only when the error gap exceeds the summed
truth-inclusive envelopes. Against the critical fixed degree that happens almost
always; against the rule's own work-matched fixed degree it rarely does. Before
spending compute on tighter runs, this script separates the two possible causes:

  (a) the two force models really are equivalent at the resolution of a seven-day
      arc, or
  (b) the numerical envelope of this experiment is larger than it needs to be.

For every orbit of both populations it records the resolution margin

  M_res = |E_At - E_fixed| / (E_num,At + E_num,fixed),

the absolute errors next to the threshold, the composition of the envelope
(truth self-difference versus policy self-difference), a hypothetical in which
the policy self-differences are zero (an upper bound on what a perfect policy
integration could buy), and---the decisive test---whether the measured error of
each policy is stable under a one-decade tolerance refinement. A physical
truncation signal is tolerance-stable; integration noise is not.

Selection for the targeted ultra-tighter campaign follows from the margins:
every resolved case (as a validation set) plus every borderline case with
0.5 < M_res <= 1.

Usage:
    python rev13_resolution_diagnosis.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

METRICS = Path(__file__).resolve().parents[1] / "metrics"
OUTPUT = METRICS / "r13_resolution_diagnosis.json"
TABLE = METRICS / "r13_resolution_diagnosis_table.tex"
SELECTION = METRICS / "r13_ultratight_selection.json"

CAMPAIGN = {"A": METRICS / "r12_atallah_campaign.json",
            "B": METRICS / "r12_atallah_campaign_designB.json"}
BORDERLINE = (0.5, 1.0)
COMPARISON = "atallah_vs_fixed_work_atallah"
CONTROL = "atallah_vs_fixed_critical"


def stat(v):
    a = np.asarray([x for x in v if x is not None and np.isfinite(x)], dtype=float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p25": float(np.percentile(a, 25)),
            "p75": float(np.percentile(a, 75)), "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def tolerance_stability(policy: dict) -> float:
    """|E(tighter) - E(tight)| / E(tight) for one policy of one orbit."""
    a = policy["error_tight"]["pos_rms_m"]
    b = policy["error_tighter"]["pos_rms_m"]
    return abs(b - a) / a if a > 0 else None


def per_design(design: str) -> dict:
    payload = json.loads(CAMPAIGN[design].read_text())
    rows = []
    for r in payload["rows"]:
        c = r["comparisons"][COMPARISON]
        pa, pf = r["policies"]["atallah"], r["policies"]["fixed_work_atallah"]
        pc = r["policies"]["fixed_critical"]
        # the truth self-difference is the part of each envelope that is common
        truth_self = pa["truth_inclusive_envelope_m"] - pa["self_difference_rms_m"]
        m_res = c["absolute_error_difference_m"] / c["resolution_threshold_m"]
        # gap at each tolerance level: a physical difference keeps its size and
        # sign, integration noise does not
        gap_tight = (pa["error_tight"]["pos_rms_m"] - pf["error_tight"]["pos_rms_m"])
        gap_tighter = (pa["error_tighter"]["pos_rms_m"]
                       - pf["error_tighter"]["pos_rms_m"])
        rows.append({
            "sobol_index": r["sobol_index"], "hp_km": r["design_point"]["hp_km"],
            "incl_deg": r["design_point"]["incl_deg"],
            "error_atallah_m": c["atallah_error_m"],
            "error_fixed_work_m": c["comparator_error_m"],
            "gap_m": c["absolute_error_difference_m"],
            "threshold_m": c["resolution_threshold_m"],
            "m_res": m_res, "resolved": c["resolved"],
            "self_difference_atallah_m": pa["self_difference_rms_m"],
            "self_difference_fixed_work_m": pf["self_difference_rms_m"],
            "truth_self_difference_m": truth_self,
            "truth_share_of_threshold": 2.0 * truth_self / c["resolution_threshold_m"],
            "m_res_zero_policy_self": c["absolute_error_difference_m"] / (2.0 * truth_self),
            "gap_tight_m": gap_tight, "gap_tighter_m": gap_tighter,
            "gap_sign_stable": bool(np.sign(gap_tight) == np.sign(gap_tighter)),
            "stability_atallah": tolerance_stability(pa),
            "stability_fixed_work": tolerance_stability(pf),
            "stability_fixed_critical": tolerance_stability(pc),
            "m_res_critical": (r["comparisons"][CONTROL]["absolute_error_difference_m"]
                               / r["comparisons"][CONTROL]["resolution_threshold_m"]),
        })
    m = np.array([x["m_res"] for x in rows])
    summary = {
        "orbits": len(rows),
        "m_res": stat(m),
        "counts": {
            "resolved_m_gt_1": int((m > 1).sum()),
            "borderline_0p5_to_1": int(((m > BORDERLINE[0]) & (m <= BORDERLINE[1])).sum()),
            "0p2_to_0p5": int(((m > 0.2) & (m <= 0.5)).sum()),
            "below_0p2": int((m <= 0.2).sum())},
        "errors": {
            "atallah_m": stat([x["error_atallah_m"] for x in rows]),
            "fixed_work_m": stat([x["error_fixed_work_m"] for x in rows]),
            "threshold_m": stat([x["threshold_m"] for x in rows]),
            "both_errors_below_threshold": int(sum(
                x["error_atallah_m"] < x["threshold_m"]
                and x["error_fixed_work_m"] < x["threshold_m"] for x in rows))},
        "envelope": {
            "truth_self_difference_m": stat([x["truth_self_difference_m"] for x in rows]),
            "self_difference_atallah_m": stat(
                [x["self_difference_atallah_m"] for x in rows]),
            "self_difference_fixed_work_m": stat(
                [x["self_difference_fixed_work_m"] for x in rows]),
            "truth_share_of_threshold": stat(
                [x["truth_share_of_threshold"] for x in rows]),
            "resolved_if_policy_self_zero": int(sum(
                x["m_res_zero_policy_self"] > 1 for x in rows))},
        "tolerance_stability": {
            "atallah": stat([x["stability_atallah"] for x in rows]),
            "fixed_work": stat([x["stability_fixed_work"] for x in rows]),
            "fixed_critical": stat([x["stability_fixed_critical"] for x in rows])},
        "gap_behaviour": {
            "median_abs_gap_tight_m": float(np.median(
                [abs(x["gap_tight_m"]) for x in rows])),
            "median_abs_gap_tighter_m": float(np.median(
                [abs(x["gap_tighter_m"]) for x in rows])),
            "sign_stable": int(sum(x["gap_sign_stable"] for x in rows))},
        "critical_control": {
            "m_res": stat([x["m_res_critical"] for x in rows])},
    }
    return {"rows": rows, "summary": summary}


def selection(designs: dict) -> dict:
    """Resolved (validation) plus borderline orbits, per design."""
    out = {}
    for d, payload in designs.items():
        res = [r["sobol_index"] for r in payload["rows"] if r["resolved"]]
        bor = [r["sobol_index"] for r in payload["rows"]
               if BORDERLINE[0] < r["m_res"] <= BORDERLINE[1]]
        out[d] = {"resolved": sorted(res), "borderline": sorted(bor),
                  "all": sorted(set(res) | set(bor))}
    out["rule"] = ("every resolved matched-work comparison (validation set) plus "
                   "every comparison with 0.5 < M_res <= 1 (reachable by one more "
                   "decade of tolerance if the gap is a physical signal)")
    return out


def build_table(designs: dict) -> str:
    def row(d):
        s = designs[d]["summary"]
        m, c, e, env, ts, gb = (s["m_res"], s["counts"], s["errors"],
                                s["envelope"], s["tolerance_stability"],
                                s["gap_behaviour"])
        return (f"    {d} & {m['median']:.2f} & {c['resolved_m_gt_1']} & "
                f"{c['borderline_0p5_to_1']} & {c['below_0p2']} & "
                f"{e['atallah_m']['median']:.3f} & {e['fixed_work_m']['median']:.3f} & "
                f"{e['threshold_m']['median']:.2f} & "
                f"{env['resolved_if_policy_self_zero']} & "
                f"{ts['fixed_critical']['median'] * 100:.1f} & "
                f"{ts['atallah']['median'] * 100:.0f} & "
                f"{ts['fixed_work']['median'] * 100:.0f} & "
                f"{gb['sign_stable']}\\\\")

    body = "\n".join(row(d) for d in ("A", "B"))
    return f"""% auto-generated by rev13_resolution_diagnosis.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\scriptsize
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Diagnosis of the unresolved matched-work comparisons. $M_{{\\mathrm{{res}}}}
  = |E_{{\\mathrm{{At}}}}-E_{{\\mathrm{{fix}}}}|/(E_{{\\mathrm{{num,At}}}}+E_{{\\mathrm{{num,fix}}}})$
  is the resolution margin, and a comparison resolves at
  $M_{{\\mathrm{{res}}}}>1$. Column ``$\\to0$ self'' counts the comparisons that
  would resolve if both policy self-differences were driven to zero, leaving only
  the truth envelope; it bounds what a more accurate policy integration alone
  could buy. The stability columns give the median relative change of each
  policy's measured seven-day error under a one-decade tolerance refinement, a
  physical truncation signal being tolerance-stable and integration noise not.
  The last column counts, out of 64, the orbits in which the sign of
  $E_{{\\mathrm{{At}}}}-E_{{\\mathrm{{fix}}}}$ survives that refinement.}}
  \\label{{tab:resolution-diagnosis}}
  \\begin{{tabular}}{{l r r r r r r r r r r r r}}
    \\toprule
    & \\multicolumn{{4}}{{c}}{{$M_{{\\mathrm{{res}}}}$}} & \\multicolumn{{3}}{{c}}{{medians [m]}} & &
      \\multicolumn{{3}}{{c}}{{stability [\\%]}} & \\\\
    \\cmidrule(lr){{2-5}}\\cmidrule(lr){{6-8}}\\cmidrule(lr){{10-12}}
    Design & med & $>1$ & $0.5$--$1$ & $\\le0.2$ & $E_{{\\mathrm{{At}}}}$ &
      $E_{{\\mathrm{{fix}}}}$ & thr. & $\\to0$ self & crit & At & fix & sign\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""


def main() -> int:
    designs = {d: per_design(d) for d in ("A", "B")}
    payload = {"schema": "r13_resolution_diagnosis_v1",
               "comparison": COMPARISON,
               "designs": designs}
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    sel = selection(designs)
    SELECTION.write_text(json.dumps(sel, indent=2), encoding="utf-8")
    TABLE.write_text(build_table(designs), encoding="utf-8")
    for d in ("A", "B"):
        s = designs[d]["summary"]
        print(f"[design {d}] M_res median {s['m_res']['median']:.3f}; "
              f"resolved {s['counts']['resolved_m_gt_1']}, "
              f"borderline {s['counts']['borderline_0p5_to_1']}, "
              f"<=0.2 {s['counts']['below_0p2']}")
        print(f"           errors med At {s['errors']['atallah_m']['median']:.4f} m, "
              f"fix {s['errors']['fixed_work_m']['median']:.4f} m, "
              f"threshold {s['errors']['threshold_m']['median']:.4f} m")
        print(f"           stability med: critical {s['tolerance_stability']['fixed_critical']['median']:.3f}, "
              f"atallah {s['tolerance_stability']['atallah']['median']:.3f}, "
              f"fixed_work {s['tolerance_stability']['fixed_work']['median']:.3f}")
        print(f"           gap median {s['gap_behaviour']['median_abs_gap_tight_m']:.4f} -> "
              f"{s['gap_behaviour']['median_abs_gap_tighter_m']:.4f} m, "
              f"sign stable {s['gap_behaviour']['sign_stable']}/64")
        print(f"           selection: {sel[d]['all']}")
    print(f"[written] {OUTPUT.name}, {SELECTION.name}, {TABLE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
