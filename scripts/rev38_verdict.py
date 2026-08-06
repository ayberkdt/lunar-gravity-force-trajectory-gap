"""R38 verdict: score the cap-lifted control against its own registration.

rev30_verdict.py cannot score this campaign, and it does not fail when asked to.
It scores the interior-versus-constant ladder, which is a different comparison
from the one R38 exists to check, and it maps outcomes by taking
``sorted(prereg["outcomes"])`` and reading the first three keys as
[holds, disagrees, undecided]. R38 registers four outcomes, so that mapping
lands on P by alphabetical position and prints a sentence about a collapse that
the numbers do not show. The generic record it wrote is kept, renamed, and
pointed at from here rather than deleted.

What this scores instead is the comparison the registration names: the radial
rule against the work-matched constant degree, from
``r14_trajectory_OEU_beta_*.json``, tallied by the archive's resolution rule,
next to the same tally on the capped parent and decomposed by whether the
parent's schedule reached the ceiling on that orbit.

One honest note on the registration itself. Outcomes P and S overlap -- a
constant-degree win is also a case where the tally no longer favours the radial
rule -- and Q and R share a clause about the ceiling flattering the magnitude.
The scoring order below resolves the overlap by specificity: S before P, and R
before Q whenever the median ratio falls by a decade or more. That ordering was
fixed here, after the numbers, and is declared as post-hoc. It changes no
threshold: the decade in R and the direction in P and S are the registered
words.

Usage:
    python rev38_verdict.py
    python rev38_verdict.py --beta 0.75
"""

from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

KEY, PARENT_KEY = "OEU", "OE"
DECADE = 10.0


def rows_of(path: Path) -> dict:
    return {r["sobol_index"]: r for r in
            json.loads(path.read_text(encoding="utf-8"))["rows"]}


def ceiling_orbits(beta: float) -> set:
    """Orbits whose calibrated schedule reached the ceiling in the capped run."""
    p = METRICS / "r31_budget_pareto_operational_elliptical.json"
    return {r["sobol_index"] for r in json.loads(p.read_text(encoding="utf-8")
                                                 )["designs"][PARENT_KEY]["rows"]
            if r["budgets"][f"beta_{beta:.2f}"]["atallah"]["allocation"][
                "fraction_at_cap"] > 0.0}


def tally(rows: dict, idxs) -> dict:
    idxs = sorted(idxs)
    a = f = u = 0
    rho_all, rho_res, undefined = [], [], 0
    for i in idxs:
        c = rows[i]["comparison"]
        rho = c["rho_budget"]
        if rho is None:
            undefined += 1
        else:
            rho_all.append(rho)
        if not c["resolved"]:
            u += 1
            continue
        if c["resolved_winner"] == "atallah":
            a += 1
        else:
            f += 1
        if rho is not None:
            rho_res.append(rho)
    if a + f == 0:
        verdict = "undecided"
    elif a > f:
        verdict = "radial"
    elif f > a:
        verdict = "constant"
    else:
        verdict = "split"
    work = [rows[i]["cost"]["total_work_ratio_atallah_over_fixed"] for i in idxs]
    return {
        "verdict": verdict, "orbits": len(idxs),
        "resolved": a + f, "unresolved": u,
        "resolved_radial_wins": a, "resolved_constant_wins": f,
        "median_rho_all_orbits": st.median(rho_all) if rho_all else None,
        "median_rho_resolved": st.median(rho_res) if rho_res else None,
        "rho_undefined_orbits": undefined,
        "median_realized_work_ratio_radial_over_constant": st.median(work),
    }


def net(t: dict) -> int:
    return t["resolved_radial_wins"] - t["resolved_constant_wins"]


def score(unc: dict, cap: dict, sub_unc: dict) -> tuple[str, str, str | None]:
    """The registered outcomes as a partition on the registered words only.

    The classification uses nothing but the registered criteria -- direction,
    and the decade in the median ratio -- so it is well defined and mutually
    exclusive. What it cannot see is that an outcome's *descriptive* clause can
    be false while its criterion is true: Q is selected on the ratio alone, but
    Q's sentence also claims the tally holds "by a margin of the same order",
    and a tally can keep its sign while losing almost all of its margin. That
    is the same failure mode as reporting "verdicts changed: 0" while the
    resolution collapses underneath it, and it is caught here rather than in
    review. The caveat is a statement of measured fact; the half-margin
    trigger for printing it was chosen after the numbers and buys no
    classification.
    """
    caveat = None
    if cap["verdict"] != "radial":
        # The registered outcomes are all sentences about a ceiling that
        # carried a radial win. At a budget where the capped run already
        # favours the constant degree there is no such win to confound, and
        # the labels describe something the numbers are not about.
        caveat = (
            f"the registered outcomes were written about the budget that "
            f"carries the published claim, where the capped run favours the "
            f"radial rule. Here the capped run already favours the constant "
            f"degree, {cap['resolved_radial_wins']}-"
            f"{cap['resolved_constant_wins']}. The label below is the "
            f"criterion applied mechanically, not a finding about the ceiling.")
    if unc["verdict"] == "constant":
        return ("S_reversal", "the uncapped tally favours the constant degree",
                caveat)
    if unc["verdict"] != "radial":
        return ("P_ceiling_carried_it",
                "the uncapped tally no longer resolves in favour of the "
                "radial rule", caveat)
    a, b = cap["median_rho_all_orbits"], unc["median_rho_all_orbits"]
    drop = a / b if (a and b) else None
    if net(unc) * 2 < net(cap):
        caveat = (
            f"the label is selected on the median ratio, and on the ratio it "
            f"is correct, but the tally margin does not survive with it: "
            f"{cap['resolved_radial_wins']}-{cap['resolved_constant_wins']} "
            f"(net {net(cap):+d}) becomes "
            f"{unc['resolved_radial_wins']}-{unc['resolved_constant_wins']} "
            f"(net {net(unc):+d}), and on the orbits whose schedule never "
            f"reached the ceiling the uncapped tally is "
            f"{sub_unc['resolved_radial_wins']}-"
            f"{sub_unc['resolved_constant_wins']}. The outcome label must not "
            f"be quoted without these numbers.")
    if drop is not None and drop >= DECADE:
        return ("R_direction_holds_magnitude_falls",
                f"the radial rule still wins the tally, "
                f"{unc['resolved_radial_wins']}-{unc['resolved_constant_wins']}"
                f" against {cap['resolved_radial_wins']}-"
                f"{cap['resolved_constant_wins']}, while the median error "
                f"ratio falls {drop:.1f}-fold, from {a:.2f} to {b:.2f} -- more "
                f"than the registered decade", caveat)
    return ("Q_survives_intact",
            f"the median error ratio falls {drop:.1f}-fold, less than the "
            f"registered decade, and the tally keeps its sign", caveat)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=1.00)
    a = ap.parse_args()
    b = a.beta

    unc_p = METRICS / f"r14_trajectory_{KEY}_beta_{b:.2f}.json"
    cap_p = METRICS / f"r14_trajectory_{PARENT_KEY}_beta_{b:.2f}.json"
    if not unc_p.exists():
        print(f"[abort] {unc_p.name} missing; beta {b:.2f} has not run")
        return 2
    unc, cap = rows_of(unc_p), rows_of(cap_p)
    touched = ceiling_orbits(b)
    every = set(cap)

    blocks = {
        "uncapped_all": tally(unc, every),
        "uncapped_never_capped_subset": tally(unc, every - touched),
        "uncapped_formerly_capped_subset": tally(unc, touched),
        "capped_all": tally(cap, every),
        "capped_never_capped_subset": tally(cap, every - touched),
        "capped_formerly_capped_subset": tally(cap, touched),
    }
    outcome, why, caveat = score(blocks["uncapped_all"], blocks["capped_all"],
                                 blocks["uncapped_never_capped_subset"])

    prereg = json.loads((METRICS / "r38_preregistration.json").read_text(
        encoding="utf-8"))

    generic = METRICS / "r38_verdict.json"
    superseded = METRICS / "r38_verdict_rev30_generic.json"
    if generic.exists() and not superseded.exists():
        superseded.write_text(generic.read_text(encoding="utf-8"),
                              encoding="utf-8")

    payload = {
        "schema": "r38_verdict_v1",
        "beta": b,
        "comparison": ("the radial rule against the work-matched constant "
                       "degree at equal realized total quadratic work; the "
                       "comparison the registration names"),
        "verdict_rule": prereg["verdict_rule"],
        "orbits_touching_ceiling_in_capped_run": len(touched),
        "blocks": blocks,
        "outcome": outcome,
        "outcome_text": prereg["outcomes"][outcome],
        "why": why,
        "caveat": caveat,
        "label_is_safe_to_quote_alone": caveat is None,
        "scoring_order_is_post_hoc": (
            "outcomes P and S overlap and Q and R share a clause; the order "
            "S, P, then R before Q was fixed after the numbers. No threshold "
            "was chosen after the numbers: the decade in R and the direction "
            "in P and S are the registered words."),
        "superseded_generic_record": {
            "file": superseded.name,
            "written_by": "rev30_verdict.py",
            "why_it_does_not_apply": (
                "it tallies the interior member against the work-matched "
                "constant degree, which is a different comparison, and it maps "
                "outcomes by sorted(outcomes)[:3]; this registration declares "
                "four, so the mapping lands on P by alphabetical position"),
        },
        "preregistration_sha256": prereg["preregistration_sha256"],
    }
    # One record per budget. A single r38_verdict.json would mean the last
    # budget scored silently erased the previous one, which is how a campaign
    # ends up quoting a number no file on disk supports.
    out = METRICS / f"r38_verdict_beta_{b:.2f}.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    index = {
        "schema": "r38_verdict_index_v1",
        "note": ("one record per budget; this file indexes them and does not "
                 "carry a verdict of its own"),
        "records": {f"beta_{x:.2f}": p.name for x in (1.00, 0.75, 0.50)
                    if (p := METRICS / f"r38_verdict_beta_{x:.2f}.json").exists()},
        "superseded_generic_record": payload["superseded_generic_record"],
    }
    (METRICS / "r38_verdict.json").write_text(json.dumps(index, indent=2),
                                              encoding="utf-8")

    print(f"radial vs work-matched constant degree, beta = {b:.2f}\n")
    hdr = f"{'block':<34}{'tally':>9}{'unres':>7}{'n':>5}{'rho_all':>10}{'rho_res':>10}{'work':>8}"
    print(hdr)
    for name in ("capped_all", "capped_never_capped_subset",
                 "capped_formerly_capped_subset", "uncapped_all",
                 "uncapped_never_capped_subset",
                 "uncapped_formerly_capped_subset"):
        t = blocks[name]
        ra = f"{t['median_rho_all_orbits']:.2f}" if t['median_rho_all_orbits'] else "n/a"
        rr = f"{t['median_rho_resolved']:.2f}" if t['median_rho_resolved'] else "n/a"
        print(f"{name:<34}"
              f"{t['resolved_radial_wins']:>5}-{t['resolved_constant_wins']:<3}"
              f"{t['unresolved']:>6}{t['orbits']:>5}{ra:>10}{rr:>10}"
              f"{t['median_realized_work_ratio_radial_over_constant']:>8.3f}")
    print(f"\nOUTCOME: {outcome}\n  {why}")
    if caveat:
        print(f"\n  !! CAVEAT: {caveat}")
    print(f"[written] r38_verdict.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
