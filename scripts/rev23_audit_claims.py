"""Mechanical audit of the Tur-5 review's arithmetic claims against the archive.

HISTORICAL. This script checks sentences the manuscript carried at the time of
the fifth review round. Most of them have since been rewritten or removed, so
its "contradicts the manuscript as written" lines now report on text that no
longer exists and must not be read as live findings: the current instruments
are claims_ledger.py and submission_gate.py. It is kept because the round it
audits is part of the record, not because its verdicts still hold.

The review verified four of its own findings independently and left the rest
resting on member reports. Before any of those are edited into the manuscript
they are checked here against the records, because a member report is exactly
as fallible as the text it is auditing, and an edit made on a wrong correction
is worse than the original error.

Each check prints the manuscript's wording, the archive's value, and a verdict.
Nothing is written; this only reads.

Usage:  python rev23_audit_claims.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PASS, FAIL, INFO = "PASS", "FAIL", "----"
_results: list[tuple[str, str, str]] = []


def load(name: str):
    p = METRICS / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def check(tag: str, where: str, says, archive, ok: bool | None,
          note: str = "") -> None:
    verdict = INFO if ok is None else (PASS if ok else FAIL)
    _results.append((verdict, tag, where))
    print(f"[{verdict}] {tag}  ({where})")
    print(f"        text    : {says}")
    print(f"        archive : {archive}")
    if note:
        print(f"        note    : {note}")


def span_rows(design: str, beta: float = 1.0):
    d = load(f"r18_span_sweep_{design}_beta_{beta:.2f}.json")
    return d["rows"] if d else []


# ---------------------------------------------------------------- finding 9
def finding_9() -> None:
    """'sixfold' is a ratio of medians quoted against a median of ratios."""
    for design in ("A", "B"):
        e_fix, e_int, ratios = [], [], []
        for r in span_rows(design):
            e0 = r["entries"]["0.00"].get("error_m")
            ek = r["entries"].get("0.50", {}).get("error_m")
            if not e0 or not ek:
                continue
            e_fix.append(e0)
            e_int.append(ek)
            ratios.append(e0 / ek)
        mor = float(np.median(ratios))
        rom = float(np.median(e_fix) / np.median(e_int))
        check(
            f"F9 per-call ratio, design {design}",
            "08_budget.tex:558,577; 09_conclusion.tex:68",
            "'roughly sixfold' / 'five- to sixfold'",
            f"median-of-ratios {mor:.2f}, ratio-of-medians {rom:.2f}",
            abs(mor - 6.0) < 0.5,
            "the realized-work figure it is compared against (2.56/2.58) is a "
            "median of ratios, so the comparable per-call number is "
            f"{mor:.2f}, not {rom:.2f}")


# --------------------------------------------------------------- finding 20
def finding_20() -> None:
    d = load("r21_gradient_sensitivity.json")
    if not d:
        return
    ep = sorted({int(r["epochs"]) for r in d["rows"]})
    check("F20 epoch count", "supp_budget_pareto.tex:174",
          "sampled at 241 epochs", f"epochs = {ep}",
          ep == [241],
          "241 was the target; step = len(t)//241 realizes 253")


# --------------------------------------------------------------- finding 19
def finding_19() -> None:
    d = load("r21_gradient_sensitivity.json")
    if not d:
        return
    s = d["summary"]
    v = s["max_neglected_over_forcing_perilune_ge_50km"]
    check("F19a ratio bound", "08_budget.tex (gradient paragraph)",
          "keeps the ratio below 0.61", f"{v:.6f}",
          v < 0.61, "0.62 is the smallest bound its own record supports")

    low = [(f"{r['design']}{r['sobol_index']:03d}", r["hp_km"],
            r["policies"]["atallah_budget"]["neglected_over_forcing_max"])
           for r in d["rows"] if r["hp_km"] < 50]
    detail = ", ".join(f"{n} (hp {h:.1f} km) -> {v:.3g}" for n, h, v in low)
    check("F19b the two 31 km orbits", "08_budget.tex (gradient paragraph)",
          "on the two 31 km-perilune orbits it reaches 40", detail,
          all(v > 30 for _, _, v in low),
          "only one of the two reaches 40; the other is over an order of "
          "magnitude smaller")


# ---------------------------------------------------------------- finding 7
def finding_7() -> None:
    d = load("r20_span_longarc.json")
    if not d:
        return
    bad = []
    for r in d["rows"]:
        e = r["entries"]
        if "1.00" not in e:
            continue
        worst = max(e, key=lambda k: e[k].get("error_m") or 0.0)
        if worst != "1.00":
            bad.append((r.get("name", r.get("sobol_index")), r["hp_km"],
                        {k: round(e[k]["error_m"], 1) for k in sorted(e)},
                        worst))
    lines = "; ".join(f"{n} hp={h:.0f}km worst={w} {vals}"
                      for n, h, vals, w in bad)
    check("F7 radial endpoint worst on the 60-day panel",
          "08_discussion.tex / long-arc paragraph",
          "the radial endpoint remains the worst by a wide margin",
          f"{len(bad)} of {len(d['rows'])} orbits contradict it: "
          f"{lines or 'none'}",
          not bad,
          "restrict the sentence to the median and name the exceptions")


# --------------------------------------------------------------- finding 25
def finding_25() -> None:
    for design in ("A", "B"):
        rows = span_rows(design)
        if not rows:
            continue
        lowest = sum(1 for r in rows
                     if r.get("best_k") not in (None, "0.00", "1.00"))
        beats_both = sum(1 for r in rows
                         if r.get("best_beats_both_endpoints_resolved"))
        check(f"F25 interior counts, design {design}",
              "main.tex:43; 09_conclusion.tex:58",
              "one of the three interior members has the lowest measured "
              "error on 51 of 64 orbits",
              f"interior has the lowest measured error on {lowest} of "
              f"{len(rows)}; beats both endpoints past the envelope on "
              f"{beats_both}",
              lowest == 51,
              "the body separates three counts and warns they answer "
              "different questions")


# --------------------------------------------------------------- finding 21
def finding_21() -> None:
    """The empirical schedule versus its work-matched fixed degree."""
    d = load("r11_manuscript_descriptives.json")
    if not d:
        return
    for label, key in (("A", "full64"), ("B", "designB")):
        block = (d.get(key) or {}).get("decisions", {})
        cmp_ = block.get("schedule_empirical_vs_fixed_work")
        if not cmp_:
            continue
        pairs = cmp_.get("pairs")
        res = cmp_.get("resolved")
        sched = cmp_.get("resolved_schedule_wins")
        fixed = (res - sched) if (res is not None and sched is not None) \
            else None
        says = ("63 resolve --- 54 for the fixed degree and 9 for the schedule"
                if label == "A" else
                "design B gives '55 of 64 ... the same way', with no split")
        check(f"F21 empirical schedule accounting, design {label}",
              "07_results.tex:688--690", says,
              f"pairs {pairs}, resolved {res} = fixed {fixed} + "
              f"schedule {sched}, unresolved "
              f"{pairs - res if (pairs and res is not None) else '?'}",
              None,
              "design A is given as a full split and design B is not; the "
              "review reads B's '55' as a resolve count when it is the "
              "fixed-win count")


# --------------------------------------------------------------- finding 28
def finding_28() -> None:
    d = load("r19_manuscript_descriptives.json")
    if not d:
        return
    for key, v in d.items():
        if not key.endswith("beta_1.00") or not v.get(
                "comparator_degree_shift"):
            continue
        s = v["comparator_degree_shift"]
        check(f"F28a comparator degree gain, {key}", "08_budget.tex:571",
              "the comparator gains a median 8 degrees",
              f"median {s['median']:.1f} [{s['min']}, {s['max']}]",
              abs(s["median"] - 8.0) < 0.25)

    # "at beta = 0.5 and beta = 1 the best member is k = 0.5 in all four
    # design-budget pairs, at a median span near 2--3"; all four pairs, not
    # just the two at beta = 1.
    med = {}
    for design in ("A", "B"):
        for beta in (0.5, 1.0):
            spans = [r["entries"]["0.50"]["span"]
                     for r in span_rows(design, beta)
                     if r["entries"].get("0.50", {}).get("span")]
            if spans:
                med[f"{design}@{beta:.2f}"] = float(np.median(spans))
    if med:
        vals = list(med.values())
        check("F28b interior span, all four design-budget pairs",
              "08_budget.tex:615",
              "at a median span near 2--3",
              ", ".join(f"{k} {v:.2f}" for k, v in med.items()),
              all(2.0 <= v <= 3.0 for v in vals),
              "the two beta = 1 pairs are inside 2--3 and the two beta = 0.5 "
              "pairs are not; the review's 3.19 and 3.09 are exactly those "
              "two. The sentence covers all four, so it needs the range "
              "widened rather than the numbers swapped")


def main() -> int:
    print("=" * 72)
    print("Tur-5 review claims checked against the archive")
    print("=" * 72)
    for fn in (finding_9, finding_20, finding_19, finding_7,
               finding_25, finding_21, finding_28):
        print()
        try:
            fn()
        except Exception as exc:                                # noqa: BLE001
            print(f"[{FAIL}] {fn.__name__} raised {type(exc).__name__}: {exc}")
    print()
    print("=" * 72)
    n_fail = sum(1 for v, _, _ in _results if v == FAIL)
    for verdict, tag, where in _results:
        print(f"  [{verdict}] {tag}")
    print(f"\n{n_fail} of {len(_results)} checks contradict the manuscript "
          f"as written.")
    print("A FAIL here means the review was right and the text needs the edit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
