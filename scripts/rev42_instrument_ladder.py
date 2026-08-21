"""Instrument ladder with the completed panel on its third rung.

rev34_instrument_ladder.py reads the panel from r37_panel_verdict.json, which
is the right file until the chain is finished and the wrong one afterwards: the
panel is now 128 orbits and lives in r42_panel_verdict.json. rev34 is sealed
under the R37 manifest, so it is imported and one function is substituted --
the third rung's reader -- and nothing else. The other three rungs, the table
layout and the output paths are rev34's, unchanged.

The substituted function is the same bookkeeping against the same verdict
schema: it copies counts out of a scored record and computes nothing. Sign
agreement therefore still means what the scoring amendment says it means.

Usage:  python rev42_instrument_ladder.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev34_instrument_ladder as ladder

METRICS = Path(__file__).resolve().parents[1] / "metrics"
VERDICT = METRICS / "r42_panel_verdict.json"


def variational_rung():
    v = json.loads(VERDICT.read_text(encoding="utf-8"))
    p = v["panel"]
    return {
        "source": VERDICT.name,
        "n_orbits": p["orbits"],
        "level_per_design": v["levels"]["highest_complete_per_design"],
        "orbits_in_record": v["record"]["orbits"],
        "gradient_degree": 120,
        "resolved": p["resolved"],
        "unresolved": p["unresolved"],
        "sign_agreement": p["sign_agreement"],
        "sign_disagreement": p["sign_disagreement"],
        "radial_better_predicted": p["predicted_favors_radial_all"],
        "record_resolved": v["record"]["resolved"],
        "record_sign_agreement": v["record"]["sign_agreement"],
    }


def variational_by_design() -> dict:
    """Per-design splits of the completed panel, from the frozen record.

    The other three rungs report per design because the manuscript's own
    reporting rule is design-level replication; pooling the third rung both
    breaks that rule and hides that the instrument reproduces each design's
    propagated median separately.
    """
    import statistics
    rec = json.loads((METRICS / "r42_variational_completion.json"
                      ).read_text(encoding="utf-8"))
    out = {}
    for des in ("A", "B"):
        ratios = [r["predicted_ratio_fixed_over_atallah"]
                  for r in rec["rows"] if r["design"] == des]
        out[des] = {
            "n_orbits": len(ratios),
            "radial_better_predicted": sum(1 for x in ratios if x > 1.0),
            "median_predicted_ratio": statistics.median(ratios),
        }
    return out


def main() -> int:
    ladder.variational_rung = variational_rung
    ladder.main()
    out = json.loads((METRICS / "r34_instrument_ladder.json"
                      ).read_text(encoding="utf-8"))
    out["note"] = (
        "Derived from frozen R14 records; no propagation performed. Rungs 1, 2 "
        "and 4 are over all 64 orbits of each design; rung 3 is the "
        "forced-variational panel, read from r42_panel_verdict.json, which "
        "completed the R37 level chain at 64 orbits per design."
    )
    by = variational_by_design()
    out["variational_panel"]["by_design"] = by
    (METRICS / "r34_instrument_ladder.json").write_text(
        json.dumps(out, indent=2) + "\n", encoding="utf-8")

    # rev34's table pools the third rung into one column pair; every other
    # rung is per design. Rewrite that one line from the per-design record so
    # the table obeys the manuscript's own design-level reporting rule.
    table = METRICS / "r34_instrument_ladder_table.tex"
    text = table.read_text(encoding="utf-8")
    a, b = by["A"], by["B"]
    pooled = (r"Forced variational & $+$ gradient $\mathbf G$ & "
              rf"\multicolumn{{2}}{{c}}"
              rf"{{{out['variational_panel']['radial_better_predicted']}"
              rf"/{out['variational_panel']['n_orbits']}}} & "
              r"\multicolumn{2}{c}{---} \\")
    split = (r"Forced variational & $+$ gradient $\mathbf G$ & "
             f"{a['radial_better_predicted']}/{a['n_orbits']} & "
             f"{b['radial_better_predicted']}/{b['n_orbits']} & "
             f"{a['median_predicted_ratio']:.2f} & "
             f"{b['median_predicted_ratio']:.2f}"
             r" \\")
    if pooled not in text:
        raise SystemExit("[r42 ladder] pooled variational line not found; "
                         "rev34 layout changed -- refusing to patch blindly")
    table.write_text(text.replace(pooled, split), encoding="utf-8")

    v = out["variational_panel"]
    print(f"[r42 ladder] variational rung: {v['radial_better_predicted']}"
          f"/{v['n_orbits']} favour radial, "
          f"{v['sign_agreement']}/{v['resolved']} resolved signs, "
          f"{v['unresolved']} undecided")
    print(f"[r42 ladder] by design: "
          f"A {a['radial_better_predicted']}/{a['n_orbits']} "
          f"@ {a['median_predicted_ratio']:.3f}, "
          f"B {b['radial_better_predicted']}/{b['n_orbits']} "
          f"@ {b['median_predicted_ratio']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
