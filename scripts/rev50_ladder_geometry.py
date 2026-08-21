"""What raising the apolune actually changes, level by level.

The main text says the paired ladder does not isolate one variable: across its
four apolune levels the eccentricity, the orbital period and the number of
revolutions inside the seven-day horizon all move together, and that is the
paper's own answer to a reader who reads the ladder as a clean radial-span
experiment. Those three numbers were quoted from the frozen designs without a
table behind them, which put them in the one class the manuscript otherwise
does not have: a printed number with no carrier a reader can check.

This builds the carrier. It reads the two blocks' frozen designs, which contain
every orbit's semimajor axis and eccentricity as drawn, and reports the median
of each quantity per apolune level over the two blocks together. Pooling the
blocks is the one pooling the ladder's registration licenses, and it is what
the main text's numbers already do; nothing here is a verdict, so no tally is
combined.

Period and revolution count are computed rather than stored: the period from
the two-body relation at the archived semimajor axis, and the revolutions from
the seven-day horizon divided by it. Neither is a propagated quantity and
neither is used as one; they describe the geometry the ladder draws.

Usage:  python rev50_ladder_geometry.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

# the value the whole paper reads from the GL1800F header, in m^3/s^2
MU = 4902.8003063302e9
HORIZON_S = 7 * 86400.0


def orbits() -> list[dict]:
    out = []
    for blk in ("a", "b"):
        p = METRICS / f"r50_span_ladder_{blk}_design_frozen.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        for o in d["orbits"]:
            out.append({"block": blk.upper(),
                        "level_km": float(o["apolune_level_km"]),
                        "ecc": float(o["eccentricity"]),
                        "a_m": float(o["semimajor_axis_m"])})
    return out


def main() -> int:
    rows = orbits()
    levels = sorted({r["level_km"] for r in rows})
    table = []
    for lvl in levels:
        sel = [r for r in rows if r["level_km"] == lvl]
        periods = [2.0 * math.pi * math.sqrt(r["a_m"] ** 3 / MU) for r in sel]
        table.append({
            "apolune_level_km": lvl,
            "orbits": len(sel),
            "median_eccentricity": median(r["ecc"] for r in sel),
            "median_period_h": median(periods) / 3600.0,
            "median_revolutions_7d": median(HORIZON_S / p for p in periods),
        })

    rec = {"schema": "r50_ladder_geometry_v1",
           "what_this_is": (
               "median eccentricity, orbital period and seven-day revolution "
               "count per apolune level of the paired ladder, over the two "
               "identity blocks together. Period and revolution count are "
               "computed from the archived semimajor axis by the two-body "
               "relation; nothing is propagated."),
           "mu_m3_s2": MU, "horizon_s": HORIZON_S,
           "sources": ["r50_span_ladder_a_design_frozen.json",
                       "r50_span_ladder_b_design_frozen.json"],
           "levels": table}
    (METRICS / "r50_ladder_geometry.json").write_text(
        json.dumps(rec, indent=1), encoding="utf-8")

    lines = ["\\begin{tabular}{@{}r r r r r@{}}", "\\toprule",
             "apolune level [km] & orbits & median $e$ & median period [h] & "
             "median revs.\\ in 7~d \\\\", "\\midrule"]
    for t in table:
        lines.append(f"    {t['apolune_level_km']:.0f} & {t['orbits']} & "
                     f"{t['median_eccentricity']:.2f} & "
                     f"{t['median_period_h']:.2f} & "
                     f"{t['median_revolutions_7d']:.0f} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    (METRICS / "r50_ladder_geometry_table.tex").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")

    print("[written] r50_ladder_geometry.json, r50_ladder_geometry_table.tex")
    for t in table:
        print(f"  {t['apolune_level_km']:6.0f} km  n={t['orbits']:3d}  "
              f"e={t['median_eccentricity']:.4f}  "
              f"T={t['median_period_h']:.3f} h  "
              f"revs={t['median_revolutions_7d']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
