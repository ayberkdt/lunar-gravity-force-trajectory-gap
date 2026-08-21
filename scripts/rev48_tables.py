"""Per-orbit table for the R48 (O48) measured-time panel.

Emits metrics/r48_interior_timing_table.tex from the frozen campaign record,
so no number is transcribed by hand. One row per panel orbit: the
serial-time-matched comparator degree, the achieved total-kernel-time ratio,
the tighter-level error ratio and the envelope verdict.

Usage:  python rev48_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

VERDICT = {"interior": "interior", "fixed": "fixed", None: "--"}


def main() -> int:
    d = json.loads((METRICS / "r48_interior_timing.json").read_text(
        encoding="utf-8"))
    lines = ["\\begin{tabular}{@{}l r r r r l@{}}", "\\toprule",
             "orbit & $h_p$ (km) & $N_{\\mathrm{time}}$ & time match & "
             "$\\rho$ & verdict \\\\",
             "\\midrule"]
    for r in d["rows"]:
        lines.append(
            f"{r['design']}{r['sobol_index']:03d} & {r['hp_km']:.1f} & "
            f"{r['comparator_degree']} & {r['achieved_time_ratio']:.3f} & "
            f"{r['rho_fixed_over_member']:.2f} & "
            f"{VERDICT[r['winner']]} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    (METRICS / "r48_interior_timing_table.tex").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    s = d["summary"]
    a = s["achieved_time_ratio"]
    print(f"[written] r48_interior_timing_table.tex  "
          f"({s['orbits']} orbits, {s['resolved']} resolved "
          f"{s['interior_wins']}--{s['fixed_wins']}, time match "
          f"{a['median']:.3f} [{a['min']:.3f},{a['max']:.3f}])")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
