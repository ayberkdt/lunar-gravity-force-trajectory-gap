"""Generate the manuscript table from the completed Rev-9 long-arc JSON."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "metrics" / "r9_potential_blend_longarc.json"
OUTPUT = ROOT / "metrics" / "r9_potential_blend_longarc_table.tex"


def main() -> int:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    if not payload.get("complete"):
        raise RuntimeError("Rev-9 long-arc result is not complete")
    rows = {row["name"]: row for row in payload["rows"]}
    ordered = (
        ("c2_50x300_polar", r"$50\times300$ polar"),
        ("c6_lro_30x216", r"LRO-like $30\times216$"),
    )
    lines = [
        r"\begin{tabular}{@{}lrrrr@{}}",
        r"  \toprule",
        r"  Geometry & fixed $N=120$ & discrete & accel. blend & corrected blend \\",
        r"  \midrule",
    ]
    for key, label in ordered:
        policies = rows[key]["policies"]
        values = [
            policies["fixed_N120"]["pos_rms_m"],
            policies["switch_N30_N120"]["pos_rms_m"],
            policies["blend_acceleration"]["pos_rms_m"],
            policies["blend_potential_corrected"]["pos_rms_m"],
        ]
        lines.append(
            f"  {label} & " + " & ".join(f"{value:.1f}" for value in values) + r" \\"
        )
    lines.extend([r"  \bottomrule", r"\end{tabular}", ""])
    OUTPUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"[written] {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
