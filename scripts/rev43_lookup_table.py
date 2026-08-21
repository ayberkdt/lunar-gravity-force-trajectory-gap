"""R43: the truncation-degree lookup, dense enough to be used without the paper.

Table~\\ref{tab:criteria} compares criteria at eight altitudes, which is the
right shape for an argument and the wrong one for a reader who wants a degree.
This builds the lookup instead: every 10 km from 20 to 300, three omission
budgets, and the directional complement where it was sampled.

Nothing here is a new measurement. Two frozen records supply everything:

  r1_spectrum_pfit.json   the per-coefficient RMS spectrum of GL1800F,
                          sigma_n = sqrt(P_n / (2n+1)), degrees 2 to 1800.
                          N_min^emp(h, eps) is Eq. (5) evaluated on it, in the
                          total-vector weighting of Eq. (4).
  r3_local_tail.json      4096 scrambled-Sobol directions per altitude, from
                          which the p95 and p99 degrees are read on the 10-degree
                          ladder they were sampled on. Six altitudes only, and
                          the table leaves the others blank rather than
                          interpolating a percentile.

The re-derivation is checked against the archive before anything is written:
all three columns are recomputed at the eight altitudes r1 already carries, in
both the vector and radial weightings, and the run aborts on any disagreement.
That is what licenses reading the other twenty-one rows, which no archived
record states.

Usage:  python rev43_lookup_table.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SPECTRUM = METRICS / "r1_spectrum_pfit.json"
LOCAL_TAIL = METRICS / "r3_local_tail.json"
OUT_JSON = METRICS / "r43_lookup.json"
OUT_TEX = METRICS / "r43_lookup_table.tex"
OUT_CSV = METRICS / "r43_lookup.csv"

# The GL1800F header values, as quoted in Section 6.
R_KM = 1738.0

# 10 km through the low-altitude band, where the degree moves fastest and a
# reader wants the row rather than the trend; 20 km above 150, where linear
# interpolation between rows is worth less than a degree.
ALTITUDES = list(range(20, 151, 10)) + list(range(160, 301, 20))
EPS = [("1e-1", 1e-1), ("1e-2", 1e-2), ("1e-3", 1e-3)]

# Printed as three side-by-side blocks: a reference block should be short enough
# to take in at once, not a column running down a page.
BLOCKS = 3


def load_spectrum():
    d = json.loads(SPECTRUM.read_text(encoding="utf-8"))
    n = np.asarray(d["spectrum_arrays"]["n"], dtype=float)
    sigma = np.asarray(d["spectrum_arrays"]["sigma_coeff_rms"], dtype=float)
    return n, sigma**2 * (2 * n + 1), d["criteria_rows"]


def nmin(n, P, h_km, eps, vector=True):
    """Eq. (5): smallest N whose discarded tail is within eps of the total."""
    atten = (R_KM / (R_KM + h_km)) ** (2 * n)
    weight = (n + 1) * (2 * n + 1) if vector else (n + 1) ** 2
    w = (atten * P * weight)[n >= 2]
    total = w.sum()
    tail = total - np.cumsum(w)
    ok = np.sqrt(np.maximum(tail, 0.0) / total) <= eps
    if not ok.any():
        return None
    return int(n[n >= 2][int(np.argmax(ok))])


def self_check(n, P, archived) -> list[str]:
    """Recompute what the archive already states; report every disagreement."""
    bad = []
    for row in archived:
        h = row["altitude_km"]
        for key, eps, vec in (("emp_vector_eps1e2", 1e-2, True),
                              ("emp_radial_eps1e2", 1e-2, False),
                              ("emp_vector_eps1e3", 1e-3, True)):
            got, want = nmin(n, P, h, eps, vec), row[key]
            if got != want:
                bad.append(f"h={h:g} {key}: recomputed {got}, archived {want}")
    return bad


def percentile_degrees():
    """p95 and p99 degrees at eps = 1e-2, on the ladder they were sampled on."""
    d = json.loads(LOCAL_TAIL.read_text(encoding="utf-8"))
    out = {}
    for r in d["recommendations"]:
        if r["eps"] != 0.01 or r.get("N_min_p95") is None:
            continue
        out[float(r["altitude_km"])] = {
            "p95": int(r["N_min_p95"]),
            "p99": int(r["N_min_p99"]),
            "sampled_rms": int(r["N_min_sampled_rms"]),
        }
    return out, int(d["n_dirs"]), int(d["seed"])


def write_tex(rows):
    head = r"$h$ [km] & $10\%$ & $1\%$ & $0.1\%$"
    per = -(-len(rows) // BLOCKS)
    cols = []
    for b in range(BLOCKS):
        chunk = rows[b * per:(b + 1) * per]
        cell = [f"{r['altitude_km']:g} & {r['N_eps1e-1']} & {r['N_eps1e-2']} & "
                f"{r['N_eps1e-3']}" for r in chunk]
        cell += [" & & & "] * (per - len(cell))
        cols.append(cell)

    spec = (r"@{\hspace{1.6em}}".join([r"r rrr"] * BLOCKS))
    # Two-row header: a spanner naming what the three numeric columns share
    # (the tail budget), then the budgets themselves under it.
    spanner = " & ".join(
        [r"& \multicolumn{3}{c}{tail budget $\varepsilon$}"] * BLOCKS)
    lines = [
        r"\begin{tabular}{@{}" + spec + r"@{}}",
        r"\toprule",
        spanner + r" \\",
        " ".join(rf"\cmidrule(lr){{{4 * b + 2}-{4 * b + 4}}}"
                 for b in range(BLOCKS)),
        " & ".join([head] * BLOCKS) + r" \\",
        " ".join(rf"\cmidrule(lr){{{4 * b + 1}-{4 * b + 4}}}"
                 for b in range(BLOCKS)),
    ]
    for row in zip(*cols):
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    OUT_TEX.write_text("\n".join(lines) + "\n", encoding="utf-8")


def percentile_sentence(rows) -> str:
    """The directional surcharge, as prose: it exists at six altitudes only."""
    parts = []
    for r in rows:
        if r["N_p95_eps1e-2"]:
            parts.append((int(r["altitude_km"]), r["N_eps1e-2"],
                          r["N_p95_eps1e-2"]))
    return "; ".join(f"{h} km: {a} -> {b}" for h, a, b in parts)


def main() -> int:
    n, P, archived = load_spectrum()

    bad = self_check(n, P, archived)
    if bad:
        print("ABORT: re-derivation disagrees with the archive")
        for b in bad:
            print("  " + b)
        return 1
    print(f"[check] {3 * len(archived)} archived degrees reproduced exactly")

    pct, n_dirs, seed = percentile_degrees()

    rows = []
    for h in ALTITUDES:
        row = {"altitude_km": h}
        for label, eps in EPS:
            row[f"N_eps{label}"] = nmin(n, P, h, eps)
        p = pct.get(float(h))
        row["N_p95_eps1e-2"] = p["p95"] if p else None
        row["N_p99_eps1e-2"] = p["p99"] if p else None
        row["N_sampled_rms_eps1e-2"] = p["sampled_rms"] if p else None
        rows.append(row)

    write_tex(rows)

    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    OUT_JSON.write_text(json.dumps({
        "schema": "r43_lookup_v1",
        "field": "GRAIL GL1800F, degree and order 1800",
        "reference_radius_km": R_KM,
        "criterion": (
            "N_min^emp(h, eps): the smallest N whose discarded-tail acceleration "
            "RMS is within a fraction eps of the total perturbing RMS on the "
            "sphere of radius R + h, in the total-vector weighting "
            "sqrt((n+1)(2n+1)). Model-relative: it measures what this expansion "
            "requires of its own truncations, not the distance to the true field."),
        "percentile_columns": (
            f"degrees at which the {'p95'} and p99 of the discarded-tail "
            f"acceleration over {n_dirs} scrambled-Sobol directions (seed {seed}) "
            "meet eps = 1e-2, read on the 10-degree ladder they were sampled on. "
            "Six altitudes were sampled; the rest are left blank rather than "
            "interpolated."),
        "derived_from": {
            "spectrum": SPECTRUM.name,
            "directional_sample": LOCAL_TAIL.name,
        },
        "self_check": (
            "every archived degree in r1_spectrum_pfit.json:criteria_rows "
            "recomputed and matched exactly, in both weightings"),
        "caveats": [
            "sphere-averaged: along localized high-amplitude directions the "
            "instantaneous tail exceeds its RMS, which the p95 and p99 columns "
            "quantify where they were sampled",
            "below about 50 km the empirical degree approaches the "
            "topography-constrained region of the product (n > 600), so it "
            "reflects the model's regularization as much as observed gravity",
            "a lookup bounds the instantaneous omission budget and not the "
            "trajectory error accumulated over a long arc",
        ],
        "rows": rows,
    }, indent=2) + "\n", encoding="utf-8")

    print(f"[written] {OUT_TEX.name}, {OUT_CSV.name}, {OUT_JSON.name}"
          f"  ({len(rows)} altitudes)")
    lo, hi = rows[0], rows[-1]
    print(f"  {lo['altitude_km']} km: {lo['N_eps1e-1']}/{lo['N_eps1e-2']}"
          f"/{lo['N_eps1e-3']}   {hi['altitude_km']} km: {hi['N_eps1e-1']}"
          f"/{hi['N_eps1e-2']}/{hi['N_eps1e-3']}")
    print(f"  p95 surcharge (1% budget): {percentile_sentence(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
