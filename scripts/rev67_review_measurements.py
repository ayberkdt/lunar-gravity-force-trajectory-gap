"""Two measurements a review round asked for, taken from sealed inputs only.

Neither is a new experiment. Both re-read records that already exist and answer
a question the manuscript stated qualitatively, so that the sentence carrying
the answer has a record behind it like every other number in the paper.

**Rank association (R14).** The R14 mechanism block records Pearson
correlations of log10 rho_budget against the per-orbit degree span and the
binned switch count. Pearson answers a narrower question than the sentence it
was used to support: rho_budget runs over eight decades per orbit, so a few
orbits at the ends of that range can set the sign and the size. This script
reproduces the Pearson values from the same archived inputs --- and refuses to
write if they disagree with the R14 record, because a disagreement would mean
the rows were joined differently --- then adds Spearman's rho with its
two-sided p-value.

**Undecided orbits by perilune (R42).** The forced-variational panel leaves 28
of its 128 comparisons undecided. The supplement described them as spread
across the perilune range rather than concentrated at its low end. This counts
them in perilune tertiles of the panel itself, so the description can be
replaced by the distribution.

Nothing is propagated and nothing is re-scored.

Usage:  python rev67_review_measurements.py
"""

from __future__ import annotations

import datetime
import io
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

BUDGETS = ["0.50", "0.75", "1.00", "1.25", "1.50", "2.00", "3.00"]


def load(name: str):
    return json.load(io.open(METRICS / name, encoding="utf-8"))


def columns(pareto: dict, design: str, beta: str):
    """The three columns the contrast is made of, joined on the Sobol index."""
    path = METRICS / f"r14_trajectory_{design}_beta_{beta}.json"
    if not path.exists():
        return None
    traj = json.load(io.open(path, encoding="utf-8"))
    alloc = {r["sobol_index"]: r["budgets"].get(f"beta_{beta}")
             for r in pareto["designs"][design]["rows"]}
    y, span, switches = [], [], []
    for r in traj["rows"]:
        i, c = r["sobol_index"], r["comparison"]
        b = alloc.get(i)
        if b is None or not c["rho_budget"] or c["rho_budget"] <= 0:
            continue
        a = b["atallah"]["allocation"]
        y.append(math.log10(c["rho_budget"]))
        span.append(a["max_degree"] / max(a["min_degree"], 1))
        switches.append(a["switch_count_binned"])
    if len(y) < 3:
        return None
    return np.array(y), np.array(span, float), np.array(switches, float)


def rank_association() -> tuple[dict, float]:
    pareto = load("r14_budget_pareto.json")
    record = load("r14_descriptives.json")["mechanism"]
    out: dict = {}
    worst_gap = 0.0
    for design in ("A", "B"):
        out[design] = {}
        for beta in BUDGETS:
            cols = columns(pareto, design, beta)
            if cols is None:
                continue
            y, span, switches = cols
            entry = {"orbits": int(len(y))}
            for name, x in (("degree_span", span), ("switches", switches)):
                r = float(np.corrcoef(y, x)[0, 1])
                rho = stats.spearmanr(x, y)
                entry[name] = {"pearson": r,
                               "spearman": float(rho.statistic),
                               "spearman_p": float(rho.pvalue)}
                ref = (record.get(design, {}).get(f"beta_{beta}", {})
                       .get("correlations_with_log_rho", {}).get(name))
                if ref is not None:
                    entry[name]["r14_record_pearson"] = ref
                    worst_gap = max(worst_gap, abs(ref - r))
            out[design][f"beta_{beta}"] = entry
    return out, worst_gap


def undecided_by_perilune() -> dict:
    verdict = load("r42_panel_verdict.json")
    undecided = [o["hp_km"] for o in verdict["unresolved"]["orbits"]]
    panel = []
    for design in ("A", "B"):
        traj = load(f"r14_trajectory_{design}_beta_1.00.json")
        panel += [r["design_point"]["hp_km"] for r in traj["rows"]]
    panel_a, und = np.array(panel), np.array(undecided)
    lo, hi = np.quantile(panel_a, [1 / 3, 2 / 3])
    bands = []
    for a, b, name in ((-np.inf, lo, "low"), (lo, hi, "mid"), (hi, np.inf, "high")):
        total = int(((panel_a > a) & (panel_a <= b)).sum())
        n = int(((und > a) & (und <= b)).sum())
        bands.append({"tertile": name, "hp_km_upper": None if np.isinf(b) else float(b),
                      "orbits": total, "undecided": n,
                      "rate": n / total if total else None})
    return {
        "panel_orbits": int(len(panel_a)),
        "undecided": int(len(und)),
        "perilune_km": {"min": float(panel_a.min()), "max": float(panel_a.max()),
                        "median": float(np.median(panel_a))},
        "tertile_edges_km": [float(lo), float(hi)],
        "bands": bands,
        "undecided_below_panel_median": int((und < np.median(panel_a)).sum()),
    }


def main() -> int:
    ranks, gap = rank_association()
    if gap > 5e-9:
        raise SystemExit(f"Pearson disagrees with the R14 record by {gap:.2e}; "
                         "the join differs and the rank values cannot be trusted")
    peri = undecided_by_perilune()

    payload = {
        "schema": "r67_review_measurements_v1",
        "created_utc": datetime.datetime.now(datetime.timezone.utc)
                        .isoformat().replace("+00:00", "Z"),
        "scope": ("Two re-readings of sealed records asked for by a review "
                  "round: the rank form of the R14 span/switch association, and "
                  "the perilune distribution of the R42 panel's undecided "
                  "comparisons."),
        "propagation": "none",
        "inputs": ["r14_budget_pareto.json", "r14_descriptives.json",
                   "r14_trajectory_{A,B}_beta_*.json", "r42_panel_verdict.json"],
        "pearson_agrees_with_r14_record_to": gap,
        "rank_association": ranks,
        "undecided_by_perilune": peri,
    }
    path = METRICS / "r67_review_measurements.json"
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] metrics/{path.name}  "
          f"(Pearson reproduces the R14 record to {gap:.1e})")
    for design in ("A", "B"):
        e = ranks[design]["beta_1.00"]
        print(f"  design {design} at beta = 1: span r = "
              f"{e['degree_span']['pearson']:+.3f}, rank "
              f"{e['degree_span']['spearman']:+.3f} (p = "
              f"{e['degree_span']['spearman_p']:.2f}); switches r = "
              f"{e['switches']['pearson']:+.3f}, rank "
              f"{e['switches']['spearman']:+.3f} (p = "
              f"{e['switches']['spearman_p']:.2f})")
    for b in peri["bands"]:
        print(f"  undecided, {b['tertile']:>4} perilune tertile: "
              f"{b['undecided']} of {b['orbits']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
