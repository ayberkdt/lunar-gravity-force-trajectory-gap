"""Where does the reversal turn on, and do both solutions agree on where?

J1 finds the reversal at the declared budget and finds it absent at half that
budget: at beta = 0.5 the radial allocation loses the force metric as well as
the trajectory metric. That is not a contradiction, it is a statement about
where the reversal lives -- and it is only worth making if the primary solution
says the same thing under the same statistic.

The archived budget records already sweep the grid, but they summarize the
defect as an RMS. The campaigns use a time average. Comparing 25/64 from one
statistic against 3/32 from another would be comparing nothing, so the primary
field is swept again here with the campaign's statistic, on the archived
reference trajectories.

No propagation: the defect is a deterministic function of the reference
trajectory and the degree history, so the whole grid runs on trajectories that
already exist. Nothing here speaks to trajectory error, which is exactly why
the campaigns propagate the declared budget separately.

Usage:
    python revJ_budget_grid.py --workers 11
"""

from __future__ import annotations

import argparse
import json
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import revJ_common as J

J.select_field("JGGRX_1800F")
J.install_field()

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "rJ_budget_grid_primary.json"
LOG = Path(__file__).resolve().parent / "rJ_budget_grid.log"

PARETO = METRICS / "r14_budget_pareto.json"
ARCHIVED = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "conv": METRICS / "r11_raw" / "convergence"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "conv": METRICS / "r11_raw" / "designB_convergence"},
}
BUDGET_GRID = (0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00)


def log(msg: str) -> None:
    J.log_line(LOG, f"budget-grid {msg}")


def task(payload: dict) -> dict:
    design, index = payload["design"], int(payload["index"])
    try:
        g = payload["geom"]
        adopted = int(g["adopted_truth_degree"])
        raw = (ARCHIVED[design]["conv"] / f"sobolA_{index:03d}"
               / "truth_tighter.npz")
        if not raw.exists():
            return {"ok": False, "design": design, "index": index,
                    "error": f"missing archived reference {raw}"}
        t, Y = J.load_times(raw), J.load_states(raw)
        model, args = J.model_for(adopted)
        h_km = (np.linalg.norm(Y[:3], axis=0) - model.r_ref) / 1e3
        n_crit = int(g["n_critical"])

        budgets, degree_sets = {}, {}
        for beta in BUDGET_GRID:
            cal = J.calibrate_radial(adopted, g["hp_km"], g["ha_km"], adopted,
                                     h_km, beta * n_crit ** 2)
            table = {float(k): int(v) for k, v in cal["table"].items()}
            deg_r = J.degrees_from_table(table, h_km)
            n_f, capped = J.fixed_degree_for(beta, n_crit, adopted)
            degree_sets[f"radial_{beta:.2f}"] = deg_r
            degree_sets[f"constant_{beta:.2f}"] = np.full(len(h_km), n_f,
                                                          dtype=int)
            budgets[f"beta_{beta:.2f}"] = {
                # N = N_ref is censored as well as N > N_ref: a policy that has
                # reached the reference has no truncation error to measure.
                "beta_requested": beta,
                "censored": bool(capped or int(deg_r.max()) >= adopted
                                 or n_f >= adopted),
                "work_mismatch": cal["mismatch"],
                "constant_degree": int(n_f),
                "radial_degree_span": [int(deg_r.min()), int(deg_r.max())]}
        defects = J.force_defects(t, Y[:3], degree_sets, adopted, args)
        for beta in BUDGET_GRID:
            key = f"beta_{beta:.2f}"
            dr = defects[f"radial_{beta:.2f}"]
            dc = defects[f"constant_{beta:.2f}"]
            den_mean = dc["J_force_mean_m_s2"]
            den_rms = dc["J_force_rms_m_s2"]
            budgets[key]["rho_force_mean"] = (dr["J_force_mean_m_s2"] / den_mean
                                              if den_mean > 0.0 else None)
            budgets[key]["rho_force_rms"] = (dr["J_force_rms_m_s2"] / den_rms
                                             if den_rms > 0.0 else None)
            if den_mean <= 0.0:
                budgets[key]["censored"] = True
        return {"ok": True, "design": design, "index": index,
                "hp_km": g["hp_km"], "n_critical": n_crit,
                "adopted_truth_degree": adopted, "budgets": budgets}
    except Exception:
        return {"ok": False, "design": design, "index": index,
                "error": traceback.format_exc()}


def summarize(rows: list[dict]) -> dict:
    out = {}
    for beta in BUDGET_GRID:
        key = f"beta_{beta:.2f}"
        live = [r["budgets"][key] for r in rows
                if not r["budgets"][key]["censored"]
                and r["budgets"][key]["rho_force_mean"] is not None]
        if not live:
            continue
        out[key] = {
            "orbits_uncensored": len(live),
            "orbits_censored": len(rows) - len(live),
            "radial_wins_force": sum(1 for b in live
                                     if b["rho_force_mean"] < 1.0),
            "median_rho_force_mean": float(np.median(
                [b["rho_force_mean"] for b in live])),
            "median_rho_force_rms": float(np.median(
                [b["rho_force_rms"] for b in live])),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()

    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    tasks = []
    for design in ("A", "B"):
        src = json.loads(ARCHIVED[design]["rows"].read_text(encoding="utf-8"))
        crit = {int(r["sobol_index"]): int(r["n_critical"])
                for r in src["rows"]}
        for row in src["rows"]:
            index = int(row["sobol_index"])
            g = row.get("design_point", row)
            tasks.append({"design": design, "index": index, "geom": {
                "hp_km": float(g["hp_km"]), "ha_km": float(g["ha_km"]),
                "n_critical": crit[index],
                "adopted_truth_degree": int(row["adopted_truth_degree"])}})
    log(f"{len(tasks)} archived orbits over {len(BUDGET_GRID)} budgets")

    results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fut in as_completed([pool.submit(task, t) for t in tasks]):
            results.append(fut.result())
    good = [r for r in results if r.get("ok")]
    bad = [r for r in results if not r.get("ok")]

    payload = {
        "schema": "rJ_budget_grid_primary_v1", "created_utc": J.utc_now(),
        "field": J.field_key(),
        "purpose": ("sweep the primary solution over the budget grid under the "
                    "campaign's own force statistic, so the cross-solution "
                    "budget dependence is compared like with like"),
        "note": "force metric only; no propagation is involved",
        "budget_grid": list(BUDGET_GRID),
        "definition": "rho_force = J_force(radial)/J_force(constant), "
                      "J_force being the time average of the defect magnitude",
        "censoring_rule": ("a budget is censored for an orbit when the policy "
                           "reaches the adopted reference degree"),
        "counts": {"orbits": len(good), "failed": len(bad)},
        "summary_by_design": {d: summarize([r for r in good
                                            if r["design"] == d])
                              for d in ("A", "B")},
        "summary_pooled": summarize(good),
        "complete": bool(good and not bad),
        "rows": sorted(good, key=lambda r: (r["design"], r["index"])),
        "failures": bad,
        "provenance": J.provenance(),
    }
    J.atomic_json(OUT, payload)
    for key, s in payload["summary_pooled"].items():
        log(f"{key}: radial wins force {s['radial_wins_force']}/"
            f"{s['orbits_uncensored']}, median rho_F(mean) "
            f"{s['median_rho_force_mean']:.3f}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
