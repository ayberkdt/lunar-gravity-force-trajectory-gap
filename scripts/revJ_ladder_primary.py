"""The two metrics, budget by budget, on the primary solution.

The cross-solution campaign propagates its population at several budgets and
finds the sharpest form of this paper's point at the crossing: where the force
metric can no longer tell the two policies apart, the trajectory metric
separates them by more than an order of magnitude. That statement should not
rest on the second gravity solution alone, and it does not have to: the primary
field's trajectories at those budgets were propagated long ago and are in the
archive. Nothing new is integrated here.

What is new is the scoring. The archived campaign reports its force metric as an
RMS and its ratio in the reciprocal convention, so the archived summaries cannot
be laid beside the campaign's without inverting one of them. Both are recomputed
here with the campaign's definitions, from the archived trajectories:

    rho_force = J_force(radial)/J_force(constant),   J_force = time average
    rho_traj  = J_traj(radial)/J_traj(constant),     J_traj  = position RMS

and the two-level resolution rule is applied exactly as elsewhere.

Usage:
    python revJ_ladder_primary.py --workers 11
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

import rev12_atallah as at                                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUT = METRICS / "rJ_ladder_primary.json"
LOG = Path(__file__).resolve().parent / "rJ_ladder_primary.log"

PARETO = METRICS / "r14_budget_pareto.json"
BUDGETS = (0.50, 0.75, 1.00, 1.25)
ARCHIVED = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "conv": METRICS / "r11_raw" / "convergence"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "conv": METRICS / "r11_raw" / "designB_convergence"},
}


def log(msg: str) -> None:
    J.log_line(LOG, f"ladder {msg}")


def _raw(design: str, beta: float, index: int, policy: str,
         level: str) -> Path | None:
    name = f"sobolA_{index:03d}"
    conv = ARCHIVED[design]["conv"]
    budget = METRICS / "r14_raw" / f"{design}_beta_{beta:.2f}" / name
    if policy == "reference":
        p = conv / name / f"truth_{level}.npz"
    elif policy == "radial":
        p = budget / f"atallah_budget_{level}.npz"
    else:
        p = budget / f"fixed_budget_{level}.npz"
        if not p.exists():
            # At beta = 1 the comparator is the critical degree itself and the
            # archived campaign reuses that run rather than repropagating it.
            p = conv / name / f"fixed_critical_{level}.npz"
    return p if p.exists() else None


def task(payload: dict) -> dict:
    design, index, beta = payload["design"], payload["index"], payload["beta"]
    try:
        g, budget = payload["geom"], payload["budget"]
        adopted = int(g["adopted_truth_degree"])
        paths = {(p, lv): _raw(design, beta, index, p, lv)
                 for p in ("reference", "constant", "radial")
                 for lv in ("tight", "tighter")}
        if any(v is None for v in paths.values()):
            return {"ok": False, "design": design, "index": index, "beta": beta,
                    "reason": "incomplete archived record at this budget"}
        Y = {k: J.load_states(v) for k, v in paths.items()}
        t = J.load_times(paths[("reference", "tighter")])
        n = min(v.shape[1] for v in Y.values())

        model, args = J.model_for(adopted)
        R = Y[("reference", "tighter")][:3, :n]
        h_km = (np.linalg.norm(R, axis=0) - model.r_ref) / 1e3
        _, table = at.atallah_binned_schedule(
            model, J.atallah_g(adopted),
            float(budget["atallah"]["tol_accel_m_s2"]),
            g["hp_km"], g["ha_km"], floor=J.FLOOR, cap=adopted,
            bin_km=J.BIN_KM)
        table = {float(k): int(v) for k, v in table.items()}
        n_f = int(budget["fixed"]["degree"])
        deg_r = J.degrees_from_table(table, h_km)
        if n_f >= adopted or int(deg_r.max()) >= adopted:
            return {"ok": False, "design": design, "index": index, "beta": beta,
                    "reason": "policy reaches the reference degree; censored"}
        defect = J.force_defects(
            t[:n], R, {"radial": deg_r,
                       "constant": np.full(n, n_f, dtype=int)}, adopted, args)

        self_ref = J.self_difference(Y[("reference", "tight")][:, :n],
                                     Y[("reference", "tighter")][:, :n])
        err, env = {}, {}
        for policy in ("constant", "radial"):
            err[policy] = J.trajectory_error(
                Y[(policy, "tighter")][:, :n],
                Y[("reference", "tighter")][:, :n])["J_traj_rms_m"]
            env[policy] = self_ref + J.self_difference(
                Y[(policy, "tight")][:, :n], Y[(policy, "tighter")][:, :n])
        rho_f = (defect["radial"]["J_force_mean_m_s2"]
                 / defect["constant"]["J_force_mean_m_s2"])
        rho_x = err["radial"] / err["constant"]
        return {"ok": True, "design": design, "index": index, "beta": beta,
                "hp_km": g["hp_km"], "rho_force": float(rho_f),
                "rho_traj": float(rho_x),
                "resolved": bool(J.resolved(err["constant"], err["radial"],
                                            env["constant"], env["radial"])),
                "reversal": bool(rho_f < 1.0 < rho_x)}
    except Exception:
        return {"ok": False, "design": design, "index": index, "beta": beta,
                "reason": traceback.format_exc()}


def summarize(rows: list[dict]) -> dict:
    res = [r for r in rows if r["resolved"]]
    if not res:
        return {"orbits": len(rows), "resolved": 0}
    return {
        "orbits": len(rows), "resolved": len(res),
        "radial_wins_force": sum(1 for r in res if r["rho_force"] < 1.0),
        "radial_loses_trajectory": sum(1 for r in res if r["rho_traj"] > 1.0),
        "reversal": sum(1 for r in res if r["reversal"]),
        "median_rho_force": float(np.median([r["rho_force"] for r in res])),
        "median_rho_traj": float(np.median([r["rho_traj"] for r in res])),
    }


def write_table(summary: dict) -> None:
    lines = [r"\begin{tabular}{lrrrrrrr}", r"\hline",
             r"& & & \multicolumn{2}{c}{force metric} & "
             r"\multicolumn{2}{c}{trajectory metric} & \\",
             r"$\beta$ & $n$ & resolved & radial wins & $\tilde{\rho}_F$ & "
             r"radial loses & $\tilde{\rho}_X$ & reversal \\", r"\hline"]
    for beta in BUDGETS:
        s = summary.get(f"beta_{beta:.2f}")
        if not s or not s.get("resolved"):
            continue
        lines.append(
            f"{beta:.2f} & {s['orbits']} & {s['resolved']} & "
            f"{s['radial_wins_force']}/{s['resolved']} & "
            f"{s['median_rho_force']:.2f} & "
            f"{s['radial_loses_trajectory']}/{s['resolved']} & "
            f"{s['median_rho_traj']:.1f} & {s['reversal']} \\\\")
    lines += [r"\hline", r"\end{tabular}"]
    (METRICS / "rJ_table_budget_ladder_primary.tex").write_text(
        "\n".join(["% generated by revJ_ladder_primary.py -- do not edit",
                   *lines]) + "\n", encoding="utf-8")
    print("[written] rJ_table_budget_ladder_primary.tex")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()

    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    tasks = []
    for design in ("A", "B"):
        src = json.loads(ARCHIVED[design]["rows"].read_text(encoding="utf-8"))
        crit = {int(r["sobol_index"]): r for r in src["rows"]}
        par = {int(r["sobol_index"]): r
               for r in pareto["designs"][design]["rows"]}
        for index, row in crit.items():
            g = row.get("design_point", row)
            geom = {"hp_km": float(g["hp_km"]), "ha_km": float(g["ha_km"]),
                    "adopted_truth_degree": int(row["adopted_truth_degree"])}
            for beta in BUDGETS:
                budget = par[index]["budgets"].get(f"beta_{beta:.2f}")
                if budget is None or budget.get("censored"):
                    continue
                tasks.append({"design": design, "index": index, "beta": beta,
                              "geom": geom, "budget": budget})
    log(f"{len(tasks)} archived orbit-budget pairs to re-score")

    results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fut in as_completed([pool.submit(task, t) for t in tasks]):
            results.append(fut.result())
    good = [r for r in results if r.get("ok")]
    skipped = [r for r in results if not r.get("ok")]

    summary = {f"beta_{b:.2f}": summarize([r for r in good if r["beta"] == b])
               for b in BUDGETS}
    by_design = {d: {f"beta_{b:.2f}": summarize(
        [r for r in good if r["beta"] == b and r["design"] == d])
        for b in BUDGETS} for d in ("A", "B")}

    J.atomic_json(OUT, {
        "schema": "rJ_ladder_primary_v1", "created_utc": J.utc_now(),
        "field": J.field_key(),
        "purpose": ("the same budget ladder the cross-solution campaign runs, "
                    "on the primary solution, from archived trajectories and "
                    "with no new integration"),
        "definitions": {
            "rho_force": "J_force(radial)/J_force(constant), time-averaged "
                         "truncation defect along the archived reference",
            "rho_traj": "J_traj(radial)/J_traj(constant), position RMS",
            "resolution": "the two-level envelope rule, unchanged"},
        "budgets": list(BUDGETS),
        "counts": {"scored": len(good), "skipped": len(skipped)},
        "summary_pooled": summary, "summary_by_design": by_design,
        "complete": bool(good),
        "rows": sorted(good, key=lambda r: (r["beta"], r["design"], r["index"])),
        # every skip, with its reason -- an earlier revision stored only the
        # first 50, which left 56 skips unaccounted for; the full derivation of
        # that run's skip list is archived in rJ_ladder_skip_audit.json
        "skipped": skipped,
        "provenance": J.provenance()})
    write_table(summary)
    for beta in BUDGETS:
        s = summary[f"beta_{beta:.2f}"]
        if s.get("resolved"):
            log(f"beta={beta:.2f}: force {s['radial_wins_force']}/"
                f"{s['resolved']} rho_F={s['median_rho_force']:.3f} | "
                f"traj {s['radial_loses_trajectory']}/{s['resolved']} "
                f"rho_X={s['median_rho_traj']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
