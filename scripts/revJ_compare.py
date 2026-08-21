"""Primary field against cross-solution field, under one set of definitions.

J1 reports its own counts, and the manuscript reports its own, but the two are
not built from the same statistic: the archived budget records summarize the
force defect as an RMS and quote the ratio constant-over-radial, while the
campaigns quote a time-averaged defect as radial-over-constant. Those are two
sensible conventions and one of them inverted would silently turn a
confirmation into a contradiction.

So the comparison is not made by putting two published tables side by side. The
archived populations are re-scored here with the *campaign's* definitions,
from the archived trajectories, and only then compared:

    rho_force = J_force(radial) / J_force(constant),   J_force = time average
    rho_traj  = J_traj(radial)  / J_traj(constant),    J_traj  = position RMS

A reversal is rho_force < 1 < rho_traj: the radial allocation buys a smaller
force error and pays for it with a larger trajectory error at the same declared
budget. The question this file answers is whether the primary solution and the
independently produced one give the same answer under one convention, applied
to both.

Usage:
    python revJ_compare.py --workers 11
"""

from __future__ import annotations

import argparse
import json
import os
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
# A tag lets the comparison be rebuilt at the end of a campaign night without
# overwriting the record the night started from.
TAG = os.environ.get("JCAMP_COMPARE_TAG", "")
OUT = METRICS / f"rJ_field_comparison{TAG}.json"
TABLE = METRICS / f"rJ_field_comparison{TAG}_table.tex"
LOG = Path(__file__).resolve().parent / "rJ_compare.log"

PARETO = METRICS / "r14_budget_pareto.json"
BETA = 1.00
ARCHIVED = {
    "A": {"rows": METRICS / "r10_sobolA_baseline_truth_corrected.json",
          "conv": METRICS / "r11_raw" / "convergence",
          "budget": METRICS / "r14_raw" / f"A_beta_{BETA:.2f}"},
    "B": {"rows": METRICS / "r11_designB_rows.json",
          "conv": METRICS / "r11_raw" / "designB_convergence",
          "budget": METRICS / "r14_raw" / f"B_beta_{BETA:.2f}"},
}


def log(msg: str) -> None:
    J.log_line(LOG, f"compare {msg}")


def _raw(design: str, index: int, policy: str, level: str) -> Path | None:
    src = ARCHIVED[design]
    name = f"sobolA_{index:03d}"
    if policy == "reference":
        p = src["conv"] / name / f"truth_{level}.npz"
    elif policy == "radial":
        p = src["budget"] / name / f"atallah_budget_{level}.npz"
    else:
        p = src["budget"] / name / f"fixed_budget_{level}.npz"
        if not p.exists():
            p = src["conv"] / name / f"fixed_critical_{level}.npz"
    return p if p.exists() else None


def archived_task(payload: dict) -> dict:
    design, index = payload["design"], int(payload["index"])
    try:
        g, budget = payload["geom"], payload["budget"]
        adopted = int(g["adopted_truth_degree"])
        paths = {(p, lv): _raw(design, index, p, lv)
                 for p in ("reference", "constant", "radial")
                 for lv in ("tight", "tighter")}
        if any(v is None for v in paths.values()):
            return {"ok": False, "design": design, "index": index,
                    "error": "incomplete archived two-level record"}
        Y = {k: J.load_states(v) for k, v in paths.items()}
        t = J.load_times(paths[("reference", "tighter")])
        n = min(v.shape[1] for v in Y.values())

        model, args = J.model_for(adopted)
        h_km = (np.linalg.norm(Y[("reference", "tighter")][:3, :n], axis=0)
                - model.r_ref) / 1e3
        _, table = at.atallah_binned_schedule(
            model, J.atallah_g(adopted),
            float(budget["atallah"]["tol_accel_m_s2"]),
            g["hp_km"], g["ha_km"], floor=J.FLOOR, cap=adopted,
            bin_km=J.BIN_KM)
        table = {float(k): int(v) for k, v in table.items()}
        defect = J.force_defects(
            t[:n], Y[("reference", "tighter")][:3, :n],
            {"radial": J.degrees_from_table(table, h_km),
             "constant": np.full(n, int(budget["fixed"]["degree"]), dtype=int)},
            adopted, args)

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
        return {"ok": True, "design": design, "index": index,
                "hp_km": g["hp_km"], "adopted_truth_degree": adopted,
                "rho_force": float(rho_f), "rho_traj": float(rho_x),
                "resolved": bool(J.resolved(err["constant"], err["radial"],
                                            env["constant"], env["radial"])),
                "reversal": bool(rho_f < 1.0 < rho_x)}
    except Exception:
        return {"ok": False, "design": design, "index": index,
                "error": traceback.format_exc()}


def summarize(name: str, field: str, rows: list[dict]) -> dict:
    res = [r for r in rows if r["resolved"]]
    return {
        "population": name, "field": field,
        "orbits": len(rows), "resolved": len(res),
        "radial_wins_force": sum(1 for r in res if r["rho_force"] < 1.0),
        "radial_loses_trajectory": sum(1 for r in res if r["rho_traj"] > 1.0),
        "reversal": sum(1 for r in res if r["reversal"]),
        "median_rho_force": float(np.median([r["rho_force"] for r in rows])),
        "median_rho_traj": float(np.median([r["rho_traj"] for r in rows])),
        "median_rho_force_resolved":
            float(np.median([r["rho_force"] for r in res])) if res else None,
        "median_rho_traj_resolved":
            float(np.median([r["rho_traj"] for r in res])) if res else None,
    }


def write_table(summaries: list[dict]) -> None:
    lines = [
        r"% generated by revJ_compare.py -- do not edit by hand",
        r"\begin{tabular}{llrrrrrr}",
        r"\hline",
        r"Population & Gravity solution & $n$ & resolved & force & traj. & "
        r"reversal & $\tilde{\rho}_F$ / $\tilde{\rho}_X$ \\",
        r"\hline",
    ]
    for s in summaries:
        lines.append(
            f"{s['population']} & {s['field'].replace('_', r'\_')} & "
            f"{s['orbits']} & {s['resolved']} & {s['radial_wins_force']} & "
            f"{s['radial_loses_trajectory']} & {s['reversal']} & "
            f"{s['median_rho_force_resolved']:.2f} / "
            f"{s['median_rho_traj_resolved']:.1f} \\\\")
    lines += [r"\hline", r"\end{tabular}"]
    TABLE.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    ap.add_argument("--stale-check", action="store_true",
                    help="re-score only when some campaign record is newer "
                         "than the comparison, so the stage can sit at the end "
                         "of a repeating queue without re-running for hours")
    ap.add_argument("--table-only", action="store_true",
                    help="rebuild the table from the saved record, so a "
                         "formatting fix costs no recomputation")
    a = ap.parse_args()

    if a.stale_check and OUT.is_file() and OUT.stat().st_size:
        newest = max((p.stat().st_mtime for p in METRICS.glob("rJ*_score*.json")),
                     default=0.0)
        if OUT.stat().st_mtime >= newest:
            log(f"{OUT.name} is newer than every scored record; nothing to do")
            write_table(json.loads(OUT.read_text(encoding="utf-8"))["summaries"])
            return 0
        log("a scored record is newer than the comparison; re-scoring")

    if a.table_only:
        saved = json.loads(OUT.read_text(encoding="utf-8"))
        write_table(saved["summaries"])
        log(f"table rebuilt from {OUT.name}")
        return 0

    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    tasks, censored = [], 0
    for design in ("A", "B"):
        src = json.loads(ARCHIVED[design]["rows"].read_text(encoding="utf-8"))
        geom = {}
        for row in src["rows"]:
            g = row.get("design_point", row)
            geom[int(row["sobol_index"])] = {
                "hp_km": float(g["hp_km"]), "ha_km": float(g["ha_km"]),
                "adopted_truth_degree": int(row["adopted_truth_degree"])}
        for row in pareto["designs"][design]["rows"]:
            index = int(row["sobol_index"])
            budget = row["budgets"].get(f"beta_{BETA:.2f}")
            if budget is None or budget.get("censored") or index not in geom:
                censored += 1
                continue
            tasks.append({"design": design, "index": index,
                          "geom": geom[index], "budget": budget})
    log(f"{len(tasks)} archived orbits to re-score under campaign definitions")

    results = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for fut in as_completed([pool.submit(archived_task, t) for t in tasks]):
            results.append(fut.result())
    good = [r for r in results if r.get("ok")]
    dropped = [r for r in results if not r.get("ok")]

    summaries = []
    for design in ("A", "B"):
        rows = [r for r in good if r["design"] == design]
        if rows:
            summaries.append(summarize(f"design {design}", "JGGRX_1800F",
                                       rows))
    if (METRICS / "rJ1_score.json").exists():
        j1 = json.loads((METRICS / "rJ1_score.json").read_text(
            encoding="utf-8"))
        j1_rows = [{"rho_force": r["rho_force"], "rho_traj": r["rho_traj"],
                    "resolved": r["resolved"]["tighter"],
                    "reversal": r["reversal"]} for r in j1["rows"]]
        summaries.append(summarize("J1 cross-solution", "GRGM1200A", j1_rows))

    payload = {
        "schema": "rJ_field_comparison_v1", "created_utc": J.utc_now(),
        "beta": BETA,
        "definitions": {
            "rho_force": "J_force(radial)/J_force(constant), J_force being the "
                         "time average of the truncation defect magnitude "
                         "along the reference",
            "rho_traj": "J_traj(radial)/J_traj(constant), J_traj being the "
                        "position RMS against the reference",
            "reversal": "rho_force < 1 < rho_traj",
            "note": "the archived budget records quote the reciprocal of "
                    "rho_force as an RMS; the archived populations are "
                    "re-scored here rather than transcribed, so both fields "
                    "are read under one convention",
        },
        "counts": {"archived_rescored": len(good), "dropped": len(dropped),
                   "censored_or_missing": censored},
        "summaries": summaries,
        "rows": sorted(good, key=lambda r: (r["design"], r["index"])),
        "dropped_detail": dropped,
        "complete": bool(good and not dropped),
        "provenance": J.provenance(),
    }
    J.atomic_json(OUT, payload)
    write_table(summaries)
    for s in summaries:
        log(f"{s['population']:20s} {s['field']:14s} reversal "
            f"{s['reversal']}/{s['resolved']} resolved of {s['orbits']}, "
            f"median rho_F {s['median_rho_force_resolved']:.3f}, "
            f"rho_X {s['median_rho_traj_resolved']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
