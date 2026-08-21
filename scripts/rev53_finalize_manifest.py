"""SHA-256 integrity manifest for R53: the post-hoc budget column.

R53 is unlike the campaigns before it in one way that decides how this manifest
is built. It is not a population. It is a budget, run on seven populations that
already exist, so its records land inside four other registries' directory
trees and under four other campaigns' key names. Ownership therefore has to be
stated rather than inferred from a prefix, and this manifest claims exactly what
this campaign wrote: the three ladder records of every completed cell at
beta = 0.62 and the case and raw subtrees underneath them, nothing else.

The consequence is worth writing down where it will be found. The finalizers of
the parent campaigns -- rev29, rev30, rev31 and rev38 -- discover their budgets
by globbing the disk, so re-running any of them now would claim these same
records a second time and break the partition the integrity check enforces.
That is the intended failure mode: loud, at the gate, rather than a silent
double count. If a parent ever has to be re-sealed, its budget list must first
be restricted to the budgets its own registration declared.

What this campaign did not write it indexes as reused: the seven calibration
records, each of which already carried 0.62 before the first propagation, and
each population's own registration and frozen design. An edit to any of them
fails this gate as well as its own.

A cell counts only with all three ladder stages on disk, so a column stopped by
its clock records what it reached rather than what it intended.

Usage:  python rev53_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"

REG = "r53"
TAG = "beta_0.62"
K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev53_campaign.py", "rev53_verdict.py",
           "rev53_finalize_manifest.py", "rev53_resume_stages.py",
           "rev53_relaunch_watch.py", "rev53_write_departure.py",
           "rev54_night_chain.py",
           "rev30_stratum_ops.py", "rev29_designC_ladder.py",
           "population_registry.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]

# Several of these scripts postdate the first five cells, and saying which is
# cheaper than letting a reader infer a history that did not happen. None of
# them changes how a record is computed; what they change is when a stage is
# abandoned and where a stage starts from.
SCRIPTS_NOTE = (
    "The first five cells were written by rev53_campaign.py driving "
    "rev30_stratum_ops.py and rev29_designC_ladder.py at eleven workers. The "
    "last two were written on 15 August by rev54_night_chain.py driving the "
    "same supervisor at four workers, after abrupt pool-worker deaths at "
    "eleven, eight and six were diagnosed as contention on the host. Two "
    "amendments to the supervisor date from that run and neither touches a "
    "computation: a stall guard that ends a stage which has stopped producing "
    "output rather than waiting on it, and a branch that hands a cell whose "
    "trajectory record is already complete to rev53_resume_stages.py, so that "
    "the equatorial cell's archived 64-row R14 record was resumed from rather "
    "than recomputed and overwritten. rev53_resume_stages.py therefore did "
    "produce records indexed here, the equatorial cell's span and "
    "work-matched stages; rev53_relaunch_watch.py produced none. "
    "rev53_write_departure.py writes the departure record from the "
    "supervisor's own progress log and the sealed verdicts, and reads no "
    "trajectory.")

# read, never written: the calibrations that already carried this budget, and
# the registration and frozen design of every population the column runs on
REUSED = ["r28_calibration_amendment.json",
          "r29_preregistration.json", "r29_budget_pareto_designC.json",
          "r30_preregistration.json",
          "r30_budget_pareto_high_apolune.json",
          "r30_budget_pareto_polar.json",
          "r30_budget_pareto_equatorial.json",
          "r30_budget_pareto_frozen_like.json",
          "r31_preregistration.json",
          "r31_budget_pareto_operational_elliptical.json",
          "r38_preregistration.json",
          "r38_budget_pareto_operational_elliptical_uncapped.json",
          "r14_budget_pareto.json", "r14_preregistration.json"]

SCOPE = ("R53: the declared post-hoc budget beta = 0.62, already computed on "
         "designs A and B and on the two identity blocks of the paired ladder, "
         "extended to the seven populations of the regime map that lacked it. "
         "No new population, orbit, parameter or reference degree; a column "
         "added to a grid at a budget every one of those populations' "
         "calibrations already carried.")


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path, required: bool):
    out, absent = {}, []
    for n in names:
        p = base / n
        if p.exists():
            out[n] = {"sha256": sha(p), "bytes": p.stat().st_size}
        else:
            out[n] = {"absent": True}
            if required:
                absent.append(n)
    return out, absent


def index_tree(rels) -> dict:
    sidecars, roll, n_raw = {}, hashlib.sha256(), 0
    for rel in rels:
        d = METRICS / rel
        if not d.exists():
            continue
        for p in sorted(d.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
        for p in sorted(d.rglob("*.npz")):
            roll.update(sha(p).encode())
            n_raw += 1
    return {"n_sidecars": len(sidecars), "sidecars": sidecars,
            "n_raw": n_raw, "raw_rollup_sha256": roll.hexdigest()}


def ladder_trees(key: str) -> list[str]:
    rels = [f"r14_cases/{key}_{TAG}", f"r14_raw/{key}_{TAG}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{key}_{TAG}_k_{k}", f"r18_raw/{key}_{TAG}_k_{k}"]
    rels += [f"r19_cases/{key}_workmatched_{TAG}",
             f"r19_raw/{key}_workmatched_{TAG}"]
    return rels


def cell_records(key: str) -> list[str]:
    return [f"r14_trajectory_{key}_{TAG}.json",
            f"r18_span_sweep_{key}_{TAG}.json",
            f"r19_equal_total_work_{key}_{TAG}.json"]


def main() -> int:
    prereg_name = f"{REG}_preregistration.json"
    prereg_p = METRICS / prereg_name
    if not prereg_p.exists():
        print(f"[abort] {prereg_name} missing; nothing to manifest")
        return 2
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))

    verdict_name = f"{REG}_verdict.json"
    if not (METRICS / verdict_name).exists():
        print(f"[abort] {verdict_name} missing; run rev53_verdict.py first")
        return 2
    verdict = json.loads((METRICS / verdict_name).read_text(encoding="utf-8"))

    registration, a1 = index_files([prereg_name], METRICS, required=True)
    scripts, a2 = index_files(SCRIPTS, CODE, required=True)
    reused, a3 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3

    result_names = [verdict_name, f"{REG}_registration_departure.json"]
    # The verdict that was sealed at five cells, before the window was extended
    # and the last two cells ran. It is kept and owned rather than overwritten,
    # because the departure record points at it and a reader checking whether
    # the extension changed the outcome needs the version it changed from.
    superseded = f"{REG}_verdict.sealed_20260814.json"
    if (METRICS / superseded).exists():
        result_names.append(superseded)
    trees, cells, not_run = {}, {}, []
    partial_names, partial_trees = [], {}
    for spec in sorted(prereg["cells"], key=lambda c: c["order"]):
        key = spec["design_key"]
        names = cell_records(key)
        if not all((METRICS / n).exists() for n in names):
            # A cell can be incomplete and still have left something real on
            # disk. The equatorial cell finished its trajectory stage -- 64
            # rows, no failures -- and then its span sweep killed a pool worker
            # abruptly on every retry, at two, six, eight and eleven workers,
            # attached and detached. That record is not deleted, because it was
            # produced and it is evidence of how far the cell got, and it is not
            # promoted into results, because no verdict is read from a cell
            # missing two stages. It is owned here and labelled for what it is.
            present = [n for n in names if (METRICS / n).exists()]
            entry = {"design_key": key, "population": spec["population"],
                     "status": ("declared and not run" if not present
                                else "started, stages on disk, cell incomplete")}
            if present:
                entry["stage_records"] = present
                partial_names.extend(present)
                partial_trees[key] = index_tree(ladder_trees(key))
            not_run.append(entry)
            continue
        result_names += names
        trees[key] = index_tree(ladder_trees(key))
        r14 = json.loads((METRICS / names[0]).read_text(encoding="utf-8"))
        r19 = json.loads((METRICS / names[2]).read_text(encoding="utf-8"))
        cells[key] = {
            "population": spec["population"],
            "parent_registry": spec["registry"],
            "radial_endpoint_vs_constant": r14.get("summary"),
            "interior_vs_work_matched_constant": r19.get("summary"),
        }

    results, a4 = index_files(sorted(set(result_names)), METRICS, required=True)
    partial, _ = index_files(sorted(set(partial_names)), METRICS,
                             required=False)
    absent += a4

    payload = {
        "schema": f"{REG}_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "declared before the first propagation in "
            f"{prereg_name}, which fixed the budget, the cell list, the cell "
            "order, four outcomes and what the campaign is prohibited from "
            "doing. The budget itself is post hoc and is declared as such "
            "wherever it is drawn; it was fixed for other designs in "
            "r28_calibration_amendment.json and is inherited here, not "
            "reopened."),
        "ownership_note": (
            "this manifest claims the beta = 0.62 ladder records and subtrees "
            "of the cells listed under cells, and nothing else. Those files sit "
            "in the r14, r18 and r19 trees under key names belonging to four "
            "other registries, whose finalizers discover budgets by globbing; "
            "re-running one of them without restricting it to its own "
            "registered budgets would claim these records twice and fail the "
            "partition check."),
        "partition_note": (
            "every population's base, operating point and calibration remain "
            "owned by that population's own manifest and are indexed here as "
            "reused inputs. This campaign propagated trajectories and wrote "
            "nothing else."),
        "outcome": verdict["outcome"],
        "outcome_text": verdict["outcome_text"],
        "cells": cells,
        "cells_declared_and_not_run": not_run,
        "partial_cell_records": partial,
        "partial_cell_trees": partial_trees,
        "partial_cell_note": (
            "stages that completed inside a cell that did not. They are owned here so that nothing this campaign wrote sits under no manifest, and they are kept out of results and out of the verdict, which read only cells carrying all three stages."),
        "n_cells_declared": len(prereg["cells"]),
        "n_cells_complete": len(cells),
        "registration": registration,
        "scripts": scripts,
        "scripts_note": SCRIPTS_NOTE,
        "reused_inputs": reused,
        "results": results,
        "trees": trees,
        "absent_required": absent,
    }
    out = METRICS / f"{REG}_final_experiment_manifest.json"
    # the self-seal is the digest of the payload without the seal field, in the
    # canonical form check_manifest_integrity.seal_ok recomputes: sorted keys
    # and no whitespace. Hashing the indented form instead is a seal only this
    # file can verify, which is not a seal.
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    print(f"[{REG}] {len(cells)} of {len(prereg['cells'])} cells complete "
          f"({', '.join(cells) or 'none'})")
    if not_run:
        print(f"[{REG}] declared and not run: "
              f"{', '.join(c['design_key'] for c in not_run)}")
    if absent:
        print(f"[FAIL] {len(absent)} required files absent: "
              f"{', '.join(absent[:8])}")
        return 1
    print(f"[written] {out.name}: {len(results)} results, "
          f"{sum(t['n_sidecars'] for t in trees.values())} sidecars, "
          f"{sum(t['n_raw'] for t in trees.values())} raw arrays")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
