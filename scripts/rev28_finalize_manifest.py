"""SHA-256 integrity manifest for the R28 post-hoc midpoint (O34).

R28 holds the one budget in the campaign that is not pre-registered: the
midpoint beta = 0.62, calibrated and propagated on both designs after the
bracketing budgets of (O33) had been read, to test whether the crossover
interval could be halved. It could not. The result is reported in the
supplement only, and this manifest is separate from R25 for the same reason:
an index that mixed the one unregistered budget in with the registered ones
would make that harder to see rather than easier.

Naming exception, same convention as R23 and R25. Every stage reuses an
archived driver called with a budget argument it already accepts, so the
outputs keep their r14_, r18_ and r19_ prefixes and carry a beta_0.62 suffix.
Those files belong to R28 and are indexed here only. R19 owns the records with
no beta suffix, R23 owns beta_0.50 and R25 owns beta_0.75, beta_1.25 and
beta_1.50, so the partition holds without any of those manifests moving.

What this campaign does not touch. The frozen Phase-A calibration record
r14_budget_pareto.json is pinned by the R14, R18 and R21 manifests and was not
edited; the beta = 0.62 calibration lives in its own record, which is indexed
here. The R19 equal-work table and its generator are pinned by R19 and R25, so
the midpoint table is generated into its own file by rev28_tables.py rather
than by regenerating theirs. Those five pinned files are listed under reused
inputs, where a stale digest would be caught by the integrity check.

Usage:  python rev28_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r28_final_experiment_manifest.json"

BETA_TAG = "beta_0.62"

SCRIPTS = ["rev28_calibration_amendment.py",
           "rev28_budget_pareto_extension.py",
           "rev28_r14_beta062.py", "rev28_queue.py", "rev28_tables.py",
           "rev28_finalize_manifest.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]

REGISTRATION = ["r28_calibration_amendment.json"]

RESULT_JSON = [
    "r28_budget_pareto_beta_0.62.json",
    f"r14_trajectory_A_{BETA_TAG}.json",
    f"r18_span_sweep_A_{BETA_TAG}.json",
    f"r19_equal_total_work_A_{BETA_TAG}.json",
    f"r14_trajectory_B_{BETA_TAG}.json",
    f"r18_span_sweep_B_{BETA_TAG}.json",
    f"r19_equal_total_work_B_{BETA_TAG}.json",
]

TABLES = ["r28_midpoint_table.tex", "r28_equal_work_table.tex",
          "r28_manuscript_descriptives.json"]

# Read but not produced here. The first is the frozen calibration this campaign
# extends without editing; the next two are the pinned table and generator the
# midpoint table exists in order not to overwrite; the last two are the
# registrations this amendment executes and amends.
REUSED = ["r14_budget_pareto.json",
          "r19_equal_work_table.tex",
          "r14_preregistration.json",
          "r25_preregistration_amendment.json"]

PRODUCED = ("A", "B")
K_INTERIOR = ("0.25", "0.50", "0.75")


def _tree_for(design: str) -> list[str]:
    rels = [f"r14_cases/{design}_{BETA_TAG}", f"r14_raw/{design}_{BETA_TAG}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{design}_{BETA_TAG}_k_{k}",
                 f"r18_raw/{design}_{BETA_TAG}_k_{k}"]
    rels += [f"r19_cases/{design}_workmatched_{BETA_TAG}",
             f"r19_raw/{design}_workmatched_{BETA_TAG}"]
    return rels


TREES = {f"O34_design_{d}_beta_0.62": _tree_for(d) for d in PRODUCED}


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


def index_trees() -> dict:
    out = {}
    for name, rels in TREES.items():
        sidecars, roll, n_raw = {}, hashlib.sha256(), 0
        for rel in rels:
            base = METRICS / rel
            if not base.exists():
                continue
            for p in sorted(base.rglob("*.json")):
                sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
            for p in sorted(base.rglob("*.npz")):
                roll.update(sha(p).encode())
                n_raw += 1
        out[name] = {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
                     "sidecar_sha256": sidecars,
                     "raw_rollup_sha256": roll.hexdigest()}
    return out


def completeness() -> dict:
    """Read from the records, so a truncated stage cannot be quoted as full."""
    out = {}
    for design in PRODUCED:
        entry = {}
        p = METRICS / f"r14_trajectory_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["budget_campaign"] = {
                "complete": d.get("complete"),
                "stopped_for_deadline": d.get("stopped_for_deadline")}
        p = METRICS / f"r18_span_sweep_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["span_sweep_orbits"] = len(d.get("rows", []))
        p = METRICS / f"r19_equal_total_work_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["equal_work"] = d.get("summary")
        out[design] = entry or {"not_run": True}
    return out


def main() -> int:
    scripts, a1 = index_files(SCRIPTS, CODE, required=True)
    registration, a2 = index_files(REGISTRATION, METRICS, required=True)
    results, a3 = index_files(RESULT_JSON, METRICS, required=True)
    tables, a4 = index_files(TABLES, METRICS, required=True)
    reused, a5 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3 + a4 + a5

    cal = json.loads(
        (METRICS / "r28_budget_pareto_beta_0.62.json").read_text(
            encoding="utf-8"))

    payload = {
        "schema": "r28_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R28 (O34): the post-hoc midpoint beta = 0.62. One Phase-A "
                  "calibration point added to the frozen grid after the "
                  "bracketing budgets had been read, and the trajectory, span "
                  "and realized-work stages that follow from it on both "
                  "designs."),
        "registration_status": (
            "NOT pre-registered. This is the one budget in the campaign that "
            "extends a frozen calibration grid after results were seen. The "
            "departure, its constraints and its reporting commitment are "
            "recorded in r28_calibration_amendment.json and in (O34) of the "
            "experiment contract. The manuscript's crossover bracket is the "
            "pre-registered one and does not rest on this budget."),
        "outcome": (
            "the interval was not narrowed: the two designs answer in "
            "opposite directions, by one resolved comparison of 39 on design A "
            "and three of 45 on design B"),
        "naming_exception": (
            "every stage reuses an archived driver with a budget argument it "
            "already accepts, so outputs keep the r14_, r18_ and r19_ prefixes "
            "and carry a beta_0.62 suffix. Those files belong to R28 and are "
            "indexed here only. R19 owns the records with no beta suffix, R23 "
            "owns beta_0.50 and R25 owns beta_0.75, beta_1.25 and beta_1.50."),
        "archive_untouched": {
            "statement": (
                "the frozen Phase-A calibration record and the pinned R19 "
                "equal-work table were read and not written. Their digests are "
                "carried under reused_inputs, so a modification would fail the "
                "integrity check here as well as in the manifests that pin "
                "them."),
            "r14_budget_pareto_sha256": cal["parent_record_sha256"],
            "pinned_by": ["r14_final_experiment_manifest.json",
                          "r18_final_experiment_manifest.json",
                          "r21_final_experiment_manifest.json"],
        },
        "admissibility_check": cal["admissibility_check"],
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "result_json": results,
        "generated_tables": tables,
        "panel_completeness": completeness(),
        "trajectory_tree": index_trees(),
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")
                   ).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in payload["trajectory_tree"].values())
    n_raw = sum(t["n_raw_arrays"] for t in payload["trajectory_tree"].values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}  "
          f"{n_side} sidecars, {n_raw} raw arrays")
    for design, info in payload["panel_completeness"].items():
        print(f"  design {design}: {json.dumps(info)[:110]}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
