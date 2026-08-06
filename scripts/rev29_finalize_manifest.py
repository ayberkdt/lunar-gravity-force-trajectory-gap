"""SHA-256 integrity manifest for R29: the third coverage design.

R29 propagates design C, the third scrambled-Sobol draw R26 froze, and runs the
archived budget ladder on it. Everything it produces is new: a base tree no
manifest has seen, an accuracy-target operating point regenerated under the R12
rule because design C has no R12 campaign, a Phase-A calibration in its own
record because the archived one is pinned in three sealed manifests, and one
ladder per propagated budget.

Naming exception, same convention as R23, R25 and R28: the ladder reuses the
archived drivers, so its outputs keep the r14_, r18_ and r19_ prefixes and carry
the design key C together with a budget suffix. Those files belong to R29 and
are indexed here only. No existing manifest claims them: R14 and R18 own
explicit (design, budget) allow-lists that contain designs A and B alone, and
R19 owns the subtrees with no budget suffix -- which is why design C writes even
beta = 1 with an explicit suffix.

The budget list is read off the disk rather than declared, because the campaign
runs to a wall clock and a manifest that named a budget the clock did not reach
would record an intention rather than an archive. Each budget indexed here is
one whose three stages are all present.

Usage:  python rev29_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r29_final_experiment_manifest.json"

DESIGN = "C"
K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev26_designC_freeze.py", "rev26_designC_base.py",
           "rev29_preregister.py", "rev29_designC_operating_point.py",
           "rev29_designC_pareto.py", "rev29_designC_ladder.py",
           "rev29_campaign.py", "rev29_verdict.py", "rev29_tables.py",
           "rev29_finalize_manifest.py",
           "rev11_full_convergence.py", "rev11_designB_convergence.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]

REGISTRATION = ["r26_preregistration.json", "r26_sobolC_design_frozen.json",
                "r29_preregistration.json"]

BASE_JSON = ["r26_designC_rows.json", "r26_designC_convergence.json",
             "r29_designC_operating_point.json",
             "r29_budget_pareto_designC.json"]

TABLES = ["r29_designC_table.tex", "r29_three_design_table.tex",
          "r29_manuscript_descriptives.json", "r29_verdict.json"]

# Read, not written. The frozen calibration this campaign refuses to edit, the
# registration whose outcomes it reports against, and the amendment that made
# beta = 0.62 a declared post-hoc budget on every design including this one.
REUSED = ["r14_budget_pareto.json", "r14_preregistration.json",
          "r28_calibration_amendment.json"]

BASE_TREES = ["r11_cases/designC_convergence", "r11_raw/designC_convergence"]
OP_TREE = ["r29_cases/atallah_designC"]

BETA_RE = re.compile(r"^r19_equal_total_work_C_beta_(\d+\.\d+)\.json$")


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def complete_budgets() -> list[str]:
    """A budget counts as propagated only with all three stages on disk."""
    out = []
    for p in sorted(METRICS.glob("r19_equal_total_work_C_beta_*.json")):
        m = BETA_RE.match(p.name)
        if not m:
            continue
        tag = f"beta_{m.group(1)}"
        need = [METRICS / f"r14_trajectory_C_{tag}.json",
                METRICS / f"r18_span_sweep_C_{tag}.json"]
        if all(q.exists() for q in need):
            out.append(tag)
    return out


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
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def ladder_trees(tag: str) -> list[str]:
    rels = [f"r14_cases/{DESIGN}_{tag}", f"r14_raw/{DESIGN}_{tag}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{DESIGN}_{tag}_k_{k}",
                 f"r18_raw/{DESIGN}_{tag}_k_{k}"]
    rels += [f"r19_cases/{DESIGN}_workmatched_{tag}",
             f"r19_raw/{DESIGN}_workmatched_{tag}"]
    return rels


def completeness(tags) -> dict:
    out = {}
    p = METRICS / "r26_designC_convergence.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        out["base"] = {"complete": d.get("complete"),
                       "orbits": len(d.get("rows", [])),
                       "failures": len(d.get("failures", []))}
    for tag in tags:
        entry = {}
        q = METRICS / f"r14_trajectory_{DESIGN}_{tag}.json"
        if q.exists():
            d = json.loads(q.read_text(encoding="utf-8"))
            entry["budget_campaign"] = {
                "complete": d.get("complete"),
                "stopped_for_deadline": d.get("stopped_for_deadline")}
        q = METRICS / f"r18_span_sweep_{DESIGN}_{tag}.json"
        if q.exists():
            entry["span_sweep_orbits"] = len(
                json.loads(q.read_text(encoding="utf-8")).get("rows", []))
        q = METRICS / f"r19_equal_total_work_{DESIGN}_{tag}.json"
        if q.exists():
            entry["equal_work"] = json.loads(
                q.read_text(encoding="utf-8")).get("summary")
        out[tag] = entry
    return out


def main() -> int:
    tags = complete_budgets()
    results = [f"r14_trajectory_{DESIGN}_{t}.json" for t in tags]
    results += [f"r18_span_sweep_{DESIGN}_{t}.json" for t in tags]
    results += [f"r19_equal_total_work_{DESIGN}_{t}.json" for t in tags]

    scripts, a1 = index_files(SCRIPTS, CODE, required=True)
    registration, a2 = index_files(REGISTRATION, METRICS, required=True)
    base_json, a3 = index_files(BASE_JSON, METRICS, required=True)
    result_json, a4 = index_files(sorted(results), METRICS, required=True)
    tables, a5 = index_files(TABLES, METRICS, required=False)
    reused, a6 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3 + a4 + a6

    cal_path = METRICS / "r29_budget_pareto_designC.json"
    cal = (json.loads(cal_path.read_text(encoding="utf-8"))
           if cal_path.exists() else {})

    trees = {"base_design_C_convergence": index_tree(BASE_TREES),
             "accuracy_target_operating_point": index_tree(OP_TREE)}
    for tag in tags:
        trees[f"design_C_{tag}"] = index_tree(ladder_trees(tag))

    payload = {
        "schema": "r29_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R29: design C, the third scrambled-Sobol coverage design "
                  "frozen by R26. Its base, its regenerated accuracy-target "
                  "operating point, its Phase-A calibration, and the archived "
                  "budget ladder at every budget the campaign reached."),
        "registration_status": (
            "pre-registered. The design and its outcomes were frozen in "
            "r26_preregistration.json before any propagation; the steps R26 "
            "left unwritten -- the separate calibration record, the "
            "regenerated operating point, the ladder order and the verdict "
            "rule -- were fixed in r29_preregistration.json, also before any "
            "design-C budget number existed. beta = 0.62 remains a declared "
            "post-hoc budget here exactly as it is on designs A and B."),
        "naming_exception": (
            "the ladder reuses the archived drivers, so its outputs carry the "
            "r14_, r18_ and r19_ prefixes with the design key C and a budget "
            "suffix, and are indexed here only. R14 and R18 own explicit "
            "(design, budget) allow-lists over designs A and B; R19 owns the "
            "subtrees with no budget suffix, which is why design C writes "
            "beta = 1 with an explicit suffix rather than bare."),
        "operating_point_note": (
            "rev14_budget_pareto.worker reads an archived R12 accuracy-target "
            "configuration per orbit and design C has no R12 campaign. The "
            "point was regenerated under the R12 rule into r29_cases/, and the "
            "regeneration was checked against the archived R12 configurations "
            "of designs A and B; the check is recorded in "
            "r29_designC_operating_point.json."),
        "archive_untouched": {
            "statement": ("the frozen Phase-A calibration was read and not "
                          "written; design C is calibrated in a record of its "
                          "own. Its digest is carried under reused inputs, so "
                          "a modification would fail the integrity check here "
                          "as well as in the three manifests that pin it."),
            "r14_budget_pareto_sha256": cal.get("parent_record_sha256"),
            "pinned_by": ["r14_final_experiment_manifest.json",
                          "r18_final_experiment_manifest.json",
                          "r21_final_experiment_manifest.json"],
        },
        "admissibility_check": cal.get("admissibility_check"),
        "budgets_indexed": tags,
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "base_json": base_json,
        "result_json": result_json,
        "generated_tables": tables,
        "panel_completeness": completeness(tags),
        "trajectory_tree": trees,
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in trees.values())
    n_raw = sum(t["n_raw_arrays"] for t in trees.values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}  "
          f"{n_side} sidecars, {n_raw} raw arrays")
    print(f"  budgets indexed: {', '.join(tags) if tags else 'none yet'}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
