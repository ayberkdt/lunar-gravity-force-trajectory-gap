"""SHA-256 integrity manifest for R38: the cap-lifted control.

R38 is not a new population. It is the R31 operational elliptical population
propagated a second time with one number changed -- the adopted reference
degree, 300 to 600 -- so that the calibrated radial schedule is never clamped
to the reference model. Everything else, including the initial states and the
comparator degrees, is the parent's.

That makes the parent an input rather than a neighbour, and the manifest treats
it as one. The parent's frozen design, its prepass rows, its registration and
its Phase-A calibration are indexed here under reused inputs, so an edit to any
of them fails this gate as well as R31's. The derived-from block records the
same relation in the form the reader needs: which fields were copied, which one
was changed, and the digest each was copied from.

What is new and claimed here: the derived design and rows, the R38 registration,
the base propagated in r11_cases/stratum_operational_elliptical_uncapped_*, the
regenerated operating point in r38_cases/, the Phase-A calibration, the three
ladders under design key OEU, and the per-budget verdict records. No earlier
manifest claims any of it -- R14 and R18 own allow-lists over designs A and B,
R19 owns the subtrees carrying no budget suffix, and OEU appears in no
allow-list -- so this manifest and its parent's partition rather than overlap.

Budgets are read off the disk, and a budget counts only with all three ladder
stages present. A manifest that named a budget the campaign did not reach would
record an intention rather than an archive.

Usage:  python rev38_finalize_manifest.py
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

REG = "r38"
NAME = "operational_elliptical_uncapped"
KEY = "OEU"
PARENT_REG, PARENT_NAME, PARENT_KEY = "r31", "operational_elliptical", "OE"

K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev38_uncapped_freeze.py", "rev38_campaign.py",
           "rev38_verdict.py", "rev38_finalize_manifest.py",
           "rev30_stratum_base.py", "rev30_stratum_ops.py",
           "rev32_reindex_two_policy_base.py", "population_registry.py",
           "rev29_designC_operating_point.py",
           "rev11_full_convergence.py", "rev11_designB_convergence.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]

# Read, never written. The parent supplies the orbits, the prepass and the
# calibration this control is paired against; the R12-rule check validated the
# regenerated operating point once, for every population that regenerates one.
REUSED = [f"{PARENT_REG}_preregistration.json",
          f"{PARENT_REG}_{PARENT_NAME}_design_frozen.json",
          f"{PARENT_REG}_{PARENT_NAME}_rows.json",
          f"{PARENT_REG}_budget_pareto_{PARENT_NAME}.json",
          f"{PARENT_REG}_{PARENT_NAME}_convergence.json",
          "r14_budget_pareto.json", "r14_preregistration.json",
          "r29_designC_operating_point.json"]

SCOPE = ("R38: the R31 operational elliptical population re-propagated with "
         "its adopted reference degree raised from 300 to 600, so that the "
         "budget-calibrated radial schedule is never clamped to the reference "
         "model. A paired control on the one population where the paper's "
         "negative result reverses, run because the confound it removes was "
         "found in the degree-ceiling audit and could not be settled by "
         "subsetting the capped run.")


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
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def ladder_trees(tag: str) -> list[str]:
    rels = [f"r14_cases/{KEY}_{tag}", f"r14_raw/{KEY}_{tag}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{KEY}_{tag}_k_{k}", f"r18_raw/{KEY}_{tag}_k_{k}"]
    rels += [f"r19_cases/{KEY}_workmatched_{tag}",
             f"r19_raw/{KEY}_workmatched_{tag}"]
    return rels


def complete_budgets() -> list[str]:
    out = []
    rx = re.compile(rf"^r19_equal_total_work_{KEY}_beta_(\d+\.\d+)\.json$")
    for p in sorted(METRICS.glob(f"r19_equal_total_work_{KEY}_beta_*.json")):
        m = rx.match(p.name)
        if not m:
            continue
        tag = f"beta_{m.group(1)}"
        if all((METRICS / f"{stem}_{KEY}_{tag}.json").exists()
               for stem in ("r14_trajectory", "r18_span_sweep")):
            out.append(tag)
    return out


def main() -> int:
    OUT = METRICS / f"{REG}_final_experiment_manifest.json"
    prereg = json.loads((METRICS / f"{REG}_preregistration.json").read_text(
        encoding="utf-8"))
    design = json.loads(
        (METRICS / f"{REG}_{NAME}_design_frozen.json").read_text(
            encoding="utf-8"))
    rows = json.loads((METRICS / f"{REG}_{NAME}_rows.json").read_text(
        encoding="utf-8"))

    registration, a1 = index_files(
        [f"{REG}_preregistration.json", f"{REG}_{NAME}_design_frozen.json"],
        METRICS, required=True)
    scripts, a2 = index_files(SCRIPTS, CODE, required=True)
    reused, a3 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3

    conv_p = METRICS / f"{REG}_{NAME}_convergence.json"
    if not conv_p.exists():
        print(f"[abort] {conv_p.name} missing; nothing to manifest")
        return 2
    conv = json.loads(conv_p.read_text(encoding="utf-8"))

    result_names = [f"{REG}_{NAME}_rows.json", f"{REG}_{NAME}_convergence.json",
                    f"{REG}_{NAME}_operating_point.json",
                    f"{REG}_budget_pareto_{NAME}.json",
                    f"{REG}_campaign_progress.json",
                    f"{REG}_manuscript_descriptives.json"]
    tree_map = {
        "base": index_tree([f"r11_cases/stratum_{NAME}_convergence",
                            f"r11_raw/stratum_{NAME}_convergence"]),
        "operating_point": index_tree([f"{REG}_cases/atallah_{NAME}"]),
    }

    budgets = {}
    for tag in complete_budgets():
        result_names += [f"r14_trajectory_{KEY}_{tag}.json",
                         f"r18_span_sweep_{KEY}_{tag}.json",
                         f"r19_equal_total_work_{KEY}_{tag}.json"]
        v = METRICS / f"{REG}_verdict_{tag}.json"
        if v.exists():
            result_names.append(v.name)
        tree_map[tag] = index_tree(ladder_trees(tag))
        q = json.loads((METRICS / f"r19_equal_total_work_{KEY}_{tag}.json"
                        ).read_text(encoding="utf-8"))
        r = json.loads((METRICS / f"r14_trajectory_{KEY}_{tag}.json"
                        ).read_text(encoding="utf-8"))
        budgets[tag] = {
            "interior_vs_work_matched_constant": q.get("summary"),
            "radial_endpoint_vs_constant": r.get("summary"),
        }

    results, a4 = index_files(sorted(set(result_names)), METRICS, required=True)
    tables, _ = index_files([f"{REG}_verdict.json",
                             f"{REG}_verdict_rev30_generic.json"],
                            METRICS, required=False)
    absent += a4

    payload = {
        "schema": f"{REG}_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The reference degree, the budget order, the "
            "verdict rule and four outcomes were fixed in "
            "r38_preregistration.json before the first trajectory. The "
            "sufficiency of degree 600 was established before propagation by a "
            "demand probe on the archived altitude histories, and that probe is "
            "recorded in the registration. The order in which the four "
            "outcomes are tested was fixed after the numbers and is declared "
            "post-hoc in each verdict record."),
        "population": NAME,
        "design_key": KEY,
        "derived_from": {
            **design["derived_from"],
            "rows_file": rows["derived_from_rows"]["file"],
            "rows_sha256": rows["derived_from_rows"]["sha256"],
            "fields_changed": rows["derived_from_rows"]["fields_changed"],
            "fields_copied_unchanged":
                rows["derived_from_rows"]["fields_copied_unchanged"],
            "prepass_rerun": rows["prepass_rerun"],
            "why": ("the comparator degrees are the parent's, so the two "
                    "campaigns are paired orbit by orbit and any difference is "
                    "attributable to the ceiling"),
        },
        "ceiling_probe": prereg["ceiling_probe"],
        "base_scope": (
            "inherited from the parent unchanged: two policies, the truth and "
            "the critical-degree comparator, which is what the ladder reads. "
            "The four schedule policies of the convergence study are not "
            "propagated and nothing is quoted from them. The index is "
            "assembled by rev32_reindex_two_policy_base.py, because the pinned "
            "convergence driver builds its summary from those four."),
        "pooling": (
            "never pooled: not with the strata, not with the coverage designs, "
            "and not with its own capped parent. The parent and this control "
            "are reported side by side, which is what the registration "
            "commits to."),
        "naming_exception": (
            "the ladder reuses the archived drivers, so its outputs carry the "
            "r14_, r18_ and r19_ prefixes with design key OEU and an explicit "
            "budget suffix on every budget including beta = 1. No earlier "
            "manifest claims them: R14 and R18 own allow-lists over designs A "
            "and B, and R19 owns the subtrees with no budget suffix."),
        "operating_point_note": (
            "the accuracy-target operating point is regenerated under the R12 "
            "rule, as for every population without an R12 campaign, and is not "
            "the parent's: with the ceiling lifted the rule allocates "
            "differently, reaching degree 465 where the parent was clamped at "
            "300. The regeneration rule itself was validated against the "
            "archived R12 configurations of designs A and B in R29, indexed "
            "here under reused inputs."),
        "scoring_note": (
            "rev30_verdict.py cannot score this campaign and does not fail "
            "when asked to: it tallies the interior member rather than the "
            "radial endpoint, and maps outcomes by sorted(outcomes)[:3] over a "
            "registration that declares four. The record it produced is kept "
            "as r38_verdict_rev30_generic.json and is superseded by "
            "rev38_verdict.py, one record per budget."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "result_json": results,
        "generated_tables": tables,
        "panel_completeness": {
            NAME: {
                "design_key": KEY,
                "base": {"complete": conv.get("complete"),
                         "orbits": len(conv.get("rows", [])),
                         "failures": len(conv.get("failures", []))},
                "budgets": budgets,
            }
        },
        "trajectory_tree": tree_map,
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in tree_map.values())
    n_raw = sum(t["n_raw_arrays"] for t in tree_map.values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}  "
          f"{n_side} sidecars, {n_raw} raw arrays")
    print(f"  base {len(conv.get('rows', []))}/64 "
          f"complete={conv.get('complete')}  budgets={sorted(budgets)}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
