"""SHA-256 integrity manifest for R52: the cap-lifted control on block B.

R52 is not a new population. It is block B of the R50 paired ladder propagated
a second time with one number changed, the adopted reference degree, 300 to
600, so that the budget-calibrated radial schedule is never clamped to the
reference model. The identities, the four apolune levels, the initial states
and the comparator degrees are the parent's.

That makes the parent an input rather than a neighbour, and the manifest treats
it as one: R50's registration, its frozen design for block B, its prepass rows
and the ceiling probe that fixed this control's reference degree are indexed
under reused inputs, so an edit to any of them fails this gate as well as
R50's.

The probe is worth naming twice, because it is what makes this campaign
answerable rather than merely expensive. It measured, before anything
propagated, that the ceiling binds at neither level inside the factor box, at
some orbits of the 1200 km level and at every orbit of the 2400 km level. The
confound this control removes is therefore aligned with the very axis the
ladder varies, which is the strongest reason to remove it by propagation
rather than by subsetting.

Budgets are read off the disk, and a budget counts only with all three ladder
stages present, so a control stopped by its clock records what it reached
rather than what it intended.

Usage:  python rev52_finalize_manifest.py
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

REG = "r52"
NAME = "span_ladder_b_uncapped"
KEY = "RS2U"
PARENT_REG, PARENT_NAME, PARENT_KEY = "r50", "span_ladder_b", "RS2"
K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev52_uncapped_freeze.py", "rev52_campaign.py",
           "rev51_ceiling_probe.py", "rev52_finalize_manifest.py",
           "rev30_stratum_base.py", "rev30_stratum_ops.py",
           "rev32_reindex_two_policy_base.py", "population_registry.py",
           "rev29_designC_operating_point.py",
           "rev11_full_convergence.py", "rev11_designB_convergence.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]

REUSED = [f"{PARENT_REG}_preregistration.json",
          f"{PARENT_REG}_{PARENT_NAME}_design_frozen.json",
          f"{PARENT_REG}_{PARENT_NAME}_rows.json",
          f"{PARENT_REG}_budget_pareto_{PARENT_NAME}.json",
          f"{PARENT_REG}_{PARENT_NAME}_convergence.json",
          "r52_ceiling_probe.json",
          "r14_budget_pareto.json", "r14_preregistration.json",
          "r29_designC_operating_point.json", "r30_preregistration.json"]

SCOPE = ("R52: block B of the R50 paired radial-span ladder re-propagated with "
         "its adopted reference degree raised from 300 to 600, so that the "
         "budget-calibrated radial schedule is never clamped to the reference "
         "model. A paired control on the axis the ladder measures, run because "
         "the ceiling was found to bind at the wide levels and nowhere else, "
         "which aligns it with that axis. The replication partner of R51.")


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
    prereg_p = METRICS / f"{REG}_preregistration.json"
    if not prereg_p.exists():
        print(f"[abort] {prereg_p.name} missing; nothing to manifest")
        return 2
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))
    design = json.loads(
        (METRICS / f"{REG}_{NAME}_design_frozen.json").read_text(
            encoding="utf-8"))
    rows = json.loads((METRICS / f"{REG}_{NAME}_rows.json").read_text(
        encoding="utf-8"))

    conv_p = METRICS / f"{REG}_{NAME}_convergence.json"
    if not conv_p.exists():
        print(f"[abort] {conv_p.name} missing; the base has not been indexed. "
              f"A manifest sealed now would claim a campaign still running.")
        return 2
    conv = json.loads(conv_p.read_text(encoding="utf-8"))
    if not conv.get("complete"):
        print(f"[abort] the base is incomplete; nothing is manifested from a "
              f"partial base")
        return 2

    registration, a1 = index_files(
        [f"{REG}_preregistration.json", f"{REG}_{NAME}_design_frozen.json"],
        METRICS, required=True)
    scripts, a2 = index_files(SCRIPTS, CODE, required=True)
    reused, a3 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3

    result_names = [f"{REG}_{NAME}_rows.json", f"{REG}_{NAME}_convergence.json",
                    f"{REG}_{NAME}_operating_point.json",
                    f"{REG}_budget_pareto_{NAME}.json",
                    f"{REG}_campaign_progress.json",
                    f"{REG}_verdict.json",
                    # R51 printed its block against its own capped parent;
                    # this campaign's table is the two blocks side by side, so
                    # the generated names differ from the ones inherited here.
                    f"{REG}_block_replication_table.tex",
                    f"{REG}_block_replication_macros.tex",
                    # The departure from the registered consequence, and the
                    # wider probe sweeps the supplement quotes. A deviation
                    # argued only in prose is the one a reader cannot check,
                    # and a probe number quoted from an unsealed file is the
                    # one that can drift.
                    f"{REG}_registration_departure.json"]
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
    absent += a4
    registered = [f"beta_{b:.2f}" for b in prereg["budgets"]]
    not_run = [t for t in registered if t not in budgets]

    payload = {
        "schema": f"{REG}_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The reference degree, the budget order, the "
            "verdict rule and four outcomes were fixed in "
            f"{REG}_preregistration.json before the first trajectory, and "
            "block A's ceiling-free result is pinned there by digest so this "
            "registration cannot be read as having been written against a "
            "different one. The sufficiency of degree 600 was established "
            "before propagation by a demand probe run on block B rather than "
            "inherited from block A, and the freezing script refuses to "
            "register this control at all if that probe reports either that "
            "600 fails to clear the demand or that the ceiling never binds."),
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
            "why": ("the comparator degrees and the identities are the "
                    "parent's, so the two campaigns are paired orbit by orbit "
                    "and level by level, and anything that moves is "
                    "attributable to the ceiling"),
        },
        "ceiling_probe": prereg["ceiling_probe"],
        "levels_km": design["apolune_levels_km"],
        "budgets_complete": sorted(budgets),
        "budgets_registered_not_run": not_run,
        "budgets": budgets,
        "base_scope": (
            "inherited from the parent unchanged: the truth and the "
            "critical-degree comparator, which is what the ladder reads. The "
            "index is assembled by rev32_reindex_two_policy_base.py, because "
            "the pinned convergence driver builds its summary from the four "
            "schedule policies this scope switches off."),
        "pooling": (
            "never pooled: not with its capped parent, not with the second "
            "block of that parent, not with the strata or the coverage "
            "designs. The capped and uncapped ladders are read level by level "
            "against each other and printed together wherever a wide-level "
            "ratio is quoted."),
        "registration": registration,
        "scripts": scripts,
        "reused_inputs": reused,
        "results": results,
        "trees": tree_map,
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in tree_map.values())
    n_raw = sum(t["n_raw_arrays"] for t in tree_map.values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  base {len(conv.get('rows', []))}/64 complete={conv.get('complete')}")
    print(f"  budgets complete: {', '.join(sorted(budgets)) or 'none'}")
    if not_run:
        print(f"  registered but not run: {', '.join(not_run)}")
    print(f"  {len(results)} result files, {n_side} sidecars, {n_raw} raw arrays")
    if absent:
        print(f"[FAIL] {len(absent)} required files absent: {absent[:6]}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
