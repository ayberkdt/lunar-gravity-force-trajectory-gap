"""SHA-256 integrity manifest for R50: the paired radial-span ladder.

R50 is one population read at four prescribed apolunes, drawn twice. Sixteen
orbit identities per block, each emitted at 300, 600, 1200 and 2400 km with
every other element of the identity held fixed to the bit, propagated through
four registered budgets and two amended ones.

Two things about the shape of this campaign decide what the manifest has to
carry. First, a level is not a file: a block's records are one population of
sixty-four orbits carrying all four levels, so every level-wise number in the
manuscript is an aggregation over per-orbit rows and the tables that perform it
are indexed here beside the records they read. Second, the amendment that added
the two highest budgets is itself evidence, because it claims to have been
written before any ladder existed; its record is indexed so that the claim can
be checked against the digests rather than taken on trust.

The ceiling probe is indexed here rather than left for the control it was run
for. It reads only R50 records and it exists whether or not that control is ever
propagated, so leaving it unowned would put an unclaimed measurement in the
metrics tree.

Budgets are read off the disk, and a budget counts for a block only when all
three ladder stages are present for it. A manifest that named a budget the
campaign did not reach would record an intention rather than an archive.

Usage:  python rev50_finalize_manifest.py
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

REG = "r50"
BLOCKS = {"span_ladder_a": "RS1", "span_ladder_b": "RS2"}
K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev50_span_ladder_freeze.py", "rev50_campaign.py",
           "rev50_budget_extension_amendment.py", "rev50_tables.py",
           "rev50_ladder_geometry.py",
           "rev50_verdict.py", "rev50_finish_chain.py",
           "rev50_finalize_manifest.py", "rev51_ceiling_probe.py",
           "rev30_stratum_base.py", "rev30_stratum_ops.py",
           "rev32_reindex_two_policy_base.py", "population_registry.py",
           "rev29_designC_operating_point.py",
           "rev11_full_convergence.py", "rev11_designB_convergence.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py",
           "make_figures_r36_regime.py"]

# Read, never written by this campaign.
REUSED = ["r14_budget_pareto.json", "r14_preregistration.json",
          "r29_designC_operating_point.json", "r30_preregistration.json",
          "r28_calibration_amendment.json"]

SCOPE = ("R50: a paired apolune ladder. Sixteen orbit identities per block, "
         "each flown at 300, 600, 1200 and 2400 km with the perilune, the "
         "inclination, the argument of perilune, the requested perilune "
         "longitude, the derived right ascension and the true anomaly at epoch "
         "identical across the four members, so that radial span is the only "
         "geometric quantity that differs. Run to replace a direction with a "
         "measured dependence on one factor, on the axis the paper names as "
         "its clearest geometric discriminator.")


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


def ladder_trees(key: str, tag: str) -> list[str]:
    rels = [f"r14_cases/{key}_{tag}", f"r14_raw/{key}_{tag}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{key}_{tag}_k_{k}", f"r18_raw/{key}_{tag}_k_{k}"]
    rels += [f"r19_cases/{key}_workmatched_{tag}",
             f"r19_raw/{key}_workmatched_{tag}"]
    return rels


def complete_budgets(key: str) -> list[str]:
    out = []
    rx = re.compile(rf"^r19_equal_total_work_{key}_beta_(\d+\.\d+)\.json$")
    for p in sorted(METRICS.glob(f"r19_equal_total_work_{key}_beta_*.json")):
        m = rx.match(p.name)
        if not m:
            continue
        tag = f"beta_{m.group(1)}"
        if all((METRICS / f"{stem}_{key}_{tag}.json").exists()
               for stem in ("r14_trajectory", "r18_span_sweep")):
            out.append(tag)
    return out


def pairing_check(name: str) -> dict:
    """Re-verify the claim the design rests on, from the frozen file.

    The freezing script checks this before writing, but a manifest that only
    repeated the claim would be indexing an assertion. Recomputing it here means
    the gate fails if the design file is ever edited into one whose identities
    do not hold.
    """
    d = json.loads((METRICS / f"{REG}_{name}_design_frozen.json").read_text(
        encoding="utf-8"))
    held = d["held_fixed_within_identity"]
    groups = {}
    for o in d["orbits"]:
        groups.setdefault(o["identity_index"], []).append(o)
    bad = []
    for identity, members in groups.items():
        for field in held:
            if len({m[field] for m in members}) != 1:
                bad.append({"identity": identity, "field": field})
    return {"identities": len(groups),
            "levels_per_identity": sorted({len(v) for v in groups.values()}),
            "fields_held_fixed": held,
            "violations": bad,
            "verified": not bad}


def main() -> int:
    OUT = METRICS / f"{REG}_final_experiment_manifest.json"
    prereg_p = METRICS / f"{REG}_preregistration.json"
    if not prereg_p.exists():
        print(f"[abort] {prereg_p.name} missing; nothing to manifest")
        return 2
    prereg = json.loads(prereg_p.read_text(encoding="utf-8"))

    registration_names = [f"{REG}_preregistration.json",
                          f"{REG}_budget_extension_amendment.json"]
    registration_names += [f"{REG}_{n}_design_frozen.json" for n in BLOCKS]
    registration, a1 = index_files(registration_names, METRICS, required=True)
    scripts, a2 = index_files(SCRIPTS, CODE, required=True)
    reused, a3 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3

    result_names = [f"{REG}_campaign_progress.json", f"{REG}_verdict.json",
                    "r51_ceiling_probe.json",
                    f"{REG}_span_ladder_primary_table.tex",
                    f"{REG}_span_ladder_secondary_table.tex",
                    f"{REG}_span_ladder_macros.tex",
                    # the geometry the ladder's levels carry with them, read
                    # from the two frozen designs and propagated from nothing
                    f"{REG}_ladder_geometry.json",
                    f"{REG}_ladder_geometry_table.tex"]
    # r36_regime_map.json is not claimed here. The regime map reads this
    # campaign's records and draws four of its levels as rows, but the figure
    # and its record belong to R35, which owns the generator; claiming it in
    # two manifests would make the owner a matter of which one was sealed last.
    tree_map, blocks_out, pairing = {}, {}, {}

    for name, key in BLOCKS.items():
        conv_p = METRICS / f"{REG}_{name}_convergence.json"
        if not conv_p.exists():
            print(f"[abort] {conv_p.name} missing; the block is not complete")
            return 2
        conv = json.loads(conv_p.read_text(encoding="utf-8"))
        if not conv.get("complete"):
            print(f"[abort] {name} base is incomplete; nothing is manifested "
                  f"from a partial base")
            return 2
        result_names += [f"{REG}_{name}_rows.json",
                         f"{REG}_{name}_convergence.json",
                         f"{REG}_{name}_operating_point.json",
                         f"{REG}_budget_pareto_{name}.json"]
        tree_map[f"{name}/base"] = index_tree(
            [f"r11_cases/stratum_{name}_convergence",
             f"r11_raw/stratum_{name}_convergence"])
        tree_map[f"{name}/operating_point"] = index_tree(
            [f"{REG}_cases/atallah_{name}"])
        pairing[name] = pairing_check(name)

        budgets = {}
        for tag in complete_budgets(key):
            result_names += [f"r14_trajectory_{key}_{tag}.json",
                             f"r18_span_sweep_{key}_{tag}.json",
                             f"r19_equal_total_work_{key}_{tag}.json"]
            tree_map[f"{key}/{tag}"] = index_tree(ladder_trees(key, tag))
            q = json.loads((METRICS / f"r19_equal_total_work_{key}_{tag}.json"
                            ).read_text(encoding="utf-8"))
            r = json.loads((METRICS / f"r14_trajectory_{key}_{tag}.json"
                            ).read_text(encoding="utf-8"))
            budgets[tag] = {
                "interior_vs_work_matched_constant": q.get("summary"),
                "radial_endpoint_vs_constant": r.get("summary"),
            }
        blocks_out[name] = {"design_key": key,
                            "orbits": len(conv.get("rows", [])),
                            "base_complete": conv.get("complete"),
                            "budgets": budgets}

    results, a4 = index_files(sorted(set(result_names)), METRICS, required=True)
    absent += a4

    amendment = json.loads(
        (METRICS / f"{REG}_budget_extension_amendment.json").read_text(
            encoding="utf-8"))
    verdict = json.loads((METRICS / f"{REG}_verdict.json").read_text(
        encoding="utf-8"))

    payload = {
        "schema": f"{REG}_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The identities, the four apolune levels, the four "
            "budgets, the readouts, the verdict rule and four outcome classes "
            "were fixed in r50_preregistration.json before the first "
            "trajectory, one outcome class contradicting a sentence the "
            "manuscript prints. Budgets 1.25 and 1.50 were added by the "
            "amendment indexed here, which records the state of the campaign "
            "at the moment it was written and refuses to write itself if any "
            "ladder record, trajectory record or calibration exists."),
        "blocks": blocks_out,
        "pairing_verified": pairing,
        "amendment": {
            "file": f"{REG}_budget_extension_amendment.json",
            "added_budgets": amendment["added_budgets"],
            "registered_budgets": amendment["registered_budgets"],
            "state_when_written": amendment["campaign_state_at_amendment"],
            "why_not_post_hoc": amendment["why_this_is_not_post_hoc"],
        },
        "outcome": {
            "class": verdict["outcome"],
            "text": verdict["outcome_text"],
            "evidence": verdict["outcome_evidence"],
            "note": ("scored by rev50_verdict.py under the registered rule "
                     "rather than by inspection; the rule was not relaxed "
                     "after the numbers were known"),
        },
        "ceiling_probe": {
            "file": "r51_ceiling_probe.json",
            "why_here": ("it reads only R50 records and propagates nothing, "
                         "and it exists whether or not the cap-lifted control "
                         "is ever run, so it is claimed here rather than left "
                         "for a campaign that may not exist"),
        },
        "base_scope": (
            "two policies, the truth and the critical-degree comparator, which "
            "is what the ladder reads. The four schedule policies of the "
            "convergence study are not propagated and nothing is quoted from "
            "them; the reduction is declared in r30_preregistration.json and "
            "inherited. The index is assembled by "
            "rev32_reindex_two_policy_base.py, because the pinned convergence "
            "driver builds its summary from those four."),
        "pooling": (
            "the two blocks are pooled level by level, which the registration "
            "licenses because they are the same design at the same levels, and "
            "the unpooled reading is reported beside the pooled one. Levels are "
            "never pooled with each other, and this population is never pooled "
            "with the coverage designs, the strata or the operational "
            "elliptical population, whose boxes it overlaps."),
        "level_aggregation": (
            "no per-level record exists on disk: a block's records are one "
            "population of sixty-four orbits carrying all four levels. Every "
            "level-wise number in the manuscript is aggregated from per-orbit "
            "rows by rev50_tables.py, rev50_verdict.py and "
            "make_figures_r36_regime.py, all three indexed here."),
        "registration": registration,
        "scripts": scripts,
        "reused_inputs": reused,
        "results": results,
        "trees": tree_map,
    }
    # The separators matter: the auditor recomputes this seal with the
    # archive's convention, and a manifest hashed with the default spacing
    # reports as hand-edited after sealing.
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in tree_map.values())
    n_raw = sum(t["n_raw_arrays"] for t in tree_map.values())
    print(f"[written] {OUT.name}")
    print(f"  blocks: {', '.join(f'{n} ({k})' for n, k in BLOCKS.items())}")
    for name, b in blocks_out.items():
        print(f"  {name}: {len(b['budgets'])} budgets "
              f"{', '.join(sorted(b['budgets']))}")
    for name, p in pairing.items():
        print(f"  pairing {name}: {p['identities']} identities x "
              f"{p['levels_per_identity']} levels, verified={p['verified']}")
    print(f"  {len(results)} result files, {n_side} sidecars, "
          f"{n_raw} raw arrays")
    print(f"  outcome: {verdict['outcome']}")
    if absent:
        print(f"[FAIL] {len(absent)} required files absent: {absent[:6]}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
