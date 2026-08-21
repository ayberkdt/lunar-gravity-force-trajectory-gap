"""SHA-256 integrity manifest for R30: the geometry strata.

R30 propagates sub-populations of the same frozen factor box: four 64-orbit
strata drawn with the pinned orbit map on a restricted sub-box, to locate the
allocation result rather than to replicate it. Each stratum has its own frozen
design, its own base, its own regenerated accuracy-target operating point, its
own Phase-A calibration, and one ladder per budget the campaign reached.

Everything indexed here is new. The strata write into trees no earlier campaign
names (r11_cases/stratum_*, r30_cases/*) and their ladders carry design keys
SP, SE, SH and SF, which appear in no manifest's allow-list. As with design C,
even beta = 1 is written with an explicit budget suffix, because the R19
manifest claims the subtrees that carry none.

Two scope facts are carried in the manifest rather than left to the reader.
A stratum base holds two policies, not six -- the truth and the critical-degree
comparator, which is what the ladder reads -- and that reduction is declared in
the registration this manifest hashes. And the operating-point regeneration was
validated once, against the archived R12 configurations of designs A and B, in
R29; the record of that check is indexed here as a reused input so a later edit
to it would fail this gate too.

Strata and budgets are read off the disk. A campaign that runs to a wall clock
cannot know in advance which of them it will reach, and a manifest that named
one it did not reach would record an intention rather than an archive.

Usage:  python rev30_finalize_manifest.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import population_registry as registry

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"

K_INTERIOR = ("0.25", "0.50", "0.75")

SCRIPTS = ["rev30_stratum_base.py", "rev30_stratum_ops.py",
           "rev30_campaign.py", "rev30_verdict.py", "rev30_tables.py",
           "rev30_finalize_manifest.py", "population_registry.py",
           "rev29_designC_operating_point.py",
           "rev11_full_convergence.py", "rev11_designB_convergence.py",
           "rev14_budget_pareto.py", "rev14_budget_trajectory.py",
           "rev18_span_sweep.py", "rev19_equal_total_work.py"]
FREEZER = {"r30": "rev30_strata_freeze.py",
           "r31": "rev31_operational_freeze.py"}

SCOPE = {
    "r30": ("R30: geometry strata of the frozen factor box. Each stratum is 64 "
            "orbits drawn with the pinned orbit map on a restricted sub-box, "
            "propagated through the same base, calibration and ladder as the "
            "coverage designs."),
    "r31": ("R31: one 64-orbit population of operational elliptical lunar "
            "orbits, perilune 80-120 km with apolune 700-2500 km, bracketing "
            "the two elliptical subsatellite orbits Kaguya flew. This is the "
            "one population that leaves the frozen factor box, and its result "
            "is a scope extension rather than a population-level statement."),
}

# Read, not written: the frozen calibration the strata never merge into, the
# registration whose grid they inherit, the post-hoc declaration beta = 0.62
# carries on every population, and the check that validated the regenerated
# operating point once for all of them.
REUSED = ["r14_budget_pareto.json", "r14_preregistration.json",
          "r28_calibration_amendment.json",
          "r29_designC_operating_point.json"]


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


def claimed_elsewhere(reg: str) -> set[str]:
    """Every metrics-relative path some other campaign manifest already names.

    Reading budgets off the disk was right while this campaign was the only
    one writing into these trees. It stopped being right when a later campaign
    added a budget column to the same strata: re-running this finalizer then
    claims that column too, and two manifests indexing one trajectory breaks
    the partition the integrity check enforces. A budget whose records another
    manifest already owns is therefore skipped and reported, not claimed.
    """
    owned: set[str] = set()
    for mp in sorted(METRICS.glob("r*_final_experiment_manifest.json")):
        if mp.name.startswith(f"{reg}_"):
            continue
        for m in re.finditer(r'"((?:r\d+_(?:cases|raw))/[^"]+)"',
                             mp.read_text(encoding="utf-8")):
            owned.add(m.group(1))
    return owned


def complete_budgets(key: str, owned: set[str]) -> tuple[list[str], list[str]]:
    """A budget counts only with all three ladder stages on disk, and only if
    no other manifest already owns its records."""
    out, skipped = [], []
    rx = re.compile(rf"^r19_equal_total_work_{key}_beta_(\d+\.\d+)\.json$")
    for p in sorted(METRICS.glob(f"r19_equal_total_work_{key}_beta_*.json")):
        m = rx.match(p.name)
        if not m:
            continue
        tag = f"beta_{m.group(1)}"
        if not all((METRICS / f"{stem}_{key}_{tag}.json").exists()
                   for stem in ("r14_trajectory", "r18_span_sweep")):
            continue
        # owned holds file paths; ladder_trees returns the directories above
        # them, so the test is a prefix test and not membership
        prefixes = tuple(f"{rel}/" for rel in ladder_trees(key, tag))
        if any(o.startswith(prefixes) for o in owned):
            skipped.append(tag)
            continue
        out.append(tag)
    return out, skipped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="r30")
    a = ap.parse_args()
    REG = a.registry
    OUT = METRICS / f"{REG}_final_experiment_manifest.json"
    prereg = registry.registration(REG)
    strata = registry.populations(REG)
    scripts_wanted = SCRIPTS + [FREEZER[REG]]

    registration_names = [f"{REG}_preregistration.json"] + [
        s["file"] for s in strata.values()]
    registration, a1 = index_files(registration_names, METRICS, required=True)
    scripts, a2 = index_files(scripts_wanted, CODE, required=True)
    reused, a3 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a3

    result_names, tree_map, completeness, propagated = [], {}, {}, []
    owned_elsewhere = claimed_elsewhere(REG)
    skipped_report: list[str] = []
    for name, spec in strata.items():
        key = spec["design_key"]
        conv = METRICS / f"{REG}_{name}_convergence.json"
        if not conv.exists():
            completeness[name] = {"not_run": True}
            continue
        propagated.append(name)
        result_names += [f"{REG}_{name}_rows.json",
                         f"{REG}_{name}_convergence.json"]
        for extra in (f"{REG}_{name}_operating_point.json",
                      f"{REG}_budget_pareto_{name}.json"):
            if (METRICS / extra).exists():
                result_names.append(extra)
        tree_map[f"{name}_base"] = index_tree(
            [f"r11_cases/stratum_{name}_convergence",
             f"r11_raw/stratum_{name}_convergence"])
        if (METRICS / f"{REG}_cases" / f"atallah_{name}").exists():
            tree_map[f"{name}_operating_point"] = index_tree(
                [f"{REG}_cases/atallah_{name}"])

        d = json.loads(conv.read_text(encoding="utf-8"))
        entry = {"design_key": key,
                 "base": {"complete": d.get("complete"),
                          "orbits": len(d.get("rows", [])),
                          "failures": len(d.get("failures", []))},
                 "budgets": {}}
        tags, skipped = complete_budgets(key, owned_elsewhere)
        if skipped:
            skipped_report.append(f"{name}: {', '.join(skipped)}")
        for tag in tags:
            result_names += [f"r14_trajectory_{key}_{tag}.json",
                             f"r18_span_sweep_{key}_{tag}.json",
                             f"r19_equal_total_work_{key}_{tag}.json"]
            tree_map[f"{name}_{tag}"] = index_tree(ladder_trees(key, tag))
            q = METRICS / f"r19_equal_total_work_{key}_{tag}.json"
            entry["budgets"][tag] = json.loads(
                q.read_text(encoding="utf-8")).get("summary")
        completeness[name] = entry

    results, a4 = index_files(sorted(set(result_names)), METRICS, required=True)
    # One verdict record per budget, plus the single-file form the campaign
    # wrote before it was split. The split exists because the one-file form
    # held whichever budget ran last, so the sealed record could carry a
    # different outcome class from the one the manuscript quotes.
    verdicts = sorted(p.name for p in METRICS.glob(f"{REG}_verdict*.json"))
    tables, _ = index_files([f"{REG}_population_table.tex", *verdicts,
                             f"{REG}_manuscript_descriptives.json"],
                            METRICS, required=False)
    absent += a4

    payload = {
        "schema": f"{REG}_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE[REG],
        "registration_status": (
            "pre-registered. Designs, sub-boxes, seed rule, outcomes, the "
            "no-pooling rule and the reduced base scope were all fixed in "
            "r30_preregistration.json before any stratum propagated."),
        "populations_propagated": propagated,
        "populations_declared": sorted(strata),
        "base_scope": (prereg.get("base_scope")
                       or registry.registration("r30")["base_scope"]),
        "pooling": (prereg.get("pooling")
                    or prereg.get("why_separate_from_R30")),
        "naming_exception": (
            "the ladder reuses the archived drivers, so its outputs carry the "
            "r14_, r18_ and r19_ prefixes with the stratum design keys SP, SE, "
            "SH and SF and a budget suffix. No earlier manifest claims them: "
            "R14 and R18 own explicit allow-lists over designs A and B, and R19 "
            "owns the subtrees with no budget suffix, which is why the strata "
            "write beta = 1 suffixed."),
        "operating_point_note": (
            "each stratum's accuracy-target operating point is regenerated "
            "under the R12 rule, because a stratum has no R12 campaign. The "
            "regeneration was validated against the archived R12 "
            "configurations of designs A and B once, in R29; that record is "
            "indexed here under reused inputs."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "result_json": results,
        "generated_tables": tables,
        "panel_completeness": completeness,
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
    for name in sorted(strata):
        info = completeness[name]
        if info.get("not_run"):
            print(f"  {name:<14} not run")
            continue
        print(f"  {name:<14} base {info['base']['orbits']}/64 "
              f"complete={info['base']['complete']} "
              f"budgets={sorted(info['budgets']) or 'none'}")
    for line in skipped_report:
        print(f"  [skipped, owned by another manifest] {line}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
