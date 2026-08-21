"""SHA-256 integrity manifest for R61 (O42-ext): the scoring-tolerance
rematch carried to the third coverage design and the five geometry strata.

R61 changes nothing about R44's method and everything about its population
set. It reuses each population's R18 span members and their tighter-level
telemetry unchanged and propagates one new constant-degree trajectory per
orbit per tolerance level, matched on realized total quadratic work at the
tighter level, in eighteen registered population-budget cells.

Partition. R44's eight A/B cells and its sealed manifest are untouched, and
no r44_* file is written: the driver redirects the case root, the raw root
and the record path to r61_*. The R19 records this campaign pairs against are
read-only inputs and are indexed here as such, under their own manifests
elsewhere.

Usage:  python rev61_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r61_final_experiment_manifest.json"

KEYS = ("C", "SL", "SP", "SE", "SF", "SH")
BETAS = ("0.50", "0.75", "1.00")
CELLS = [(k, b) for b in BETAS for k in KEYS]

SCRIPTS = ["rev61_equal_work_tighter_ext.py", "rev61_preregister.py",
           "rev61_campaign.py", "rev61_tables.py",
           "rev61_finalize_manifest.py", "rev44_equal_work_tighter.py",
           "population_registry.py"]
RESULT_JSON = ([f"r61_equal_work_tighter_{k}_beta_{b}.json"
                for k, b in CELLS]
               + ["r61_manuscript_descriptives.json",
                  "r61_campaign_progress.json"])
REGISTRATION = ["r61_preregistration.json"]
TABLES = ["r61_equal_work_table.tex", "r61_bracket_shift_table.tex"]
REUSED = ([f"r18_span_sweep_{k}_beta_{b}.json" for k, b in CELLS]
          + [f"r19_equal_total_work_{k}_beta_{b}.json" for k, b in CELLS]
          + ["r26_designC_rows.json"]
          + [f"r30_{s}_rows.json" for s in
             ("low_perilune", "polar", "equatorial", "frozen_like",
              "high_apolune")])


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path) -> dict:
    out = {}
    for n in names:
        p = base / n
        out[n] = ({"sha256": sha(p), "bytes": p.stat().st_size}
                  if p.exists() else {"missing": True})
    return out


def index_tree() -> dict:
    """R61 owns its whole tree: every case directory carries both the
    population key and the workmatched_tighter tag, so no suffix partition
    against another campaign is needed."""
    sidecars = {}
    root = METRICS / "r61_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r61_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome_note() -> dict:
    """What the campaign returned, read from the descriptives rather than
    restated by hand."""
    p = METRICS / "r61_manuscript_descriptives.json"
    if not p.exists():
        return {"missing": True}
    d = json.loads(p.read_text(encoding="utf-8"))
    moved = {k: v["population"] for k, v in d.items() if v["bracket_moved"]}
    boundary = {k: v["boundary_cells"] for k, v in d.items()
                if v["boundary_cells"]}
    return {
        "brackets": {k: {"tight_level": v["bracket_tight_level"],
                         "scoring_tolerance": v["bracket_level_consistent"],
                         "moved": v["bracket_moved"]}
                     for k, v in d.items()},
        "moved_count": len(moved),
        "moved": sorted(moved.values()),
        "boundary_cells": boundary,
    }


def main() -> int:
    payload = {
        "schema": "r61_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R61 (O42-ext): the interior span member (k = 0.5) against "
                  "a constant degree matched on realized total quadratic work "
                  "at the tighter tolerance, the level every error in the "
                  "comparison is scored at, in eighteen registered "
                  "population-budget cells covering the third coverage design "
                  "and the five geometry strata."),
        "why": ("R44 established this match on designs A and B only, and one "
                "of the two moved its tally crossing above 0.75. Section IX.B "
                "nevertheless reads the crossing bracket across seven "
                "populations and attributes it to the matching convention as "
                "much as to the dynamics. This campaign measures the "
                "convention on the populations the sentence generalises "
                "over."),
        "relationship_to_r44": (
            "method identical, population set disjoint, records disjoint. The "
            "R44 cells, tree and manifest are untouched and no r44_* file is "
            "written; r18_cases is read read-only."),
        "relationship_to_r18_r19": (
            "additive. The member trajectories are not re-propagated; their "
            "archived tighter-level telemetry supplies the work target. The "
            "R19 records are the paired tight-level baseline and are inputs, "
            "not outputs: the tight-level counterpart of an R61 cell is the "
            "R19 cell, never the R30 ladder, which is a different comparator."
        ),
        "comparator_rule": (
            "N* = round(N_0 * sqrt(W_k/W_0)) with both works at the tighter "
            "level; where the beta-specific fixed_budget sidecar is not "
            "addressable through r14.reuse_paths, n_RHS for W_0 comes from "
            "the orbit's fixed_critical run and the source policy is recorded "
            "in the case config; the achieved ratio is measured from the "
            "propagated runs at both levels. Comparators at or above the "
            "adopted reference degree are censored, not clamped; no orbit "
            "required this in any cell."),
        "completion": (
            "all eighteen registered cells completed, 64 orbits each, 0 "
            "failed, 0 censored; the window ran 2026-08-19 02:14 to 09:19 "
            "local at 10 workers. Per-cell wall times are in "
            "r61_campaign_progress.json."),
        "declared_outcome_returned": (
            "B (see r61_preregistration.json): further populations move the "
            "crossing above 0.75. Two of the three bracket moves are stable "
            "across the resolution cuts and one is a boundary cell; the "
            "per-cell re-tallies are in r61_manuscript_descriptives.json and "
            "the boundary cells are marked in both generated tables."),
        "outcome_as_measured": outcome_note(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R44"},
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": index_tree(),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    t = payload["trajectory_tree"]
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  {t['n_sidecars']} sidecars, {t['n_raw_arrays']} raw arrays")
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "reused_inputs", "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
