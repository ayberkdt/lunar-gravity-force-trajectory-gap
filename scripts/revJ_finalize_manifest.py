"""SHA-256 integrity manifest for the JGCD replication campaigns (J1--J3).

One manifest for the three campaigns because they were registered, run and
reported as one block by one harness: J1 (cross-solution replication on
GRGM1200A), J2 (expanded-dynamics population replication) and J3 (four-level
tolerance control), together with the harness fidelity check and the archived
re-scorings (budget grids, primary-solution ladder and its skip audit) that
read frozen R-series records and write no trajectory of their own.

Three things a reader would otherwise take on trust are indexed explicitly:

  * the registrations, including the staging. J2's population size was
    enlarged in stages after earlier stages' verdicts were read; every plan
    file of every stage is hashed here, so the timeline the supplement
    discloses is checkable against the sealed records rather than asserted.

  * the raw state arrays. They live outside the repository
    (D:\\makale_raw_offload\\jgcd) because the working drive lacked space, and
    they are not shipped in the package; their per-file digests are recorded
    here in the same way the gravity coefficient products are, so
    regeneration can be verified file by file.

  * the partition. The J sidecar trees (metrics/rJ1_cases, rJ2_cases,
    rJ3_cases) are claimed by this manifest and by no other, and the ladder
    and budget-grid re-scorings claim no trajectory records at all: the
    archived trajectories they read stay sealed under R10/R11/R14.

Usage:  python revJ_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
RAW = Path(r"D:\makale_raw_offload\jgcd")

OUT = METRICS / "rJ_final_experiment_manifest.json"

SCOPE = (
    "J1-J3: pre-submission replication block for the JGCD retarget. J1 "
    "recalibrates the entire recipe on GSFC GRGM1200A and repeats the budget "
    "comparison on a fresh scrambled-Sobol population (seed 20260808, 32 "
    "orbits extended to 64 under the continuation declared in the campaign "
    "plan); J2 repropagates the frozen confirmatory designs under "
    "DE440/MOON_PA with Earth/Sun third-body gravity and eclipsed cannonball "
    "SRP, in a carried-over and a recalibrated arm; J3 propagates a "
    "margin-stratified sample at four tolerance levels. The block also holds "
    "the harness fidelity check (128/128 archived defects recomputed "
    "bit-identically) and the integration-free re-scorings of archived "
    "populations: the budget grids, the primary-solution ladder and its "
    "complete skip audit."
)

REGISTRATION_PREFIXES = (
    "rJ1_preregistration", "rJ1_design", "rJ2_plan", "rJ3_plan",
)


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_named(paths, base: Path) -> dict:
    out = {}
    for p in sorted(paths):
        out[p.relative_to(base).as_posix()] = {
            "sha256": sha(p), "bytes": p.stat().st_size}
    return out


def main() -> int:
    records = sorted(METRICS.glob("rJ*.json"))
    records = [p for p in records if p.name != OUT.name]
    registration = [p for p in records
                    if p.name.startswith(REGISTRATION_PREFIXES)]
    results = [p for p in records if p not in registration]
    tables = sorted(METRICS.glob("rJ_*.tex"))
    scripts = sorted(CODE.glob("revJ*.py"))
    # scripts the harness imports from the sealed campaigns; hashed here at
    # the state this block ran them, indexed under their own manifests too
    # (a script, unlike a trajectory record, may be held by more than one)
    reused_scripts = [CODE / n for n in (
        "rev3_common.py", "rev10_sobol_confirmatory.py", "rev12_atallah.py",
        "rev4_robustness_controls.py", "rev14_budget_pareto.py")
        if (CODE / n).exists()]

    sidecars = {}
    for tree in ("rJ1_cases", "rJ2_cases", "rJ3_cases"):
        base = METRICS / tree
        if base.is_dir():
            sidecars.update(index_named(base.rglob("*.json"), METRICS))

    raw_files = sorted(RAW.rglob("*")) if RAW.is_dir() else []
    raw_files = [p for p in raw_files if p.is_file()]
    raws, roll = {}, hashlib.sha256()
    for p in raw_files:
        digest = sha(p)
        rel = p.relative_to(RAW).as_posix()
        raws[rel] = {"sha256": digest, "bytes": p.stat().st_size,
                     "shipped_in_package": False}
        roll.update(f"{rel}:{digest}\n".encode())

    payload = {
        "schema": "rJ_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "registered before propagation, with one disclosed staging: "
            "population identity was fixed by rule throughout, and J2's "
            "population size was enlarged 24 -> 48 -> 64 after earlier "
            "stages' verdicts were read. Every stage's plan and score are "
            "hashed here and the supplement states the timeline; the later "
            "plan files repeat the first plan's frozen-before-propagation "
            "status line, which is true of the orbits each stage adds and "
            "not of the campaign."),
        "partition_note": (
            "this manifest claims the rJ1/rJ2/rJ3 sidecar trees and no "
            "other trajectory record. The re-scorings (rJ_budget_grid*, "
            "rJ_ladder_primary, rJ_ladder_skip_audit, rJ_fidelity_check, "
            "rJ_field_comparison*) integrate nothing and read archived "
            "trajectories that stay sealed under the R10/R11/R14 manifests."),
        "skip_accounting_note": (
            "rJ_ladder_primary.json stored only the first 50 of its 106 "
            "skip entries; rJ_ladder_skip_audit.json re-derives all 106 "
            "integration-free, verifies the stored 50 identically, and "
            "records what the one-sided censor removed at each budget. Both "
            "are sealed here and the driver now stores the full list."),
        "raw_location_note": (
            "raw state arrays are under D:\\makale_raw_offload\\jgcd, "
            "outside the repository, because the working drive lacked "
            "space; they are regenerable from the sealed drivers and their "
            "per-file digests are recorded below as unshipped inputs."),
        "registration": index_named(registration, METRICS),
        "scripts": index_named(scripts, CODE),
        "reused_scripts": index_named(reused_scripts, CODE),
        "result_json": index_named(results, METRICS),
        "generated_tables": index_named(tables, METRICS),
        "sidecars": sidecars,
        "input_products": raws,
        "raw_rollup": {"files": len(raws),
                       "bytes": sum(v["bytes"] for v in raws.values()),
                       "raw_rollup_sha256": roll.hexdigest()},
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  registration {len(registration)}, results {len(results)}, "
          f"tables {len(tables)}, scripts {len(scripts)}"
          f"+{len(reused_scripts)} reused")
    print(f"  sidecars {len(sidecars)}, raw files {len(raws)} "
          f"({payload['raw_rollup']['bytes'] / 1e6:.0f} MB, not shipped)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
