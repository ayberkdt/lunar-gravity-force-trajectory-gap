"""SHA-256 integrity manifest for the R15 audit-response items.

R15 does not re-run R14 to obtain a different answer. Each item either
strengthens a comparator, removes an information advantage given to one policy,
or tests a sampling assumption the campaign rested on. The R14 protocol stands as
issued; R15 adds a separately hashed amendment covering the rules these items
needed.

Contents: the frozen amendment; the output-cadence convergence of the
force-defect statistic (R15-D); the deployable budget calibrations, two per
population (R15-B); the budget-saturating and best-under-budget fixed comparators
(R15-A); and the grid convergence of the tolerance sample (R15-E).

Usage:  python rev15_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r15_final_experiment_manifest.json"

SCRIPTS = [
    "rev15_preregister_amendment.py",
    "rev15_cadence_check.py",
    "rev15_deployable_calibration.py",
    "rev15_fixed_oracle.py",
    "rev15_atallah_grid.py",
    "rev15_tables.py",
    "wait_for_idle.py",
    "check_labels.py",
    "check_assets.py",
    "run_r15_queue.sh",
    "rev15_finalize_manifest.py",
]

RESULT_JSON = [
    "r15_preregistration_amendment.json",
    "r15_cadence_check.json",
    "r15_deployable_calibration_A.json",
    "r15_deployable_calibration_B.json",
    "r15_fixed_oracle.json",
    "r15_atallah_grid.json",
]

TABLES = [
    "r15_cadence_check_table.tex",
    "r15_deployable_table.tex",
    "r15_fixed_oracle_table.tex",
    "r15_atallah_grid_table.tex",
]


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


def index_tree(case_dir: Path, raw_dir: Path) -> dict:
    sidecars = {}
    for p in sorted(case_dir.rglob("*.json")):
        sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n_raw = 0
    for p in sorted(raw_dir.rglob("*.npz")):
        roll.update(sha(p).encode())
        n_raw += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    payload = {
        "schema": "r15_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "scope": ("R15 audit response: output-cadence convergence of the "
                  "force-defect statistic (R15-D), deployable truth-free budget "
                  "calibration (R15-B), budget-saturating and best-under-budget "
                  "fixed comparators (R15-A), and grid convergence of the "
                  "tolerance sample (R15-E)."),
        "relationship_to_r14": (
            "additive. No R14 result is recomputed to a different answer; each "
            "R15 item strengthens a comparator, removes an information advantage, "
            "or tests a sampling assumption. Three moved the outcome slightly "
            "further against the radial rule; R15-E failed its convergence "
            "criterion and its consequence is bounded and reported."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R14"},
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_trees": {},
    }
    prereg = METRICS / "r14_preregistration.json"
    amend = METRICS / "r15_preregistration_amendment.json"
    if prereg.exists() and amend.exists():
        payload["preregistration"] = {
            "parent": json.loads(prereg.read_text(encoding="utf-8"))["protocol_sha256"],
            "amendment": json.loads(amend.read_text(encoding="utf-8"))["amendment_sha256"],
            "note": "the R14 protocol is unchanged; the amendment is additive"}
    case_root, raw_root = METRICS / "r15_cases", METRICS / "r15_raw"
    if case_root.exists():
        for case_dir in sorted(case_root.iterdir()):
            if not case_dir.is_dir():
                continue
            t = index_tree(case_dir, raw_root / case_dir.name)
            payload["trajectory_trees"][case_dir.name] = t
            print(f"[tree] {case_dir.name}: {t['n_sidecars']} sidecars, "
                  f"{t['n_raw_arrays']} raw arrays")
    # the fixed-comparator ladder writes raw arrays without sidecars
    fo = raw_root / "fixed_oracle"
    if fo.exists():
        roll = hashlib.sha256()
        n = 0
        for p in sorted(fo.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
        payload["trajectory_trees"]["fixed_oracle_ladder"] = {
            "n_sidecars": 0, "n_raw_arrays": n,
            "raw_rollup_sha256": roll.hexdigest(),
            "note": "degree-ladder arrays; per-orbit results in r15_fixed_oracle.json"}
        print(f"[tree] fixed_oracle_ladder: {n} raw arrays")
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256={payload['manifest_sha256'][:16]}")
    missing = [k for sec in ("scripts", "result_json", "generated_tables")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[note] not produced: " + ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
