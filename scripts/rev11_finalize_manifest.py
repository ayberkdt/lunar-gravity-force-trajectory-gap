"""SHA-256 integrity manifest for the R11 extension.

The R10 manifest (``r10_final_experiment_manifest.json``) indexes the
confirmatory campaign only. The R11 extension adds the full-population
vector-tolerance rerun, the independent design-B replication, the 24-phase
dispersion, the mission-arc geometric/symplectic control, and the
corrected-blend vector-tolerance rerun. This script indexes those artifacts so
every number added to the manuscript is traceable.

The Lunaris numerical kernel is unchanged: every R11 trajectory records tag
``paper-truncation-v1.0`` at commit ``27e9ab86...`` with the same kernel and
gravity-file hashes as R10. What is new here is the driver scripts and their
outputs, which live in the manuscript repository rather than in Lunaris.

Usage:  python rev11_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r11_final_experiment_manifest.json"

SCRIPTS = [
    "rev11_full_convergence.py",
    "rev11_designB_convergence.py",
    "rev11_geometric_verification.py",
    "rev11_phase_sweep.py",
    "rev11_blend_lro_vector.py",
    "rev11_manuscript_tables.py",
    "rev11_finalize_manifest.py",
    "run_r11_overnight.sh",
    "run_r11_designB.sh",
    "run_r11_finish.sh",
    "run_r11_postprocess.sh",
    "wait_designB.sh",
]

RESULT_JSON = [
    "r11_full_convergence.json",
    "r11_designB_convergence.json",
    "r11_designB_rows.json",
    "r11_geometric_verification.json",
    "r11_phase_sweep.json",
    "r11_blend_lro_vector.json",
    "r11_manuscript_descriptives.json",
    # supporting provenance for the compact-rule exponent quoted in the text:
    # review_postprocess.json holds p_fit(eps=1e-2)=1.759 and p_fit(eps=1e-3)=1.800
    # for the 50-300 km grid, previously not indexed in any manifest.
    "review_postprocess.json",
]

TABLES = [
    "r11_full64_convergence_table.tex",
    "r11_full64_per_orbit_table.tex",
    "r11_designB_convergence_table.tex",
    "r11_designB_per_orbit_table.tex",
    "r11_phase_sweep_table.tex",
    "r11_blend_vector_table.tex",
]

TREES = [
    ("design_A_vector_convergence", METRICS / "r11_cases" / "convergence",
     METRICS / "r11_raw" / "convergence"),
    ("design_B_vector_convergence", METRICS / "r11_cases" / "designB_convergence",
     METRICS / "r11_raw" / "designB_convergence"),
    ("phase_sweep", METRICS / "r11_cases" / "phase_sweep",
     METRICS / "r11_raw" / "phase_sweep"),
    ("blend_lro_vector", METRICS / "r11_cases" / "blend_lro_vector",
     METRICS / "r11_raw" / "blend_lro_vector"),
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
        out[n] = {"sha256": sha(p), "bytes": p.stat().st_size} if p.exists() \
            else {"missing": True}
    return out


def index_tree(case_dir: Path, raw_dir: Path) -> dict:
    """Per-file hashes for sidecars, rolled-up digest for the raw arrays."""
    sidecars = {}
    for p in sorted(case_dir.rglob("*.json")):
        if "invalid" in p.name or "smoke" in p.name:
            continue
        sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n_raw = 0
    for p in sorted(raw_dir.rglob("*.npz")):
        if "invalid" in p.name or "smoke" in p.name:
            continue
        roll.update(sha(p).encode())
        n_raw += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars,
            "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    payload = {
        "schema": "r11_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat()
                       .replace("+00:00", "Z"),
        "scope": ("R11 extension: full-population design-A vector-tolerance "
                  "rerun, independent design-B replication, 24-phase "
                  "dispersion, mission-arc geometric/symplectic control, and "
                  "corrected-blend vector-tolerance rerun"),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": ("unchanged from R10; every R11 trajectory sidecar records "
                     "this tag and the same kernel and gravity-file hashes"),
        },
        "supersedes_note": ("complements r10_final_experiment_manifest.json, "
                            "which indexes the R10 confirmatory campaign"),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_trees": {},
    }
    for name, case_dir, raw_dir in TREES:
        if case_dir.exists():
            payload["trajectory_trees"][name] = index_tree(case_dir, raw_dir)
            t = payload["trajectory_trees"][name]
            print(f"[tree] {name}: {t['n_sidecars']} sidecars, "
                  f"{t['n_raw_arrays']} raw arrays")
        else:
            payload["trajectory_trees"][name] = {"missing": True}

    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256={payload['manifest_sha256'][:16]}")
    missing = [k for sec in ("scripts", "result_json", "generated_tables")
               for k, v in payload[sec].items() if v.get("missing")]
    missing += [f"tree:{k}" for k, v in payload["trajectory_trees"].items()
                if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
