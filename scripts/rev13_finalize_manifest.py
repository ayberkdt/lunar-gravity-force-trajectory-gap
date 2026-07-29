"""SHA-256 integrity manifest for the R13 extension.

R13 answers one question about the R12 benchmark: what the unresolved
matched-work verdicts mean. It contains the resolution diagnosis, the targeted
third-tolerance-level retest of the contested orbits, the integration-noise-free
force-defect comparison over both populations, its forced-variational
calibration, the extended kernel cost curve, and the measured-time-matched
comparator (including the serial Atallah re-runs that make the timing
comparable).

The numerical kernel is unchanged: every R13 trajectory sidecar records Lunaris
tag ``paper-truncation-v1.0`` at commit ``27e9ab86...`` with the same kernel and
gravity-file hashes as R10--R12. The truth and Atallah trajectories that R13
compares against are indexed in the R11 and R12 manifests.

Usage:  python rev13_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r13_final_experiment_manifest.json"

SCRIPTS = [
    "rev13_resolution_diagnosis.py",
    "rev13_ultratight.py",
    "rev13_force_defect.py",
    "rev13_variational_check.py",
    "rev13_cost_curve_high.py",
    "rev13_timing_match.py",
    "rev13_timing_repair.py",
    "make_figures_r13.py",
    "rev13_finalize_manifest.py",
]

RESULT_JSON = [
    "r13_resolution_diagnosis.json",
    "r13_ultratight_selection.json",
    "r13_ultratight.json",
    "r13_force_defect.json",
    "r13_variational_check.json",
    "r13_kernel_cost_curve_high.json",
    "r13_timing_match_selection.json",
    "r13_timing_match.json",
    "r13_timing_repair.json",
]

TABLES = [
    "r13_resolution_diagnosis_table.tex",
    "r13_force_defect_table.tex",
    "r13_variational_check_table.tex",
    "r13_timing_match_table.tex",
]

FIGURES = ["fig_atallah_benchmark.pdf"]

TREES = [
    ("ultratight_third_level", METRICS / "r13_cases" / "ultratight",
     METRICS / "r13_raw" / "ultratight"),
    ("timing_match", METRICS / "r13_cases" / "timing_match",
     METRICS / "r13_raw" / "timing_match"),
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
        "schema": "r13_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat()
                       .replace("+00:00", "Z"),
        "scope": ("R13 extension: diagnosis of the unresolved matched-work "
                  "Atallah comparisons, targeted third-tolerance-level retest, "
                  "integration-noise-free force-defect comparison with a "
                  "forced-variational calibration, extended kernel cost curve, "
                  "and measured-time-matched comparator with serial timing "
                  "re-runs"),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R12"},
        "reused_evidence": ("truth and critical-degree trajectories from the R11 "
                            "trees; Atallah and work-matched trajectories from "
                            "the R12 trees; both indexed in their own manifests"),
        "timing_note": ("the campaign's Atallah kernel times were recorded under "
                        "five concurrent workers and are not comparable; the "
                        "serial re-runs indexed here are the timing reference"),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "generated_figures": index_files(FIGURES, ROOT / "figures"),
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
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256={payload['manifest_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
