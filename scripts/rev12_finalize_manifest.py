"""SHA-256 integrity manifest for the R12 extension.

R12 is the direct benchmark of the published Atallah radial-adaptive rule: the
implementation and its verification, the two 64-orbit campaigns and their
tolerance-selection records and degree tables, the per-orbit matching record,
the bin-resolution control, the tolerance sweep, and the kernel cost curve. The
campaigns reuse the R11 truth and critical-degree trajectory trees, which are
indexed in ``r11_final_experiment_manifest.json``; only the new trajectories are
indexed here.

The Lunaris numerical kernel is unchanged: every R12 trajectory sidecar records
tag ``paper-truncation-v1.0`` at commit ``27e9ab86...`` with the same kernel and
gravity-file hashes as R10 and R11.

Usage:  python rev12_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r12_final_experiment_manifest.json"

SCRIPTS = [
    "rev12_atallah.py",
    "rev12_atallah_campaign.py",
    "rev12_atallah_tables.py",
    "rev12_atallah_verification.py",
    "rev12_atallah_transparency.py",
    "rev12_atallah_bincontrol.py",
    "rev12_atallah_sweep.py",
    "rev12_kernel_cost_curve.py",
    "rev12_finalize_manifest.py",
    "run_r12_overnight.sh",
    "run_r12_sweep_after.sh",
]

RESULT_JSON = [
    "r12_atallah_campaign.json",
    "r12_atallah_campaign_designB.json",
    "r12_atallah_descriptives.json",
    "r12_atallah_descriptives_designB.json",
    "r12_atallah_verification.json",
    "r12_atallah_transparency.json",
    "r12_atallah_bincontrol.json",
    "r12_atallah_sweep.json",
    "r12_kernel_cost_curve.json",
]

TABLES = [
    "r12_atallah_combined_table.tex",
    "r12_atallah_benchmark_table.tex",
    "r12_atallah_benchmark_table_designB.tex",
    "r12_atallah_verification_table.tex",
    "r12_atallah_matching_table_A.tex",
    "r12_atallah_matching_table_B.tex",
    "r12_atallah_bincontrol_table.tex",
]

LOGS = [
    "r12_atallah_phase0_verification.txt",
]

TREES = [
    ("design_A_atallah", METRICS / "r12_cases" / "atallah",
     METRICS / "r12_raw" / "atallah"),
    ("design_B_atallah", METRICS / "r12_cases" / "atallah_designB",
     METRICS / "r12_raw" / "atallah_designB"),
    ("atallah_tolerance_sweep", METRICS / "r12_cases" / "atallah_sweep",
     METRICS / "r12_raw" / "atallah_sweep"),
    ("atallah_bin_control", METRICS / "r12_cases" / "atallah_bincontrol",
     METRICS / "r12_raw" / "atallah_bincontrol"),
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
    """Per-file hashes for sidecars, rolled-up digest for the raw arrays.

    Each Atallah sidecar carries its own tolerance-selection record
    (``atallah_tol_accel_m_s2``), the full binned degree table, and the realized
    degree histogram, so hashing the sidecars pins those inputs as well.
    """
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
        "schema": "r12_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat()
                       .replace("+00:00", "Z"),
        "scope": ("R12 extension: faithful implementation and verification of "
                  "the Atallah (2022) analytical radial-adaptive rule, its "
                  "benchmark on both 64-orbit scrambled-Sobol populations, the "
                  "per-orbit matching record, the bin-resolution control, the "
                  "tolerance sweep, and the measured kernel cost curve"),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": ("unchanged from R10 and R11; every R12 trajectory sidecar "
                     "records this tag and the same kernel and gravity-file "
                     "hashes"),
        },
        "reused_evidence": (
            "the same-tolerance truth and critical-degree trajectories are "
            "reused from the R11 design-A and design-B convergence trees and "
            "are indexed in r11_final_experiment_manifest.json"),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "verification_logs": index_files(LOGS, METRICS),
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
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "verification_logs")
               for k, v in payload[sec].items() if v.get("missing")]
    missing += [f"tree:{k}" for k, v in payload["trajectory_trees"].items()
                if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
