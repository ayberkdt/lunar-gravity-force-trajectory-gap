"""Create the machine-readable R10 overnight aggregate after the queue ends."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
OUTPUT = METRICS / "r10_aggregate_summary.json"


def read(name: str):
    path = METRICS / name
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def ratio_summary(rows: list[dict], key: str) -> dict:
    values = np.asarray([float(row[key]) for row in rows])
    return {
        "n": int(len(values)),
        "median": float(np.median(values)),
        "p10": float(np.percentile(values, 10)),
        "p90": float(np.percentile(values, 90)),
        "raw_schedule_wins": int(np.sum(values > 1.0)),
    }


def convergence_counts(rows: list[dict], comparator: str) -> dict:
    raw = resolved_schedule = resolved_fixed = unresolved = 0
    for row in rows:
        comparison = row["comparisons"][comparator]
        raw += comparison["rho_tight"] > 1.0
        if not comparison["resolved"]:
            unresolved += 1
        elif comparison["winner_if_resolved"] == "schedule_empirical":
            resolved_schedule += 1
        else:
            resolved_fixed += 1
    return {
        "n_completed_convergence_orbits": len(rows),
        "raw_schedule_wins": int(raw),
        "resolved_schedule_wins": resolved_schedule,
        "resolved_fixed_wins": resolved_fixed,
        "unresolved": unresolved,
    }


def main() -> int:
    protocol = base.protocol_payload()
    baseline = read("r10_sobolA_baseline.json")
    audit = read("r10_sobolA_truth_audit.json")
    corrected = read("r10_sobolA_baseline_truth_corrected.json")
    selection = read("r10_sobolA_convergence_selection.json")
    convergence = read("r10_sobolA_convergence.json")
    blend = read("r10_blend_lro_convergence.json")
    queue = read("r10_overnight_queue_manifest.json")
    if baseline is None or audit is None or corrected is None or selection is None:
        raise RuntimeError("required baseline/audit/selection artifacts are missing")
    corrected_rows = corrected["rows"]
    convergence_rows = [] if convergence is None else convergence.get("rows", [])
    payload = {
        "schema": "r10_aggregate_summary_v1",
        "created_utc": base.utc_now(),
        "protocol_sha256": protocol["protocol_sha256"],
        "artifact_hashes": {
            "baseline": base.file_hash(METRICS / "r10_sobolA_baseline.json"),
            "truth_audit": base.file_hash(METRICS / "r10_sobolA_truth_audit.json"),
            "truth_corrected_baseline": base.file_hash(
                METRICS / "r10_sobolA_baseline_truth_corrected.json"
            ),
            "selection": base.file_hash(
                METRICS / "r10_sobolA_convergence_selection.json"
            ),
        },
        "baseline": {
            "complete": baseline["complete"],
            "designed_orbits": 64,
            "truth_surviving": baseline["summary"]["truth_surviving"],
            "truth_impacts": baseline["summary"]["truth_impacts"],
            "numerical_failures": baseline["summary"]["orbit_numerical_failures"],
            "timing_comparable": baseline["timing_comparable"],
            "truth_corrected_rho_work": ratio_summary(corrected_rows, "rho_work"),
            "truth_corrected_rho_crit": ratio_summary(corrected_rows, "rho_crit"),
        },
        "truth_audit": {
            "complete": audit["complete"],
            "cases": len(audit["rows"]),
            "passes": audit["summary"]["passes"],
            "failures": audit["summary"]["failures"],
            "adopted_N900_indices": audit["summary"]["adopted_N900_indices"],
        },
        "selection": {
            "selected_count": selection["selected_count"],
            "selected_indices": selection["selected_indices"],
            "selection_sha256": selection["selection_sha256"],
        },
        "convergence": {
            "artifact_present": convergence is not None,
            "complete": False if convergence is None else convergence["complete"],
            "completed_orbits": len(convergence_rows),
            "selected_orbits": selection["selected_count"],
            "vs_fixed_work": convergence_counts(convergence_rows, "fixed_work"),
            "vs_fixed_critical": convergence_counts(
                convergence_rows, "fixed_critical"
            ),
        },
        "blend_lro_convergence": {
            "artifact_present": blend is not None,
            "complete": False if blend is None else blend["complete"],
            "record_count": 0 if blend is None else len(blend["records"]),
            "summary": None if blend is None else blend["summary"],
        },
        "overnight_queue": queue,
        "reporting_ready": bool(
            convergence is not None and convergence["complete"]
            and blend is not None and blend["complete"]
        ),
    }
    if convergence is not None:
        payload["artifact_hashes"]["convergence"] = base.file_hash(
            METRICS / "r10_sobolA_convergence.json"
        )
    if blend is not None:
        payload["artifact_hashes"]["blend_lro_convergence"] = base.file_hash(
            METRICS / "r10_blend_lro_convergence.json"
        )
    payload["aggregate_sha256"] = base.object_hash(payload)
    base.atomic_json(OUTPUT, payload)
    print(json.dumps({
        "output": str(OUTPUT),
        "reporting_ready": payload["reporting_ready"],
        "blend_complete": payload["blend_lro_convergence"]["complete"],
        "convergence_completed": payload["convergence"]["completed_orbits"],
        "convergence_selected": payload["convergence"]["selected_orbits"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
