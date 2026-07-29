"""Extended N=600/N=900 truth-degree audit for the remaining sub-50 km Sobol A cases.

The frozen primary audit (``rev10_truth_audit.py``) covered the four
lowest-perilune truth survivors and every sub-50-km raw empirical-schedule win,
i.e. indices (4, 20, 27, 28, 36, 59). Reviewer feedback asked why the remaining
sub-50 km survivors were left on the N=600 truth degree when five of the six
audited cases failed the N=600/N=900 acceptance test. This extension audits
exactly those remaining cases -- the sub-50-km truth survivors that are neither
in the frozen primary set nor raw wins -- so that every orbit with
h_p < 50 km carries an audited adopted truth degree.

The derivation ``sub-50 km survivors minus frozen primary set`` yields
(11, 35, 43, 52); the script asserts this before propagating anything. As in the
primary audit, existing validated N=600 and policy trajectories are reused; only
the four N=900 truth trajectories are newly propagated, and every policy error
ratio is recomputed against the newly adopted audited truth degree.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import rev10_sobol_confirmatory as base
import rev10_truth_audit as primary


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
BASELINE_PATH = METRICS / "r10_sobolA_baseline.json"
OUTPUT_PATH = METRICS / "r10_sobolA_truth_audit_extended.json"
SMOKE_PATH = METRICS / "r10_sobolA_truth_audit_extended_smoke.json"
ACTIVE_PATH = METRICS / "r10_sobolA_truth_audit_extended_active.json"

EXTENDED_INDICES = (11, 35, 43, 52)
SCHEMA = "r10_sobolA_truth_audit_extended_v1"


def derive_extended_indices(baseline: dict) -> tuple[int, ...]:
    """Remaining sub-50 km truth survivors not covered by the frozen primary set."""
    survivors = [row for row in baseline["rows"]
                 if row.get("truth_survives_full_arc")]
    sub50 = {
        row["sobol_index"]
        for row in survivors
        if row["design_point"]["hp_km"] < 50.0
    }
    derived = tuple(sorted(sub50 - set(primary.AUDIT_INDICES)))
    if derived != EXTENDED_INDICES:
        raise RuntimeError(
            f"frozen extended indices {EXTENDED_INDICES} differ from derived {derived}"
        )
    return derived


def write_index(path: Path, rows: list[dict], complete: bool,
                timing_comparable: bool, started_utc: str,
                wall_start: int, cpu_start: int) -> None:
    payload = {
        "schema": SCHEMA,
        "extends": "r10_sobolA_truth_audit_v1",
        "primary_audit_indices": list(primary.AUDIT_INDICES),
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "design_sha256": base.load_design_a()["design_sha256"],
        "baseline_sha256": base.file_hash(BASELINE_PATH),
        "primary_audit_script_sha256": base.file_hash(
            Path(primary.__file__).resolve()
        ),
        "audit_script_sha256": base.file_hash(Path(__file__).resolve()),
        "audit_indices": list(EXTENDED_INDICES),
        "audit_rationale": (
            "sub-50 km truth survivors not covered by the frozen primary audit"
        ),
        "acceptance_rule": (
            "E600_900 < min(5 m, 0.05*smallest interpreted policy RMS vs N900)"
        ),
        "started_utc": started_utc,
        "updated_utc": base.utc_now(),
        "complete": complete,
        "timing_comparable": timing_comparable,
        "rows": rows,
        "summary": {
            "completed_cases": len(rows),
            "passes": sum(row.get("passes") is True for row in rows),
            "failures": sum(row.get("passes") is False for row in rows),
            "numerical_failures": sum(
                row.get("batch_status") == "numerical_failure" for row in rows
            ),
            "adopted_N900_indices": [
                row["sobol_index"] for row in rows
                if row.get("adopted_truth_degree") == 900
            ],
        },
        "current_session_wall_ns": time.perf_counter_ns() - wall_start,
        "current_session_cpu_ns": time.process_time_ns() - cpu_start,
    }
    if complete:
        payload["ended_utc"] = base.utc_now()
    base.atomic_json(path, payload)


def run(smoke: bool) -> int:
    baseline = primary.load_baseline()
    indices = derive_extended_indices(baseline)
    other = base.other_python_processes()
    if other:
        raise RuntimeError(
            "other Python processes are active; refusing timing-comparable audit: "
            + json.dumps(other)
        )
    timing_comparable = True
    output_path = SMOKE_PATH if smoke else OUTPUT_PATH
    run_kind = "truth_audit_extended_smoke" if smoke else "truth_audit_extended"
    duration = 2.0 * 3600.0 if smoke else base.DURATION_S
    run_indices = indices[:1] if smoke else indices
    existing: dict[int, dict] = {}
    if output_path.exists():
        previous = json.loads(output_path.read_text(encoding="utf-8"))
        existing = {row["sobol_index"]: row for row in previous.get("rows", [])}
    design = base.load_design_a()
    orbit_by_index = {orbit["sobol_index"]: orbit for orbit in design["orbits"]}
    model = base.load_model(900)
    args = base.kernel_args(model)
    base.warmup(model, args)
    started_utc = base.utc_now()
    wall_start = time.perf_counter_ns()
    cpu_start = time.process_time_ns()
    print(
        f"[{run_kind}] {len(run_indices)} cases, timing_comparable=True",
        flush=True,
    )
    for position, index in enumerate(run_indices, start=1):
        current_other = base.other_python_processes()
        if current_other:
            ordered = [existing[key] for key in sorted(existing)]
            write_index(output_path, ordered, False, False, started_utc,
                        wall_start, cpu_start)
            raise RuntimeError("another Python process appeared during audit")
        orbit = orbit_by_index[index]
        print(
            f"[{position}/{len(run_indices)}] {orbit['name']} "
            f"hp={orbit['hp_km']:.3f} km",
            flush=True,
        )
        try:
            row = primary.audit_orbit(
                run_kind, orbit, model, args, duration, timing_comparable
            )
            row["batch_status"] = "completed"
            print(
                f"  E600-900={row['N600_to_N900']['pos_rms_m']:.6f} m, "
                f"limit={row['acceptance_threshold_m']:.6f} m, "
                f"passes={row['passes']}",
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001 -- preserved as a batch failure
            row = {
                "name": orbit["name"],
                "sobol_index": index,
                "batch_status": "numerical_failure",
                "failure_message": f"{type(exc).__name__}: {exc}",
                "passes": None,
                "adopted_truth_degree": None,
            }
            print(f"  [failure preserved] {row['failure_message']}", flush=True)
        existing[index] = row
        ordered = [existing[key] for key in sorted(existing)]
        write_index(output_path, ordered, False, timing_comparable,
                    started_utc, wall_start, cpu_start)
    ordered = [existing[key] for key in sorted(existing)]
    write_index(output_path, ordered, True, timing_comparable,
                started_utc, wall_start, cpu_start)
    if ACTIVE_PATH.exists():
        ACTIVE_PATH.unlink()
    print(f"[{run_kind}] complete", flush=True)
    return 0


def status() -> int:
    if not OUTPUT_PATH.exists():
        print("no formal extended truth-audit index yet")
    else:
        payload = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
        print(json.dumps({
            "complete": payload.get("complete"),
            "completed_cases": len(payload.get("rows", [])),
            "planned_cases": len(payload.get("audit_indices", [])),
            "timing_comparable": payload.get("timing_comparable"),
            "summary": payload.get("summary"),
            "updated_utc": payload.get("updated_utc"),
        }, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("smoke", "run", "status"))
    args = parser.parse_args()
    if args.command == "smoke":
        return run(True)
    if args.command == "run":
        return run(False)
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
