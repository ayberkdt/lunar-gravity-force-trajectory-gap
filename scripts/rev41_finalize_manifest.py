"""SHA-256 integrity manifest for R41: the reference-degree control (O41).

Every truncation error in the paper is measured against an adopted reference
degree, 300 on most coverage-design orbits. R41 repeats the equal-budget
comparison of (O25) on a perilune-stratified subset of sixteen of them with the
reference raised to 600, and asks whether any verdict moves.

What the cap audit licenses is what makes the control cheap: at beta = 1 the
calibrated radial schedule reaches the reference on 9 of 64 orbits in design A
and 6 of 64 in design B, so on the rest raising the reference changes neither
policy. Only the reference is re-propagated. Both archived policy trajectories
are reused byte for byte.

Three things are indexed here that a reader would otherwise take on trust:

  * the registration, hashed. The eligibility test, the perilune-stratified
    selection, the reuse protocol and the outcomes were fixed before any
    reference was propagated at the raised degree.
  * the trajectories this campaign owns. Unlike R37 and R39, which return
    predictions and spool nothing, R41 propagates: thirty-two reference arrays,
    sixteen orbits at two tolerance levels, each with a sidecar carrying the
    configuration that determines it and the digest of the array.
  * the archived records it reads and must never overwrite, listed as reused
    inputs with the digests they were read at.

Usage:  python rev41_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"

OUT = METRICS / "r41_final_experiment_manifest.json"
RECORD = "r41_reference_degree_control.json"
PANEL = "r14_trajectory_A_beta_1.00.json"

SCOPE = (
    "R41: the reference-degree control (O41). Every truncation error in the "
    "paper is measured against an adopted reference degree, 300 on most "
    "coverage-design orbits. This repeats the equal-budget comparison of (O25) "
    "on a perilune-stratified subset of sixteen orbits with the reference "
    "raised to 600. The cap audit licenses reusing both policy trajectories "
    "unchanged, so only the reference is re-propagated, at both tolerance "
    "levels, and the truth-inclusive resolution rule is re-applied."
)

REGISTRATION = ["r41_preregistration.json"]

SCRIPTS = ["rev41_reference_degree_control.py",
           "rev41_sidecars.py",
           "rev41_finalize_manifest.py",
           "rev14_budget_trajectory.py", "rev12_atallah.py",
           "rev10_sobol_confirmatory.py"]

# Read, never written. Each stays indexed under the manifest that produced it.
REUSED = [PANEL, "r14_trajectory_B_beta_1.00.json",
          "r14_budget_pareto.json",
          "r10_sobolA_baseline_truth_corrected.json",
          "r11_designB_rows.json"]

RESULTS = [RECORD, "r41_trajectory_index.json"]
TABLES = []


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path):
    out, absent = {}, []
    for n in names:
        p = base / n
        if p.exists():
            out[n] = {"sha256": sha(p), "bytes": p.stat().st_size}
        else:
            out[n] = {"absent": True}
            absent.append(n)
    return out, absent


def _sidecars() -> dict:
    """The reference trajectories this campaign owns, by rolled-up digest.
    The arrays are excluded from the repository as regenerable; the sidecars
    that carry their digests are not."""
    idx = METRICS / "r41_trajectory_index.json"
    if not idx.exists():
        return {"absent": True}
    d = json.loads(idx.read_text(encoding="utf-8"))
    return {"trajectories": d["trajectories"],
            "rolled_up_raw_digest": d["rolled_up_raw_digest"],
            "index": "r41_trajectory_index.json",
            "arrays_shipped": False,
            "note": d["note"]}


def completion() -> dict:
    """What the control found, lifted from its own record."""
    p = METRICS / RECORD
    if not p.exists():
        return {"absent": True}
    r = json.loads(p.read_text(encoding="utf-8"))
    s = r["summary"]
    moved = []
    for row in r["rows"]:
        a = (row.get("archived") or {}).get("rho_budget")
        n = (row.get("comparison") or {}).get("rho_budget")
        if a and n and (n / a > 1.02 or n / a < 0.98):
            moved.append({"orbit": f"{row['design']}{row['index']:03d}",
                          "rho_at_300": a, "rho_at_600": n,
                          "factor": n / a})
    return {
        "orbits": s["orbits"],
        "verdict_changes": s["verdict_changes"],
        "decidability_changes": s["decidability_changes"],
        "resolved_now": s["resolved_now"],
        "base_reference": r["base_reference"],
        "new_reference": r["new_reference"],
        "ratios_moving_more_than_2_percent": moved,
        "reuse_note": r["protocol_note"],
        "registration_sha256": r["registration_sha256"],
        "trajectory_sidecars": _sidecars(),
    }


def main() -> int:
    registration, a1 = index_files(REGISTRATION, METRICS)
    scripts, a2 = index_files(SCRIPTS, CODE)
    reused, a3 = index_files(REUSED, METRICS)
    results, a4 = index_files(RESULTS, METRICS)
    tables, a5 = index_files(TABLES, METRICS)
    absent = a1 + a2 + a3 + a4 + a5

    payload = {
        "schema": "r41_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The eligibility test, the perilune-stratified "
            "selection of eight orbits per design, the reuse protocol and three "
            "outcomes were fixed in r41_preregistration.json before any "
            "reference was propagated at the raised degree. Eligibility is read "
            "from the archived record rather than chosen: an orbit qualifies if "
            "its adopted reference is 300 and its archived per-call mean "
            "squared degree is below 300^2, which is the cap test. "),
        "partition_note": (
            "R41 owns the thirty-two raised-reference trajectory arrays and "
            "their sidecars: sixteen orbits at two tolerance levels. It owns no "
            "policy trajectory, because both archived policy trajectories are "
            "reused unchanged. The raw arrays are not redistributed; their "
            "rolled-up digest and per-trajectory sidecars are indexed by this "
            "manifest. "),
        "panel_note": (
            "the subset is drawn from the archived beta = 1 records, "
            "r14_trajectory_A_beta_1.00.json and its design-B counterpart, "
            "which are sealed under the R14 manifest and are read, never "
            "written. Their digests as read are recorded under reused_inputs; a "
            "mismatch against the R14 manifest would mean they had been "
            "modified. "),
        "reporting_limit": (
            "the campaign's reported errors remain those measured against the "
            "adopted reference of 300. This control reports only whether "
            "raising the reference to 600 would move a verdict, and the "
            "registration forbids quoting a reference-600 error as a campaign "
            "number. "),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "scripts": scripts,
        "reused_inputs": reused,
        "result_json": results,
        "generated_tables": tables,
        "completion": completion(),
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    c = payload["completion"]
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
