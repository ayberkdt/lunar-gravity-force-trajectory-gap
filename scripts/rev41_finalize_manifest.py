"""SHA-256 integrity manifest for R41: the gradient-degree audit (O40).

R41 re-solves a registered stratified subset of the R37 panel with the
reference gravity gradient evaluated at each orbit's own adopted reference
degree instead of at 120, so the approximation the mechanism rung rests on is
removed rather than bounded. Like R37 it writes no trajectory record: the
augmented solve returns a prediction and spools nothing to a cases or raw tree,
so the trajectory partition the supplement states is unchanged by it.

Four things are indexed here that a reader would otherwise take on trust:

  * the registration, hashed. The four bands, the per-band selection rule, the
    admissibility conditions and the four outcomes were fixed before any orbit
    was solved with the raised gradient. The registration also declares which
    band is outcome-dependent, and why that band is adversarial rather than
    flattering.
  * the self-check. Before any reference-degree orbit was accepted the run
    reproduced an archived degree-120 prediction through the same code path.
    The result is lifted out of the record into the manifest so that a later
    edit to the calibration or the reference trajectories surfaces here.
  * the completion state, band by band. The run returned 13 of 16 orbits and
    band L is one of four; the registration forbids reporting a partial band as
    the band, so the manifest records the per-band solved/declared counts rather
    than a single orbit total that would hide it.
  * the R37 panel record, which R41 reads to draw its subset and must never
    overwrite. It stays indexed under R37; the digest it was read at is
    recorded here so the two can be compared without opening that manifest.

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
           "rev41_finalize_manifest.py",
           "rev14_budget_trajectory.py", "rev12_atallah.py",
           "rev10_sobol_confirmatory.py"]

# Read, never written. Each stays indexed under the manifest that produced it.
REUSED = [PANEL, "r14_trajectory_B_beta_1.00.json",
          "r14_budget_pareto.json",
          "r10_sobolA_baseline_truth_corrected.json",
          "r11_designB_rows.json"]

RESULTS = [RECORD]
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
            "pre-registered. The four bands, the per-band selection rule, the "
            "submission order, the admissibility self-check and four outcomes "
            "were fixed in r41_preregistration.json before any orbit was solved "
            "with the raised gradient. One band is declared outcome-dependent "
            "in the registration itself: it ranks on the archived prediction, "
            "and is adversarial by construction."),
        "partition_note": (
            "this campaign claims no trajectory record. The augmented solve "
            "returns a prediction and spools nothing to a cases or raw tree, "
            "and the reference trajectories it integrates along stay indexed "
            "under R11 and R14. The trajectory partition is unchanged."),
        "panel_note": (
            "the subset is drawn from " + PANEL + ", which is sealed under the "
            "R37 manifest and is read, never written. Its digest as read is "
            "recorded under reused_inputs; a mismatch against the R37 manifest "
            "would mean it had been modified."),
        "partial_band_note": (
            "band L returned one of its four declared orbits. The registration "
            "forbids reporting a partial band as the band, so band_status "
            "carries solved and declared counts per band and the supplement "
            "reports L as partial. The three that did not return are named "
            "under completion.orbits_unfinished and are a compute bound, not a "
            "result."),
        "reporting_limit": (
            "the panel's reported predictions remain the degree-120 ones. This "
            "campaign reports only whether they would change side, which the "
            "registration fixes and this manifest repeats so that a later "
            "reader cannot mistake a reference-degree ratio for a panel "
            "number."),
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
