"""SHA-256 integrity manifest for R39: the gradient-degree audit (O40).

R39 re-solves a registered stratified subset of the R37 panel with the
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
  * the R37 panel record, which R39 reads to draw its subset and must never
    overwrite. It stays indexed under R37; the digest it was read at is
    recorded here so the two can be compared without opening that manifest.

Usage:  python rev39_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"

OUT = METRICS / "r39_final_experiment_manifest.json"
RECORD = "r39_gradient_degree_panel.json"
PANEL = "r37_variational_extension.json"

SCOPE = (
    "R39: the gradient-degree audit of the enlarged forced-variational panel "
    "(O40). The forced solve of (O25) and (O38) evaluates the reference gravity "
    "gradient at degree 120; (O21) bounds what that neglects on eight orbits "
    "and the bound is loose at the lowest perilunes. This campaign removes the "
    "approximation instead of bounding it, on a registered stratified subset of "
    "the eighty-orbit panel, and reports only whether the panel's predictions "
    "would change side -- never a reference-degree prediction as the panel's."
)

REGISTRATION = ["r39_preregistration.json"]

SCRIPTS = ["rev39_gradient_degree_panel.py", "rev39_table.py",
           "rev39_finalize_manifest.py", "rev14_variational_budget.py",
           "rev13_variational_check.py", "rev14_budget_pareto.py",
           "rev10_sobol_confirmatory.py", "rev12_atallah.py"]

# Read, never written. Each stays indexed under the manifest that produced it.
REUSED = [PANEL, "r14_variational_budget.json", "r14_budget_pareto.json",
          "r21_gradient_sensitivity.json",
          "r10_sobolA_baseline_truth_corrected.json", "r11_designB_rows.json",
          "r14_trajectory_A_beta_1.00.json", "r14_trajectory_B_beta_1.00.json"]

RESULTS = [RECORD]
TABLES = ["r39_gradient_degree_table.tex"]


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
    """What the compute actually left, lifted from the run's own record."""
    p = METRICS / RECORD
    if not p.exists():
        return {"absent": True}
    r = json.loads(p.read_text(encoding="utf-8"))
    s = r["summary"]
    declared = len(json.loads(
        (METRICS / REGISTRATION[0]).read_text(encoding="utf-8")
    )["selection_rule"]["orbits"])
    hp = [c["hp_km"] for c in r["comparison"]]
    return {
        "orbits_declared": declared,
        "orbits_solved": s["orbits"],
        "orbits_unfinished": [f"{u['design']}{u['sobol_index']:03d}"
                              for u in r["unfinished"]],
        "band_status": s["band_status"],
        "bands_fully_complete": s["bands_fully_complete"],
        "audited_perilune_km": [min(hp), max(hp)] if hp else None,
        "resolved": s["resolved"],
        "side_changes_resolved": s["side_changes_resolved"],
        "side_changes_unresolved": s["side_changes_unresolved"],
        "abs_relative_change": s["abs_relative_change"],
        "abs_relative_change_by_band": s["abs_relative_change_by_band"],
        "self_check": r.get("self_check"),
        "gradient_degree": r.get("gradient_degree"),
        "panel_gradient_degree": r.get("panel_gradient_degree"),
    }


def main() -> int:
    registration, a1 = index_files(REGISTRATION, METRICS)
    scripts, a2 = index_files(SCRIPTS, CODE)
    reused, a3 = index_files(REUSED, METRICS)
    results, a4 = index_files(RESULTS, METRICS)
    tables, a5 = index_files(TABLES, METRICS)
    absent = a1 + a2 + a3 + a4 + a5

    payload = {
        "schema": "r39_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The four bands, the per-band selection rule, the "
            "submission order, the admissibility self-check and four outcomes "
            "were fixed in r39_preregistration.json before any orbit was solved "
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
    if not c.get("absent"):
        print(f"  {c['orbits_solved']} of {c['orbits_declared']} orbits; "
              f"unfinished {', '.join(c['orbits_unfinished']) or 'none'}")
        for b, v in c["band_status"].items():
            print(f"    {b:<10} {v['solved']}/{v['declared']}"
                  f"{'  COMPLETE' if v['complete'] else '  partial'}")
        print(f"  side changes {c['side_changes_resolved']} of "
              f"{c['resolved']} resolved, {c['side_changes_unresolved']} "
              f"among unresolved")
        lo, hi = c["audited_perilune_km"]
        print(f"  audited perilune {lo:.1f}-{hi:.1f} km")
        sc = c.get("self_check") or {}
        if sc:
            print(f"  self-check {sc.get('orbit')}: rel {sc.get('rel'):.1e} "
                  f"({'passed' if sc.get('passed') else 'FAILED'})")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
