"""SHA-256 integrity manifest for R37: the enlarged forced-variational panel.

R37 is a propagating campaign: each orbit is an augmented integration of the
forced variational system along that orbit's archived reference. What it does
not do is write trajectory records of its own --- the augmented solve returns
its prediction and nothing is spooled to a cases or raw tree --- so R37 claims
no trajectory record and the partition the supplement states over the
propagating campaigns is unchanged by it.

Three things are indexed here that a reader would otherwise have to take on
trust:

  * the registration and its scoring amendment, both hashed. The level chain,
    the wall clock and the commitment to report whatever the clock left
    incomplete were fixed before any new orbit was solved; the amendment
    records why sign agreement is scored on resolved comparisons only, and was
    written after the first undecided orbit appeared.
  * the self-check. Before any new orbit was accepted the run recomputed one
    archived orbit and required its prediction to reproduce the archived value.
    The result is lifted out of the record and into the manifest, so a later
    edit to the calibration or the reference trajectories surfaces here.
  * the archived eight-orbit panel, which R37 reads and carries over but must
    never overwrite. It stays indexed under R14; this manifest lists it as a
    reused input and additionally records the digest it was read at, so the
    two can be compared without opening the R14 manifest.

Usage:  python rev37_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
FIGURES = ROOT / "figures"

OUT = METRICS / "r37_final_experiment_manifest.json"
ARCHIVED = "r14_variational_budget.json"

SCOPE = (
    "R37: extension of the forced-variational panel of Section 7.2 from the "
    "archived eight orbits toward both full designs, along the same "
    "perilune-stratified rule refined into a nested level chain "
    "(8/14/26/50/80/112/128). Each orbit is one augmented integration of the "
    "forced system along that orbit's archived reference, at beta = 1. The run "
    "was bounded by a wall clock fixed in advance; the highest completed level "
    "is the panel and the orbits solved beyond it are reported separately."
)

REGISTRATION = ["r37_preregistration.json", "r37_scoring_amendment.json"]

SCRIPTS = ["rev37_variational_extend.py", "rev37_score.py", "rev37_table.py",
           "rev37_finalize_manifest.py", "rev14_variational_budget.py",
           "rev13_variational_check.py", "rev14_budget_pareto.py",
           "rev10_sobol_confirmatory.py", "rev12_atallah.py",
           "rev34_instrument_ladder.py", "make_figures_r19.py",
           "paper_style.py"]

# Read, never written. Each stays indexed under the manifest that produced it.
REUSED = [ARCHIVED, "r14_budget_pareto.json", "r14_preregistration.json",
          "r10_sobolA_baseline_truth_corrected.json", "r11_designB_rows.json",
          "r14_trajectory_A_beta_1.00.json", "r14_trajectory_B_beta_1.00.json"]

# The ladder table and the parity figure were last written by R42, which
# completed this campaign's level chain from 80 orbits to 128. A file is indexed
# under the campaign that last wrote it, so both are held in the R42 manifest and
# named there with the manifest they came from. R37's own record, verdict and
# panel table stay here.
HANDED_TO_R42 = ["r34_instrument_ladder_table.tex", "fig_variational_parity.pdf"]

RESULTS = ["r37_variational_extension.json", "r37_panel_verdict.json"]
TABLES = ["r37_panel_extension_table.tex"]
FIGS = []


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
    """What the wall clock actually left, lifted from the scored record."""
    v = METRICS / "r37_panel_verdict.json"
    r = METRICS / "r37_variational_extension.json"
    if not (v.exists() and r.exists()):
        return {"absent": True}
    vd = json.loads(v.read_text(encoding="utf-8"))
    rd = json.loads(r.read_text(encoding="utf-8"))
    solved = sum(1 for x in rd["rows"] if not x.get("reused_from_r14"))
    return {
        "levels_per_design": rd.get("levels_per_design"),
        "deadline_local": rd.get("deadline_local"),
        "orbits_in_record": vd["record"]["orbits"],
        "orbits_solved_here": solved,
        "orbits_reused_from_r14": vd["record"]["orbits"] - solved,
        "panel_orbits": vd["panel"]["orbits"],
        "panel_level_per_design": vd["levels"]["highest_complete_per_design"],
        "next_level_orbits": vd["levels"]["next_level_orbits"],
        "next_level_present": vd["levels"]["next_level_present"],
        "panel_resolved": vd["panel"]["resolved"],
        "panel_sign_agreement": vd["panel"]["sign_agreement"],
        "panel_sign_disagreement": vd["panel"]["sign_disagreement"],
        "record_resolved": vd["record"]["resolved"],
        "record_sign_agreement": vd["record"]["sign_agreement"],
        "self_check": rd.get("self_check"),
    }


def main() -> int:
    registration, a1 = index_files(REGISTRATION, METRICS)
    scripts, a2 = index_files(SCRIPTS, CODE)
    reused, a3 = index_files(REUSED, METRICS)
    results, a4 = index_files(RESULTS, METRICS)
    tables, a5 = index_files(TABLES, METRICS)
    figures, a6 = index_files(FIGS, FIGURES)
    absent = a1 + a2 + a3 + a4 + a5 + a6

    payload = {
        "schema": "r37_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered. The selection rule, the nested level chain, the "
            "wall clock and the rule that a truncated level is never reported "
            "as a panel were all fixed in r37_preregistration.json before any "
            "orbit outside the archived eight was solved. The scoring "
            "amendment is post hoc and says so."),
        "partition_note": (
            "this campaign claims no trajectory record. The augmented solve "
            "returns a prediction and spools nothing to a cases or raw tree, "
            "and the reference trajectories it integrates along stay indexed "
            "under R11 and R14. The trajectory partition is unchanged."),
        "archived_panel_note": (
            "the eight archived orbits are carried over from " + ARCHIVED +
            ", which is sealed under the R14 manifest and is read, never "
            "written. Its digest as read is recorded under reused_inputs; a "
            "mismatch against the R14 manifest would mean it had been "
            "modified."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "scripts": scripts,
        "reused_inputs": reused,
        "result_json": results,
        "generated_tables": tables,
        "generated_figures": figures,
        "handed_to_r42": HANDED_TO_R42,
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
        print(f"  panel {c['panel_orbits']} orbits "
              f"(level {c['panel_level_per_design']} per design), "
              f"{c['panel_sign_agreement']}/{c['panel_resolved']} resolved signs")
        print(f"  record {c['orbits_in_record']} orbits "
              f"({c['orbits_solved_here']} solved here, "
              f"{c['orbits_reused_from_r14']} reused), "
              f"{c['record_sign_agreement']}/{c['record_resolved']} resolved signs")
        print(f"  next level {c['next_level_orbits']} orbits: "
              f"{c['next_level_present']} present (truncated by the clock)")
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
