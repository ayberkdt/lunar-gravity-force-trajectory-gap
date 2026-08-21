"""SHA-256 integrity manifest for R42: the completed forced-variational panel.

R42 is R37 finished rather than R37 repeated, and the manifest has to make that
checkable. Three things are indexed here that a reader would otherwise take on
trust:

  * what was carried. 102 of the 128 rows were computed by R37 and are
    reproduced unchanged; only 26 were solved here. The carried file is listed
    as a reused input with the digest it was read at, and the record itself
    stores that digest twice --- once as read and once re-read after the last
    orbit --- so a concurrent edit would show as a mismatch rather than as a
    silent mixture.

  * what licensed carrying it. Before the first new orbit the run recomputed
    R37's own admissibility orbit and required its prediction to reproduce the
    archived value under the same threshold. That result is lifted out of the
    record into the manifest, where a later edit to the calibration or the
    reference trajectories would surface.

  * how the run was scheduled. It shared its machine with the R30 stratum
    ladders and was stopped once and restarted at a lower worker count. Every
    process is listed with its worker count and deadline. None of it touches
    which orbits were attempted or in what order, which is why it can be
    recorded plainly instead of explained away.

  * what it rewrote. Four artifacts carry an earlier campaign's name because
    that campaign first produced them, and are regenerated here because the
    panel their third rung reads grew from 80 orbits to 128: the ladder record
    and its table, the ladder figure and the parity figure. A file is indexed
    under the campaign that last wrote it, so these are held here and leave the
    R35 and R37 indexes, which is recorded in both directions --- named here
    with the manifest they came from, and dropped there.

The R37 records themselves are read and never written: this campaign writes its
own record, its own verdict and its own table, and R37's stay sealed under the
R37 manifest.

Usage:  python rev42_finalize_manifest.py
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

OUT = METRICS / "r42_final_experiment_manifest.json"
CARRIED = "r37_variational_extension.json"

SCOPE = (
    "R42: completion of the R37 forced-variational level chain to its last two "
    "levels, 56 and 64 per design (112 and 128 orbits). The 26 orbits R37's "
    "wall clock never submitted are solved here with the worker R37 imported "
    "from R14; the 102 rows R37 had already written are carried forward "
    "unchanged. Level 64 completed, so the panel is both full designs and the "
    "chain has no next level."
)

REGISTRATION = ["r42_preregistration.json"]

SCRIPTS = ["rev42_variational_complete.py", "rev42_preregister.py",
           "rev42_score.py", "rev42_table.py", "rev42_instrument_ladder.py",
           "rev42_finalize_manifest.py", "make_figures_r42.py",
           "make_figures_r42_ladder.py",
           "rev37_variational_extend.py", "rev37_score.py", "rev37_table.py",
           "rev14_variational_budget.py", "rev13_variational_check.py",
           "rev14_budget_pareto.py", "rev10_sobol_confirmatory.py",
           "rev12_atallah.py", "rev34_instrument_ladder.py",
           "make_figures_r19.py", "make_figures_r34.py", "paper_style.py"]

# Read, never written. Each stays indexed under the manifest that produced it.
REUSED = [CARRIED, "r37_preregistration.json", "r37_scoring_amendment.json",
          "r37_final_experiment_manifest.json", "r14_variational_budget.json",
          "r14_budget_pareto.json", "r14_preregistration.json",
          "r10_sobolA_baseline_truth_corrected.json", "r11_designB_rows.json",
          "r14_trajectory_A_beta_1.00.json", "r14_trajectory_B_beta_1.00.json"]

RESULTS = ["r42_variational_completion.json", "r42_panel_verdict.json",
           "r34_instrument_ladder.json"]
TABLES = ["r42_panel_extension_table.tex", "r34_instrument_ladder_table.tex"]
FIGS = ["fig_variational_parity.pdf", "fig_instrument_ladder.pdf"]

# Four artifacts carry an earlier campaign's name and are regenerated here
# because the panel they read grew: the ladder record and its table (R34/R35),
# the ladder figure (R35) and the parity figure (R37). A file is indexed under
# the campaign that last wrote it, so these leave those manifests and are held
# here. Nothing else moves: R35's and R37's own records stay where they are.
REGENERATED_FROM = {
    "r34_instrument_ladder.json": "r35",
    "r34_instrument_ladder_table.tex": "r37",
    "fig_instrument_ladder.pdf": "r35",
    "fig_variational_parity.pdf": "r37",
}


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
    """What the chain actually reached, lifted from the scored record."""
    v = METRICS / "r42_panel_verdict.json"
    r = METRICS / "r42_variational_completion.json"
    if not (v.exists() and r.exists()):
        return {"absent": True}
    vd = json.loads(v.read_text(encoding="utf-8"))
    rd = json.loads(r.read_text(encoding="utf-8"))
    solved = sum(1 for x in rd["rows"] if not x.get("carried_from_r37"))
    return {
        "levels_per_design": rd.get("levels_per_design"),
        "orbits_in_record": vd["record"]["orbits"],
        "orbits_solved_here": solved,
        "orbits_carried_from_r37": vd["record"]["orbits"] - solved,
        "panel_orbits": vd["panel"]["orbits"],
        "panel_level_per_design": vd["levels"]["highest_complete_per_design"],
        "next_level_orbits": vd["levels"]["next_level_orbits"],
        "chain_exhausted": vd["levels"]["next_level_orbits"] is None,
        "panel_resolved": vd["panel"]["resolved"],
        "panel_unresolved": vd["panel"]["unresolved"],
        "panel_sign_agreement": vd["panel"]["sign_agreement"],
        "panel_sign_disagreement": vd["panel"]["sign_disagreement"],
        "panel_sign_agreement_raw": vd["panel"]["sign_agreement_raw"],
        "panel_sign_disagreement_raw": vd["panel"]["sign_disagreement_raw"],
        "panel_calibrated_within_band": vd["panel"]["calibrated_within_band"],
        "panel_calibration_median": vd["panel"]["calibration_median"],
        "predicted_favors_radial_all": vd["panel"]["predicted_favors_radial_all"],
        "carried": rd.get("carried"),
        "runs": rd.get("runs"),
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

    c = completion()
    carried_ok = None
    if not c.get("absent"):
        car = c.get("carried") or {}
        carried_ok = (car.get("sha256") == car.get("sha256_after_run")
                      == reused.get(CARRIED, {}).get("sha256"))

    payload = {
        "schema": "r42_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": SCOPE,
        "registration_status": (
            "pre-registered, and inheriting a registration rather than "
            "restating one. The selection rule, the nested level chain and the "
            "rule that a truncated level is never reported as a panel are "
            "R37's, pinned by digest in r42_preregistration.json, which was "
            "written before any of the 26 remaining orbits was submitted. That "
            "registration is not blind and says so: the outcome of every level "
            "up to 40 per design was already reported when it was written."),
        "continuation_note": (
            "the 102 carried rows were not recomputed. Two checks stand in "
            "their place: the carried file matched the digest sealed in the "
            "R37 manifest before the run and again after it, and R37's "
            "admissibility orbit B005 was recomputed and reproduced its "
            "archived prediction exactly."),
        "carried_record_unchanged": carried_ok,
        "partition_note": (
            "this campaign claims no trajectory record. The augmented solve "
            "returns a prediction and spools nothing to a cases or raw tree, "
            "and the reference trajectories it integrates along stay indexed "
            "under R11 and R14. The trajectory partition is unchanged."),
        "scoring_note": (
            "scored by rev37_score.py, imported verbatim and pointed at the "
            "R42 record, so the scoring amendment applies unchanged: sign "
            "agreement is counted among resolved comparisons only and the raw "
            "tally is archived beside it."),
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
        "regenerated_from_earlier_campaign": REGENERATED_FROM,
        "completion": c,
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    if not c.get("absent"):
        print(f"  panel {c['panel_orbits']} orbits "
              f"(level {c['panel_level_per_design']} per design), "
              f"chain exhausted: {c['chain_exhausted']}")
        print(f"  {c['panel_sign_agreement']}/{c['panel_resolved']} resolved "
              f"signs, {c['panel_unresolved']} undecided; raw tally "
              f"{c['panel_sign_agreement_raw']}/{c['panel_orbits']}")
        print(f"  record {c['orbits_in_record']} orbits "
              f"({c['orbits_solved_here']} solved here, "
              f"{c['orbits_carried_from_r37']} carried from R37)")
        print(f"  carried record unchanged: {carried_ok}")
        sc = c.get("self_check") or {}
        if sc:
            print(f"  self-check {sc.get('orbit')}: rel {sc.get('rel'):.1e} "
                  f"({'passed' if sc.get('passed') else 'FAILED'})")
        for r in (c.get("runs") or []):
            print(f"  run {r.get('started_local')}: {r.get('workers')} workers, "
                  f"deadline {r.get('deadline_local')}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
