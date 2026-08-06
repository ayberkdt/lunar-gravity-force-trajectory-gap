"""R39: the gradient-degree audit, carried from the archived eight to the panel.

The forced variational reconstruction evaluates the reference gravity gradient
G at degree 120 while the reference trajectory and the differenced accelerations
are at the orbit's adopted truth degree. R21 measures what that approximation
neglects, and R33 was written to remove it rather than bound it -- but both
address the original eight-orbit panel, and R37 enlarged the panel to eighty.
The sign reconstruction is therefore supported at eighty orbits while the
approximation underneath it is audited at eight, which is the asymmetry this
campaign closes.

It closes it on a stratified subset rather than on all eighty, because the
audit is a doubling-to-2.5x cost on top of a solve that already costs a
propagation, and the quantity it interrogates -- whether the neglected part of
G can move a prediction across ratio one -- is a property of the regime, not of
the orbit count. Sixteen orbits chosen to span the regimes bound it; eighty
would only bound it again.

  selection   four perilune bands over the eighty-orbit panel. Bands L, M and H
              are hp < 50, 50 <= hp < 100 and hp >= 100 km, and within each the
              four orbits are taken at linspace(0, len-1, 4) of the
              perilune-sorted band -- the same rule the panel itself was drawn
              with. Band X takes the four the audit has most to lose on: two
              orbits whose propagated comparison is unresolved and whose
              predicted ratio is closest to one, where a small perturbation of
              the prediction is most able to move its sign, and the two orbits
              of most extreme predicted ratio, where the amplitude claim is
              largest. Bands L, M and H are drawn on geometry alone; band X is
              drawn on the archived prediction and is declared as such.

  gradient    vb.GRADIENT_DEGREE is raised past every adopted truth degree, so
              rev14's own min(GRADIENT_DEGREE, adopted) resolves the gradient to
              the reference degree on every orbit. The patch is applied at
              module import, because ProcessPoolExecutor on Windows re-imports
              this module in each worker to unpickle the callable below, and a
              patch applied in main() would leave the children computing the
              degree-120 gradient while the parent reported otherwise. Each
              record carries the module value the child actually saw, and a run
              in which any child disagrees with the parent is not written.

  clock       orbits are submitted in a fixed order and the run stops
              submitting at a deadline fixed before it started. What the clock
              leaves incomplete is reported as incomplete; a partial band is
              never reported as the band.

The registration is metrics/r39_preregistration.json. This script writes only
metrics/r39_gradient_degree_panel.json. The R37 and R14 records are read and
never touched.

Usage:
    python rev39_gradient_degree_panel.py preregister
    python rev39_gradient_degree_panel.py run --workers 8 --deadline-h 20
"""

from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_variational_budget as vb

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

PANEL = METRICS / "r37_variational_extension.json"
PREREG = METRICS / "r39_preregistration.json"
OUTPUT = METRICS / "r39_gradient_degree_panel.json"

# Above every adopted truth degree in the panel (300, 600, 900), so that
# rev14's own cap resolves the gradient to the reference degree on every orbit.
# Module level on purpose: see the "gradient" note in the docstring.
REFERENCE_SENTINEL = 10_000
ARCHIVED_GRADIENT_DEGREE = 120
vb.GRADIENT_DEGREE = REFERENCE_SENTINEL

# The panel is the highest complete level of R37: level 40 per design and the
# eight orbits carried over from R14, which the record stores without a level.
PANEL_LEVEL = 40
BANDS = (("L", 0.0, 50.0), ("M", 50.0, 100.0), ("H", 100.0, math.inf))
PER_BAND = 4
SELF_CHECK_ORBIT = "A002"


def worker(task: dict) -> dict:
    """rev14's worker, reached through this module so the child applies the
    gradient-degree patch above before running it.

    A task may carry ``gradient_degree_override``; the self-check uses it to
    reproduce an archived degree-120 prediction through this same path.
    """
    override = task.get("gradient_degree_override")
    if override is not None:
        vb.GRADIENT_DEGREE = int(override)
    rec = vb.worker(task)
    rec["gradient_degree_module"] = vb.GRADIENT_DEGREE
    if override is not None:
        vb.GRADIENT_DEGREE = REFERENCE_SENTINEL
    return rec


def panel_rows() -> list[dict]:
    """The eighty orbits R37 reports as its panel."""
    rec = json.loads(PANEL.read_text(encoding="utf-8"))
    if int(rec["gradient_degree"]) != ARCHIVED_GRADIENT_DEGREE:
        raise SystemExit(f"[abort] {PANEL.name} is not the degree-120 record")
    rows = [r for r in rec["rows"]
            if r.get("status") == "complete" and r.get("level", 0) <= PANEL_LEVEL]
    if len(rows) != 2 * PANEL_LEVEL:
        raise SystemExit(f"[abort] panel is {len(rows)} orbits, expected "
                         f"{2 * PANEL_LEVEL}")
    return rows


def log_ratio(row: dict) -> float:
    p = row.get("predicted_ratio_fixed_over_atallah")
    return abs(math.log10(p)) if p and p > 0 else math.inf


def key(row: dict) -> tuple:
    return (row["design"], int(row["sobol_index"]))


def select(rows: list[dict]) -> list[dict]:
    """The declared stratified subset. Deterministic; no result reordering
    beyond the two keys band X is declared to use."""
    chosen, taken = [], set()
    for name, lo, hi in BANDS:
        band = sorted((r for r in rows if lo <= r["hp_km"] < hi),
                      key=lambda r: r["hp_km"])
        idx = [int(i) for i in np.linspace(0, len(band) - 1, PER_BAND).round()]
        for i in idx:
            r = band[i]
            if key(r) in taken:
                continue
            taken.add(key(r))
            chosen.append(dict(r, band=name, band_size=len(band),
                               band_position=int(i)))
    rest = [r for r in rows if key(r) not in taken]

    fragile = sorted((r for r in rest if not r["measured"]["resolved"]),
                     key=lambda r: (log_ratio(r), r["hp_km"]))[:2]
    for r in fragile:
        taken.add(key(r))
        chosen.append(dict(r, band="X-fragile", band_size=len(rest),
                           band_position=None))
    extreme = sorted((r for r in rest if key(r) not in taken),
                     key=lambda r: (-log_ratio(r), r["hp_km"]))[:2]
    for r in extreme:
        taken.add(key(r))
        chosen.append(dict(r, band="X-extreme", band_size=len(rest),
                           band_position=None))
    return chosen


def make_task(entry: dict, rows_by_design: dict) -> dict:
    row = rows_by_design[entry["design"]][int(entry["sobol_index"])]
    return {"design": entry["design"], "row": row, "spec": entry["spec"]}


def load_specs() -> tuple[dict, dict]:
    """The beta = 1 schedule specs and design rows the panel was solved with."""
    beta_key = f"beta_{vb.BETA:.2f}"
    pareto = json.loads(vb.PARETO.read_text(encoding="utf-8"))
    specs = {}
    for d in ("A", "B"):
        for r in pareto["designs"][d]["rows"]:
            specs[(d, int(r["sobol_index"]))] = r["budgets"][beta_key]
    rows_by_design = {
        d: {int(r["sobol_index"]): r
            for r in json.loads(vb.ROWS[d].read_text(encoding="utf-8"))["rows"]}
        for d in ("A", "B")}
    return specs, rows_by_design


def preregister(args) -> int:
    rows = panel_rows()
    chosen = select(rows)
    payload = {
        "schema": "r39_preregistration_v1",
        "campaign": ("R39 -- gradient-degree audit of the enlarged "
                     "forced-variational panel"),
        "question": (
            "The mechanism rung is now an eighty-orbit result (R37), but the "
            "audit of the approximation it rests on -- the reference gravity "
            "gradient evaluated at degree 120 rather than at the orbit's "
            "adopted reference degree -- is an eight-orbit result (R21), and "
            "it is a bound rather than a measurement. Does removing the "
            "approximation change the sign of any prediction on a subset "
            "chosen to span the panel's regimes?"),
        "written_before": (
            "any orbit was solved with the gradient at the reference degree. "
            "No result of this audit had been seen when this file was written; "
            "the subset below is drawn from the archived degree-120 record "
            "alone."),
        "not_blind": (
            "the archived degree-120 predictions are already reported, and "
            "R21's bound is known: below 0.62 on every panel orbit with "
            "perilune at or above 50 km, and not small on the two 31-km "
            "orbits. What is not known is what the reference-degree gradient "
            "does to any prediction."),
        "selection_rule": {
            "statement": (
                "Four perilune bands over the eighty-orbit R37 panel. Bands L, "
                "M and H are hp < 50, 50 <= hp < 100 and hp >= 100 km; within "
                "each, four orbits are taken at linspace(0, len-1, 4).round() "
                "of the perilune-sorted band, which is the rule the panel "
                "itself was drawn with. Band X takes two orbits whose "
                "propagated comparison is unresolved and whose predicted ratio "
                "is closest to one, and the two orbits of most extreme "
                "predicted ratio, each from what the first three bands leave."),
            "geometry_bands": [b[0] for b in BANDS],
            "per_band": PER_BAND,
            "band_X_is_outcome_dependent": (
                "band X ranks on the archived prediction, not on geometry. It "
                "is declared as outcome-dependent because it is: it "
                "deliberately selects the orbits where a perturbation of the "
                "prediction has the most to change. It cannot flatter the "
                "audit -- it is the adversarial half of the subset -- and its "
                "four orbits are reported separately from the twelve "
                "geometry-drawn ones so either half can be read alone."),
            "why_sixteen": (
                "the audit costs a factor of two to two and a half on top of a "
                "solve that already costs a propagation, and the two "
                "900-degree low-perilune orbits dominate the wall clock. "
                "Sixteen spanning orbits bound the neglected term in every "
                "regime the panel contains; eighty would bound it again at "
                "five times the cost."),
            "orbits": [
                {"design": r["design"], "sobol_index": int(r["sobol_index"]),
                 "band": r["band"], "hp_km": r["hp_km"],
                 "adopted_truth_degree": r["adopted_truth_degree"],
                 "archived_predicted_ratio":
                     r["predicted_ratio_fixed_over_atallah"],
                 "archived_measured_resolved": bool(r["measured"]["resolved"])}
                for r in chosen],
        },
        "stopping_rule": {
            "statement": (
                "orbits are submitted in the fixed band order L, M, H, "
                "X-fragile, X-extreme, perilune ascending within a band, and "
                "the run stops submitting at a deadline fixed before it "
                "started."),
            "reporting": (
                "orbits the clock leaves unfinished are reported as "
                "unfinished, with their band. A partially completed band is "
                "never reported as the band."),
        },
        "admissibility_self_check": {
            "statement": (
                "the audit imports rev14's worker verbatim rather than "
                "reimplementing it, and before any reference-degree orbit is "
                "accepted it recomputes one panel orbit at degree 120 through "
                "this same path and compares the predicted ratio with the "
                "archived value."),
            "recomputed_orbit": SELF_CHECK_ORBIT,
            "abort_threshold_rel": 0.001,
            "patch_check": (
                "every record carries the gradient-degree module value the "
                "child process actually saw. If any child disagrees with the "
                "parent, the run is not written -- the failure mode of a patch "
                "that does not survive re-import is the one that would be "
                "hardest to notice and worst to publish."),
        },
        "outcomes": {
            "E1": (
                "if no resolved orbit in the subset changes the side of ratio "
                "one, the degree-120 gradient is reported as audited across "
                "the panel's regimes rather than on the original eight, and "
                "the main text's limitation sentence is narrowed to the "
                "subset actually run -- not removed."),
            "E2": (
                "if a resolved orbit changes side, it is reported with its "
                "perilune, both ratios and its band, the panel's sign tally is "
                "qualified in the main text, and the subset is neither "
                "enlarged nor redrawn in search of a different outcome."),
            "E3": (
                "the relative change in predicted ratio is reported as a "
                "distribution over the subset whatever it is. If it is not "
                "small at the lowest perilunes, the manuscript's existing "
                "refusal to interpret amplitudes there stands and is restated "
                "at this subset's size."),
            "E4": (
                "an unresolved orbit that changes side is reported, and "
                "carries no verdict either way, exactly as it does in the "
                "panel's own scoring."),
        },
        "prohibited": [
            "dropping an orbit because it is slow, or because its result is "
            "inconvenient",
            "redrawing the bands, moving a band edge, or changing the linspace "
            "positions after any orbit has returned",
            "reporting the subset as though it were the eighty-orbit panel",
            "overwriting metrics/r37_variational_extension.json or "
            "metrics/r14_variational_budget.json, both sealed under their own "
            "manifests; this campaign writes its own record",
            "quoting a reference-degree prediction as the panel's prediction: "
            "the panel's reported numbers stay the archived degree-120 ones, "
            "and this campaign reports only whether they would change",
        ],
        "writes": [OUTPUT.name],
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    PREREG.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r39 prereg] {PREREG.name} "
          f"sha256={payload['preregistration_sha256'][:16]} "
          f"{len(chosen)} orbits")
    for r in chosen:
        print(f"  {r['band']:<9} {r['design']}{int(r['sobol_index']):03d} "
              f"hp={r['hp_km']:6.1f} N_ref={r['adopted_truth_degree']:3d} "
              f"ratio={r['predicted_ratio_fixed_over_atallah']:.4g} "
              f"resolved={int(r['measured']['resolved'])}")
    return 0


def compare(new: dict, old: dict) -> dict:
    n = new["predicted_ratio_fixed_over_atallah"]
    o = old["predicted_ratio_fixed_over_atallah"]
    return {
        "design": new["design"], "sobol_index": int(new["sobol_index"]),
        "band": old["band"], "hp_km": old["hp_km"],
        "adopted_truth_degree": new["adopted_truth_degree"],
        "gradient_degree_used": min(new["gradient_degree_module"],
                                    new["adopted_truth_degree"]),
        "measured_resolved": bool(old["measured"]["resolved"]),
        "measured_ratio": (old["measured"]["fixed_budget"]
                           / old["measured"]["atallah_budget"]),
        "ratio_gradient_120": o,
        "ratio_gradient_reference": n,
        "relative_change": (n / o - 1.0) if o else None,
        "side_changed": bool((n > 1.0) != (o > 1.0)),
    }


def band_status(comparison: list[dict], declared: list[dict]) -> dict:
    """Per band, solved against declared. A band is complete only when every
    orbit the registration put in it returned; the registration forbids
    reporting a partial band as the band."""
    want, got = {}, {}
    for o in declared:
        want[o["band"]] = want.get(o["band"], 0) + 1
    for c in comparison:
        got[c["band"]] = got.get(c["band"], 0) + 1
    return {b: {"solved": got.get(b, 0), "declared": want[b],
                "complete": got.get(b, 0) == want[b]}
            for b in sorted(want)}


def summarize(comparison: list[dict], declared: list[dict] = None) -> dict:
    res = [c for c in comparison if c["measured_resolved"]]
    ch = [abs(c["relative_change"]) for c in comparison
          if c["relative_change"] is not None]
    bs = band_status(comparison, declared) if declared else None
    return {
        "band_status": bs,
        "bands_fully_complete": ([b for b, v in bs.items() if v["complete"]]
                                 if bs else None),
        "orbits": len(comparison),
        "resolved": len(res),
        "side_changes_resolved": sum(1 for c in res if c["side_changed"]),
        "side_changes_unresolved": sum(1 for c in comparison
                                       if c["side_changed"]
                                       and not c["measured_resolved"]),
        "side_changed_detail": [c for c in comparison if c["side_changed"]],
        "abs_relative_change": ({
            "median": float(np.median(ch)), "max": float(np.max(ch)),
            "p90": float(np.percentile(ch, 90))} if ch else None),
        "abs_relative_change_by_band": {
            b: float(np.median([abs(c["relative_change"]) for c in comparison
                                if c["band"] == b
                                and c["relative_change"] is not None]))
            for b in sorted({c["band"] for c in comparison})},
        "bands_complete": sorted({c["band"] for c in comparison}),
    }


def run(args) -> int:
    if not PREREG.exists():
        print(f"[abort] {PREREG.name} missing; preregister first")
        return 2
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))

    rows = panel_rows()
    chosen = select(rows)
    declared = [(o["design"], int(o["sobol_index"]))
                for o in prereg["selection_rule"]["orbits"]]
    if [key(r) for r in chosen] != declared:
        print("[abort] the selection no longer reproduces the registered one")
        return 2

    specs, rows_by_design = load_specs()
    for r in chosen:
        r["spec"] = specs[key(r)]

    order = {b: i for i, b in enumerate(
        ["L", "M", "H", "X-fragile", "X-extreme"])}
    chosen.sort(key=lambda r: (order[r["band"]], r["hp_km"]))
    by_key = {key(r): r for r in chosen}

    deadline = time.time() + args.deadline_h * 3600.0
    print(f"[r39] {len(chosen)} orbits, gradient at the reference degree "
          f"(panel used {ARCHIVED_GRADIENT_DEGREE}); deadline in "
          f"{args.deadline_h:g} h", flush=True)

    # ---- admissibility self-check: one panel orbit reproduced at degree 120
    cd, ci = SELF_CHECK_ORBIT[0], int(SELF_CHECK_ORBIT[1:])
    ref = next(r for r in rows if key(r) == (cd, ci))
    task = make_task({"design": cd, "sobol_index": ci,
                      "spec": specs[(cd, ci)]}, rows_by_design)
    task["gradient_degree_override"] = ARCHIVED_GRADIENT_DEGREE
    print(f"[r39] self-check: recomputing {SELF_CHECK_ORBIT} at degree "
          f"{ARCHIVED_GRADIENT_DEGREE}", flush=True)
    t0 = time.time()
    got = worker(task)
    if got["status"] != "complete":
        print(f"[r39] ABORT: self-check orbit failed: {got.get('message')}")
        return 1
    old = ref["predicted_ratio_fixed_over_atallah"]
    new = got["predicted_ratio_fixed_over_atallah"]
    rel = abs(new - old) / abs(old)
    thr = float(prereg["admissibility_self_check"]["abort_threshold_rel"])
    print(f"[r39] self-check {SELF_CHECK_ORBIT}: archived {old:.10g}, "
          f"recomputed {new:.10g}, rel {rel:.2e} (threshold {thr:g}), "
          f"{time.time() - t0:.0f}s", flush=True)
    if rel > thr:
        print("[r39] ABORT: archived record and current source disagree; "
              "nothing written.")
        return 1

    meta = {
        "registration_sha256": prereg["preregistration_sha256"],
        "panel_record": PANEL.name,
        "panel_gradient_degree": ARCHIVED_GRADIENT_DEGREE,
        "gradient_degree": "adopted reference degree of each orbit",
        "deadline_h": args.deadline_h,
        "self_check": {"orbit": SELF_CHECK_ORBIT, "archived": old,
                       "recomputed": new, "rel": rel, "threshold": thr,
                       "passed": True},
        "scope_note": (
            "a stratified sixteen-orbit subset of the eighty-orbit panel, not "
            "the panel. The panel's reported predictions remain the archived "
            "degree-120 ones; this record says only whether raising the "
            "gradient to the reference degree would move them."),
    }

    def write(results, failures, unfinished):
        comparison = [compare(r, by_key[key(r)]) for r in results]
        comparison.sort(key=lambda c: (order[c["band"]], c["hp_km"]))
        payload = {
            "schema": "r39_gradient_degree_panel_v1",
            "created_utc": base.utc_now(), "beta": vb.BETA, **meta,
            "rows": sorted(results, key=lambda r: (r["design"],
                                                   r["sobol_index"])),
            "comparison": comparison,
            "failures": failures,
            "unfinished": unfinished,
            "summary": summarize(comparison),
            "source": base.provenance(),
        }
        base.atomic_json(OUTPUT, payload)
        return payload

    results, failures = [], []
    submitted = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs, it = {}, iter(chosen)
        for _ in range(min(args.workers, len(chosen))):
            e = next(it, None)
            if e is None:
                break
            futs[pool.submit(worker, make_task(e, rows_by_design))] = e
            submitted += 1
        while futs:
            for fut in as_completed(list(futs)):
                e = futs.pop(fut)
                rec = fut.result()
                if rec["status"] != "complete":
                    failures.append({"design": e["design"],
                                     "sobol_index": int(e["sobol_index"]),
                                     "band": e["band"],
                                     "message": rec.get("message")})
                    print(f"  !! {e['design']}{int(e['sobol_index']):03d} "
                          f"{rec.get('message')}", flush=True)
                elif rec["gradient_degree_module"] != REFERENCE_SENTINEL:
                    print(f"[r39] ABORT: child solved "
                          f"{e['design']}{int(e['sobol_index']):03d} with "
                          f"gradient module {rec['gradient_degree_module']}, "
                          f"parent has {REFERENCE_SENTINEL}; nothing written.")
                    return 1
                else:
                    rec["band"] = e["band"]
                    results.append(rec)
                    c = compare(rec, e)
                    print(f"  [{len(results)}/{len(chosen)}] {e['band']:<9} "
                          f"{e['design']}{int(e['sobol_index']):03d} "
                          f"hp={e['hp_km']:6.1f} G={c['gradient_degree_used']} "
                          f"{c['ratio_gradient_120']:.4g} -> "
                          f"{c['ratio_gradient_reference']:.4g} "
                          f"({c['relative_change']:+.2%})"
                          f"{'  SIDE CHANGE' if c['side_changed'] else ''} "
                          f"elapsed={(time.time() - t0) / 60:.1f}min",
                          flush=True)
                done = {key(r) for r in results} | {(f["design"],
                                                     f["sobol_index"])
                                                    for f in failures}
                write(results, failures,
                      [{"design": r["design"],
                        "sobol_index": int(r["sobol_index"]),
                        "band": r["band"]}
                       for r in chosen if key(r) not in done])
                if time.time() < deadline:
                    nxt = next(it, None)
                    if nxt is not None:
                        futs[pool.submit(worker,
                                         make_task(nxt, rows_by_design))] = nxt
                        submitted += 1
                break

    done = {key(r) for r in results} | {(f["design"], f["sobol_index"])
                                        for f in failures}
    payload = write(results, failures,
                    [{"design": r["design"],
                      "sobol_index": int(r["sobol_index"]), "band": r["band"]}
                     for r in chosen if key(r) not in done])
    s = payload["summary"]
    print(f"\n[r39] {s['orbits']} of {len(chosen)} orbits solved "
          f"({len(payload['unfinished'])} never submitted, "
          f"{len(failures)} failed)")
    print(f"[r39] side changes: {s['side_changes_resolved']} of "
          f"{s['resolved']} resolved, {s['side_changes_unresolved']} among "
          f"unresolved")
    if s["abs_relative_change"]:
        a = s["abs_relative_change"]
        print(f"[r39] |relative change| median {a['median']:.2e}, "
              f"p90 {a['p90']:.2e}, max {a['max']:.2e}")
        for b, v in s["abs_relative_change_by_band"].items():
            print(f"    {b:<9} median {v:.2e}")
    for c in s["side_changed_detail"]:
        print(f"[r39] SIDE CHANGE {c['design']}{c['sobol_index']:03d} "
              f"hp={c['hp_km']:.1f} band={c['band']} "
              f"resolved={int(c['measured_resolved'])} "
              f"{c['ratio_gradient_120']:.4g} -> "
              f"{c['ratio_gradient_reference']:.4g}")
    return 1 if failures else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preregister").set_defaults(fn=preregister)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--deadline-h", type=float, default=20.0)
    r.set_defaults(fn=run)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
