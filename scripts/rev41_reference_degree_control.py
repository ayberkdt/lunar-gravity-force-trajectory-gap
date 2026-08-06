"""R41 (O41): does the beta = 1 reversal survive a higher reference degree?

Every truncation error in this paper is measured against an adopted reference
degree, 300 on most coverage-design orbits and 600 or 900 at the lowest
perilunes, while the coefficient product runs to 1800. The reference-degree
audit establishes that 300 is adequate by its own acceptance rule. It does not
answer the different question a referee is entitled to ask: would the reported
reversal still be there if the reference were higher?

R38 answered it for the operational elliptical population, where raising the
reference from 300 to 600 left the direction standing and moved only the
magnitude. This asks it where the confirmatory claim actually lives.

What makes the control cheap is the cap audit. At beta = 1 the calibrated
radial schedule reaches the adopted reference degree on 9 of 64 orbits in
design A and 6 of 64 in design B. On the remaining 113 the schedule never asks
for a degree the cap withholds, so raising the reference cannot change either
policy: the two policy trajectories are the same trajectories. Only the
reference against which their error is measured changes. The control therefore
re-propagates the reference and nothing else, and reuses the archived policy
trajectories byte for byte.

That is the whole design. It is not a re-run of the campaign, and it is not
comparable to R38, which had to recalibrate because there the cap did bind.

Usage:  python rev41_reference_degree_control.py preregister
        python rev41_reference_degree_control.py plan
"""

from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PREREG = METRICS / "r41_preregistration.json"
OUTROOT = METRICS / "r41_raw"
CASEROOT = METRICS / "r41_cases"

NEW_REFERENCE = 600
BASE_REFERENCE = 300
N_PER_DESIGN = 8


def eligible(design: str) -> list[dict]:
    """Orbits whose reference is the base degree and whose radial schedule
    never reaches the cap at beta = 1. Both conditions come from the archived
    record; neither is a choice made here."""
    p = METRICS / f"r14_trajectory_{design}_beta_1.00.json"
    out = []
    for r in json.loads(p.read_text(encoding="utf-8"))["rows"]:
        if int(r.get("adopted_truth_degree", 0)) != BASE_REFERENCE:
            continue
        c = r.get("cost") or {}
        # a schedule that touches the cap would have its mean degree pinned at
        # the reference; the archived per-call mean is the test
        mean_sq = c.get("per_call_mean_degree_sq_atallah")
        if mean_sq is None or mean_sq >= BASE_REFERENCE ** 2:
            continue
        cmp_ = r.get("comparison") or {}
        out.append({
            "design": design,
            "sobol_index": int(r["sobol_index"]),
            "hp_km": float(r["design_point"]["hp_km"]),
            "n_critical": int(r["n_critical"]),
            "archived_resolved": bool(cmp_.get("resolved")),
            "archived_winner": cmp_.get("resolved_winner"),
            "archived_rho_budget": cmp_.get("rho_budget"),
        })
    return out


def select() -> list[dict]:
    """Perilune-stratified, the rule the panels of this paper are drawn with:
    sort the eligible orbits by perilune and take evenly spaced ranks."""
    chosen = []
    for d in ("A", "B"):
        el = sorted(eligible(d), key=lambda r: r["hp_km"])
        idx = [int(i) for i in np.linspace(0, len(el) - 1, N_PER_DESIGN).round()]
        for i in idx:
            if el[i] not in chosen:
                chosen.append(el[i])
    return chosen


def preregister(_) -> int:
    chosen = select()
    res = [c for c in chosen if c["archived_resolved"]]
    payload = {
        "schema": "r41_preregistration_v1",
        "campaign": "R41 -- reference-degree control on the coverage designs",
        "question": (
            "The confirmatory reversal at beta = 1 is measured against an "
            "adopted reference degree of 300 on most coverage-design orbits. "
            "Would it still be there if the reference were 600?"),
        "written_before": (
            "any reference trajectory was propagated at the raised degree. No "
            "result of this control had been seen when this file was written."),
        "not_blind": (
            "the archived verdicts at reference 300 are reported, and R38 has "
            "already shown that raising the reference on the operational "
            "elliptical population left the direction standing. What is not "
            "known is what it does on the coverage designs, where the "
            "confirmatory claim lives."),
        "why_this_is_cheap": (
            "at beta = 1 the calibrated radial schedule reaches the cap on 9 "
            "of 64 orbits in design A and 6 of 64 in design B. On the rest it "
            "never asks for a degree the cap withholds, so raising the "
            "reference changes neither policy. The archived policy "
            "trajectories are reused unchanged and only the reference is "
            "re-propagated. This is why the control does not need the "
            "recalibration R38 needed."),
        "selection_rule": {
            "statement": (
                "eligible orbits are those whose adopted reference is 300 and "
                "whose archived per-call mean squared degree is below 300^2, "
                "which is the cap test. Eligible orbits are sorted by perilune "
                "and eight per design taken at evenly spaced ranks, the rule "
                "every panel in this paper is drawn with. No orbit is chosen "
                "or dropped on its verdict."),
            "base_reference": BASE_REFERENCE,
            "new_reference": NEW_REFERENCE,
            "per_design": N_PER_DESIGN,
            "orbits": chosen,
        },
        "protocol": {
            "reused_unchanged": (
                "both policy trajectories of each orbit, byte for byte from "
                "the archived beta = 1 records, together with their degree "
                "tables and calibrated tolerance"),
            "re_propagated": (
                "the reference trajectory only, at degree 600, at both vector "
                "tolerance levels so that the reference-inclusive envelope is "
                "rebuilt rather than carried over"),
            "scoring": (
                "errors are recomputed against the new reference and the "
                "truth-inclusive resolution rule is applied exactly as "
                "archived. A verdict is the resolved winner; an orbit that "
                "resolves at one reference and not at the other is reported as "
                "a change in decidability, not as a flip."),
        },
        "outcomes": {
            "E1": (
                "no resolved verdict changes winner. The reversal is reported "
                "as insensitive to the reference degree over 300 to 600 on "
                "this subset, and the manuscript says so in one sentence."),
            "E2": (
                "a resolved verdict changes winner. It is reported with its "
                "perilune, both error pairs and both envelopes, the "
                "confirmatory claim is qualified in the main text, and the "
                "subset is neither enlarged nor redrawn."),
            "E3": (
                "decidability moves without any winner changing. Reported as "
                "such, with the counts at both references."),
        },
        "prohibited": [
            "dropping an orbit because its result is inconvenient",
            "redrawing the subset or changing the rank rule after any orbit "
            "has returned",
            "quoting a reference-600 error as a campaign number: the campaign "
            "reports reference-300 errors and this control reports only "
            "whether they would change verdict",
            "overwriting the archived beta = 1 records",
        ],
        "writes": ["r41_reference_degree_control.json"],
        "archived_state_of_the_subset": {
            "orbits": len(chosen),
            "resolved": len(res),
            "resolved_for_fixed": sum(1 for c in res
                                      if c["archived_winner"] == "fixed"),
            "resolved_for_radial": sum(1 for c in res
                                       if c["archived_winner"] != "fixed"),
            "median_rho_budget": st.median(
                [c["archived_rho_budget"] for c in res
                 if c["archived_rho_budget"]]) if res else None,
        },
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    PREREG.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r41 prereg] {PREREG.name} "
          f"sha256={payload['preregistration_sha256'][:16]}")
    a = payload["archived_state_of_the_subset"]
    print(f"  {a['orbits']} orbits, {a['resolved']} resolved at reference 300: "
          f"{a['resolved_for_fixed']} fixed, {a['resolved_for_radial']} radial")
    for c in chosen:
        print(f"    {c['design']}{c['sobol_index']:03d} hp={c['hp_km']:6.1f} "
              f"N_crit={c['n_critical']:3d} "
              f"{'resolved:' + str(c['archived_winner']) if c['archived_resolved'] else 'unresolved'}")
    return 0


def plan(_) -> int:
    """Cost, from the cost curve and the archived call counts. Nothing is run."""
    curve = json.loads((METRICS / "r12_kernel_cost_curve.json").read_text(
        encoding="utf-8"))
    f = curve["fit_quadratic_full"]
    def t(n):
        return f["a_N2"] * n * n + f["b_N"] * n + f["c_0"]
    chosen = select()
    rhs = []
    for d in ("A", "B"):
        p = METRICS / f"r14_trajectory_{d}_beta_1.00.json"
        by = {int(r["sobol_index"]): r for r in
              json.loads(p.read_text(encoding="utf-8"))["rows"]}
        for c in chosen:
            if c["design"] == d:
                rhs.append((by[c["sobol_index"]].get("cost") or {})
                           .get("rhs_fixed", 0))
    med = st.median([r for r in rhs if r])
    one = med * t(NEW_REFERENCE)
    print(f"[r41 plan] {len(chosen)} orbits, reference {BASE_REFERENCE} "
          f"-> {NEW_REFERENCE}")
    print(f"  per-call at N={NEW_REFERENCE}: {t(NEW_REFERENCE)*1e6:.0f} us; "
          f"median archived RHS: {med:.0f}")
    print(f"  one reference trajectory: ~{one/60:.0f} min")
    print(f"  two tolerance levels, tighter costing about twice: "
          f"~{3*one/60:.0f} min per orbit")
    print(f"  {len(chosen)} orbits serial: ~{len(chosen)*3*one/3600:.1f} h; "
          f"on 8 workers: ~{len(chosen)*3*one/3600/8:.1f} h")
    return 0


def worker(task: dict) -> dict:
    """One orbit: propagate the raised reference at both tolerance levels and
    rescore the archived policy trajectories against it.

    Nothing is reimplemented. The propagator, the tolerance levels, the arc,
    the output grid, the maximum step and the error metric are rev14's, reached
    through the module rather than copied.
    """
    import numpy as _np
    import rev14_budget_trajectory as bt

    row, design = task["row"], task["design"]
    index = int(row["sobol_index"])
    try:
        y0 = _np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = bt._model(NEW_REFERENCE)
        deg = lambda t, h: NEW_REFERENCE                          # noqa: E731

        new_ref = {}
        for lv in bt.LEVELS:
            t, y, st, ev, fail, tel = bt._propagate(model, args, y0, deg, lv)
            if st == "numerical_failure":
                return {"index": index, "design": design,
                        "status": "numerical_failure",
                        "where": f"reference_{NEW_REFERENCE}/{lv}",
                        "message": fail}
            new_ref[lv] = (t, y)
            out = OUTROOT / f"{design}_ref{NEW_REFERENCE}" / f"sobolA_{index:03d}"
            out.mkdir(parents=True, exist_ok=True)
            npz = out / f"reference_{lv}.npz"
            base.atomic_npz(npz, t_s=t, state_si=y)
            # a raw array without a sidecar is not auditable; every other
            # campaign writes one and so does this
            cfg = {
                "campaign": "R41",
                "purpose": ("reference trajectory at the raised degree, for "
                            "the reference-degree control (O41)"),
                "design": design, "sobol_index": index,
                "reference_degree": NEW_REFERENCE,
                "base_reference_degree": BASE_REFERENCE,
                "tolerance_level": lv,
                "rtol": bt.LEVELS[lv]["rtol"],
                "atol_position_m": bt.LEVELS[lv]["atol_position_m"],
                "atol_velocity_m_s": bt.LEVELS[lv]["atol_velocity_m_s"],
                "max_step_s": bt.MAX_STEP, "duration_s": bt.DURATION,
                "output_step_s": bt.OUTPUT_STEP,
                "initial_state_si": list(map(float, y0)),
                "frame": ("inertial, Moon rotating uniformly at its sidereal "
                          "rate about the polar axis, gravity only"),
                "integrator": "DOP853",
                "propagator_source": ("rev14_budget_trajectory._propagate, "
                                      "called through this module"),
            }
            side = CASEROOT / npz.relative_to(OUTROOT).with_suffix(".json")
            side.parent.mkdir(parents=True, exist_ok=True)
            base.atomic_json(side, {
                "schema": "r41_reference_trajectory_v1",
                "created_utc": base.utc_now(),
                "config": cfg, "config_sha256": base.object_hash(cfg),
                "status": st, "event": ev, "telemetry": tel,
                "raw_path": str(npz.relative_to(ROOT)).replace("\\", "/"),
                "raw_sha256": base.file_hash(npz),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1]),
            })

        # the archived policies, byte for byte
        spec_reuse = bool(task["reuse_fixed_critical"])
        pol = {}
        pol["atallah_budget"] = {
            lv: bt._load(bt.paths(design, 1.0, index, "atallah_budget", lv)[1])
            for lv in bt.LEVELS}
        key = "fixed_critical" if spec_reuse else "fixed_budget"
        getter = bt.reuse_paths if spec_reuse else (
            lambda d, i, p, l: bt.paths(d, 1.0, i, p, l))
        pol["fixed_budget"] = {lv: bt._load(getter(design, index, key, lv)[1])
                               for lv in bt.LEVELS}

        def err(p, lv, ref):
            return base.common_error(p[lv][0], p[lv][1], ref[lv][0], ref[lv][1])

        ref_self = base.common_error(
            new_ref["tight"][0], new_ref["tight"][1],
            new_ref["tighter"][0], new_ref["tighter"][1])["pos_rms_m"]

        res = {"index": index, "design": design, "status": "complete",
               "reference_degree": NEW_REFERENCE,
               "reference_self_difference_rms_m": ref_self, "policies": {}}
        for name, p in pol.items():
            sd = base.common_error(p["tight"][0], p["tight"][1],
                                   p["tighter"][0], p["tighter"][1])["pos_rms_m"]
            res["policies"][name] = {
                "error_tight_m": err(p, "tight", new_ref)["pos_rms_m"],
                "self_difference_rms_m": sd,
                "truth_inclusive_envelope_m": sd + ref_self}
        a = res["policies"]["atallah_budget"]
        f = res["policies"]["fixed_budget"]
        diff = abs(a["error_tight_m"] - f["error_tight_m"])
        thr = a["truth_inclusive_envelope_m"] + f["truth_inclusive_envelope_m"]
        res["comparison"] = {
            "atallah_error_m": a["error_tight_m"],
            "fixed_error_m": f["error_tight_m"],
            "rho_budget": (f["error_tight_m"] / a["error_tight_m"]
                           if a["error_tight_m"] > 0 else None),
            "absolute_error_difference_m": diff,
            "resolution_threshold_m": thr,
            "resolution_margin": (diff / thr) if thr > 0 else None,
            "resolved": bool(diff > thr),
            "resolved_winner": (None if diff <= thr else
                                ("atallah" if a["error_tight_m"]
                                 < f["error_tight_m"] else "fixed")),
        }
        return res
    except Exception as exc:                                   # noqa: BLE001
        import traceback
        return {"index": index, "design": design, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def run(args) -> int:
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time
    import rev14_budget_trajectory as bt

    if not PREREG.exists():
        print(f"[abort] {PREREG.name} missing; preregister first")
        return 2
    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    declared = [(o["design"], int(o["sobol_index"]))
                for o in prereg["selection_rule"]["orbits"]]
    if [(c["design"], c["sobol_index"]) for c in select()] != declared:
        print("[abort] the selection no longer reproduces the registered one")
        return 2

    rows = {}
    for d in ("A", "B"):
        p = METRICS / f"r14_trajectory_{d}_beta_1.00.json"
        rows[d] = {int(r["sobol_index"]): r
                   for r in json.loads(p.read_text(encoding="utf-8"))["rows"]}
    src = {d: {int(r["sobol_index"]): r for r in
               json.loads(bt.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]}
           for d in ("A", "B")}

    tasks = []
    for dsg, idx in declared:
        arch = rows[dsg][idx]
        tasks.append({"design": dsg, "row": src[dsg][idx],
                      "reuse_fixed_critical": bool(arch.get(
                          "reuse_fixed_critical", True)),
                      "archived": arch.get("comparison") or {}})

    print(f"[r41] {len(tasks)} orbits, reference {BASE_REFERENCE} -> "
          f"{NEW_REFERENCE}, {args.workers} workers", flush=True)
    results, fails, t0 = [], [], time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for fut in as_completed(futs):
            t = futs[fut]
            r = fut.result()
            if r["status"] != "complete":
                fails.append(r)
                print(f"  !! {r['design']}{r['index']:03d} {r.get('message')}",
                      flush=True)
                continue
            arch = t["archived"]
            r["archived"] = {k: arch.get(k) for k in
                             ("resolved", "resolved_winner", "rho_budget",
                              "resolution_margin")}
            n = r["comparison"]
            r["verdict_change"] = (bool(arch.get("resolved"))
                                   and n["resolved"]
                                   and arch.get("resolved_winner")
                                   != n["resolved_winner"])
            r["decidability_change"] = (bool(arch.get("resolved"))
                                        != n["resolved"])
            results.append(r)
            print(f"  {r['design']}{r['index']:03d}  "
                  f"rho {arch.get('rho_budget'):.3g} -> "
                  f"{n['rho_budget']:.3g}   "
                  f"{arch.get('resolved_winner') or 'undecided'} -> "
                  f"{n['resolved_winner'] or 'undecided'}"
                  f"{'   VERDICT CHANGE' if r['verdict_change'] else ''}"
                  f"   [{len(results)}/{len(tasks)}] "
                  f"{(time.time()-t0)/60:.0f}min", flush=True)
            write(results, fails, prereg)
    write(results, fails, prereg)
    vc = sum(1 for r in results if r["verdict_change"])
    dc = sum(1 for r in results if r["decidability_change"])
    print(f"\n[r41] {len(results)} solved, {len(fails)} failed")
    print(f"[r41] verdict changes: {vc}   decidability changes: {dc}")
    return 1 if fails else 0


def write(results, fails, prereg):
    base.atomic_json(METRICS / "r41_reference_degree_control.json", {
        "schema": "r41_reference_degree_control_v1",
        "created_utc": base.utc_now(),
        "registration_sha256": prereg["preregistration_sha256"],
        "base_reference": BASE_REFERENCE, "new_reference": NEW_REFERENCE,
        "protocol_note": ("policy trajectories reused unchanged; only the "
                          "reference was re-propagated, at both tolerance "
                          "levels, and the truth-inclusive resolution rule "
                          "re-applied"),
        "rows": sorted(results, key=lambda r: (r["design"], r["index"])),
        "failures": fails,
        "summary": {
            "orbits": len(results),
            "verdict_changes": sum(1 for r in results if r["verdict_change"]),
            "decidability_changes": sum(1 for r in results
                                        if r["decidability_change"]),
            "resolved_now": sum(1 for r in results
                                if r["comparison"]["resolved"]),
        },
        "source": base.provenance(),
    })


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preregister").set_defaults(fn=preregister)
    sub.add_parser("plan").set_defaults(fn=plan)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=8)
    r.set_defaults(fn=run)
    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    raise SystemExit(main())
