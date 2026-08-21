"""R24-A: the oracle control at a third tolerance level.

R23-B applied the oracle control to the interior member and the interior member
survived it, but thinly: five resolved comparisons to three, with eight of the
sixteen left undecided. That tally is the weakest evidence the constructive
claim rests on, and the reason it is weak is visible in the record. In six of
the eight undecided cases the interior member has the *smaller* raw error, in
four of them by more than an order of magnitude; what blocks the verdict is the
interior member's own numerical envelope, which is one to nine meters where the
oracle's is centimeters. The envelopes, not the physics, are deciding.

R23-C already established that this is fixable and that fixing it does not
manufacture winners: refining tight-to-tighter envelopes into tighter-to-ultra
self-differences resolved nineteen of twenty-three borderline span comparisons
with zero verdict flips and the gap sign stable in eighty-two of eighty-three.
This applies the same refinement to the oracle panel.

What is reused and what is recomputed, following R23-C's convention exactly:
errors stay measured at the tighter level against the tighter truth, so the
numbers the paper already reports do not move. Only the envelopes are rebuilt,
from tighter-to-ultra self-differences instead of tight-to-tighter ones. A
comparison can therefore become resolved, but the error that decides it is the
same error as before.

The panel is frozen before the run: all sixteen comparisons in the R23-B
record, both comparators, no subsetting. Running only the undecided eight would
resolve them against an unexamined background -- the previously resolved eight
are what test whether the refinement flips anything.

Usage:
    python rev24_oracle_ultra.py preregister
    python rev24_oracle_ultra.py run --workers 11 --deadline-min 150
    python rev24_oracle_ultra.py summarize
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import scipy

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev18_span_sweep as r18
import rev23_oracle_vs_interior as r23b

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r24_cases" / "oracle_ultra"
RAW_ROOT = METRICS / "r24_raw" / "oracle_ultra"
SPAN_CASE = METRICS / "r23_cases" / "ultra_span"
SPAN_RAW = METRICS / "r23_raw" / "ultra_span"
OUTPUT = METRICS / "r24_oracle_ultra.json"
PREREG = METRICS / "r24_preregistration.json"

ULTRA = {"rtol": 3.0e-14,
         "atol": np.array([1.0e-7] * 3 + [1.0e-10] * 3),
         "atol_position_m": 1.0e-7, "atol_velocity_m_s": 1.0e-10}
BETA = 1.00
K = 0.50
K_TAG = "0.50"
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP


# --------------------------------------------------------------- paths
def paths(design: str, index: int, stem: str):
    return (CASE_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.json",
            RAW_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.npz")


def span_ultra(design: str, index: int, policy: str):
    """Where R23-C put its ultra runs, if it ran this orbit at all."""
    return (SPAN_CASE / design / f"sobolA_{index:03d}" / f"{policy}_ultra.json",
            SPAN_RAW / design / f"sobolA_{index:03d}" / f"{policy}_ultra.npz")


def _degree_fn(table: dict):
    tab = {float(a): int(b) for a, b in table.items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, r18.BIN_KM * math.floor(h_m / 1e3 / r18.BIN_KM)))
        return tab[hb]

    return degree_of


# ----------------------------------------------------------------- panel
def build_panel() -> list[dict]:
    """Every comparison in the R23-B record, ordered by design then index.
    Neither ordering key correlates with cost, so an early stop shortens the
    panel without tilting it."""
    rows = json.loads(r23b.OUTPUT.read_text(encoding="utf-8"))["rows"]
    panel = []
    for r in rows:
        panel.append({
            "design": r["design"], "index": int(r["sobol_index"]),
            "n_sat": int(r["n_sat"]), "n_oracle": int(r["n_oracle"]),
            "previous_resolved_vs_oracle": bool(r["resolved_vs_oracle"]),
            "previous_winner_vs_oracle": r["winner_vs_oracle"],
            "previous_resolved_vs_sat": bool(r["resolved_vs_sat"]),
            "previous_winner_vs_sat": r["winner_vs_sat"],
        })
    panel.sort(key=lambda d: (d["design"], d["index"]))
    return panel


# ------------------------------------------------------------- reuse gate
def _reusable(donor_sidecar: Path, donor_raw: Path, expect: dict) -> Path | None:
    """A donor run is the same object only if every field that can move the
    trajectory matches. Anything less than an exact match is refused, because a
    near-match would silently compare two different trajectories."""
    if not (donor_sidecar.exists() and donor_raw.exists()):
        return None
    try:
        rec = json.loads(donor_sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if rec.get("status") != "ok":
        return None
    cfg = rec.get("config", {})
    for key, want in expect.items():
        got = cfg.get(key)
        if isinstance(want, float):
            if got is None or float(got) != want:
                return None
        elif isinstance(want, dict):
            if {str(a): int(b) for a, b in (got or {}).items()} \
                    != {str(a): int(b) for a, b in want.items()}:
                return None
        elif got != want:
            return None
    return donor_raw


# ---------------------------------------------------------------- worker
def worker(task: dict) -> dict:
    design, index = task["design"], int(task["index"])
    t0 = time.time()
    try:
        r18_side, _ = r18.paths(design, BETA, index, K, "tight")
        cfg18 = json.loads(r18_side.read_text(encoding="utf-8"))["config"]
        adopted = int(cfg18["adopted_truth_degree"])
        y0 = np.asarray(cfg18["initial_state_si"], dtype=float)
        table_k = cfg18["degree_table"]

        jobs = {
            "truth": {"spec": {"kind": "fixed_truth", "degree": adopted},
                      "sel": (lambda t, h, n=adopted: n),
                      "donor": "truth"},
            "interior": {"spec": {"kind": "binned_span_member", "k": K,
                                  "degree_table": table_k},
                         "sel": _degree_fn(table_k),
                         "donor": "interior"},
            "sat": {"spec": {"kind": "fixed_budget_saturating",
                             "degree": int(task["n_sat"])},
                    "sel": (lambda t, h, n=int(task["n_sat"]): n),
                    "donor": None},
        }
        if int(task["n_oracle"]) != int(task["n_sat"]):
            jobs["oracle"] = {
                "spec": {"kind": "fixed_ladder_oracle",
                         "degree": int(task["n_oracle"])},
                "sel": (lambda t, h, n=int(task["n_oracle"]): n),
                "donor": None}

        model, args = r14._model(adopted)
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        reused, ran = [], []

        for policy, job in jobs.items():
            sidecar, raw = paths(design, index, policy)
            cfg = {
                "design": design, "sobol_index": index, "beta": BETA,
                "policy": policy, "level": "ultra",
                "adopted_truth_degree": adopted, "interior_k": K,
                "policy_spec": job["spec"],
                "initial_state_si": [float(v) for v in y0],
                "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
                "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
                "atol_kind": "vector", "rtol": ULTRA["rtol"],
                "atol_position_m": ULTRA["atol_position_m"],
                "atol_velocity_m_s": ULTRA["atol_velocity_m_s"],
                "timing_comparable": False,
                "purpose": ("third tolerance level for the oracle control on "
                            "the interior member at beta = 1"),
                "source": task.get("provenance", {})}
            digest = base.object_hash(cfg)

            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == digest \
                        and prev.get("status") == "ok":
                    continue

            # R23-C already propagated the truth and the interior member at
            # ultra for the orbits in its own panel. The integrator is
            # deterministic, so an exact field match means the same trajectory.
            donor_raw = None
            if job["donor"]:
                expect = {"policy": job["donor"], "level": "ultra",
                          "adopted_truth_degree": adopted,
                          "rtol": ULTRA["rtol"],
                          "atol_position_m": ULTRA["atol_position_m"],
                          "atol_velocity_m_s": ULTRA["atol_velocity_m_s"],
                          "max_step_s": float(MAX_STEP),
                          "duration_s": float(DURATION),
                          "output_step_s": float(OUTPUT_STEP)}
                if job["donor"] == "interior":
                    expect["policy_spec"] = None  # compared below instead
                donor_sidecar, dr = span_ultra(design, index, job["donor"])
                cand = _reusable(donor_sidecar, dr,
                                 {k: v for k, v in expect.items()
                                  if v is not None})
                if cand is not None and job["donor"] == "interior":
                    # the degree table has to be the same member, not merely
                    # the same k label
                    dcfg = json.loads(
                        donor_sidecar.read_text(encoding="utf-8"))["config"]
                    dtab = (dcfg.get("policy_spec") or {}).get("degree_table")
                    if {str(a): int(b) for a, b in (dtab or {}).items()} \
                            != {str(a): int(b) for a, b in table_k.items()}:
                        cand = None
                donor_raw = cand

            if donor_raw is not None:
                d = np.load(donor_raw)
                base.atomic_npz(raw, t_s=d["t_s"], state_si=d["state_si"])
                base.atomic_json(sidecar, {
                    "schema": "r24_oracle_ultra_v1",
                    "created_utc": base.utc_now(), "config": cfg,
                    "config_sha256": digest, "status": "ok",
                    "event": None, "telemetry": None,
                    "reused_from": str(donor_raw.relative_to(ROOT)),
                    "reuse_basis": ("identical policy, degrees, tolerances, "
                                    "max step, duration and output grid; the "
                                    "integrator is deterministic"),
                    "raw_path": str(raw.relative_to(ROOT)),
                    "raw_sha256": base.file_hash(raw),
                    "n_output_epochs": int(len(d["t_s"])),
                    "last_output_epoch_s": float(d["t_s"][-1])})
                reused.append(policy)
                continue

            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, job["sel"], args,
                ULTRA["rtol"], ULTRA["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"design": design, "index": index, "policy": policy,
                        "status": "numerical_failure", "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r24_oracle_ultra_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": digest, "status": "ok",
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
            ran.append(policy)

        return {"design": design, "index": index, "status": "ok",
                "reused": reused, "propagated": ran,
                "minutes": (time.time() - t0) / 60.0}
    except Exception as exc:                                    # noqa: BLE001
        return {"design": design, "index": index, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ------------------------------------------------------------ preregister
OUTCOMES = {
    "A_refinement_resolves_and_sign_holds": (
        "Envelopes shrink, previously undecided comparisons become resolved, "
        "and no previously resolved verdict flips. The oracle control is "
        "reported at its refined tally and the constructive claim keeps the "
        "oracle check with a firmer margin."),
    "B_a_previously_resolved_verdict_flips": (
        "Some comparison resolved at the tighter level reverses at the ultra "
        "level. The tight-to-tighter envelopes were then not conservative, "
        "the oracle comparison is reported as unstable under tolerance "
        "refinement, and the oracle check is withdrawn rather than restated "
        "at the new tally."),
    "C_envelopes_do_not_shrink": (
        "The tighter-to-ultra self-differences are comparable to the "
        "tight-to-tighter ones, so the interior member's error is not "
        "tolerance-limited and the undecided comparisons stay undecided. "
        "Reported as a limit of the measurement, not as evidence either way."),
    "D_refinement_resolves_against_the_interior_member": (
        "Comparisons resolve but the added verdicts favour the oracle. "
        "Reported as found; the oracle margin narrows or reverses and the "
        "constructive claim is narrowed accordingly."),
}


def preregister(args) -> int:
    panel = build_panel()
    prev_or = sum(1 for p in panel if p["previous_resolved_vs_oracle"])
    payload = {
        "schema": "r24_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": "R24-A: oracle control at the ultra tolerance level",
        "why": ("The R23-B oracle tally is five to three with eight of sixteen "
                "undecided, and in six of the eight undecided cases the "
                "interior member holds the smaller raw error while its own "
                "tight-to-tighter envelope blocks the verdict. R23-C showed "
                "the same refinement resolves such cases without flipping "
                "any."),
        "frozen_panel": {
            "rule": ("all sixteen comparisons in r23_oracle_vs_interior.json, "
                     "both comparators, no subsetting; ordered by design then "
                     "Sobol index"),
            "n_comparisons": len(panel),
            "previously_resolved_vs_oracle": prev_or,
            "previously_undecided_vs_oracle": len(panel) - prev_or,
            "orbits": [{"design": p["design"], "index": p["index"],
                        "n_sat": p["n_sat"], "n_oracle": p["n_oracle"]}
                       for p in panel],
        },
        "measurement_convention": {
            "errors": ("unchanged: the tighter-level error against the "
                       "tighter truth, as already recorded by R23-B, so no "
                       "number the paper reports moves"),
            "envelopes": ("rebuilt as the policy's tighter-to-ultra "
                          "self-difference plus the truth's"),
            "level": {"rtol": ULTRA["rtol"],
                      "atol_position_m": ULTRA["atol_position_m"],
                      "atol_velocity_m_s": ULTRA["atol_velocity_m_s"],
                      "max_step_s": MAX_STEP},
            "note": ("the level is R13's frozen ultra setting, unchanged, and "
                     "max step is held so the refinement isolates tolerance"),
        },
        "stability_check": (
            "verdict_flips counts only a resolved verdict reversing its "
            "winner. A comparison moving from undecided to resolved is not a "
            "flip; miscounting that was a reporting bug caught in R23-C and it "
            "is not repeated here."),
        "outcomes": OUTCOMES,
        "not_blind": (
            "The direction is guessable from the R23-B record: the six "
            "undecided cases where the interior member leads on raw error "
            "will tend to resolve in its favour if the envelopes shrink. The "
            "value of this registration is that it fixes what will be written "
            "for every outcome, including the ones that weaken the claim, "
            "before any fresh number exists."),
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    PREREG.parent.mkdir(parents=True, exist_ok=True)
    PREREG.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r24 prereg] {PREREG.name} "
          f"sha256={payload['preregistration_sha256'][:16]} "
          f"panel={len(panel)} comparisons, "
          f"{len(panel) - prev_or} currently undecided vs the oracle")
    return 0


# -------------------------------------------------------------------- run
def run(args) -> int:
    if not PREREG.exists():
        print("[r24] no preregistration on disk; run `preregister` first")
        return 1
    panel = build_panel()
    print(f"[r24] {len(panel)} orbits, workers={args.workers}, "
          f"deadline={args.deadline_min:.0f} min")
    t0 = time.time()
    done = fail = 0
    reused_total = propagated_total = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(worker, {**p, "provenance": {
            "driver": "rev24_oracle_ultra.py",
            "reuses": ("r23_oracle_vs_interior.json panel and degrees, "
                       "r23_cases/ultra_span ultra truths and interior "
                       "members, r18 span sidecars for states and tables"),
            "preregistration": "r24_preregistration.json"}}): p
            for p in panel}
        for fut in as_completed(futures):
            res = fut.result()
            if res["status"] == "ok":
                done += 1
                reused_total += len(res["reused"])
                propagated_total += len(res["propagated"])
                print(f"  [{res['design']}{res['index']:03d}] ok "
                      f"reused={res['reused']} ran={res['propagated']} "
                      f"{res['minutes']:.1f} min")
            else:
                fail += 1
                print(f"  [{res['design']}{res['index']:03d}] "
                      f"{res['status']}: {res.get('detail')}")
            if (time.time() - t0) / 60.0 > args.deadline_min:
                print("[r24] deadline reached; remaining orbits left for a "
                      "later invocation (every stage is resumable)")
                break
    print(f"[r24] {done} orbits ok, {fail} failed, "
          f"{reused_total} reused, {propagated_total} propagated, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace())


# -------------------------------------------------------------- summarize
def environment() -> dict:
    """R18 through R23 recorded no interpreter or library versions, which is
    why placing those campaigns in an environment needed the bytecode cache
    rather than the records. This campaign records them. Aggregation runs in
    the interpreter that propagated it, so these are that interpreter's."""
    return {"python": platform.python_version(),
            "python_full": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "captured": "at aggregation, in the interpreter that ran the "
                        "propagations"}


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def _err(a, b) -> float:
    return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]


def summarize(args) -> int:
    prev_rows = {(r["design"], int(r["sobol_index"])): r for r in
                 json.loads(r23b.OUTPUT.read_text(encoding="utf-8"))["rows"]}
    rows, missing = [], []

    for p in build_panel():
        design, index = p["design"], p["index"]
        prev = prev_rows[(design, index)]

        _, truth_tighter = r14.reuse_paths(design, index, "truth", "tighter")
        _, truth_ultra = paths(design, index, "truth")
        if not (truth_tighter.exists() and truth_ultra.exists()):
            missing.append(f"{design}{index:03d}:truth")
            continue
        truth_self = _err(_load(truth_tighter), _load(truth_ultra))

        # interior at tighter comes from R18, at ultra from this campaign
        _, int_tighter = r18.paths(design, BETA, index, K, "tighter")
        _, int_ultra = paths(design, index, "interior")
        if not (int_tighter.exists() and int_ultra.exists()):
            missing.append(f"{design}{index:03d}:interior")
            continue
        int_self = _err(_load(int_tighter), _load(int_ultra))

        comp = {}
        ok = True
        for kind in ("sat", "oracle"):
            degree = int(p[f"n_{kind}"])
            if kind == "oracle" and degree == int(p["n_sat"]):
                comp["oracle"] = dict(comp["sat"])
                comp["oracle"]["is_sat"] = True
                continue
            _, tighter = r23b.paths(design, index, kind, degree)
            _, ultra = paths(design, index, kind)
            if not (tighter.exists() and ultra.exists()):
                missing.append(f"{design}{index:03d}:{kind}")
                ok = False
                break
            self_diff = _err(_load(tighter), _load(ultra))
            comp[kind] = {"degree": degree,
                          "error_m": prev[f"{kind}_error_m"],
                          "self_difference_m": self_diff,
                          "envelope_m": self_diff + truth_self,
                          "is_sat": kind == "sat"}
        if not ok:
            continue

        int_side = {"error_m": prev["interior_error_m"],
                    "envelope_m": int_self + truth_self,
                    "self_difference_m": int_self}

        def verdict(a: dict, b: dict):
            diff = b["error_m"] - a["error_m"]
            thr = a["envelope_m"] + b["envelope_m"]
            if abs(diff) <= thr:
                return None, False
            return ("a" if diff > 0 else "b"), True

        w_or, r_or = verdict(int_side, comp["oracle"])
        w_sat, r_sat = verdict(int_side, comp["sat"])
        name = {"a": "interior", "b": None}

        def flip(prev_res, prev_win, now_res, now_win) -> bool:
            """Only a resolved verdict reversing its winner is a flip."""
            return bool(prev_res and now_res and prev_win and now_win
                        and prev_win != now_win)

        rows.append({
            "design": design, "sobol_index": index, "hp_km": prev["hp_km"],
            "n_sat": comp["sat"]["degree"], "n_oracle": comp["oracle"]["degree"],
            "interior_error_m": int_side["error_m"],
            "oracle_error_m": comp["oracle"]["error_m"],
            "sat_error_m": comp["sat"]["error_m"],
            "interior_envelope_before_m": prev["interior_envelope_m"],
            "interior_envelope_after_m": int_side["envelope_m"],
            "oracle_envelope_before_m": prev["oracle_envelope_m"],
            "oracle_envelope_after_m": comp["oracle"]["envelope_m"],
            "truth_self_difference_tighter_to_ultra_m": truth_self,
            "interior_self_difference_tighter_to_ultra_m": int_self,
            "resolved_vs_oracle_before": prev["resolved_vs_oracle"],
            "resolved_vs_oracle_after": r_or,
            "winner_vs_oracle_before": prev["winner_vs_oracle"],
            "winner_vs_oracle_after": ({"a": "interior",
                                        "b": "oracle"}.get(w_or)),
            "resolved_vs_sat_before": prev["resolved_vs_sat"],
            "resolved_vs_sat_after": r_sat,
            "winner_vs_sat_before": prev["winner_vs_sat"],
            "winner_vs_sat_after": ({"a": "interior", "b": "sat"}.get(w_sat)),
            "verdict_flip_vs_oracle": flip(
                prev["resolved_vs_oracle"], prev["winner_vs_oracle"], r_or,
                {"a": "interior", "b": "oracle"}.get(w_or)),
            "verdict_flip_vs_sat": flip(
                prev["resolved_vs_sat"], prev["winner_vs_sat"], r_sat,
                {"a": "interior", "b": "sat"}.get(w_sat)),
        })

    if not rows:
        print("[r24] nothing to summarize yet")
        return 1

    def tally(rs, suffix):
        res = [r for r in rs if r[f"resolved_vs_{suffix}_after"]]
        return {
            "comparisons": len(rs),
            "resolved_before": sum(1 for r in rs
                                   if r[f"resolved_vs_{suffix}_before"]),
            "resolved_after": len(res),
            "interior_wins": sum(1 for r in res
                                 if r[f"winner_vs_{suffix}_after"] == "interior"),
            "comparator_wins": sum(1 for r in res
                                   if r[f"winner_vs_{suffix}_after"] not in
                                   (None, "interior")),
            "verdict_flips": sum(1 for r in rs if r[f"verdict_flip_vs_{suffix}"]),
            "newly_resolved": sum(1 for r in rs
                                  if r[f"resolved_vs_{suffix}_after"]
                                  and not r[f"resolved_vs_{suffix}_before"]),
            "became_undecided": sum(1 for r in rs
                                    if r[f"resolved_vs_{suffix}_before"]
                                    and not r[f"resolved_vs_{suffix}_after"]),
        }

    shrink = [r["interior_envelope_before_m"] / r["interior_envelope_after_m"]
              for r in rows if r["interior_envelope_after_m"] > 0]
    payload = {
        "schema": "r24_oracle_ultra_v1", "created_utc": base.utc_now(),
        "environment": environment(),
        "beta": BETA, "interior_member_k": K_TAG, "level": "ultra",
        "preregistration": json.loads(
            PREREG.read_text(encoding="utf-8"))["preregistration_sha256"]
        if PREREG.exists() else None,
        "what_this_tests": (
            "whether the oracle control's undecided comparisons are undecided "
            "because of the interior member's numerical envelope rather than "
            "its error. Errors are the tighter-level errors R23-B already "
            "recorded and are not recomputed; only envelopes are rebuilt from "
            "tighter-to-ultra self-differences."),
        "panel_completeness": {
            "aggregated": len(rows),
            "missing": len(missing),
            "missing_detail": sorted(set(missing)),
            "note": ("the aggregate covers the completed comparisons only and "
                     "the incomplete ones are named; it is never computed over "
                     "a silently truncated panel"),
        },
        "envelope_shrink_factor": {
            "median": float(np.median(shrink)) if shrink else None,
            "min": float(np.min(shrink)) if shrink else None,
            "max": float(np.max(shrink)) if shrink else None,
        },
        "summary": {"vs_oracle": tally(rows, "oracle"),
                    "vs_sat": tally(rows, "sat")},
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    s = payload["summary"]
    print(f"[r24] written {OUTPUT.name}: {len(rows)} comparisons, "
          f"{len(missing)} missing")
    for tag in ("vs_oracle", "vs_sat"):
        t = s[tag]
        print(f"  {tag}: resolved {t['resolved_before']} -> "
              f"{t['resolved_after']} ({t['interior_wins']} interior, "
              f"{t['comparator_wins']} comparator), "
              f"flips={t['verdict_flips']}, new={t['newly_resolved']}")
    if shrink:
        print(f"  interior envelope shrank by median "
              f"{np.median(shrink):.1f}x (range {min(shrink):.1f}-"
              f"{max(shrink):.1f})")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preregister").set_defaults(fn=preregister)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=150.0)
    r.set_defaults(fn=run)
    sub.add_parser("summarize").set_defaults(fn=summarize)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
