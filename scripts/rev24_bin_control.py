"""R24-B: the bin-resolution control at the budget-calibrated point.

The paper's schedules quantize degree onto 10-km altitude bins. R12 tested
whether that quantization is load-bearing by re-running the published rule at
the exact instantaneous radius on every right-hand-side call, but it did so at
the *accuracy-targeted* operating point. The council's Finding 8 asked for the
same control at the *budget-calibrated* point, where the paper's constructive
claim lives, and that run was never made: the response narrowed the control's
declared scope instead. This makes the run.

Only the binning changes. Every parameter of the interior member is read from
its archived R18 sidecar and held: the span exponent k, the multiplicative
scale the work calibration settled on, the constant-degree endpoint, the
acceleration tolerance, the cap and the floor. The exact-radius variant is

    N(r) = clip( round( s * N_0^(1-k) * N_A(r)^k ), floor, adopted )

with N_A evaluated at the current radius rather than at the floor of the
enclosing 10-km bin. Nothing is recalibrated, because recalibrating would let
the budget move and the control would no longer isolate quantization.

Panel: the ten orbits of the R19 beta = 1 panel whose archived interior member
runs the widest degree span. Quantization can only matter where the degree
actually varies, so this is where a binning artifact would be largest -- which
makes a null result here stronger than a null result on a random sample, and is
the reason the rule is stated this way rather than as a random draw.

Usage:
    python rev24_bin_control.py preregister
    python rev24_bin_control.py run --workers 10 --deadline-min 120
    python rev24_bin_control.py summarize
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from rev24_oracle_ultra import environment

import rev10_sobol_confirmatory as base
import rev12_atallah as at
import rev14_budget_trajectory as r14
import rev18_span_sweep as r18
import rev19_equal_total_work as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r24_cases" / "bin_control"
RAW_ROOT = METRICS / "r24_raw" / "bin_control"
OUTPUT = METRICS / "r24_bin_control.json"
PREREG = METRICS / "r24_bin_control_preregistration.json"

BETA = 1.00
K = 0.50
K_TAG = "0.50"
N_PANEL = 10
LEVELS = ("tight", "tighter")


def paths(design: str, index: int, level: str):
    return (CASE_ROOT / design / f"sobolA_{index:03d}" / f"exact_{level}.json",
            RAW_ROOT / design / f"sobolA_{index:03d}" / f"exact_{level}.npz")


def _archived_interior_cfg(design: str, index: int) -> dict:
    sidecar, _ = r18.paths(design, BETA, index, K, "tight")
    return json.loads(sidecar.read_text(encoding="utf-8"))["config"]


def build_panel() -> list[dict]:
    """The widest-span *resolved* comparisons.

    Two filters, and the first one matters more than it looks. A comparison
    that the archive leaves undecided has no verdict for quantization to
    change, so asking whether the verdict moves is unanswerable there; the
    panel is drawn from the sixty resolved comparisons only. Among those the
    ordering is by the archived interior degree span, because quantization can
    only act where the degree varies, and span is a property of the schedule
    rather than of the comparison's outcome -- ordering by margin would be
    selecting on the thing being measured.
    """
    cand = []
    for design in ("A", "B"):
        rec = json.loads(r19.out_path(design, BETA).read_text(encoding="utf-8"))
        for row in rec["rows"]:
            index = int(row["sobol_index"])
            if not row.get("resolved"):
                continue
            try:
                cfg = _archived_interior_cfg(design, index)
            except (OSError, KeyError, json.JSONDecodeError):
                continue
            span = cfg.get("degree_span")
            if span is None or row.get("work_matched_degree") is None:
                continue
            cand.append({
                "design": design, "index": index,
                "degree_span": float(span),
                "work_matched_degree": int(row["work_matched_degree"]),
                "archived_interior_error_m": row.get("interior_error_m"),
                "archived_fixed_error_m": row.get("work_matched_error_m"),
                "archived_resolved": bool(row.get("resolved")),
                "archived_winner": row.get("winner"),
                "comparator_source": row.get("comparator_source"),
            })
    cand.sort(key=lambda d: (-d["degree_span"], d["design"], d["index"]))
    panel = cand[:N_PANEL]
    panel.sort(key=lambda d: (d["design"], d["index"]))
    return panel


def exact_degree_fn(cfg: dict, model, g):
    """The archived interior member with the binning removed."""
    k = float(cfg["k"])
    scale = float(cfg["scale"])
    n0 = float(cfg["constant_degree_endpoint"])
    cap = int(cfg["adopted_truth_degree"])
    tol = float(cfg["atallah_tol_accel_m_s2"])
    base_fn = at.atallah_degree_fn(model, g, tol, floor=r18.FLOOR, cap=cap)

    def degree_of(t, h_m):
        na = base_fn(t, h_m)
        v = scale * (n0 ** (1.0 - k)) * (float(na) ** k)
        return int(min(cap, max(r18.FLOOR, round(v))))

    return degree_of


def worker(task: dict) -> dict:
    design, index = task["design"], int(task["index"])
    t0 = time.time()
    try:
        cfg18 = _archived_interior_cfg(design, index)
        adopted = int(cfg18["adopted_truth_degree"])
        y0 = np.asarray(cfg18["initial_state_si"], dtype=float)
        model, args = r14._model(adopted)
        g = r14._g(adopted)
        sel = exact_degree_fn(cfg18, model, g)
        grid = np.arange(0.0, r14.DURATION + 0.5 * r14.OUTPUT_STEP,
                         r14.OUTPUT_STEP)
        done = []

        for level in LEVELS:
            sidecar, raw = paths(design, index, level)
            cfg = {
                "design": design, "sobol_index": index, "beta": BETA,
                "policy": "interior_member_exact_radius", "interior_k": K,
                "level": level,
                "policy_spec": {
                    "kind": "binned_span_member_without_binning",
                    "k": float(cfg18["k"]),
                    "scale": float(cfg18["scale"]),
                    "constant_degree_endpoint":
                        int(cfg18["constant_degree_endpoint"]),
                    "atallah_tol_accel_m_s2":
                        float(cfg18["atallah_tol_accel_m_s2"]),
                    "floor": int(r18.FLOOR), "cap": adopted,
                    "evaluated_at": ("exact instantaneous radius on every "
                                     "right-hand-side call"),
                },
                "archived_degree_table": cfg18["degree_table"],
                "archived_degree_span": cfg18.get("degree_span"),
                "adopted_truth_degree": adopted,
                "initial_state_si": [float(v) for v in y0],
                "duration_s": r14.DURATION,
                "output_step_s": r14.OUTPUT_STEP,
                "integrator": "InstrumentedDOP853",
                "max_step_s": r14.MAX_STEP, "atol_kind": "vector",
                "timing_comparable": False,
                "purpose": ("bin-resolution control at the budget-calibrated "
                            "point: the archived interior member with its "
                            "quantization removed and nothing else changed"),
                "source": task.get("provenance", {})}
            digest = base.object_hash(cfg)

            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == digest \
                        and prev.get("status") == "ok":
                    continue

            tol = r14.LEVELS[level]
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, r14.DURATION, grid, sel, args,
                tol["rtol"], tol["atol"], max_step=r14.MAX_STEP)
            if st == "numerical_failure":
                return {"design": design, "index": index, "level": level,
                        "status": "numerical_failure", "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r24_bin_control_v1", "created_utc": base.utc_now(),
                "config": cfg, "config_sha256": digest, "status": "ok",
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
            done.append(level)

        return {"design": design, "index": index, "status": "ok",
                "propagated": done, "minutes": (time.time() - t0) / 60.0}
    except Exception as exc:                                    # noqa: BLE001
        return {"design": design, "index": index, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


OUTCOMES = {
    "A_binning_is_not_load_bearing": (
        "Every verdict is unchanged and the exact-radius variant's realized "
        "quadratic work is within a few percent of the binned one. The control "
        "passes at the budget-calibrated point and the narrowed scope "
        "declared for the R12 control can be dropped."),
    "B_verdicts_hold_but_the_work_moves": (
        "Verdicts are unchanged while the exact-radius variant burns "
        "materially different work, meaning quantization was flattering or "
        "penalizing the schedule's budget. The work difference is reported and "
        "the comparison is restated on realized work."),
    "C_a_verdict_changes": (
        "At least one comparison changes its resolved winner. The binned "
        "comparison is then a quantization artifact at the budget-calibrated "
        "point; this is reported and the constructive claim is narrowed to the "
        "bin convention it was measured under."),
    "D_unresolvable": (
        "The exact-radius variant's envelope is large enough that the "
        "comparisons go undecided. Reported as a limit of the control, not as "
        "support either way."),
}


def preregister(args) -> int:
    panel = build_panel()
    payload = {
        "schema": "r24_bin_control_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R24-B: bin-resolution control at the budget-calibrated "
                     "point (council Finding 8)"),
        "why": ("R12 ran this control at the accuracy-targeted point only. "
                "The paper's constructive claim lives at the budget-calibrated "
                "point, where the control was never run and the response "
                "narrowed the declared scope instead."),
        "what_is_held": (
            "k, the calibrated multiplicative scale, the constant-degree "
            "endpoint, the acceleration tolerance, the cap and the floor are "
            "all read from the archived R18 sidecar and held. The budget is "
            "not recalibrated, so the only difference is quantization."),
        "frozen_panel": {
            "rule": (f"the {N_PANEL} *resolved* comparisons of the R19 "
                     "beta = 1 panel with the widest archived interior degree "
                     "span, ties by design then Sobol index"),
            "rationale": ("quantization can only matter where the degree "
                          "varies, so the widest spans are where an artifact "
                          "would be largest and a null result is strongest; "
                          "span is a property of the schedule, not of the "
                          "comparison's outcome, so ordering by it is not "
                          "selecting on the measured quantity"),
            "revision_note": (
                "the first draft of this rule ranked the whole panel by span "
                "and was replaced before any propagation ran, because nine of "
                "the ten orbits it selected are undecided in the archive and "
                "an undecided comparison has no verdict for quantization to "
                "change. No fresh number existed when the rule was changed."),
            "orbits": [{"design": p["design"], "index": p["index"],
                        "degree_span": p["degree_span"],
                        "work_matched_degree": p["work_matched_degree"],
                        "archived_resolved": p["archived_resolved"],
                        "archived_winner": p["archived_winner"]}
                       for p in panel],
        },
        "measurement_convention": (
            "the exact-radius variant is propagated at both the tight and the "
            "tighter level; its error is the tighter-level error against the "
            "tighter truth and its envelope is its tight-to-tighter "
            "self-difference plus the truth's, which is the convention the "
            "R19 beta = 1 comparison already uses"),
        "outcomes": OUTCOMES,
        "not_blind": (
            "R12's control at the accuracy-targeted point found no verdict "
            "change, so outcome A is the expected one. The registration fixes "
            "what will be written if the budget-calibrated point behaves "
            "differently, which is the case that would matter."),
    }
    payload["preregistration_sha256"] = base.object_hash(payload)
    PREREG.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[r24b prereg] {PREREG.name} "
          f"sha256={payload['preregistration_sha256'][:16]} "
          f"panel={len(panel)} orbits, spans "
          f"{panel[0]['degree_span']:.1f}-{panel[-1]['degree_span']:.1f}")
    return 0


def run(args) -> int:
    if not PREREG.exists():
        print("[r24b] no preregistration on disk; run `preregister` first")
        return 1
    panel = build_panel()
    print(f"[r24b] {len(panel)} orbits x {len(LEVELS)} levels, "
          f"workers={args.workers}")
    t0 = time.time()
    done = fail = 0
    prov = {"driver": "rev24_bin_control.py",
            "reuses": ("r18 span sidecars for every held parameter, r19 "
                       "beta = 1 record for the comparator and the archived "
                       "verdict, r14 truths"),
            "preregistration": "r24_bin_control_preregistration.json"}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(worker, {**p, "provenance": prov}): p
                   for p in panel}
        for fut in as_completed(futures):
            res = fut.result()
            if res["status"] == "ok":
                done += 1
                print(f"  [{res['design']}{res['index']:03d}] ok "
                      f"ran={res['propagated']} {res['minutes']:.1f} min")
            else:
                fail += 1
                print(f"  [{res['design']}{res['index']:03d}] "
                      f"{res['status']}: {res.get('detail')}")
            if (time.time() - t0) / 60.0 > args.deadline_min:
                print("[r24b] deadline reached; resumable")
                break
    print(f"[r24b] {done} ok, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min")
    return summarize(argparse.Namespace())


def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def _err(a, b) -> float:
    return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]


def _work(sidecar: Path) -> float | None:
    try:
        tel = json.loads(sidecar.read_text(encoding="utf-8"))["telemetry"]
    except (OSError, KeyError, json.JSONDecodeError):
        return None
    if not tel:
        return None
    n = tel.get("n_rhs")
    m = tel.get("mean_degree_sq")
    return float(n) * float(m) if n and m else None


def summarize(args) -> int:
    rows, missing = [], []
    for p in build_panel():
        design, index = p["design"], p["index"]
        s_tight, r_tight = paths(design, index, "tight")
        s_tighter, r_tighter = paths(design, index, "tighter")
        if not (r_tight.exists() and r_tighter.exists()):
            missing.append(f"{design}{index:03d}")
            continue
        _, truth_tighter = r14.reuse_paths(design, index, "truth", "tighter")
        _, truth_tight = r14.reuse_paths(design, index, "truth", "tight")
        if not (truth_tighter.exists() and truth_tight.exists()):
            missing.append(f"{design}{index:03d}:truth")
            continue

        tt, tg = _load(truth_tight), _load(truth_tighter)
        truth_self = _err(tt, tg)
        ex_t, ex_g = _load(r_tight), _load(r_tighter)
        ex_self = _err(ex_t, ex_g)
        ex_err = _err(ex_g, tg)
        ex_env = ex_self + truth_self

        # The archived record keeps only the summed threshold, so every
        # envelope here is rebuilt the same way from the raw arrays -- the
        # exact-radius one and the two archived ones alike. Mixing a stored
        # threshold with a freshly computed one would compare two conventions.
        rec = json.loads(r19.out_path(design, BETA).read_text(encoding="utf-8"))
        arch = next(r for r in rec["rows"] if int(r["sobol_index"]) == index)
        fixed_err = arch["work_matched_error_m"]
        binned_err = arch["interior_error_m"]

        _, fx_t = r19.paths(design, BETA, index, "tight")
        _, fx_g = r19.paths(design, BETA, index, "tighter")
        _, bn_t = r18.paths(design, BETA, index, K, "tight")
        _, bn_g = r18.paths(design, BETA, index, K, "tighter")
        if not all(p.exists() for p in (fx_t, fx_g, bn_t, bn_g)):
            missing.append(f"{design}{index:03d}:archived-levels")
            continue
        fixed_env = _err(_load(fx_t), _load(fx_g)) + truth_self
        binned_env = _err(_load(bn_t), _load(bn_g)) + truth_self

        # The rebuilt binned threshold must reproduce the archived one. If it
        # does not, this pipeline is not on the archive's footing and the
        # comparison is reported as unreproduced rather than patched.
        archived_thr = arch.get("resolution_threshold_m")
        thr_rel = (abs((binned_env + fixed_env) - archived_thr) / archived_thr
                   if archived_thr else None)

        def verdict(ea, va, eb, vb):
            diff = eb - ea
            thr = va + vb
            if abs(diff) <= thr:
                return None, False
            return ("a" if diff > 0 else "b"), True

        w_ex, r_ex = verdict(ex_err, ex_env, fixed_err, fixed_env)
        w_bin, r_bin = verdict(binned_err, binned_env, fixed_err, fixed_env)

        w_exact = _work(s_tight)
        w_binned = None
        s18, _ = r18.paths(design, BETA, index, K, "tight")
        w_binned = _work(s18)
        rows.append({
            "design": design, "sobol_index": index,
            "degree_span": p["degree_span"],
            "work_matched_degree": p["work_matched_degree"],
            "exact_error_m": ex_err, "exact_envelope_m": ex_env,
            "exact_self_difference_m": ex_self,
            "binned_error_m": binned_err, "binned_envelope_m": binned_env,
            "fixed_error_m": fixed_err, "fixed_envelope_m": fixed_env,
            "truth_self_difference_m": truth_self,
            "archived_resolution_threshold_m": archived_thr,
            "rebuilt_threshold_rel_diff": thr_rel,
            "error_ratio_exact_over_binned": (ex_err / binned_err
                                              if binned_err else None),
            "realized_work_exact": w_exact,
            "realized_work_binned": w_binned,
            "work_ratio_exact_over_binned": (w_exact / w_binned
                                             if w_exact and w_binned else None),
            "resolved_binned": r_bin,
            "winner_binned": {"a": "interior", "b": "fixed"}.get(w_bin),
            "resolved_exact": r_ex,
            "winner_exact": {"a": "interior", "b": "fixed"}.get(w_ex),
            "verdict_changed": bool(
                r_bin and r_ex
                and {"a": "interior", "b": "fixed"}.get(w_bin)
                != {"a": "interior", "b": "fixed"}.get(w_ex)),
        })

    if not rows:
        print("[r24b] nothing to summarize yet")
        return 1

    wr = [r["work_ratio_exact_over_binned"] for r in rows
          if r["work_ratio_exact_over_binned"]]
    er = [r["error_ratio_exact_over_binned"] for r in rows
          if r["error_ratio_exact_over_binned"]]
    # Reporting "verdicts changed: 0" on its own would be true and misleading
    # here, because the dominant effect is not reversal but loss of resolution,
    # and the loss is driven by the exact-radius variant's own envelope rather
    # than by the two errors converging. The inflation is therefore a reported
    # quantity, not something a reader has to derive from the rows.
    inflate = [r["exact_envelope_m"] / r["binned_envelope_m"] for r in rows
               if r["binned_envelope_m"]]
    lost = [r for r in rows
            if r["resolved_binned"] and not r["resolved_exact"]]
    payload = {
        "schema": "r24_bin_control_v1", "created_utc": base.utc_now(),
        "environment": environment(),
        "beta": BETA, "interior_member_k": K_TAG,
        "preregistration": json.loads(
            PREREG.read_text(encoding="utf-8"))["preregistration_sha256"]
        if PREREG.exists() else None,
        "what_this_tests": (
            "whether the interior member's advantage at the budget-calibrated "
            "point survives removing the 10-km altitude quantization, with "
            "every other parameter held at its archived value"),
        "panel_completeness": {"aggregated": len(rows), "missing": len(missing),
                               "missing_detail": sorted(set(missing))},
        "summary": {
            "comparisons": len(rows),
            "verdicts_changed": sum(1 for r in rows if r["verdict_changed"]),
            "resolved_binned": sum(1 for r in rows if r["resolved_binned"]),
            "resolved_exact": sum(1 for r in rows if r["resolved_exact"]),
            "interior_wins_binned": sum(1 for r in rows
                                        if r["winner_binned"] == "interior"),
            "interior_wins_exact": sum(1 for r in rows
                                       if r["winner_exact"] == "interior"),
            "work_ratio_median": float(np.median(wr)) if wr else None,
            "work_ratio_min": float(np.min(wr)) if wr else None,
            "work_ratio_max": float(np.max(wr)) if wr else None,
            "error_ratio_median": float(np.median(er)) if er else None,
            "error_ratio_min": float(np.min(er)) if er else None,
            "error_ratio_max": float(np.max(er)) if er else None,
            "worst_threshold_reproduction_rel_diff": max(
                (r["rebuilt_threshold_rel_diff"] for r in rows
                 if r["rebuilt_threshold_rel_diff"] is not None), default=None),
            "lost_resolution": len(lost),
            "lost_resolution_detail": [
                {"orbit": f"{r['design']}{r['sobol_index']:03d}",
                 "was": r["winner_binned"],
                 "m_res_after": (abs(r["fixed_error_m"] - r["exact_error_m"])
                                 / (r["exact_envelope_m"]
                                    + r["fixed_envelope_m"]))}
                for r in lost],
            "envelope_inflation_exact_over_binned": {
                "median": float(np.median(inflate)) if inflate else None,
                "min": float(np.min(inflate)) if inflate else None,
                "max": float(np.max(inflate)) if inflate else None,
            },
            "errors_worse_under_exact": sum(
                1 for r in rows
                if (r["error_ratio_exact_over_binned"] or 0) > 1.0),
        },
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    s = payload["summary"]
    print(f"[r24b] written {OUTPUT.name}: {len(rows)} comparisons, "
          f"{len(missing)} missing")
    print(f"  resolved {s['resolved_binned']} binned -> "
          f"{s['resolved_exact']} exact; interior wins "
          f"{s['interior_wins_binned']} -> {s['interior_wins_exact']}")
    print(f"  verdicts reversed: {s['verdicts_changed']}; "
          f"resolution lost: {s['lost_resolution']}")
    inf = s["envelope_inflation_exact_over_binned"]
    if inf["median"]:
        print(f"  exact-radius envelope inflates by median "
              f"{inf['median']:.1f}x (range {inf['min']:.2f}-{inf['max']:.1f})")
    print(f"  errors worse under exact radius: "
          f"{s['errors_worse_under_exact']} of {len(rows)}")
    print(f"  archived-threshold reproduction, worst rel diff: "
          f"{s['worst_threshold_reproduction_rel_diff']:.2e}")
    if wr:
        print(f"  realized work exact/binned: median {s['work_ratio_median']:.3f} "
              f"({s['work_ratio_min']:.3f}-{s['work_ratio_max']:.3f})")
    if er:
        print(f"  error exact/binned: median {s['error_ratio_median']:.3f} "
              f"({s['error_ratio_min']:.3f}-{s['error_ratio_max']:.3f})")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preregister").set_defaults(fn=preregister)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--deadline-min", type=float, default=120.0)
    r.set_defaults(fn=run)
    sub.add_parser("summarize").set_defaults(fn=summarize)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
