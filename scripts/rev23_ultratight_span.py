"""R23-C: a third tolerance level for the surviving constructive comparison.

Why this exists
---------------
R23-A restricted the constructive claim to beta near 1: at beta = 0.5 the
interior member loses to a constant degree once realized work is equalized. The
evidence that remains is therefore the beta = 1 realized-work comparison, and
that comparison sits closer to the integration floor than any other the paper
relies on -- its median resolution threshold is a fraction of a metre. A claim
that now rests on one operating point has to be shown to be a truncation signal
rather than integrator noise, which is what a third tolerance level tests.

The design is the R13 retest applied to the constructive comparison instead of
the Atallah one. For each panel orbit the truth, the interior member k = 0.5 and
its realized-work-matched constant comparator are propagated at

    ultra:  rtol 3e-14, atol [1e-7 m] * 3 + [1e-10 m/s] * 3

one decade below ``tighter``, with the campaign's 60 s maximum step kept so the
tighter-to-ultra difference isolates the tolerance change. Errors are then taken
at the tighter level against the tighter truth, exactly as archived, while each
policy's envelope becomes its own tighter-to-ultra self-difference plus the
truth's. If a separation is physical it survives a shrinking envelope; if it is
noise it dissolves.

Nothing is recalibrated. The interior degree table and the comparator degree are
read out of the archived sidecars, so the objects propagated here are the ones
the archive already scored, not new ones that happen to share a label.

Panel and priority are fixed in r23_preregistration_amendment.json before any
result is inspected. The run is resume-safe and deadline-aware; if it stops
early it reports exactly which orbits completed and the aggregate is computed
over those alone, never over a silently truncated panel.

Usage:
    python rev23_ultratight_span.py run --workers 11 --deadline-min 540
    python rev23_ultratight_span.py aggregate
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev18_span_sweep as r18
import rev19_equal_total_work as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r23_cases" / "ultra_span"
RAW_ROOT = METRICS / "r23_raw" / "ultra_span"
OUTPUT = METRICS / "r23_ultratight_span.json"
PANEL_FILE = METRICS / "r23_ultra_panel.json"
R13_ULTRA_RAW = METRICS / "r13_raw" / "ultratight"
R13_ULTRA_CASE = METRICS / "r13_cases" / "ultratight"

ULTRA = {"rtol": 3.0e-14,
         "atol": np.array([1.0e-7] * 3 + [1.0e-10] * 3),
         "atol_position_m": 1.0e-7, "atol_velocity_m_s": 1.0e-10}
BETA = 1.00
K = 0.50
MAX_STEP = r14.MAX_STEP
DURATION = r14.DURATION
OUTPUT_STEP = r14.OUTPUT_STEP
POLICIES = ("truth", "interior", "fixed")

RESOLVED_CUT = 1.0
BORDERLINE_CUT = 0.5


def paths(design: str, index: int, policy: str):
    stem = f"{policy}_ultra"
    return (CASE_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.json",
            RAW_ROOT / design / f"sobolA_{index:03d}" / f"{stem}.npz")


def r13_truth_ultra(design: str, index: int) -> tuple[Path, Path]:
    return (R13_ULTRA_CASE / design / f"sobolA_{index:03d}" / "truth_ultra.json",
            R13_ULTRA_RAW / design / f"sobolA_{index:03d}" / "truth_ultra.npz")


def _degree_fn(table: dict):
    tab = {float(a): int(b) for a, b in table.items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, r18.BIN_KM * math.floor(h_m / 1e3 / r18.BIN_KM)))
        return tab[hb]

    return degree_of


# ------------------------------------------------------------------- panel
def build_panel() -> dict:
    """Pre-declared: every resolved R19 beta=1 comparison, then every
    borderline one. Ordering is group, then design, then Sobol index -- none
    of which correlates with cost, so an early stop truncates the panel
    without tilting it."""
    groups = {"resolved": [], "borderline": []}
    for design in ("A", "B"):
        rec = json.loads(r19.out_path(design, BETA).read_text(encoding="utf-8"))
        for r in rec["rows"]:
            thr = r.get("resolution_threshold_m")
            if not thr:
                continue
            gap = abs(r["work_matched_error_m"] - r["interior_error_m"])
            m = gap / thr
            item = {"design": design, "index": int(r["sobol_index"]),
                    "m_res": m,
                    "work_matched_degree": int(r["work_matched_degree"]),
                    "comparator_source": r["comparator_source"],
                    "previous_resolved": bool(r["resolved"]),
                    "previous_winner": r["winner"]}
            if m > RESOLVED_CUT:
                groups["resolved"].append(item)
            elif m > BORDERLINE_CUT:
                groups["borderline"].append(item)
    for g in groups.values():
        g.sort(key=lambda d: (d["design"], d["index"]))
    return groups


def ordered_tasks() -> list[dict]:
    panel = build_panel()
    tasks = []
    for group in ("resolved", "borderline"):
        for item in panel[group]:
            tasks.append({**item, "group": group})
    return tasks


# ------------------------------------------------------------------ worker
def _reusable_r13_truth(design: str, index: int, adopted: int) -> Path | None:
    """The R13 truth at ultra is the same object if every field that can move
    the trajectory matches. Anything less than an exact match is refused."""
    sidecar, raw = r13_truth_ultra(design, index)
    if not (sidecar.exists() and raw.exists()):
        return None
    try:
        cfg = json.loads(sidecar.read_text(encoding="utf-8"))["config"]
    except (OSError, KeyError, json.JSONDecodeError):
        return None
    same = (int(cfg.get("adopted_truth_degree", -1)) == adopted
            and cfg.get("policy") == "truth"
            and cfg.get("level") == "ultra"
            and float(cfg.get("rtol", 0)) == ULTRA["rtol"]
            and float(cfg.get("atol_position_m", 0)) == ULTRA["atol_position_m"]
            and float(cfg.get("atol_velocity_m_s", 0)) == ULTRA["atol_velocity_m_s"]
            and float(cfg.get("max_step_s", 0)) == MAX_STEP
            and float(cfg.get("duration_s", 0)) == DURATION
            and float(cfg.get("output_step_s", 0)) == OUTPUT_STEP)
    return raw if same else None


def worker(task: dict) -> dict:
    design, index = task["design"], int(task["index"])
    t_start = time.time()
    try:
        r18_side, _ = r18.paths(design, BETA, index, K, "tight")
        r18_cfg = json.loads(r18_side.read_text(encoding="utf-8"))["config"]
        adopted = int(r18_cfg["adopted_truth_degree"])
        y0 = np.asarray(r18_cfg["initial_state_si"], dtype=float)
        table_k = r18_cfg["degree_table"]

        n_fixed = int(task["work_matched_degree"])
        selectors = {
            "truth": (lambda t, h, n=adopted: n),
            "interior": _degree_fn(table_k),
            "fixed": (lambda t, h, n=n_fixed: n),
        }
        model, args = r14._model(adopted)
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        reused = []

        for policy in POLICIES:
            sidecar, raw = paths(design, index, policy)
            cfg = {
                "design": design, "sobol_index": index, "beta": BETA,
                "policy": policy, "level": "ultra",
                "adopted_truth_degree": adopted,
                "interior_k": K,
                "policy_spec": (
                    {"kind": "fixed_truth", "degree": adopted} if policy == "truth"
                    else {"kind": "fixed_work_matched", "degree": n_fixed}
                    if policy == "fixed"
                    else {"kind": "binned_span_member", "k": K,
                          "degree_table": table_k}),
                "initial_state_si": [float(v) for v in y0],
                "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
                "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
                "atol_kind": "vector", "rtol": ULTRA["rtol"],
                "atol_position_m": ULTRA["atol_position_m"],
                "atol_velocity_m_s": ULTRA["atol_velocity_m_s"],
                "timing_comparable": False,
                "purpose": ("third tolerance level for the realized-work "
                            "interior-versus-constant comparison at beta = 1"),
                "source": task.get("provenance", {})}
            digest = base.object_hash(cfg)

            if sidecar.exists() and raw.exists():
                prev = json.loads(sidecar.read_text(encoding="utf-8"))
                if prev.get("config_sha256") == digest \
                        and prev.get("status") == "ok":
                    continue

            if policy == "truth":
                donor = _reusable_r13_truth(design, index, adopted)
                if donor is not None:
                    d = np.load(donor)
                    base.atomic_npz(raw, t_s=d["t_s"], state_si=d["state_si"])
                    base.atomic_json(sidecar, {
                        "schema": "r23_ultratight_span_v1",
                        "created_utc": base.utc_now(), "config": cfg,
                        "config_sha256": digest, "status": "ok",
                        "reused_from": str(donor.relative_to(ROOT)),
                        "reuse_basis": ("identical truth degree, tolerances, "
                                        "max step, duration and output grid; "
                                        "the integrator is deterministic"),
                        "raw_path": str(raw.relative_to(ROOT)),
                        "raw_sha256": base.file_hash(raw),
                        "n_output_epochs": int(len(d["t_s"])),
                        "last_output_epoch_s": float(d["t_s"][-1])})
                    reused.append(policy)
                    continue

            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, selectors[policy], args,
                ULTRA["rtol"], ULTRA["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"design": design, "index": index,
                        "status": "numerical_failure",
                        "where": policy, "detail": fail}
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r23_ultratight_span_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": digest, "status": "ok",
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
        return {"design": design, "index": index, "status": "ok",
                "adopted_truth_degree": adopted, "reused": reused,
                "minutes": (time.time() - t_start) / 60.0}
    except Exception as exc:                                    # noqa: BLE001
        return {"design": design, "index": index, "status": "error",
                "detail": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def run(args) -> int:
    tasks = ordered_tasks()
    PANEL_FILE.write_text(json.dumps(
        {"created_utc": base.utc_now(),
         "rule": ("every resolved R19 beta=1 realized-work comparison "
                  "(M_res > 1) as a validation set, plus every borderline one "
                  "(0.5 < M_res <= 1); order is group, design, Sobol index"),
         "resolved": [t for t in tasks if t["group"] == "resolved"],
         "borderline": [t for t in tasks if t["group"] == "borderline"]},
        indent=2), encoding="utf-8")
    print(f"[r23c] panel {len(tasks)} orbits "
          f"({sum(1 for t in tasks if t['group'] == 'resolved')} resolved, "
          f"{sum(1 for t in tasks if t['group'] == 'borderline')} borderline), "
          f"workers={args.workers}", flush=True)
    t0 = time.time()
    deadline = t0 + args.deadline_min * 60.0
    done = fail = 0
    completed: list[str] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        try:
            for fut in as_completed(futs):
                res = fut.result()
                done += 1
                tag = f"{res['design']}{res['index']:03d}"
                if res["status"] != "ok":
                    fail += 1
                    print(f"  [FAIL] {tag} {res['status']} "
                          f"{res.get('where', '')}: {res.get('detail')}",
                          flush=True)
                else:
                    completed.append(tag)
                    print(f"  [{done}/{len(tasks)}] {tag} ok "
                          f"N_truth={res['adopted_truth_degree']} "
                          f"reused={res['reused']} "
                          f"took={res['minutes']:.1f}min "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
                if time.time() > deadline:
                    print("[r23c] deadline; cancelling pending", flush=True)
                    for f in futs:
                        f.cancel()
                    break
        except CancelledError:
            pass
    print(f"[r23c] {done} finished, {fail} failed, "
          f"wall={(time.time()-t0)/60:.1f} min", flush=True)
    return aggregate(argparse.Namespace())


# --------------------------------------------------------------- aggregate
def _load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def _diff(a, b) -> float:
    return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]


def _degree_in(sidecar: Path) -> int | None:
    try:
        cfg = json.loads(sidecar.read_text(encoding="utf-8"))["config"]
    except (OSError, KeyError, json.JSONDecodeError):
        return None
    if "degree" in cfg:
        return int(cfg["degree"])
    spec = cfg.get("policy_spec") or {}
    return int(spec["degree"]) if "degree" in spec else None


def _archived_fixed_tighter(design: str, index: int, n_fixed: int):
    """The comparator trajectory the archive already scored, at tighter.

    R19 wrote its own sidecar on most orbits. Where the work-matched degree
    rounded back to the constant endpoint it scored the archived endpoint
    instead, and that endpoint is itself stored in one of two places depending
    on whether R14 had to propagate it or could reuse the critical-altitude run.
    Whichever candidate is found, its recorded degree must equal the degree R19
    actually scored; a mismatch returns nothing so the orbit is reported
    missing rather than compared against the wrong object.
    """
    for sidecar, raw in (r19.paths(design, BETA, index, "tighter"),
                         r14.paths(design, BETA, index, "fixed_budget",
                                   "tighter"),
                         r14.reuse_paths(design, index, "fixed_critical",
                                         "tighter")):
        if not (sidecar.exists() and raw.exists()):
            continue
        deg = _degree_in(sidecar)
        if deg is not None and deg != n_fixed:
            continue
        return raw, deg
    return None, None


def _archived_tighter(design: str, index: int, n_fixed: int, source: str):
    """The tighter-level trajectories the archive already scored."""
    _, truth = r14.reuse_paths(design, index, "truth", "tighter")
    _, interior = r18.paths(design, BETA, index, K, "tighter")
    fixed, _deg = _archived_fixed_tighter(design, index, n_fixed)
    return truth, interior, fixed


def aggregate(args) -> int:
    tasks = {(t["design"], t["index"]): t for t in ordered_tasks()}
    rows, missing = [], []
    for (design, index), t in sorted(tasks.items()):
        if not all(paths(design, index, p)[0].exists() for p in POLICIES):
            missing.append({"design": design, "index": index,
                            "group": t["group"]})
            continue
        truth_t, int_t, fix_t = _archived_tighter(
            design, index, t["work_matched_degree"], t["comparator_source"])
        if fix_t is None or not all(p.exists() for p in (truth_t, int_t)):
            missing.append({"design": design, "index": index,
                            "group": t["group"],
                            "why": "archived tighter trajectory absent"})
            continue
        tt, it, ft = _load(truth_t), _load(int_t), _load(fix_t)
        tu = _load(paths(design, index, "truth")[1])
        iu = _load(paths(design, index, "interior")[1])
        fu = _load(paths(design, index, "fixed")[1])

        truth_self = _diff(tt, tu)
        int_self, fix_self = _diff(it, iu), _diff(ft, fu)
        e_int = _diff(it, tt)
        e_fix = _diff(ft, tt)
        env_int = int_self + truth_self
        env_fix = fix_self + truth_self
        gap = abs(e_fix - e_int)
        thr = env_int + env_fix
        gap_u = abs(_diff(fu, tu) - _diff(iu, tu))
        rows.append({
            "design": design, "sobol_index": index, "group": t["group"],
            "adopted_truth_degree": None,
            "interior_error_tighter_m": e_int,
            "fixed_error_tighter_m": e_fix,
            "interior_envelope_m": env_int,
            "fixed_envelope_m": env_fix,
            "truth_self_difference_tighter_to_ultra_m": truth_self,
            "absolute_error_difference_m": gap,
            "resolution_threshold_m": thr,
            "m_res_after": gap / thr if thr else None,
            "m_res_before": t["m_res"],
            "resolved_after": bool(gap > thr),
            "resolved_before": t["previous_resolved"],
            "winner_after": (("interior" if e_int < e_fix else "fixed")
                             if gap > thr else None),
            "winner_before": t["previous_winner"],
            "gap_at_ultra_m": gap_u,
            "gap_sign_stable_tighter_to_ultra": bool(
                np.sign(e_fix - e_int) == np.sign(_diff(fu, tu) - _diff(iu, tu))),
        })

    if not rows:
        print("[r23c] nothing to aggregate yet")
        return 1

    summary = {}
    for group in ("resolved", "borderline", "all"):
        sub = [r for r in rows if group == "all" or r["group"] == group]
        if not sub:
            continue
        thr_b = np.array([tasks[(r["design"], r["sobol_index"])]["m_res"]
                          for r in sub])
        summary[group] = {
            "orbits": len(sub),
            "resolved_after": int(sum(r["resolved_after"] for r in sub)),
            "interior_wins": int(sum(r["winner_after"] == "interior"
                                     for r in sub)),
            "fixed_wins": int(sum(r["winner_after"] == "fixed" for r in sub)),
            "unresolved_after": int(sum(not r["resolved_after"] for r in sub)),
            "m_res_median_before": float(np.median(thr_b)),
            "m_res_median_after": float(np.median(
                [r["m_res_after"] for r in sub if r["m_res_after"]])),
            "envelope_median_m": float(np.median(
                [r["resolution_threshold_m"] for r in sub])),
            "gap_sign_stable": int(sum(
                r["gap_sign_stable_tighter_to_ultra"] for r in sub)),
            # A flip is a comparison that was decided both times and changed
            # sides. Counting "was undecided, now decided" as a flip would
            # inflate this into the headline number of the campaign, which is
            # the opposite of what it means.
            "verdict_flips": int(sum(
                r["winner_before"] is not None
                and r["winner_after"] is not None
                and r["winner_after"] != r["winner_before"] for r in sub)),
            "became_resolved": int(sum(
                not r["resolved_before"] and r["resolved_after"] for r in sub)),
            "became_unresolved": int(sum(
                r["resolved_before"] and not r["resolved_after"] for r in sub)),
        }
    payload = {
        "schema": "r23_ultratight_span_v1", "created_utc": base.utc_now(),
        "beta": BETA, "interior_k": K,
        "level": {k: v for k, v in ULTRA.items() if k != "atol"},
        "max_step_s": MAX_STEP,
        "preregistration": "r23_preregistration_amendment.json",
        "what_this_tests": (
            "whether the beta = 1 realized-work separations are truncation "
            "signal or integration noise: errors stay at the tighter level "
            "against the tighter truth, while every envelope is rebuilt from "
            "the tighter-to-ultra self-difference"),
        "panel_completeness": {
            "aggregated": len(rows), "missing": len(missing),
            "missing_detail": missing,
            "note": ("the aggregate covers the completed orbits only and the "
                     "incomplete ones are listed; it is never computed over a "
                     "silently truncated panel")},
        "summary": summary, "rows": rows}
    base.atomic_json(OUTPUT, payload)
    print(f"[r23c] written {OUTPUT.name}: {len(rows)} orbits aggregated, "
          f"{len(missing)} missing")
    print(json.dumps(summary, indent=2))
    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=11)
    r.add_argument("--deadline-min", type=float, default=540.0)
    r.set_defaults(func=run)
    a_ = sub.add_parser("aggregate")
    a_.set_defaults(func=aggregate)
    a = p.parse_args()
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
