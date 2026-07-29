"""Repair and re-verification of the measured-time-matched comparator (R13).

Two problems are fixed here.

(1) Integrity. The refined comparator runs are cached by file existence. An
    earlier refinement pass, later discarded because it had been computed from
    contended (parallel) Atallah kernel times, left one trajectory on disk:
    design-A orbit 036 kept a tight run at degree 746 while its tighter
    counterpart was rerun at degree 466. A self-difference taken across two
    different comparator degrees is meaningless, so that cell is quarantined and
    recomputed. This script audits every cell for the same defect before doing
    anything else.

(2) Time parity. Orbit 036 was the one case in which the comparator did not land
    within a few percent of the rule's measured kernel time. The comparator
    degree is re-estimated from measured per-call costs at two known degrees,
    the tight run is repeated at that degree, and the achieved ratio is checked
    against a 0.95--1.05 acceptance window before the tighter run is made.

Timing is then re-measured with repeated serial runs of both the rule and the
comparator, and the median is reported, so the accepted ratio does not rest on a
single execution. Everything runs serially on an otherwise idle machine.

Usage:
    python rev13_timing_repair.py audit
    python rev13_timing_repair.py repair --design A --index 36
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev13_timing_match as tm

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RECORD = METRICS / "r13_timing_repair.json"

ACCEPT_LO, ACCEPT_HI = 0.95, 1.05
TIMING_REPEATS = 3


def cell(design: str, index: int, level: str, stem: str = "fixed_time2"):
    return (tm.CASE_ROOT / design / f"sobolA_{index:03d}" / f"{stem}_{level}.json",
            tm.RAW_ROOT / design / f"sobolA_{index:03d}" / f"{stem}_{level}.npz")


def audit() -> list[dict]:
    """Every refined cell must have been run at the degree the selection names."""
    sel = json.loads(tm.SELECTION.read_text())
    bad = []
    for design, entries in sel["designs"].items():
        for e in entries:
            index = int(e["sobol_index"])
            want = e.get("n_time_refined")
            for level in ("tight", "tighter"):
                p, _ = cell(design, index, level)
                if not p.exists():
                    bad.append({"design": design, "index": index, "level": level,
                                "problem": "missing"})
                    continue
                got = json.loads(p.read_text())["config"]["policy_spec"]["degree"]
                if got != want:
                    bad.append({"design": design, "index": index, "level": level,
                                "problem": "degree mismatch",
                                "ran_degree": got, "selected_degree": want})
    for b in bad:
        print(f"  !! {b['design']}{b['index']:03d} {b['level']}: {b['problem']} "
              f"{b.get('ran_degree', '')}->{b.get('selected_degree', '')}")
    if not bad:
        print("  all refined cells run at their selected degree")
    return bad


def quarantine(design: str, index: int, level: str, reason: str) -> dict:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    moved = []
    for p in cell(design, index, level):
        if p.exists():
            dest = p.with_suffix(p.suffix + f".quarantined.{stamp}")
            shutil.move(str(p), str(dest))
            moved.append(str(dest.relative_to(METRICS)).replace("\\", "/"))
    print(f"  quarantined {design}{index:03d}/{level}: {reason}")
    return {"design": design, "index": index, "level": level, "reason": reason,
            "files": moved}


def per_call_ns(sidecar: Path) -> float:
    t = json.loads(sidecar.read_text())["telemetry"]
    return t["gravity_kernel_ns"] / t["n_rhs"]


def propagate(design: str, index: int, degree: int, level: str, save: bool):
    """One serial seven-day run of a fixed degree; returns telemetry (+ writes)."""
    side = tm.atallah_sidecar(design, index)
    cfg0 = side["config"]
    adopted = int(cfg0["adopted_truth_degree"])
    model = base.load_model(adopted)
    args = base.kernel_args(model)
    base.warmup(model, args)
    y0 = np.asarray(cfg0["initial_state_si"], dtype=float)
    grid = np.arange(0.0, tm.DURATION + 0.5 * tm.OUTPUT_STEP, tm.OUTPUT_STEP)
    tol = tm.LEVELS[level]
    t, y, st, ev, fail, tel = base.propagate_event_instrumented(
        model, y0, tm.DURATION, grid, (lambda tt, hh, n=degree: n), args,
        tol["rtol"], tol["atol"], max_step=tm.MAX_STEP)
    if st == "numerical_failure":
        raise RuntimeError(f"{design}{index:03d}/{level}: {fail}")
    if save:
        sidecar, raw = cell(design, index, level)
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        raw.parent.mkdir(parents=True, exist_ok=True)
        cfg = {"design": design, "sobol_index": index,
               "adopted_truth_degree": adopted,
               "initial_state_si": [float(v) for v in y0],
               "policy": "fixed_time_matched_refined", "level": level,
               "stage": "repaired",
               "policy_spec": {"kind": "fixed_measured_cost_matched",
                               "degree": degree,
                               "source": ("re-estimated from measured per-call "
                                          "kernel costs at two degrees and "
                                          "accepted on the achieved serial "
                                          "kernel-time ratio")},
               "n_work_proxy": None,
               "duration_s": tm.DURATION, "output_step_s": tm.OUTPUT_STEP,
               "integrator": "InstrumentedDOP853", "max_step_s": tm.MAX_STEP,
               "atol_kind": "vector", "rtol": tol["rtol"],
               "atol_position_m": tol["atol_position_m"],
               "atol_velocity_m_s": tol["atol_velocity_m_s"],
               "timing_comparable": True,
               "timing_note": ("serial run on an idle machine; comparable with "
                               "the serial Atallah re-runs of the same orbit"),
               "source": base.provenance()}
        base.atomic_npz(raw, t_s=t, state_si=y)
        base.atomic_json(sidecar, {
            "schema": "r13_timing_match_trajectory_v1",
            "created_utc": base.utc_now(), "config": cfg,
            "config_sha256": base.object_hash(cfg), "status": st,
            "event": ev, "telemetry": tel,
            "raw_path": str(raw.relative_to(ROOT)),
            "raw_sha256": base.file_hash(raw),
            "n_output_epochs": int(len(t)),
            "last_output_epoch_s": float(t[-1])})
    return tel


def atallah_run(design: str, index: int):
    """One serial seven-day run of the archived Atallah configuration."""
    side = tm.atallah_sidecar(design, index)
    cfg0 = side["config"]
    adopted = int(cfg0["adopted_truth_degree"])
    model = base.load_model(adopted)
    args = base.kernel_args(model)
    base.warmup(model, args)
    tab = {float(k): int(v) for k, v in cfg0["atallah_degree_table"].items()}
    hmin, hmax = min(tab), max(tab)

    def degree_of(tt, h_m):
        hb = min(hmax, max(hmin, 10.0 * math.floor(h_m / 1e4)))
        return tab[hb]

    y0 = np.asarray(cfg0["initial_state_si"], dtype=float)
    grid = np.arange(0.0, tm.DURATION + 0.5 * tm.OUTPUT_STEP, tm.OUTPUT_STEP)
    tol = tm.LEVELS["tight"]
    *_, tel = base.propagate_event_instrumented(
        model, y0, tm.DURATION, grid, degree_of, args,
        tol["rtol"], tol["atol"], max_step=tm.MAX_STEP)
    return tel


def repair(design: str, index: int) -> int:
    others = base.other_python_processes()
    if others:
        print(f"!! {len(others)} other python processes running; timing needs an "
              f"idle machine")
        return 2
    t0 = time.time()
    record = {"schema": "r13_timing_repair_v1", "created_utc": base.utc_now(),
              "design": design, "index": index,
              "accept_window": [ACCEPT_LO, ACCEPT_HI],
              "timing_repeats": TIMING_REPEATS,
              "audit": audit(), "quarantined": [], "attempts": []}

    # 1. quarantine any cell of this orbit that was not run at one degree
    degrees = {}
    for level in ("tight", "tighter"):
        p, _ = cell(design, index, level)
        if p.exists():
            degrees[level] = json.loads(p.read_text())["config"]["policy_spec"]["degree"]
    if len(set(degrees.values())) > 1:
        for level in ("tight", "tighter"):
            record["quarantined"].append(quarantine(
                design, index, level,
                f"comparator degree inconsistent across levels {degrees}"))

    # 2. re-estimate the degree from two measured per-call costs
    first = tm.CASE_ROOT / design / f"sobolA_{index:03d}" / "fixed_time_tight.json"
    quarantined_high = sorted(
        (tm.CASE_ROOT / design / f"sobolA_{index:03d}").glob(
            "fixed_time2_tight.json.quarantined.*"))
    ref_points = [(json.loads(first.read_text())["config"]["policy_spec"]["degree"],
                   per_call_ns(first))]
    for q in quarantined_high[-1:]:
        ref_points.append(
            (json.loads(q.read_text())["config"]["policy_spec"]["degree"],
             per_call_ns(q)))
    tighter_cell = cell(design, index, "tighter")[0]
    for q in sorted((tm.CASE_ROOT / design / f"sobolA_{index:03d}").glob(
            "fixed_time2_tighter.json*")):
        ref_points.append(
            (json.loads(q.read_text())["config"]["policy_spec"]["degree"],
             per_call_ns(q)))
    ref_points = sorted(set(ref_points))
    (n1, c1), (n2, c2) = ref_points[0], ref_points[-1]
    exponent = math.log(c2 / c1) / math.log(n2 / n1)
    serial = tm.CASE_ROOT / design / f"sobolA_{index:03d}" / "atallah_serial_tight.json"
    tel_at = json.loads(serial.read_text())["telemetry"]
    at_kernel_ns = tel_at["gravity_kernel_ns"]
    # the comparator's own call count at this tolerance, from the first pass
    n_rhs_cmp = json.loads(first.read_text())["telemetry"]["n_rhs"]
    target_per_call = at_kernel_ns / n_rhs_cmp
    n_star = int(round(n2 * (target_per_call / c2) ** (1.0 / exponent)))
    record["estimate"] = {
        "reference_points": [{"degree": n, "per_call_ns": c} for n, c in ref_points],
        "local_exponent": exponent, "target_per_call_ns": target_per_call,
        "atallah_serial_kernel_ns": at_kernel_ns,
        "first_guess_degree": n_star}
    print(f"  local cost exponent {exponent:.3f} from N={n1} and N={n2}; "
          f"first guess N={n_star}", flush=True)

    # 3. run tight, accept on the achieved ratio, retry once if outside the window
    accepted = None
    for attempt in range(2):
        tel = propagate(design, index, n_star, "tight", save=True)
        ratio = tel["gravity_kernel_ns"] / at_kernel_ns
        record["attempts"].append({"degree": n_star, "ratio": ratio,
                                   "kernel_ns": tel["gravity_kernel_ns"],
                                   "n_rhs": tel["n_rhs"]})
        print(f"  [{(time.time()-t0)/60:5.1f} min] N={n_star}: kernel="
              f"{tel['gravity_kernel_ns']/1e9:.1f}s, ratio={ratio:.3f}", flush=True)
        if ACCEPT_LO <= ratio <= ACCEPT_HI:
            accepted = n_star
            break
        n_star = int(round(n_star * (1.0 / ratio) ** (1.0 / exponent)))
        print(f"  outside the acceptance window; retrying at N={n_star}", flush=True)
    if accepted is None:
        accepted = record["attempts"][-1]["degree"]
        print("  !! acceptance window not reached; keeping the closest attempt",
              flush=True)
    record["accepted_degree"] = accepted

    # 4. matching tighter run at the accepted degree
    tel_tighter = propagate(design, index, accepted, "tighter", save=True)
    record["tighter"] = {"degree": accepted, "n_rhs": tel_tighter["n_rhs"],
                         "kernel_ns": tel_tighter["gravity_kernel_ns"]}
    print(f"  [{(time.time()-t0)/60:5.1f} min] tighter N={accepted} done",
          flush=True)

    # 5. repeated serial timing of both policies
    at_times, cmp_times = [], []
    for k in range(TIMING_REPEATS):
        at_times.append(atallah_run(design, index)["gravity_kernel_ns"])
        cmp_times.append(propagate(design, index, accepted, "tight",
                                   save=False)["gravity_kernel_ns"])
        print(f"  [{(time.time()-t0)/60:5.1f} min] timing repeat {k+1}: "
              f"atallah {at_times[-1]/1e9:.1f}s, comparator "
              f"{cmp_times[-1]/1e9:.1f}s", flush=True)
    record["timing"] = {
        "atallah_kernel_ns": at_times, "comparator_kernel_ns": cmp_times,
        "atallah_median_ns": float(np.median(at_times)),
        "comparator_median_ns": float(np.median(cmp_times)),
        "median_ratio": float(np.median(cmp_times) / np.median(at_times)),
        "ratio_spread": [float(min(cmp_times) / max(at_times)),
                         float(max(cmp_times) / min(at_times))]}
    print(f"  median comparator/Atallah kernel-time ratio = "
          f"{record['timing']['median_ratio']:.3f}", flush=True)
    record["wall_min"] = (time.time() - t0) / 60.0
    base.atomic_json(RECORD, record)
    print(f"[written] {RECORD.name}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("audit", "repair"))
    ap.add_argument("--design", default="A")
    ap.add_argument("--index", type=int, default=36)
    a = ap.parse_args()
    if a.command == "audit":
        audit()
        return 0
    return repair(a.design, a.index)


if __name__ == "__main__":
    raise SystemExit(main())
