"""Independent design-B confirmation population at vector tolerance (R11).

Purpose
-------
The main 64-orbit scheduling result was established on the design-A scrambled
Sobol population.  ``revision/FINAL_SUBMISSION_REPORT_R10.md`` lists a frozen,
*unpropagated* second population (design B, ``r10_sobolB_design_frozen.json``,
seed 20260724) as future work, explicitly noting it "does not represent the
unpropagated second Sobol scramble as evidence".  Propagating design B under
the same vector-tolerance convergence contract turns the main result from a
single-population finding into an independently replicated one.

Design B is a genuinely independent draw: a different Sobol seed, its own
maximin rejection, spanning the same physical bounds (perilune 30-150 km,
inclination 2-179 deg, eccentricity up to ~0.4).  Eleven of its orbits have
sub-50 km perilune, exactly the regime that required a high-degree truth audit
for design A.

What this script does
---------------------
1. ``prepass`` -- for each of the 64 design-B orbits it derives, without any
   manuscript claim attached:
     * ``original_truth_degree``  = 300 (perilune >= 50 km) else 600, the same
       rule design A used to build schedules and the work/critical degrees;
     * ``n_critical``  = min(250, empirical Nmin at the orbit's own perilune),
       a closed-form field quantity needing no propagation;
     * ``n_work``  = round(sqrt(mean N^2)) from one empirical-schedule
       propagation at the tight vector tolerance, exactly the design-A
       definition (the schedule's altitude choice is tolerance-independent);
     * ``adopted_truth_degree``  = 300 (perilune >= 50 km) else 900.  Using the
       highest degree design A's audit ever adopted for the sub-50 km regime is
       conservative: a higher truth degree is never less accurate.  The N600 vs
       N900 degree-adequacy gap is recorded per sub-50 km orbit as a
       self-check.
   The result is written as ``metrics/r11_designB_rows.json`` in the same
   schema as the design-A corrected baseline.

2. ``run`` -- rebinds the validated ``rev11_full_convergence`` machinery onto
   the design-B rows and output tree and runs the identical
   64-orbit x 6-policy x 2-level campaign.  Only the input population differs;
   the integrator, tolerances, error metric, and resolution rule are shared
   code, so the two campaigns are directly comparable.

Usage
-----
    python rev11_designB_convergence.py prepass --workers 5
    python rev11_designB_convergence.py run --workers 5 --deadline 2026-07-24T16:30:00+03:00
    python rev11_designB_convergence.py status
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev11_full_convergence as fc
from rev3_common import degree_power
from rev7_doe_screening import (CAP, alt_sched, emp_nmin_exact, emp_table,
                                initial_state)

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
DESIGN_B = METRICS / "r10_sobolB_design_frozen.json"
ROWS = METRICS / "r11_designB_rows.json"

TIGHT = {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3)}
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
MAX_STEP = 60.0

# design-B output tree (kept separate from the design-A r11 artifacts)
OUTPUT = METRICS / "r11_designB_convergence.json"
SMOKE_OUTPUT = METRICS / "r11_designB_convergence_smoke.json"
CASE_ROOT = METRICS / "r11_cases" / "designB_convergence"
RAW_ROOT = METRICS / "r11_raw" / "designB_convergence"

_MODELS: dict[int, tuple] = {}


def _model(degree: int):
    if degree not in _MODELS:
        model = base.load_model(degree)
        args = base.kernel_args(model)
        base.warmup(model, args)
        _MODELS[degree] = (model, args)
    return _MODELS[degree]


def design_points() -> list[dict]:
    payload = json.loads(DESIGN_B.read_text(encoding="utf-8"))
    orbits = (payload.get("orbits") or payload.get("rows")
              or payload.get("design_points"))
    if len(orbits) != 64:
        raise RuntimeError(f"expected 64 design-B orbits, got {len(orbits)}")
    return orbits


# ---------------------------------------------------------------- prepass task
def prepass_task(orbit: dict) -> dict:
    """Derive n_work / n_critical / degree adequacy for one orbit."""
    try:
        hp = float(orbit["hp_km"])
        original = 300 if hp >= 50.0 else 600
        adopted = 300 if hp >= 50.0 else 900
        model, args = _model(original)
        power = degree_power(model)
        n_critical = int(min(CAP, emp_nmin_exact(power, model.r_ref, hp * 1e3)))
        schedule = alt_sched(emp_table(model, power))
        y0 = np.asarray(orbit["initial_state_si"], dtype=float)
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        t, y, status, event, failure, telemetry = \
            base.propagate_event_instrumented(
                model, y0, DURATION, grid, schedule, args, TIGHT["rtol"],
                TIGHT["atol"], max_step=MAX_STEP)
        mean_n2 = telemetry.get("mean_degree_sq")
        if mean_n2 is None:
            raise RuntimeError("empirical schedule produced no RHS calls")
        n_work = int(round(math.sqrt(mean_n2)))

        record = {
            "sobol_index": int(orbit["sobol_index"]),
            "name": orbit["name"],
            "design_point": orbit,
            "original_truth_degree": original,
            "adopted_truth_degree": adopted,
            "n_work": n_work,
            "n_critical": n_critical,
            "empirical_prepass": {
                "status": status,
                "mean_degree_sq": mean_n2,
                "degree_range": telemetry.get("degree_range"),
                "surface_impact": bool(event is not None)},
        }
        # degree-adequacy self-check for the sub-50 km regime
        if hp < 50.0:
            m600, a600 = _model(600)
            m900, a900 = _model(900)
            _, y600, s600, _, _, _ = base.propagate_event_instrumented(
                m600, y0, DURATION, grid, lambda t, h: 600, a600,
                TIGHT["rtol"], TIGHT["atol"], max_step=MAX_STEP)
            _, y900, s900, _, _, _ = base.propagate_event_instrumented(
                m900, y0, DURATION, grid, lambda t, h: 900, a900,
                TIGHT["rtol"], TIGHT["atol"], max_step=MAX_STEP)
            n = min(y600.shape[1], y900.shape[1])
            d = np.linalg.norm(y600[:3, :n] - y900[:3, :n], axis=0)
            record["degree_adequacy"] = {
                "N600_status": s600, "N900_status": s900,
                "N600_vs_N900_pos_rms_m": float(np.sqrt(np.mean(d * d))),
                "N600_vs_N900_pos_max_m": float(np.max(d)),
                "adopted": adopted,
                "note": "adopted N900 as truth; gap bounds residual degree error"}
        return record
    except Exception as exc:
        return {"sobol_index": int(orbit.get("sobol_index", -1)),
                "name": orbit.get("name"), "prepass_error":
                f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def prepass(workers: int) -> int:
    orbits = design_points()
    print(f"[prepass] deriving n_work/n_critical/degree-adequacy for "
          f"{len(orbits)} design-B orbits, workers={workers}", flush=True)
    rows, failures = [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(prepass_task, o) for o in orbits]
        for n, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            if "prepass_error" in record:
                failures.append(record)
                print(f"  !! {record['sobol_index']}: "
                      f"{record['prepass_error']}", flush=True)
            else:
                rows.append(record)
            if n % 8 == 0 or n == len(orbits):
                el = time.time() - t0
                print(f"  [{n:2d}/{len(orbits)}] elapsed={el/60:5.1f} min "
                      f"eta={(len(orbits)-n)*el/n/60:5.1f} min", flush=True)
    rows.sort(key=lambda r: r["sobol_index"])
    payload = {
        "schema": "r11_designB_rows_v1",
        "created_utc": base.utc_now(),
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "design_frozen_sha256": base.file_hash(DESIGN_B),
        "script_sha256": base.file_hash(Path(__file__).resolve()),
        "truth_degree_rule": "300 if perilune>=50km else 600 (schedule basis); "
                             "adopted truth 300 if perilune>=50km else 900",
        "tight_level_for_prepass": {"rtol": TIGHT["rtol"],
                                    "atol_position_m": 1e-5,
                                    "atol_velocity_m_s": 1e-8},
        "rows": rows, "failures": failures}
    base.atomic_json(ROWS, payload)
    print(f"[prepass] wrote {len(rows)}/{len(orbits)} rows, "
          f"failures={len(failures)} -> {ROWS.name}", flush=True)
    return 0 if len(rows) == len(orbits) else 3


# ------------------------------------------------------------------- run (reuse)
def install() -> None:
    """Point the validated design-A machinery at the design-B population.

    The path override is done through the environment, NOT by rebinding
    ``fc`` module globals: ``fc.run`` executes trajectories in a
    ProcessPoolExecutor, and on Windows spawn each worker re-imports
    ``rev11_full_convergence`` fresh, so a parent-only monkey-patch would leave
    the workers writing to the design-A tree.  Environment variables are
    inherited by spawned children; ``rev11_full_convergence`` reads them at
    import, so both parent and workers resolve to the design-B tree.  The env
    must be set before this process imported ``fc`` for the parent to see it
    too, so we also rebind the already-imported parent globals here.
    """
    os.environ["R11_TREE"] = "designB_convergence"
    os.environ["R11_CORRECTED"] = str(ROWS)
    os.environ["R11_OUTPUT"] = str(OUTPUT)
    os.environ["R11_SMOKE_OUTPUT"] = str(SMOKE_OUTPUT)
    # parent already imported fc with design-A defaults; realign it now
    fc.CORRECTED = ROWS
    fc.OUTPUT = OUTPUT
    fc.SMOKE_OUTPUT = SMOKE_OUTPUT
    fc.CASE_ROOT = CASE_ROOT
    fc.RAW_ROOT = RAW_ROOT
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)


def status() -> int:
    if not OUTPUT.exists():
        print("no design-B convergence output yet")
        return 0
    data = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({"complete": data["complete"],
                      "orbits": len(data["rows"]),
                      "failures": len(data.get("failures", [])),
                      "summary": data["summary"]}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepass", "run", "smoke", "status"))
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--deadline")
    args = parser.parse_args()
    if args.command == "prepass":
        return prepass(args.workers)
    if args.command == "status":
        install()
        return status()
    install()
    if not ROWS.exists():
        raise RuntimeError("run prepass first: r11_designB_rows.json missing")
    return fc.run(args.command == "smoke", fc.parse_deadline(args.deadline),
                  args.workers)


if __name__ == "__main__":
    raise SystemExit(main())
