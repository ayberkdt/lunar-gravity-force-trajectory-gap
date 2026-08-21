"""J1: does the force--trajectory reversal survive a second GRAIL solution?

The manuscript calibrates its recipe on several independent coefficient
products, but every propagated result rests on one of them (JPL JGGRX_1800F).
That leaves one obvious question open: is the reversal -- the radial allocation
winning on the force metric while losing on the trajectory metric at the same
declared budget -- a property of the *policy pair*, or of the *coefficient
realization* it happened to be demonstrated on?

This campaign answers it on GSFC GRGM1200A, an independently produced solution
of the same body. Nothing is carried across from the primary field: the tail
exponent, the empirical degree table, the critical degree, the truth degree and
the Atallah accuracy parameter are all re-derived from GRGM1200A's own
spectrum. The population is a fresh scrambled-Sobol draw with its own seed on
the same design box, so the replication is of the *finding*, not of the
arithmetic.

Stages (each resumable, each refusing to start what it cannot finish):

    calibrate  field-level recalibration of the new solution
    design     freeze the 32-point population and register the outcome rule
    prepass    per-orbit critical / work / truth degrees on the new field
    base       reference and constant-degree trajectories, two tolerances
    radial     budget-calibrated radial policy, force defect, trajectories
    score      per-orbit and population verdicts

Usage:
    python revJ1_crossfield.py calibrate
    python revJ1_crossfield.py design
    python revJ1_crossfield.py prepass --workers 11
    python revJ1_crossfield.py base    --workers 11
    python revJ1_crossfield.py radial  --workers 11
    python revJ1_crossfield.py score
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import revJ_common as J

J.select_field("GRGM1200A")
J.install_field()

import rev3_common as rc                                          # noqa: E402
import rev10_sobol_confirmatory as base                           # noqa: E402
import rev12_atallah as at                                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RAW_ROOT = Path(os.environ.get("JCAMP_RAW_ROOT",
                               r"D:\makale_raw_offload\jgcd")) / "J1"
CASE_ROOT = METRICS / "rJ1_cases"
LOG = Path(__file__).resolve().parent / "rJ1_campaign.log"

SEED = 20260808            # drawn once, before any trajectory was inspected
N_POINTS = 32              # random_base2(m=5)
LEVEL_NAMES = ("tight", "tighter")

# The budget is read from the environment rather than from argv so that parent
# and spawned children agree on it: a Windows ProcessPoolExecutor child does
# not inherit the parent's argv, but it does inherit the environment.
BETA = float(os.environ.get("JCAMP_BETA", "1.00"))
SUFFIX = "" if abs(BETA - 1.0) < 1e-12 else f"_beta_{BETA:.2f}"

CALIB = METRICS / "rJ1_field_calibration.json"
DESIGN = METRICS / "rJ1_design.json"
PREREG = METRICS / "rJ1_preregistration.json"
ROWS = METRICS / "rJ1_rows.json"
DESIGN_EXT = METRICS / "rJ1_design_extension.json"
ROWS_EXT = METRICS / "rJ1_rows_extension.json"
RADIAL = METRICS / f"rJ1_radial{SUFFIX}.json"
SCORE = METRICS / f"rJ1_score{SUFFIX}.json"


def log(msg: str) -> None:
    J.log_line(LOG, f"J1 {msg}")


# ------------------------------------------------------------- 1. calibrate
def command_calibrate() -> int:
    """Re-derive the compact tail rule on GRGM1200A's own spectrum.

    This repeats the main text's field-level calibration rather than importing
    its answer. The primary field's p_fit is reported alongside only so that a
    reader can see whether the two solutions agree; it is never substituted.
    """
    import rev16_multibody_calibration as r16

    field = J.FIELDS["GRGM1200A"]
    path = Path(field["path"])
    spec = dict(key="GRGM1200A", body="Moon", center="GSFC",
                nmax_file=field["max_degree_in_file"], fmt="pds_km",
                role="cross_solution_trajectory_level", path=str(path))
    res = r16.analyze(spec)
    if res["file_sha256"] != field["sha256"]:
        raise SystemExit("GRGM1200A digest does not match the registered value")

    primary = None
    prev = METRICS / "r16_multibody_calibration.json"
    if prev.exists():
        for row in json.loads(prev.read_text(encoding="utf-8"))["fields"]:
            if row["key"] == "JGGRX_1800F":
                primary = {"p_fit": row["p_fit"], "p_safe": row["p_safe"],
                           "spectral_slope_p": row["spectral_slope"]["p"]}

    payload = {
        "schema": "rJ1_field_calibration_v1",
        "created_utc": J.utc_now(),
        "purpose": ("independent recalibration of the truncation recipe on the "
                    "cross-solution field, so that no number is inherited from "
                    "the primary solution"),
        "field": {k: res[k] for k in
                  ("key", "body", "center", "file", "file_sha256",
                   "max_degree_in_file", "mu_m3_s2", "reference_radius_m")},
        "spectral_slope": res["spectral_slope"],
        "p_fit": res["p_fit"], "p_safe": res["p_safe"],
        "sse_p_fit": res["sse_p_fit"], "sse_p2": res["sse_p2"],
        "rms_mismatch_p_fit": res["rms_mismatch_p_fit"],
        "rms_mismatch_p2": res["rms_mismatch_p2"],
        "holdout": res["holdout"],
        "emp_range": res["emp_range"],
        "cap_guard_ok": res["cap_guard_ok"],
        "cap_guard_note": res["cap_guard_note"],
        "criteria_rows": res["criteria_rows"],
        "primary_field_for_comparison_only": primary,
        "provenance": J.provenance(),
    }
    J.atomic_json(CALIB, payload)
    log(f"calibrate: p_fit={res['p_fit']:.3f} p_safe={res['p_safe']} "
        f"N_emp {res['emp_range']} cap_ok={res['cap_guard_ok']}")
    return 0


# ---------------------------------------------------------------- 2. design
def command_design() -> int:
    if DESIGN.exists() and PREREG.exists():
        log("design: already frozen")
        return 0
    if not CALIB.exists():
        raise SystemExit("run calibrate first")
    from scipy.stats import qmc

    model = rc.load_model(300)
    points = qmc.Sobol(d=5, scramble=True, seed=SEED).random_base2(m=5)
    orbits = [base.orbit_from_u(i, row, "sobolJ1", model)
              for i, row in enumerate(points)]
    for o in orbits:
        # The adopted truth degree follows the design-B convention: the highest
        # degree the design-A audit ever adopted for the sub-50 km regime.
        o["adopted_truth_degree"] = 300 if o["hp_km"] >= 50.0 else 900
    if len(np.unique(np.asarray([o["u"] for o in orbits]), axis=0)) != N_POINTS:
        raise SystemExit("duplicate design point")
    for o in orbits:
        if abs(o["longitude_roundtrip_error_deg"]) > 1.0e-10:
            raise SystemExit(f"{o['name']}: longitude roundtrip failed")

    design = {
        "schema": "rJ1_design_v1", "created_utc": J.utc_now(),
        "family": "sobolJ1", "seed": SEED, "sample_count": len(orbits),
        "dimension": 5, "generator": "scrambled Sobol random_base2(m=5)",
        "box": "identical to the manuscript's confirmatory design box",
        "field": J.field_key(),
        "n_rejected": 0, "n_replaced": 0, "optimized": False,
        "orbits": orbits,
    }
    design["design_sha256"] = J.object_hash(design)
    J.atomic_json(DESIGN, design)

    prereg = {
        "schema": "rJ1_preregistration_v1", "created_utc": J.utc_now(),
        "status": "frozen_before_any_J1_trajectory_was_propagated",
        "question": ("does the force--trajectory reversal reproduce under an "
                     "independently produced GRAIL gravity solution?"),
        "field": {"key": "GRGM1200A",
                  "sha256": J.FIELDS["GRGM1200A"]["sha256"],
                  "recalibrated": True,
                  "calibration_record": CALIB.name},
        "population": {"file": DESIGN.name,
                       "design_sha256": design["design_sha256"],
                       "sample_count": N_POINTS},
        "arc": {"duration_s": J.DURATION, "output_step_s": J.OUTPUT_STEP,
                "max_step_s": J.MAX_STEP},
        "policies": {
            "reference": "constant adopted truth degree",
            "constant": f"N_F(beta) = argmin_N |N^2 - beta N_crit^2|, beta={BETA}",
            "radial": ("budget-calibrated Atallah radial rule, accuracy "
                       "parameter bisected so the 10-km-binned degree history "
                       f"spends <N^2> = beta N_crit^2 with beta={BETA}"),
        },
        "levels": {k: {"rtol": J.LEVELS[k]["rtol"],
                       "atol_position_m": J.LEVELS[k]["atol_position_m"],
                       "atol_velocity_m_s": J.LEVELS[k]["atol_velocity_m_s"]}
                   for k in LEVEL_NAMES},
        "metrics": {
            "J_force": "time average of |a_policy(x_ref(t)) - a_ref(x_ref(t))|",
            "J_traj": "position RMS against the reference on the output grid",
        },
        "resolution_rule": ("a trajectory comparison counts only when the error "
                            "gap exceeds the sum of the two policies' numerical "
                            "envelopes, each envelope being the policy's own "
                            "between-level self-difference plus the reference's"),
        "declared_outcomes": {
            "A_reversal_reproduced": ("on the resolved orbits the radial policy "
                                      "wins the force metric in the majority and "
                                      "loses the trajectory metric in the "
                                      "majority, with the same sign as the "
                                      "primary field"),
            "B_reversal_absent": ("the trajectory-metric sign agrees with the "
                                  "force-metric sign, i.e. no reversal on this "
                                  "solution"),
            "C_unresolved": ("fewer than half the orbits resolve, and no "
                             "population statement is made"),
        },
        "reporting_commitment": ("all three outcomes are reported as found; the "
                                "per-orbit table is written whatever the "
                                "aggregate says"),
        "provenance": J.provenance(),
    }
    prereg["preregistration_sha256"] = J.object_hash(prereg)
    J.atomic_json(PREREG, prereg)
    hp = [o["hp_km"] for o in orbits]
    log(f"design: {len(orbits)} orbits frozen, perilune {min(hp):.1f}-"
        f"{max(hp):.1f} km, {sum(1 for h in hp if h < 50)} below 50 km, "
        f"sha {design['design_sha256'][:16]}")
    return 0


def all_rows() -> list[dict]:
    """The frozen population, plus the sequence extension if it has been run.

    A scrambled Sobol draw of 2^6 points begins with exactly the 2^5 points of
    the smaller draw from the same seed, so the extension is a continuation of
    the same sequence rather than a second, differently balanced design. That
    property is asserted, not assumed, when the extension is frozen.
    """
    rows = json.loads(ROWS.read_text(encoding="utf-8"))["rows"]
    if ROWS_EXT.exists():
        rows = rows + json.loads(ROWS_EXT.read_text(encoding="utf-8"))["rows"]
    return sorted(rows, key=lambda r: r["sobol_index"])


# --------------------------------------------------------------- 3. prepass
def prepass_task(orbit: dict) -> dict:
    try:
        hp = float(orbit["hp_km"])
        adopted = int(orbit["adopted_truth_degree"])
        original = 300 if hp >= 50.0 else 600
        model, args = J.model_for(original)
        power = J.power_for(original)
        n_crit = J.critical_degree(power, model.r_ref, hp)
        sched = J.table_sched(J.empirical_table(power, model.r_ref))
        y0 = np.asarray(orbit["initial_state_si"], dtype=float)
        grid = J.out_grid()
        Y, rhs, info = rc.propagate_instr(
            model, y0, J.DURATION, grid, sched, args,
            J.LEVELS["tight"]["rtol"], J.LEVELS["tight"]["atol"],
            max_step=J.MAX_STEP)
        mean_n2 = rhs.sum_deg_sq / rhs.n_calls
        return {"ok": True, "row": {
            "sobol_index": int(orbit["sobol_index"]),
            "name": orbit["name"], "hp_km": hp, "ha_km": float(orbit["ha_km"]),
            "incl_deg": float(orbit["incl_deg"]),
            "argp_deg": float(orbit["argp_deg"]),
            "raan_deg": float(orbit["raan_deg"]),
            "eccentricity": float(orbit["eccentricity"]),
            "initial_state_si": [float(v) for v in orbit["initial_state_si"]],
            "original_truth_degree": original,
            "adopted_truth_degree": adopted,
            "n_critical": int(n_crit),
            "n_work": int(round(math.sqrt(mean_n2))),
            "empirical_mean_degree_sq": float(mean_n2),
            "prepass_n_rhs": int(info["n_rhs"]),
            "prepass_wall_s": float(info["wall_s"]),
        }}
    except Exception:
        return {"ok": False, "sobol_index": int(orbit["sobol_index"]),
                "error": traceback.format_exc()}


def command_prepass(workers: int) -> int:
    design = json.loads(DESIGN.read_text(encoding="utf-8"))
    done = {}
    if ROWS.exists():
        done = {r["sobol_index"]: r
                for r in json.loads(ROWS.read_text(encoding="utf-8"))["rows"]}
    todo = [o for o in design["orbits"] if o["sobol_index"] not in done]
    log(f"prepass: {len(done)} done, {len(todo)} to run")
    if todo:
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(prepass_task, o): o["sobol_index"]
                    for o in todo}
            for fut in as_completed(futs):
                res = fut.result()
                if not res["ok"]:
                    log(f"prepass FAIL {res['sobol_index']}\n{res['error']}")
                    continue
                done[res["row"]["sobol_index"]] = res["row"]
                if len(done) % 4 == 0:
                    _write_rows(design, done)
        log(f"prepass: finished in {(time.time() - t0) / 60:.1f} min")
    _write_rows(design, done)
    return 0 if len(done) == len(design["orbits"]) else 1


def command_extend(workers: int) -> int:
    """Continue the same Sobol sequence from 32 points to 64.

    This is an extension, not a second design. The check that the first 32
    points of the 2^6 draw reproduce the frozen 2^5 draw exactly is what makes
    that claim true rather than merely intended, and it fails loudly if the
    generator's behaviour is not what is assumed here.
    """
    from scipy.stats import qmc

    frozen = json.loads(DESIGN.read_text(encoding="utf-8"))
    if not DESIGN_EXT.exists():
        model = rc.load_model(300)
        points = qmc.Sobol(d=5, scramble=True, seed=SEED).random_base2(m=6)
        head = np.asarray([o["u"] for o in frozen["orbits"]])
        if not np.allclose(points[:N_POINTS], head, rtol=0.0, atol=0.0):
            raise SystemExit("the 64-point draw does not begin with the frozen "
                             "32-point draw; the populations would not nest "
                             "and this cannot be called an extension")
        orbits = [base.orbit_from_u(i, points[i], "sobolJ1", model)
                  for i in range(N_POINTS, len(points))]
        for o in orbits:
            o["adopted_truth_degree"] = 300 if o["hp_km"] >= 50.0 else 900
        ext = {"schema": "rJ1_design_extension_v1", "created_utc": J.utc_now(),
               "family": "sobolJ1", "seed": SEED,
               "generator": "scrambled Sobol random_base2(m=6), points 32-63",
               "extends": {"file": DESIGN.name,
                           "design_sha256": frozen["design_sha256"],
                           "nesting_verified": True},
               "sample_count": len(orbits), "orbits": orbits,
               "field": J.field_key()}
        ext["design_sha256"] = J.object_hash(ext)
        J.atomic_json(DESIGN_EXT, ext)
        hp = [o["hp_km"] for o in orbits]
        log(f"extend: {len(orbits)} further orbits frozen, perilune "
            f"{min(hp):.1f}-{max(hp):.1f} km, "
            f"{sum(1 for h in hp if h < 50)} below 50 km")

    ext = json.loads(DESIGN_EXT.read_text(encoding="utf-8"))
    done = {}
    if ROWS_EXT.exists():
        done = {r["sobol_index"]: r for r in
                json.loads(ROWS_EXT.read_text(encoding="utf-8"))["rows"]}
    todo = [o for o in ext["orbits"] if o["sobol_index"] not in done]
    if todo:
        log(f"extend: prepass for {len(todo)} orbits")
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for fut in as_completed([pool.submit(prepass_task, o)
                                     for o in todo]):
                res = fut.result()
                if not res["ok"]:
                    log(f"extend prepass FAIL {res['sobol_index']}\n"
                        f"{res['error']}")
                    continue
                done[res["row"]["sobol_index"]] = res["row"]
    J.atomic_json(ROWS_EXT, {
        "schema": "rJ1_rows_extension_v1", "created_utc": J.utc_now(),
        "field": J.field_key(), "design_sha256": ext["design_sha256"],
        "complete": len(done) == len(ext["orbits"]),
        "rows": [done[k] for k in sorted(done)],
        "provenance": J.provenance()})
    if len(done) != len(ext["orbits"]):
        return 1
    rc_base = command_base(workers)
    rc_rad = command_radial(workers)
    if rc_base == 0 and rc_rad == 0:
        command_score()
    J.atomic_json(METRICS / "rJ1_extension.json", {
        "schema": "rJ1_extension_v1", "created_utc": J.utc_now(),
        "complete": rc_base == 0 and rc_rad == 0,
        "population_after_extension": len(all_rows()),
        "base_rc": rc_base, "radial_rc": rc_rad,
        "design_extension_sha256": ext["design_sha256"]})
    return 0 if rc_base == 0 and rc_rad == 0 else 1


def _write_rows(design: dict, done: dict) -> None:
    rows = [done[k] for k in sorted(done)]
    J.atomic_json(ROWS, {
        "schema": "rJ1_rows_v1", "created_utc": J.utc_now(),
        "field": J.field_key(),
        "design_sha256": design["design_sha256"],
        "complete": len(rows) == len(design["orbits"]),
        "definitions": {
            "n_critical": "min(250, unquantized empirical tail degree at the "
                          "orbit's own perilune) on this field",
            "n_work": "round(sqrt(<N^2>)) from one empirical-schedule "
                      "propagation at the tight vector tolerance",
            "adopted_truth_degree": "300 for perilune >= 50 km, else 900",
        },
        "rows": rows, "provenance": J.provenance()})


# ------------------------------------------------------------------ 4. base
def case_paths(index: int, policy: str, level: str):
    """The reference does not depend on the budget, so it is never re-run for
    a second budget; the two budget-dependent policies carry the tag."""
    tag = "" if policy == "reference" else SUFFIX
    stem = f"{policy}{tag}_{level}"
    return (CASE_ROOT / f"J1_{index:03d}" / f"{stem}.json",
            RAW_ROOT / f"J1_{index:03d}" / f"{stem}.npz")


def _run_policy(row: dict, policy: str, level: str, degree_of, spec: dict):
    adopted = int(row["adopted_truth_degree"])
    model, args = J.model_for(adopted)
    y0 = np.asarray(row["initial_state_si"], dtype=float)
    grid = J.out_grid()
    lv = J.LEVELS[level]
    Y, rhs, info = rc.propagate_instr(model, y0, J.DURATION, grid, degree_of,
                                      args, lv["rtol"], lv["atol"],
                                      max_step=J.MAX_STEP)
    cj, cn = case_paths(int(row["sobol_index"]), policy, level)
    J.atomic_npz(cn, t=grid, y=Y)
    record = {
        "schema": "rJ1_case_v1", "created_utc": J.utc_now(),
        "config": {"sobol_index": int(row["sobol_index"]), "policy": policy,
                   "policy_spec": spec, "level": level,
                   "adopted_truth_degree": adopted,
                   "duration_s": J.DURATION, "output_step_s": J.OUTPUT_STEP,
                   "max_step_s": J.MAX_STEP, "rtol": lv["rtol"],
                   "atol_kind": "vector",
                   "atol_position_m": lv["atol_position_m"],
                   "atol_velocity_m_s": lv["atol_velocity_m_s"],
                   "integrator": "InstrumentedDOP853",
                   "field": J.field_key(),
                   "timing_comparable": False,
                   "execution": "parallel_process_pool"},
        "telemetry": {"n_rhs": int(info["n_rhs"]),
                      "n_accepted_steps": int(info["n_accepted_steps"]),
                      "n_attempted_steps": int(info["n_attempted_steps"]),
                      "n_rejected_trials": int(info["n_rejected_trials"]),
                      "mean_degree_sq": float(rhs.sum_deg_sq / rhs.n_calls),
                      "degree_range": [int(min(rhs.deg_counts)),
                                       int(max(rhs.deg_counts))],
                      "gravity_wall_s": float(info["grav_s"]),
                      "total_wall_s": float(info["wall_s"])},
        "raw_path": str(cn), "raw_sha256": J.sha256_file(cn),
    }
    J.atomic_json(cj, record)
    return record


def base_task(task: dict) -> dict:
    row, policy, level = task["row"], task["policy"], task["level"]
    try:
        if policy == "reference":
            n = int(row["adopted_truth_degree"])
            spec = {"kind": "constant_reference", "degree": n}
        else:
            n = int(task["degree"])
            spec = {"kind": "constant_budget", "degree": n, "beta": BETA,
                    "rule": "argmin_N |N^2 - beta N_crit^2|",
                    "n_critical": int(row["n_critical"])}
        rec = _run_policy(row, policy, level, lambda t, h, _n=n: _n, spec)
        return {"ok": True, "index": int(row["sobol_index"]), "policy": policy,
                "level": level, "wall_s": rec["telemetry"]["total_wall_s"]}
    except Exception:
        return {"ok": False, "index": int(row["sobol_index"]), "policy": policy,
                "level": level, "error": traceback.format_exc()}


def command_base(workers: int, deadline: str | None = None) -> int:
    rows = all_rows()
    tasks = []
    for row in rows:
        n_f, capped = J.fixed_degree_for(BETA, int(row["n_critical"]),
                                         int(row["adopted_truth_degree"]))
        for level in LEVEL_NAMES:
            for policy, degree in (("reference", None), ("constant", n_f)):
                cj, cn = case_paths(int(row["sobol_index"]), policy, level)
                if cj.exists() and cn.exists():
                    continue
                tasks.append({"row": row, "policy": policy, "level": level,
                              "degree": degree, "capped": capped})
    log(f"base: {len(tasks)} trajectories to run")
    rc_base = _drive(tasks, base_task, workers, deadline, "base")
    J.atomic_json(METRICS / f"rJ1_base_complete{SUFFIX}.json", {
        "schema": "rJ1_base_complete_v1", "created_utc": J.utc_now(),
        "beta": BETA, "orbits": len(rows), "trajectories_run": len(tasks),
        "complete": rc_base == 0})
    return rc_base


# ---------------------------------------------------------------- 5. radial
def radial_task(task: dict) -> dict:
    row = task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        model, args = J.model_for(adopted)
        _, ref_raw = case_paths(index, "reference", "tighter")
        with np.load(ref_raw) as z:
            t_ref, Y_ref = z["t"], z["y"]
        h_km = (np.linalg.norm(Y_ref[:3], axis=0) - model.r_ref) / 1e3
        n_crit = int(row["n_critical"])
        target = BETA * n_crit ** 2
        cal = J.calibrate_radial(adopted, float(row["hp_km"]),
                                 float(row["ha_km"]), adopted, h_km, target)
        table = {float(k): int(v) for k, v in cal["table"].items()}
        deg_radial = J.degrees_from_table(table, h_km)
        n_f, capped = J.fixed_degree_for(BETA, n_crit, adopted)
        deg_fixed = np.full(len(h_km), n_f, dtype=int)

        defect = J.force_defects(t_ref, Y_ref[:3],
                                 {"radial": deg_radial, "constant": deg_fixed},
                                 adopted, args)
        defect_radial, defect_fixed = defect["radial"], defect["constant"]

        spec = {"kind": "budget_calibrated_radial", "beta": BETA,
                "accuracy_parameter_m_s2": cal["tol"],
                "target_mean_degree_sq": target,
                "achieved_mean_degree_sq": cal["work"],
                "work_mismatch": cal["mismatch"],
                "attainable": bool(cal["attainable"]),
                "limit": cal["limit"],
                "bin_km": J.BIN_KM, "floor": J.FLOOR, "cap": adopted,
                "n_critical": n_crit,
                "degree_span": [int(deg_radial.min()), int(deg_radial.max())]}
        degree_of, _ = at.atallah_binned_schedule(
            model, J.atallah_g(adopted), cal["tol"], float(row["hp_km"]),
            float(row["ha_km"]), floor=J.FLOOR, cap=adopted, bin_km=J.BIN_KM)
        walls = {}
        for level in LEVEL_NAMES:
            cj, cn = case_paths(index, "radial", level)
            if cj.exists() and cn.exists():
                walls[level] = json.loads(cj.read_text(
                    encoding="utf-8"))["telemetry"]["total_wall_s"]
                continue
            rec = _run_policy(row, "radial", level, degree_of, spec)
            walls[level] = rec["telemetry"]["total_wall_s"]
        return {"ok": True, "index": index, "record": {
            "sobol_index": index, "name": row["name"],
            "hp_km": row["hp_km"], "ha_km": row["ha_km"],
            "incl_deg": row["incl_deg"],
            "adopted_truth_degree": adopted, "n_critical": n_crit,
            "beta": BETA,
            "constant_degree": int(n_f), "constant_capped": bool(capped),
            "radial": spec,
            "radial_degree_table": {str(k): int(v) for k, v in
                                    sorted(table.items())},
            "force_defect_radial": defect_radial,
            "force_defect_constant": defect_fixed,
            "wall_s": walls}}
    except Exception:
        return {"ok": False, "index": index, "error": traceback.format_exc()}


def command_radial(workers: int, deadline: str | None = None) -> int:
    rows = all_rows()
    have = {}
    if RADIAL.exists():
        have = {r["sobol_index"]: r for r in
                json.loads(RADIAL.read_text(encoding="utf-8"))["rows"]}
    tasks = [{"row": r} for r in rows if r["sobol_index"] not in have]
    log(f"radial: {len(have)} done, {len(tasks)} to run")
    if tasks:
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(radial_task, t) for t in tasks]
            for fut in as_completed(futs):
                res = fut.result()
                if not res["ok"]:
                    log(f"radial FAIL {res['index']}\n{res['error']}")
                    continue
                have[res["index"]] = res["record"]
                _write_radial(have, len(rows))
        log(f"radial: finished in {(time.time() - t0) / 60:.1f} min")
    _write_radial(have, len(rows))
    return 0 if len(have) == len(rows) else 1


def _write_radial(have: dict, total: int) -> None:
    J.atomic_json(RADIAL, {
        "schema": "rJ1_radial_v1", "created_utc": J.utc_now(),
        "field": J.field_key(), "beta": BETA,
        "complete": len(have) == total,
        "rows": [have[k] for k in sorted(have)],
        "provenance": J.provenance()})


# ---------------------------------------------------- 6. force-only Pareto
BUDGET_GRID = (0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00)
PARETO = METRICS / "rJ1_budget_pareto.json"


def pareto_task(task: dict) -> dict:
    """Force-metric ranking across the budget grid, with no propagation.

    The defect is a deterministic function of the reference trajectory and the
    degree history, so the whole budget grid can be swept on trajectories that
    already exist. This is the cross-solution counterpart of the manuscript's
    Phase-A sweep, and it costs no integration: what it cannot say is anything
    about trajectory error, which is exactly why the campaign propagates the
    declared budget separately.
    """
    row = task["row"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        model, args = J.model_for(adopted)
        _, ref_raw = case_paths(index, "reference", "tighter")
        t_ref, Y_ref = J.load_times(ref_raw), J.load_states(ref_raw)
        h_km = (np.linalg.norm(Y_ref[:3], axis=0) - model.r_ref) / 1e3
        n_crit = int(row["n_critical"])

        budgets = {}
        degree_sets = {}
        for beta in BUDGET_GRID:
            cal = J.calibrate_radial(adopted, float(row["hp_km"]),
                                     float(row["ha_km"]), adopted, h_km,
                                     beta * n_crit ** 2)
            table = {float(k): int(v) for k, v in cal["table"].items()}
            n_f, capped = J.fixed_degree_for(beta, n_crit, adopted)
            deg_r = J.degrees_from_table(table, h_km)
            degree_sets[f"radial_{beta:.2f}"] = deg_r
            degree_sets[f"constant_{beta:.2f}"] = np.full(len(h_km), n_f,
                                                          dtype=int)
            budgets[f"beta_{beta:.2f}"] = {
                # A policy that has reached the reference degree has no
                # truncation error to measure: its defect is identically zero,
                # and a ratio against it is not a comparison. The endpoint case
                # N = N_ref is censored here as well as N > N_ref, which the
                # cap flag alone does not catch.
                "beta_requested": beta,
                "censored": bool(capped or int(deg_r.max()) >= adopted
                                 or n_f >= adopted),
                "radial": {"accuracy_parameter_m_s2": cal["tol"],
                           "achieved_mean_degree_sq": cal["work"],
                           "work_mismatch": cal["mismatch"],
                           "attainable": bool(cal["attainable"]),
                           "degree_span": [int(deg_r.min()), int(deg_r.max())]},
                "constant": {"degree": int(n_f), "capped": bool(capped)}}
        defects = J.force_defects(t_ref, Y_ref[:3], degree_sets, adopted, args)
        for beta in BUDGET_GRID:
            key = f"beta_{beta:.2f}"
            dr = defects[f"radial_{beta:.2f}"]
            dc = defects[f"constant_{beta:.2f}"]
            budgets[key]["radial"]["defect"] = dr
            budgets[key]["constant"]["defect"] = dc
            den_mean = dc["J_force_mean_m_s2"]
            den_rms = dc["J_force_rms_m_s2"]
            budgets[key]["rho_force_mean"] = (dr["J_force_mean_m_s2"] / den_mean
                                              if den_mean > 0.0 else None)
            budgets[key]["rho_force_rms"] = (dr["J_force_rms_m_s2"] / den_rms
                                             if den_rms > 0.0 else None)
            if den_mean <= 0.0:
                budgets[key]["censored"] = True
        return {"ok": True, "index": index, "record": {
            "sobol_index": index, "name": row["name"],
            "hp_km": row["hp_km"], "ha_km": row["ha_km"],
            "n_critical": n_crit, "adopted_truth_degree": adopted,
            "budgets": budgets}}
    except Exception:
        return {"ok": False, "index": index, "error": traceback.format_exc()}


def command_pareto(workers: int) -> int:
    rows = all_rows()
    have = {}
    if PARETO.exists():
        have = {r["sobol_index"]: r for r in
                json.loads(PARETO.read_text(encoding="utf-8"))["rows"]}
    todo = [{"row": r} for r in rows if r["sobol_index"] not in have]
    log(f"pareto: {len(have)} done, {len(todo)} to run over "
        f"{len(BUDGET_GRID)} budgets")
    if todo:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for fut in as_completed([pool.submit(pareto_task, t)
                                     for t in todo]):
                res = fut.result()
                if not res["ok"]:
                    log(f"pareto FAIL {res['index']}\n{res['error']}")
                    continue
                have[res["index"]] = res["record"]
    summary = {}
    for beta in BUDGET_GRID:
        key = f"beta_{beta:.2f}"
        live = [r["budgets"][key] for r in have.values()
                if not r["budgets"][key]["censored"]
                and r["budgets"][key]["rho_force_mean"] is not None]
        if not live:
            continue
        summary[key] = {
            "orbits_uncensored": len(live),
            "orbits_censored": len(have) - len(live),
            "radial_wins_force": sum(1 for b in live
                                     if b["rho_force_mean"] < 1.0),
            "median_rho_force_mean": float(np.median(
                [b["rho_force_mean"] for b in live])),
            "median_rho_force_rms": float(np.median(
                [b["rho_force_rms"] for b in live])),
            "max_work_mismatch": float(max(abs(b["radial"]["work_mismatch"])
                                           for b in live))}
    J.atomic_json(PARETO, {
        "schema": "rJ1_budget_pareto_v1", "created_utc": J.utc_now(),
        "field": J.field_key(),
        "note": ("force-metric only; no propagation is involved, so nothing "
                 "here speaks to trajectory error"),
        "budget_grid": list(BUDGET_GRID),
        "censoring_rule": ("a budget is censored for an orbit when the policy "
                           "reaches the adopted reference degree, because a "
                           "defect measured against a reference the policy has "
                           "reached is not a truncation error"),
        "complete": len(have) == len(rows),
        "summary": summary,
        "rows": [have[k] for k in sorted(have)],
        "provenance": J.provenance()})
    for key, s in summary.items():
        log(f"pareto {key}: radial wins force {s['radial_wins_force']}/"
            f"{s['orbits_uncensored']}, median rho_F(mean) "
            f"{s['median_rho_force_mean']:.3f}")
    return 0 if len(have) == len(rows) else 1


# ----------------------------------------------------------------- 7. score
def command_score() -> int:
    rows = {r["sobol_index"]: r for r in all_rows()}
    radial = {r["sobol_index"]: r
              for r in json.loads(RADIAL.read_text(encoding="utf-8"))["rows"]}
    out = []
    for index in sorted(rows):
        if index not in radial:
            continue
        Y = {}
        missing = False
        for policy in ("reference", "constant", "radial"):
            for level in LEVEL_NAMES:
                _, raw = case_paths(index, policy, level)
                if not raw.exists():
                    missing = True
                    break
                with np.load(raw) as z:
                    Y[(policy, level)] = z["y"]
            if missing:
                break
        if missing:
            continue
        self_ref = J.self_difference(Y[("reference", "tight")],
                                     Y[("reference", "tighter")])
        rec = {"sobol_index": index, "name": rows[index]["name"],
               "hp_km": rows[index]["hp_km"], "ha_km": rows[index]["ha_km"],
               "incl_deg": rows[index]["incl_deg"],
               "n_critical": rows[index]["n_critical"],
               "adopted_truth_degree": rows[index]["adopted_truth_degree"],
               "self_difference_reference_m": self_ref}
        for policy in ("constant", "radial"):
            err = {lv: J.trajectory_error(Y[(policy, lv)],
                                          Y[("reference", lv)])["J_traj_rms_m"]
                   for lv in LEVEL_NAMES}
            self_p = J.self_difference(Y[(policy, "tight")],
                                       Y[(policy, "tighter")])
            rec[policy] = {
                "J_traj_rms_m": err,
                "self_difference_m": self_p,
                "numerical_envelope_m": self_p + self_ref,
                "J_force_mean_m_s2":
                    radial[index][f"force_defect_"
                                  f"{'radial' if policy == 'radial' else 'constant'}"
                                  ]["J_force_mean_m_s2"],
                "J_force_rms_m_s2":
                    radial[index][f"force_defect_"
                                  f"{'radial' if policy == 'radial' else 'constant'}"
                                  ]["J_force_rms_m_s2"],
            }
        for lv in LEVEL_NAMES:
            e_c = rec["constant"]["J_traj_rms_m"][lv]
            e_r = rec["radial"]["J_traj_rms_m"][lv]
            rec.setdefault("resolved", {})[lv] = bool(J.resolved(
                e_c, e_r, rec["constant"]["numerical_envelope_m"],
                rec["radial"]["numerical_envelope_m"]))
            rec.setdefault("margin", {})[lv] = float(
                abs(e_c - e_r) / (rec["constant"]["numerical_envelope_m"]
                                  + rec["radial"]["numerical_envelope_m"]))
        rec["rho_traj"] = (rec["radial"]["J_traj_rms_m"]["tighter"]
                           / rec["constant"]["J_traj_rms_m"]["tighter"])
        rec["rho_force"] = (rec["radial"]["J_force_mean_m_s2"]
                            / rec["constant"]["J_force_mean_m_s2"])
        rec["reversal"] = bool(rec["rho_force"] < 1.0 < rec["rho_traj"])
        rec["work_mismatch"] = radial[index]["radial"]["work_mismatch"]
        out.append(rec)

    res = [r for r in out if r["resolved"]["tighter"]]
    n_rev = sum(1 for r in res if r["reversal"])
    n_force_win = sum(1 for r in res if r["rho_force"] < 1.0)
    n_traj_lose = sum(1 for r in res if r["rho_traj"] > 1.0)
    verdict = ("A_reversal_reproduced"
               if len(res) >= 0.5 * len(out) and n_force_win > len(res) / 2
               and n_traj_lose > len(res) / 2
               else ("C_unresolved" if len(res) < 0.5 * len(out)
                     else "B_reversal_absent"))
    payload = {
        "schema": "rJ1_score_v1", "created_utc": J.utc_now(),
        "field": J.field_key(), "beta": BETA,
        "preregistration_sha256": json.loads(
            PREREG.read_text(encoding="utf-8"))["preregistration_sha256"],
        "counts": {"orbits": len(out), "resolved_tighter": len(res),
                   "resolved_tight": sum(1 for r in out
                                         if r["resolved"]["tight"]),
                   "radial_wins_force": n_force_win,
                   "radial_loses_trajectory": n_traj_lose,
                   "reversal_orbits": n_rev},
        "aggregates": {
            "median_rho_force": float(np.median([r["rho_force"] for r in out])),
            "median_rho_traj": float(np.median([r["rho_traj"] for r in out])),
            "median_rho_force_resolved": float(np.median(
                [r["rho_force"] for r in res])) if res else None,
            "median_rho_traj_resolved": float(np.median(
                [r["rho_traj"] for r in res])) if res else None,
            "max_work_mismatch": float(max(r["work_mismatch"] for r in out)),
        },
        "verdict": verdict,
        "rows": out,
        "provenance": J.provenance(),
    }
    J.atomic_json(SCORE, payload)
    log(f"score: {verdict}; resolved {len(res)}/{len(out)}, "
        f"force-wins {n_force_win}, trajectory-loses {n_traj_lose}, "
        f"reversal {n_rev}")
    return 0


# ------------------------------------------------------------------- driver
def _drive(tasks, fn, workers, deadline, label) -> int:
    if not tasks:
        return 0
    stop = None
    if deadline:
        from datetime import datetime as _dt
        stop = _dt.fromisoformat(deadline)
    t0 = time.time()
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(fn, t) for t in tasks]
        for fut in as_completed(futs):
            res = fut.result()
            if res["ok"]:
                ok += 1
            else:
                fail += 1
                log(f"{label} FAIL {res['index']} {res['policy']} "
                    f"{res['level']}\n{res['error']}")
            if (ok + fail) % 8 == 0:
                log(f"{label}: {ok + fail}/{len(tasks)} "
                    f"({(time.time() - t0) / 60:.1f} min)")
            if stop is not None:
                from datetime import datetime as _dt
                if _dt.now() > stop:
                    log(f"{label}: deadline reached, cancelling remainder")
                    for f in futs:
                        f.cancel()
                    break
    log(f"{label}: {ok} ok, {fail} failed, {(time.time() - t0) / 60:.1f} min")
    return 0 if fail == 0 else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("calibrate", "design", "prepass",
                                        "base", "radial", "score", "extend",
                                        "pareto", "status"))
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    if a.command == "calibrate":
        return command_calibrate()
    if a.command == "design":
        return command_design()
    if a.command == "prepass":
        return command_prepass(a.workers)
    if a.command == "base":
        return command_base(a.workers, a.deadline)
    if a.command == "radial":
        return command_radial(a.workers, a.deadline)
    if a.command == "score":
        return command_score()
    if a.command == "extend":
        return command_extend(a.workers)
    if a.command == "pareto":
        return command_pareto(a.workers)
    for p in (CALIB, DESIGN, ROWS, RADIAL, SCORE):
        print(f"{p.name:34s} {'present' if p.exists() else 'MISSING'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
