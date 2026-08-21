"""J2: does the force--trajectory reversal survive realistic lunar dynamics?

The headline comparison is established on a gravity-only system: a uniformly
rotating body-fixed spherical-harmonic field and nothing else. That is the
right place to establish it, because it is the only setting in which the
measured difference between two policies can be attributed to truncation alone.
It also invites one obvious objection -- that the reversal is an artifact of the
isolated system, and that under the perturbations a real lunar orbiter actually
experiences the ordering would not survive.

This campaign answers that objection on a population rather than on a single
control orbit. The same three policies are propagated under

    a = a_Moon,SH(N(t)) + a_Earth + a_Sun + a_SRP,

with DE440 ephemerides, the MOON_PA lunar orientation, differential third-body
gravity for Earth and Sun, and cannonball SRP with lunar eclipse. Every one of
those additional accelerations is *identical* for the reference and for both
policies, so the quantity being measured is still the truncation policy
difference and nothing else.

Two deliberate choices make this a replication rather than a new experiment:

  * the population is a deterministic subset of the frozen confirmatory design
    -- the first 24 orbits by Sobol index, chosen by rule and not by result --
    so every J2 orbit has an archived gravity-only counterpart and the two can
    be compared pair by pair;
  * the policies are the archived ones. The critical degree, the constant
    comparator degree and the radial rule's accuracy parameter are read from
    the frozen budget records, not recalibrated, so the *policy* is literally
    the one the manuscript propagated. The work it realizes under the new
    dynamics is measured and reported rather than assumed.

Usage:
    python revJ2_fullforce.py select
    python revJ2_fullforce.py run --workers 11
    python revJ2_fullforce.py score
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

import revJ_common as J

J.select_field("JGGRX_1800F")
J.install_field()

import rev3_common as rc                                          # noqa: E402
import rev12_atallah as at                                        # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RAW_ROOT = Path(os.environ.get("JCAMP_RAW_ROOT",
                               r"D:\makale_raw_offload\jgcd")) / "J2"
CASE_ROOT = METRICS / "rJ2_cases"
LOG = Path(__file__).resolve().parent / "rJ2_campaign.log"

LEVEL_NAMES = ("tight", "tighter")
POLICIES = ("reference", "constant", "radial")

# Read from the environment, not argv: a spawned worker inherits the former.
BETA = float(os.environ.get("JCAMP_BETA", "1.00"))
# The population rule is "the first N by Sobol index", so a larger N contains
# the smaller one: extending the campaign never re-runs what is already on
# disk, and never changes which orbits the first 24 were.
N_ORBITS = int(os.environ.get("JCAMP_J2_ORBITS", "24"))
# Which archived confirmatory population is replicated under full dynamics.
DESIGN = os.environ.get("JCAMP_J2_DESIGN", "A").upper()
# Recalibrate the radial policy on the full-dynamics reference instead of
# reading the archived accuracy parameter. Off by default: the campaign's
# primary arm is the archived policy.
RECAL = os.environ.get("JCAMP_J2_RECAL", "") not in ("", "0")
# Spacecraft area-to-mass ratio in m^2/kg. The archived control uses
# 1 m^2 / 100 kg; a different value changes only the SRP magnitude.
AOVERM = float(os.environ.get("JCAMP_J2_AOVERM", "0.01"))

SUFFIX = "" if abs(BETA - 1.0) < 1e-12 else f"_beta_{BETA:.2f}"
POP_SUFFIX = "" if N_ORBITS == 24 else f"_n{N_ORBITS}"
DESIGN_SUFFIX = "" if DESIGN == "A" else f"_design{DESIGN}"
RECAL_SUFFIX = "_recal" if RECAL else ""
AM_SUFFIX = "" if abs(AOVERM - 0.01) < 1e-12 else f"_am{AOVERM:g}"
VARIANT = f"{RECAL_SUFFIX}{AM_SUFFIX}"

# Which archive a population's rows, budgets and gravity-only reference arcs
# live in. The two confirmatory designs were the whole of this table until the
# wide-elliptic population was replicated under the same dynamics; that
# population is the one the geometry result is read from, so the control that
# asks whether the geometry crossing survives realistic dynamics has to be able
# to name it. Nothing else about the campaign changes: the budget records have
# the same designs/rows shape, the orbit directories carry the same sobolA_NNN
# names, and the archived budget trajectories already sit under
# r14_cases/<key>_beta_<b> and r14_raw/<key>_beta_<b> for every key here.
POPULATIONS = {
    "A":   ("r10_sobolA_baseline_truth_corrected.json",
            "r14_budget_pareto.json", "convergence"),
    "B":   ("r11_designB_rows.json",
            "r14_budget_pareto.json", "designB_convergence"),
    "OE":  ("r31_operational_elliptical_rows.json",
            "r31_budget_pareto_operational_elliptical.json",
            "stratum_operational_elliptical_convergence"),
    "OEU": ("r38_operational_elliptical_uncapped_rows.json",
            "r38_budget_pareto_operational_elliptical_uncapped.json",
            "stratum_operational_elliptical_uncapped_convergence"),
}
if DESIGN not in POPULATIONS:
    raise SystemExit(f"JCAMP_J2_DESIGN={DESIGN!r} is not a population this "
                     f"campaign knows: {', '.join(sorted(POPULATIONS))}")
_ROWS_NAME, _PARETO_NAME, _CONV_DIR = POPULATIONS[DESIGN]
SOURCE_ROWS = METRICS / _ROWS_NAME
SOURCE_PARETO = METRICS / _PARETO_NAME
PLAN = METRICS / f"rJ2_plan{SUFFIX}{POP_SUFFIX}{DESIGN_SUFFIX}{VARIANT}.json"
SCORE = METRICS / f"rJ2_score{SUFFIX}{POP_SUFFIX}{DESIGN_SUFFIX}{VARIANT}.json"


def log(msg: str) -> None:
    J.log_line(LOG, f"J2 {msg}")


# ------------------------------------------------------------------- worker
_EPHEM = {}


def ephemeris():
    """One DE440/MOON_PA table set per process, built on first use."""
    if "e" not in _EPHEM:
        import rev4_robustness_controls as r4
        _EPHEM["e"] = r4.build_ephemeris(J.DURATION + 7200.0)
    return _EPHEM["e"]


def full_force_rhs(model, args, degree_of):
    import rev4_robustness_controls as r4
    from lunaris.common.type_defs import SpacecraftProps
    rhs = r4.ExpandedForceRhs(model, args, ephemeris(), degree_of,
                              use_third_body=True, use_srp=True,
                              track_components=True)
    if abs(AOVERM - 0.01) > 1e-12:
        # Same mass, scaled area: only the SRP magnitude changes, and the
        # reference feels it exactly as the policies do.
        rhs.sc = SpacecraftProps(mass_kg=100.0, area_m2=100.0 * AOVERM,
                                 cr=1.3)
    return rhs


def _recalibrated_tol(row: dict, adopted: int, model) -> float:
    """Bisect the accuracy parameter on the altitude history the policy will
    actually see, so that it spends the declared budget under these dynamics.

    This is the same calibration the budget campaign performs, pointed at the
    full-dynamics reference rather than the gravity-only one. It needs that
    reference to exist, which it does: the reference is propagated first and is
    shared across budgets.
    """
    index = int(row["sobol_index"])
    _, ref_raw = case_paths(index, "reference", "tighter")
    if not ref_raw.exists():
        raise RuntimeError(f"recalibration needs the full-dynamics reference "
                           f"for orbit {index}; {ref_raw} is missing")
    Y = J.load_states(ref_raw)
    h_km = (np.linalg.norm(Y[:3], axis=0) - model.r_ref) / 1e3
    cal = J.calibrate_radial(adopted, float(row["hp_km"]), float(row["ha_km"]),
                             adopted, h_km,
                             BETA * int(row["n_critical"]) ** 2)
    return float(cal["tol"])


def case_paths(index: int, policy: str, level: str):
    """The reference is budget-independent and is shared across budgets."""
    # The reference is budget-independent, but it *is* dynamics-dependent, so a
    # different area-to-mass ratio needs its own reference while a
    # recalibrated policy can share the one already on disk.
    tag = "" if policy == "reference" else SUFFIX
    if policy == "radial":
        tag += RECAL_SUFFIX
    stem = f"{policy}{tag}{AM_SUFFIX}_{level}"
    case = f"J2{DESIGN_SUFFIX}_{index:03d}"
    return (CASE_ROOT / case / f"{stem}.json",
            RAW_ROOT / case / f"{stem}.npz")


# --------------------------------------------------------------------- plan
def command_select() -> int:
    """Freeze which orbits and which policies J2 replicates, before running."""
    if PLAN.exists():
        log("select: plan already frozen")
        return 0
    src = json.loads(SOURCE_ROWS.read_text(encoding="utf-8"))
    rows = src.get("rows") or src.get("orbits")
    by_index = {int(r["sobol_index"]): r for r in rows}
    pareto = json.loads(SOURCE_PARETO.read_text(encoding="utf-8"))
    par = {int(r["sobol_index"]): r
           for r in pareto["designs"][DESIGN]["rows"]}

    plan_rows = []
    for index in sorted(by_index)[:N_ORBITS]:
        row, prow = by_index[index], par[index]
        geom = row.get("design_point", row)
        budget = prow["budgets"][f"beta_{BETA:.2f}"]
        if budget.get("censored"):
            raise SystemExit(f"orbit {index} is censored at beta={BETA}; the "
                             "deterministic subset rule cannot be applied "
                             "without a declared substitution rule")
        plan_rows.append({
            "sobol_index": index,
            "name": row.get("name", f"sobolA_{index:03d}"),
            "hp_km": float(geom["hp_km"]), "ha_km": float(geom["ha_km"]),
            "incl_deg": float(geom["incl_deg"]),
            "eccentricity": float(geom["eccentricity"]),
            "initial_state_si": [float(v) for v in geom["initial_state_si"]],
            "adopted_truth_degree": int(row["adopted_truth_degree"]),
            "n_critical": int(prow["n_critical"]),
            "constant_degree": int(budget["fixed"]["degree"]),
            "radial_tol_accel_m_s2": float(budget["atallah"]["tol_accel_m_s2"]),
            "radial_degree_archived": budget["atallah"]["degree"],
            "archived_beta_achieved_radial":
                budget["atallah"].get("beta_achieved"),
            "archived_beta_achieved_constant":
                budget["fixed"].get("beta_achieved"),
        })

    import rev4_robustness_controls as r4
    plan = {
        "schema": "rJ2_plan_v1", "created_utc": J.utc_now(),
        "status": "frozen_before_any_J2_trajectory_was_propagated",
        "question": ("does the force--trajectory reversal reproduce when the "
                     "same policies are propagated under realistic lunar "
                     "dynamics rather than the gravity-only system?"),
        "population_rule": (f"the first {N_ORBITS} orbits of the frozen "
                            "confirmatory design by Sobol index; a rule fixed "
                            "before any J2 result existed, so the subset "
                            "cannot have been chosen on outcome. Enlarging N "
                            "extends the same nested subset and never changes "
                            "which orbits the smaller one held"),
        "policy_rule": ("policies are read from the frozen budget records at "
                        f"beta={BETA:.2f} and are not recalibrated; the work "
                        "they realize under the new dynamics is measured"),
        "dynamics": {
            "gravity": "JGGRX_1800F spherical harmonics, degree set by policy",
            "orientation": "DE440 / MOON_PA, 60 s table",
            "third_body": ["Earth", "Sun"],
            "srp": "cannonball with lunar eclipse, Cr=1.3, A=1 m2, m=100 kg",
            "ephemeris_start_utc": r4.START_UTC,
            "common_to_all_policies": True,
        },
        "arc": {"duration_s": J.DURATION, "output_step_s": J.OUTPUT_STEP,
                "max_step_s": J.MAX_STEP},
        "levels": {k: {"rtol": J.LEVELS[k]["rtol"],
                       "atol_position_m": J.LEVELS[k]["atol_position_m"],
                       "atol_velocity_m_s": J.LEVELS[k]["atol_velocity_m_s"]}
                   for k in LEVEL_NAMES},
        "declared_outcomes": {
            "A_reversal_survives": ("the sign of J_traj(radial) - J_traj("
                                    "constant) under full dynamics agrees with "
                                    "the archived gravity-only sign on the "
                                    "majority of resolved orbits"),
            "B_reversal_does_not_survive": "the sign flips on the majority",
            "C_unresolved": "fewer than half the orbits resolve",
        },
        "source": {"rows": SOURCE_ROWS.name,
                   "pareto": SOURCE_PARETO.name,
                   "pareto_preregistration_sha256":
                       pareto.get("preregistration_sha256")},
        "rows": plan_rows,
        "provenance": J.provenance(),
    }
    plan["plan_sha256"] = J.object_hash(plan)
    J.atomic_json(PLAN, plan)
    deep = sum(1 for r in plan_rows if r["adopted_truth_degree"] > 300)
    log(f"select: {len(plan_rows)} orbits, {deep} above degree 300, "
        f"sha {plan['plan_sha256'][:16]}")
    return 0


# ---------------------------------------------------------------------- run
def run_task(task: dict) -> dict:
    row, policy, level = task["row"], task["policy"], task["level"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        model, args = J.model_for(adopted)
        if policy == "reference":
            n = adopted
            degree_of = lambda t, h, _n=n: _n
            spec = {"kind": "constant_reference", "degree": n}
        elif policy == "constant":
            n = int(row["constant_degree"])
            degree_of = lambda t, h, _n=n: _n
            spec = {"kind": "constant_budget", "degree": n, "beta": BETA,
                    "n_critical": int(row["n_critical"]),
                    "source": "frozen budget record, not recalibrated"}
        else:
            tol = float(row["radial_tol_accel_m_s2"])
            if RECAL:
                tol = _recalibrated_tol(row, adopted, model)
            degree_of, table = at.atallah_binned_schedule(
                model, J.atallah_g(adopted), tol, float(row["hp_km"]),
                float(row["ha_km"]), floor=J.FLOOR, cap=adopted,
                bin_km=J.BIN_KM)
            spec = {"kind": "budget_calibrated_radial", "beta": BETA,
                    "accuracy_parameter_m_s2": tol, "bin_km": J.BIN_KM,
                    "floor": J.FLOOR, "cap": adopted,
                    "n_critical": int(row["n_critical"]),
                    "source": ("recalibrated on the full-dynamics reference"
                               if RECAL else
                               "frozen budget record, not recalibrated"),
                    "table": {str(k): int(v) for k, v in sorted(table.items())}}

        y0 = np.asarray(row["initial_state_si"], dtype=float)
        grid = J.out_grid()
        lv = J.LEVELS[level]
        rhs = full_force_rhs(model, args, degree_of)
        Y, rhs, info = rc.propagate_instr(model, y0, J.DURATION, grid,
                                          degree_of, args, lv["rtol"],
                                          lv["atol"], max_step=J.MAX_STEP,
                                          rhs_obj=rhs)
        cj, cn = case_paths(index, policy, level)
        J.atomic_npz(cn, t=grid, y=Y)
        telem = rhs.info(info["wall_s"])
        record = {
            "schema": "rJ2_case_v1", "created_utc": J.utc_now(),
            "config": {"sobol_index": index, "policy": policy,
                       "policy_spec": spec, "level": level,
                       "adopted_truth_degree": adopted,
                       "dynamics": "moon_sh + earth + sun + srp_eclipse",
                       "orientation": "DE440/MOON_PA",
                       "duration_s": J.DURATION,
                       "output_step_s": J.OUTPUT_STEP,
                       "max_step_s": J.MAX_STEP, "rtol": lv["rtol"],
                       "atol_kind": "vector",
                       "atol_position_m": lv["atol_position_m"],
                       "atol_velocity_m_s": lv["atol_velocity_m_s"],
                       "integrator": "InstrumentedDOP853",
                       "field": J.field_key(),
                       "timing_comparable": False,
                       "execution": "parallel_process_pool"},
            "telemetry": {**telem,
                          "n_accepted_steps": int(info["n_accepted_steps"]),
                          "n_attempted_steps": int(info["n_attempted_steps"]),
                          "n_rejected_trials": int(info["n_rejected_trials"]),
                          "total_wall_s": float(info["wall_s"])},
            "raw_path": str(cn), "raw_sha256": J.sha256_file(cn),
        }
        J.atomic_json(cj, record)
        return {"ok": True, "index": index, "policy": policy, "level": level,
                "wall_s": float(info["wall_s"])}
    except Exception:
        return {"ok": False, "index": index, "policy": policy, "level": level,
                "error": traceback.format_exc()}


def command_run(workers: int, deadline: str | None) -> int:
    if not PLAN.exists():
        command_select()
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    tasks = []
    for row in plan["rows"]:
        for level in LEVEL_NAMES:
            for policy in POLICIES:
                cj, cn = case_paths(int(row["sobol_index"]), policy, level)
                if cj.exists() and cn.exists():
                    continue
                tasks.append({"row": row, "policy": policy, "level": level})
    log(f"run: {len(tasks)} trajectories to run")
    if not tasks:
        return 0
    stop = None
    if deadline:
        from datetime import datetime as _dt
        stop = _dt.fromisoformat(deadline)
    t0 = time.time()
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(run_task, t) for t in tasks]
        for fut in as_completed(futs):
            res = fut.result()
            if res["ok"]:
                ok += 1
            else:
                fail += 1
                log(f"run FAIL {res['index']} {res['policy']} {res['level']}\n"
                    f"{res['error']}")
            if (ok + fail) % 8 == 0:
                log(f"run: {ok + fail}/{len(tasks)} "
                    f"({(time.time() - t0) / 60:.1f} min)")
            if stop is not None:
                from datetime import datetime as _dt
                if _dt.now() > stop:
                    log("run: deadline reached, cancelling remainder")
                    for f in futs:
                        f.cancel()
                    break
    log(f"run: {ok} ok, {fail} failed, {(time.time() - t0) / 60:.1f} min")
    # A completion record, so that "did this stage finish?" is a question about
    # one file rather than about counting case files in a directory several
    # stages write into.
    J.atomic_json(
        METRICS / f"rJ2_run_complete{SUFFIX}{POP_SUFFIX}{DESIGN_SUFFIX}{VARIANT}.json", {
        "schema": "rJ2_run_complete_v1", "created_utc": J.utc_now(),
        "beta": BETA, "orbits": N_ORBITS, "design": DESIGN,
        "recalibrated_under_full_dynamics": RECAL,
        "area_to_mass_m2_per_kg": AOVERM,
        "plan_sha256": plan["plan_sha256"],
        "trajectories_run": ok, "failures": fail, "complete": fail == 0})
    return 0 if fail == 0 else 1


# -------------------------------------------------------------------- score
def _force_defect_moonpa(row: dict, level: str = "tighter") -> dict:
    """Truncation defect along the full-dynamics reference, in the MOON_PA
    frame. Third-body and SRP terms are common and cancel exactly, so only the
    spherical-harmonic difference survives and no ephemeris term enters."""
    from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba
    index = int(row["sobol_index"])
    adopted = int(row["adopted_truth_degree"])
    model, args = J.model_for(adopted)
    ephem = ephemeris()
    _, ref_raw = case_paths(index, "reference", level)
    with np.load(ref_raw) as z:
        t_ref, Y_ref = z["t"], z["y"]
    h_km = (np.linalg.norm(Y_ref[:3], axis=0) - model.r_ref) / 1e3

    tol = float(row["radial_tol_accel_m_s2"])
    if RECAL:
        tol = _recalibrated_tol(row, adopted, model)
    degree_of, table = at.atallah_binned_schedule(
        model, J.atallah_g(adopted), tol,
        float(row["hp_km"]), float(row["ha_km"]), floor=J.FLOOR, cap=adopted,
        bin_km=J.BIN_KM)
    table = {float(k): int(v) for k, v in table.items()}
    deg = {"radial": J.degrees_from_table(table, h_km),
           "constant": np.full(len(h_km), int(row["constant_degree"]),
                               dtype=int)}
    def to_fixed(t, r_i):
        rb = np.empty(3)
        ephem.transform_inertial_to_fixed(float(t), r_i, out=rb)
        return rb

    out = J.force_defects(t_ref, Y_ref[:3], deg, adopted, args,
                          to_fixed=to_fixed)
    out["realized_mean_degree_sq"] = {
        "radial": float(np.mean(deg["radial"].astype(float) ** 2)),
        "constant": float(int(row["constant_degree"]) ** 2)}
    out["realized_beta"] = {
        k: v / float(int(row["n_critical"]) ** 2)
        for k, v in out["realized_mean_degree_sq"].items()}
    return out


def _archived_gravity_only(index: int) -> dict | None:
    """The paired gravity-only result for the same orbit, as archived."""
    conv = _CONV_DIR
    base_dir = (METRICS / "r14_cases" / f"{DESIGN}_beta_{BETA:.2f}"
                / f"sobolA_{index:03d}")
    ref_raw = (METRICS / "r11_raw" / conv / f"sobolA_{index:03d}"
               / "truth_tighter.npz")
    rad_raw = (METRICS / "r14_raw" / f"{DESIGN}_beta_{BETA:.2f}"
               / f"sobolA_{index:03d}" / "atallah_budget_tighter.npz")
    # At beta = 1 the constant comparator is the critical degree itself and the
    # budget campaign reuses the archived run rather than repropagating it; at
    # any other budget it has its own file.
    fix_raw = (METRICS / "r14_raw" / f"{DESIGN}_beta_{BETA:.2f}"
               / f"sobolA_{index:03d}" / "fixed_budget_tighter.npz")
    if not fix_raw.exists():
        fix_raw = (METRICS / "r11_raw" / conv / f"sobolA_{index:03d}"
                   / "fixed_critical_tighter.npz")
    if not (ref_raw.exists() and rad_raw.exists() and fix_raw.exists()):
        return None
    T = J.load_states(ref_raw)
    R = J.load_states(rad_raw)
    F = J.load_states(fix_raw)
    n = min(T.shape[1], R.shape[1], F.shape[1])
    e_r = J.trajectory_error(R[:, :n], T[:, :n])["J_traj_rms_m"]
    e_f = J.trajectory_error(F[:, :n], T[:, :n])["J_traj_rms_m"]
    out = {"J_traj_rms_m_radial": e_r, "J_traj_rms_m_constant": e_f,
           "rho_traj": e_r / e_f, "case_dir": str(base_dir), "resolved": None}
    # Whether the *gravity-only* comparison was itself resolved decides how
    # much a sign agreement is worth. Comparing a resolved full-dynamics sign
    # against an archived tie says nothing, so the archived resolution verdict
    # is carried rather than assumed. It is read from the re-scoring record,
    # which computes it under the same rule this campaign uses.
    comp = METRICS / "rJ_field_comparison.json"
    if comp.exists():
        for row in json.loads(comp.read_text(encoding="utf-8"))["rows"]:
            if row["design"] == DESIGN and int(row["index"]) == index:
                out["resolved"] = bool(row["resolved"])
                break
    if out["resolved"] is None:
        # rJ_field_comparison covers the two confirmatory designs only. Any
        # other population carries the same verdict, computed under the same
        # rule, in its own budget-campaign record, and reading it there is what
        # keeps a missing lookup from being mistaken for an unresolved
        # comparison: without this the wide-elliptic run scored C_unresolved on
        # zero orbits resolved in both dynamics while fifteen of them were in
        # fact resolved on each side.
        own = METRICS / f"r14_trajectory_{DESIGN}_beta_{BETA:.2f}.json"
        if own.exists():
            try:
                rows = json.loads(own.read_text(encoding="utf-8"))["rows"]
            except (ValueError, OSError, KeyError):
                rows = []
            for row in rows:
                if int(row.get("sobol_index", -1)) == index:
                    c = row.get("comparison") or {}
                    if "resolved" in c:
                        out["resolved"] = bool(c["resolved"])
                    break
    return out


def _budget_parity(rows: list[dict], width: float = 0.10) -> dict:
    """Does the conclusion depend on the orbits where the archived policy drifts
    off its declared budget under the new dynamics?

    Declared here as what it is: a post-hoc subgroup, specified after the
    realized budgets were seen. It is reported because the drift is real and a
    reader is entitled to know whether the sign agreement rests on the orbits
    where the comparison is no longer at equal work -- not because the subgroup
    is the primary analysis. It is not.
    """
    keep = [r for r in rows
            if abs(r["realized_beta"]["radial"] - BETA) <= width * BETA]
    return {
        "status": "post_hoc; specified after the realized budgets were seen",
        "rule": f"|realized beta - {BETA:.2f}| <= {width:.0%} of the declared "
                "budget, on the orbits resolved under both dynamics",
        "orbits": len(keep), "of": len(rows),
        "sign_agrees": sum(1 for r in keep
                           if r["sign_agrees_with_gravity_only"]),
        "reversal": sum(1 for r in keep if r["reversal"]),
        "median_rho_traj": (float(np.median([r["rho_traj"] for r in keep]))
                            if keep else None),
        "median_rho_force": (float(np.median([r["rho_force"] for r in keep]))
                             if keep else None),
    }


def command_score() -> int:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    out = []
    for row in plan["rows"]:
        index = int(row["sobol_index"])
        Y, missing = {}, False
        for policy in POLICIES:
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
        defect = _force_defect_moonpa(row)
        self_ref = J.self_difference(Y[("reference", "tight")],
                                     Y[("reference", "tighter")])
        rec = {"sobol_index": index, "name": row["name"],
               "hp_km": row["hp_km"], "ha_km": row["ha_km"],
               "incl_deg": row["incl_deg"],
               "n_critical": row["n_critical"],
               "adopted_truth_degree": row["adopted_truth_degree"],
               "self_difference_reference_m": self_ref,
               "realized_beta": defect["realized_beta"]}
        for policy in ("constant", "radial"):
            err = {lv: J.trajectory_error(Y[(policy, lv)],
                                          Y[("reference", lv)])["J_traj_rms_m"]
                   for lv in LEVEL_NAMES}
            self_p = J.self_difference(Y[(policy, "tight")],
                                       Y[(policy, "tighter")])
            rec[policy] = {"J_traj_rms_m": err, "self_difference_m": self_p,
                           "numerical_envelope_m": self_p + self_ref,
                           **defect[policy]}
        for lv in LEVEL_NAMES:
            e_c = rec["constant"]["J_traj_rms_m"][lv]
            e_r = rec["radial"]["J_traj_rms_m"][lv]
            rec.setdefault("resolved", {})[lv] = bool(J.resolved(
                e_c, e_r, rec["constant"]["numerical_envelope_m"],
                rec["radial"]["numerical_envelope_m"]))
        rec["rho_traj"] = (rec["radial"]["J_traj_rms_m"]["tighter"]
                           / rec["constant"]["J_traj_rms_m"]["tighter"])
        rec["rho_force"] = (rec["radial"]["J_force_mean_m_s2"]
                            / rec["constant"]["J_force_mean_m_s2"])
        rec["reversal"] = bool(rec["rho_force"] < 1.0 < rec["rho_traj"])
        rec["gravity_only_archived"] = _archived_gravity_only(index)
        if rec["gravity_only_archived"]:
            rec["sign_agrees_with_gravity_only"] = bool(
                (rec["rho_traj"] > 1.0)
                == (rec["gravity_only_archived"]["rho_traj"] > 1.0))
        out.append(rec)

    res = [r for r in out if r["resolved"]["tighter"]]
    paired = [r for r in res if r.get("sign_agrees_with_gravity_only")
              is not None]
    # The pairing that carries weight is the one where *both* comparisons are
    # resolved; the rest are a resolved sign held against an archived tie.
    doubly = [r for r in paired
              if (r["gravity_only_archived"] or {}).get("resolved")]
    n_agree = sum(1 for r in paired if r["sign_agrees_with_gravity_only"])
    n_agree_doubly = sum(1 for r in doubly if r["sign_agrees_with_gravity_only"])
    n_rev = sum(1 for r in res if r["reversal"])
    verdict = ("C_unresolved" if len(res) < 0.5 * len(out) or not doubly
               else ("A_reversal_survives"
                     if n_agree_doubly > len(doubly) / 2
                     else "B_reversal_does_not_survive"))
    payload = {
        "schema": "rJ2_score_v1", "created_utc": J.utc_now(),
        "plan_sha256": plan["plan_sha256"], "beta": BETA,
        "counts": {"orbits": len(out), "resolved_tighter": len(res),
                   "paired_with_gravity_only": len(paired),
                   "sign_agrees": n_agree,
                   "resolved_in_both_dynamics": len(doubly),
                   "sign_agrees_resolved_in_both": n_agree_doubly,
                   "reversal_orbits": n_rev,
                   "radial_wins_force": sum(1 for r in res
                                            if r["rho_force"] < 1.0),
                   "radial_loses_trajectory": sum(1 for r in res
                                                  if r["rho_traj"] > 1.0)},
        "aggregates": {
            "median_rho_traj": float(np.median([r["rho_traj"] for r in out])),
            "median_rho_force": float(np.median([r["rho_force"] for r in out])),
            "median_rho_traj_gravity_only": float(np.median(
                [r["gravity_only_archived"]["rho_traj"] for r in out
                 if r.get("gravity_only_archived")])),
            "max_realized_beta_deviation_radial": float(max(
                abs(r["realized_beta"]["radial"] - BETA) for r in out)),
            "median_realized_beta_radial": float(np.median(
                [r["realized_beta"]["radial"] for r in out])),
            "orbits_within_10pct_of_declared_budget": sum(
                1 for r in out
                if abs(r["realized_beta"]["radial"] - BETA) <= 0.10 * BETA),
        },
        "budget_caveat": (
            "the policy is the archived one and is not recalibrated, so under "
            "the new dynamics it no longer spends exactly the declared budget: "
            "the altitude history it reads is not the gravity-only history it "
            "was calibrated on. The realized budget is reported per orbit "
            "rather than corrected, because recalibrating would change the "
            "policy and this campaign is a replication of a policy, not of a "
            "budget"),
        "verdict": verdict,
        "budget_parity_subgroup": _budget_parity(doubly),
        "rows": out, "provenance": J.provenance(),
    }
    J.atomic_json(SCORE, payload)
    log(f"score: {verdict}; resolved {len(res)}/{len(out)}, sign agreement "
        f"{n_agree_doubly}/{len(doubly)} where both dynamics resolve "
        f"({n_agree}/{len(paired)} counting archived ties), reversal {n_rev}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("select", "run", "score", "status"))
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    if a.command == "select":
        return command_select()
    if a.command == "run":
        return command_run(a.workers, a.deadline)
    if a.command == "score":
        return command_score()
    for p in (PLAN, SCORE):
        print(f"{p.name:26s} {'present' if p.exists() else 'MISSING'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
