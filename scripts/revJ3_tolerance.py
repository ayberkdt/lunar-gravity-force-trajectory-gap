"""J3: is the ranking a property of the orbits, or of the resolution criterion?

The manuscript decides each policy comparison with a two-level numerical
envelope: a comparison counts only when the error gap exceeds what the two
policies' own between-level self-differences can explain. That rule is defended
at length, but it is still *a* rule, and a reviewer is entitled to ask what the
counts would look like under a different one.

Answering it by re-running the whole population at three tolerances would cost
more compute than the question is worth. What the question actually needs is
the cases where the rule could plausibly change its mind, so this campaign
takes a deterministic sample built around the margin

    M_res = |E_constant - E_radial| / (envelope_constant + envelope_radial),

the quantity the rule thresholds at 1. Orbits are drawn by rank, not by
outcome, and the rule is fixed here before the new tolerances are run:

    8 orbits with M_res closest to 1 in log -- the boundary, where refinement
      has the best chance of moving a verdict;
    4 with the smallest M_res and 4 with the largest -- the two tails, so the
      sample is not only boundary cases.

Each is then propagated at two *new* tolerance levels, one looser and one
tighter than the archived pair, giving four levels spanning three decades of
rtol. The archived tight and tighter runs are reused unchanged, so the four
levels are the same trajectories the manuscript reports plus two extensions.

What is being tested is not the integrator's order. It is one question: does

    rho = E_radial / E_constant

stay on the same side of 1 at every level?

Usage:
    python revJ3_tolerance.py select
    python revJ3_tolerance.py run --workers 11
    python revJ3_tolerance.py score
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
                               r"D:\makale_raw_offload\jgcd")) / "J3"
CASE_ROOT = METRICS / "rJ3_cases"
LOG = Path(__file__).resolve().parent / "rJ3_campaign.log"

SOURCE_PARETO = METRICS / "r14_budget_pareto.json"

# Which archived population the control is drawn from. Read from the
# environment so parent and spawned workers agree; argv does not survive a
# Windows process spawn.
DESIGN = os.environ.get("JCAMP_J3_DESIGN", "A").upper()
DESIGN_SUFFIX = "" if DESIGN == "A" else f"_design{DESIGN}"
SOURCE_ROWS = (METRICS / "r10_sobolA_baseline_truth_corrected.json"
               if DESIGN == "A" else METRICS / "r11_designB_rows.json")
CONV_RAW = (METRICS / "r11_raw" / ("convergence" if DESIGN == "A"
                                   else "designB_convergence"))

PLAN = METRICS / f"rJ3_plan{DESIGN_SUFFIX}.json"
SCORE = METRICS / f"rJ3_score{DESIGN_SUFFIX}.json"

BETA = 1.00
ARCHIVED_LEVELS = ("tight", "tighter")
NEW_LEVELS = ("loose", "tightest")
ALL_LEVELS = ("loose", "tight", "tighter", "tightest")
N_BOUNDARY, N_LOW, N_HIGH = 8, 4, 4

R11_RAW = CONV_RAW
R14_RAW = METRICS / "r14_raw" / f"{DESIGN}_beta_{BETA:.2f}"


def log(msg: str) -> None:
    J.log_line(LOG, f"J3 {msg}")


# --------------------------------------------------------- archived reuse
def archived_raw(index: int, policy: str, level: str) -> Path | None:
    """Where the manuscript's own run of this case lives, if it exists."""
    name = f"sobolA_{index:03d}"
    if policy == "reference":
        p = R11_RAW / name / f"truth_{level}.npz"
    elif policy == "radial":
        p = R14_RAW / name / f"atallah_budget_{level}.npz"
    else:
        p = R14_RAW / name / f"fixed_budget_{level}.npz"
        if not p.exists():
            # At beta = 1 the constant comparator is the critical degree
            # itself, and R14 reuses the archived R11 run rather than
            # repropagating an identical trajectory.
            p = R11_RAW / name / f"fixed_critical_{level}.npz"
    return p if p.exists() else None


def case_paths(index: int, policy: str, level: str):
    stem = f"J3{DESIGN_SUFFIX}_{index:03d}"
    return (CASE_ROOT / stem / f"{policy}_{level}.json",
            RAW_ROOT / stem / f"{policy}_{level}.npz")


def load_states(index: int, policy: str, level: str) -> np.ndarray | None:
    if level in ARCHIVED_LEVELS:
        p = archived_raw(index, policy, level)
    else:
        _, p = case_paths(index, policy, level)
        p = p if p.exists() else None
    if p is None:
        return None
    return J.load_states(p)


# ------------------------------------------------------------------ select
def command_select() -> int:
    if PLAN.exists():
        log("select: plan already frozen")
        return 0
    src = json.loads(SOURCE_ROWS.read_text(encoding="utf-8"))
    by_index = {int(r["sobol_index"]): r for r in src["rows"]}
    pareto = json.loads(SOURCE_PARETO.read_text(encoding="utf-8"))
    par = {int(r["sobol_index"]): r
           for r in pareto["designs"][DESIGN]["rows"]}

    scored = []
    for index in sorted(by_index):
        budget = par[index]["budgets"].get(f"beta_{BETA:.2f}")
        if budget is None or budget.get("censored"):
            continue
        Y = {}
        ok = True
        for policy in ("reference", "constant", "radial"):
            for level in ARCHIVED_LEVELS:
                y = load_states(index, policy, level)
                if y is None:
                    ok = False
                    break
                Y[(policy, level)] = y
            if not ok:
                break
        if not ok:
            continue
        n = min(y.shape[1] for y in Y.values())
        self_ref = J.self_difference(Y[("reference", "tight")][:, :n],
                                    Y[("reference", "tighter")][:, :n])
        env, err = {}, {}
        for policy in ("constant", "radial"):
            self_p = J.self_difference(Y[(policy, "tight")][:, :n],
                                       Y[(policy, "tighter")][:, :n])
            env[policy] = self_p + self_ref
            err[policy] = J.trajectory_error(
                Y[(policy, "tighter")][:, :n],
                Y[("reference", "tighter")][:, :n])["J_traj_rms_m"]
        margin = abs(err["constant"] - err["radial"]) / (env["constant"]
                                                         + env["radial"])
        scored.append({
            "sobol_index": index,
            "name": by_index[index].get("name", f"sobolA_{index:03d}"),
            "M_res": float(margin),
            "rho_traj_archived": float(err["radial"] / err["constant"]),
            "E_constant_m": float(err["constant"]),
            "E_radial_m": float(err["radial"]),
            "envelope_constant_m": float(env["constant"]),
            "envelope_radial_m": float(env["radial"]),
        })
    if len(scored) < N_BOUNDARY + N_LOW + N_HIGH:
        raise SystemExit(f"only {len(scored)} orbits have a complete archived "
                         "two-level record; cannot draw the declared sample")

    by_margin = sorted(scored, key=lambda r: r["M_res"])
    boundary = sorted(scored,
                      key=lambda r: abs(math.log(max(r["M_res"], 1e-12))))
    chosen, roles = {}, {}
    for r in boundary[:N_BOUNDARY]:
        chosen[r["sobol_index"]] = r
        roles[r["sobol_index"]] = "boundary"
    for r in by_margin:
        if len([k for k in roles if roles[k] == "low"]) >= N_LOW:
            break
        if r["sobol_index"] not in chosen:
            chosen[r["sobol_index"]] = r
            roles[r["sobol_index"]] = "low"
    for r in reversed(by_margin):
        if len([k for k in roles if roles[k] == "high"]) >= N_HIGH:
            break
        if r["sobol_index"] not in chosen:
            chosen[r["sobol_index"]] = r
            roles[r["sobol_index"]] = "high"

    rows = []
    for index in sorted(chosen):
        row, prow = by_index[index], par[index]
        geom = row.get("design_point", row)
        budget = prow["budgets"][f"beta_{BETA:.2f}"]
        rows.append({**chosen[index], "role": roles[index],
                     "hp_km": float(geom["hp_km"]),
                     "ha_km": float(geom["ha_km"]),
                     "incl_deg": float(geom["incl_deg"]),
                     "initial_state_si": [float(v) for v in
                                          geom["initial_state_si"]],
                     "adopted_truth_degree": int(row["adopted_truth_degree"]),
                     "n_critical": int(prow["n_critical"]),
                     "constant_degree": int(budget["fixed"]["degree"]),
                     "radial_tol_accel_m_s2":
                         float(budget["atallah"]["tol_accel_m_s2"])})

    plan = {
        "schema": "rJ3_plan_v1", "created_utc": J.utc_now(),
        "status": "frozen_before_any_new_tolerance_level_was_run",
        "population": f"archived confirmatory design {DESIGN}",
        "question": ("does the sign of rho = E_radial / E_constant survive a "
                     "change of integration tolerance across four levels?"),
        "selection_rule": {
            "statistic": "M_res = |E_constant - E_radial| / (env_c + env_r), "
                         "evaluated on the archived two-level records",
            "boundary": f"{N_BOUNDARY} orbits with M_res closest to 1 in log",
            "tails": f"{N_LOW} smallest and {N_HIGH} largest M_res",
            "note": "drawn by rank on an archived statistic, so the sample is "
                    "fixed before any new trajectory exists",
        },
        "levels": {k: {"rtol": J.LEVELS[k]["rtol"],
                       "atol_position_m": J.LEVELS[k]["atol_position_m"],
                       "atol_velocity_m_s": J.LEVELS[k]["atol_velocity_m_s"],
                       "source": ("archived" if k in ARCHIVED_LEVELS
                                  else "run by this campaign")}
                   for k in ALL_LEVELS},
        "declared_outcomes": {
            "A_sign_stable": "rho stays on the same side of 1 at every level "
                             "for the large majority of orbits",
            "B_sign_moves": "rho crosses 1 with tolerance on a material "
                            "fraction, and the two-level counts are therefore "
                            "criterion-dependent",
        },
        "beta": BETA, "rows": rows, "provenance": J.provenance(),
    }
    plan["plan_sha256"] = J.object_hash(plan)
    J.atomic_json(PLAN, plan)
    log(f"select: {len(rows)} orbits "
        f"(M_res {min(r['M_res'] for r in rows):.2f}-"
        f"{max(r['M_res'] for r in rows):.2f}), sha {plan['plan_sha256'][:16]}")
    return 0


# --------------------------------------------------------------------- run
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
            spec = {"kind": "constant_budget", "degree": n, "beta": BETA}
        else:
            tol = float(row["radial_tol_accel_m_s2"])
            degree_of, table = at.atallah_binned_schedule(
                model, J.atallah_g(adopted), tol, float(row["hp_km"]),
                float(row["ha_km"]), floor=J.FLOOR, cap=adopted,
                bin_km=J.BIN_KM)
            spec = {"kind": "budget_calibrated_radial", "beta": BETA,
                    "accuracy_parameter_m_s2": tol,
                    "source": "frozen budget record, not recalibrated"}
        y0 = np.asarray(row["initial_state_si"], dtype=float)
        grid = J.out_grid()
        lv = J.LEVELS[level]
        Y, rhs, info = rc.propagate_instr(model, y0, J.DURATION, grid,
                                          degree_of, args, lv["rtol"],
                                          lv["atol"], max_step=J.MAX_STEP)
        cj, cn = case_paths(index, policy, level)
        J.atomic_npz(cn, t=grid, y=Y)
        J.atomic_json(cj, {
            "schema": "rJ3_case_v1", "created_utc": J.utc_now(),
            "config": {"sobol_index": index, "policy": policy,
                       "policy_spec": spec, "level": level,
                       "adopted_truth_degree": adopted,
                       "duration_s": J.DURATION,
                       "output_step_s": J.OUTPUT_STEP,
                       "max_step_s": J.MAX_STEP, "rtol": lv["rtol"],
                       "atol_kind": "vector",
                       "atol_position_m": lv["atol_position_m"],
                       "atol_velocity_m_s": lv["atol_velocity_m_s"],
                       "integrator": "InstrumentedDOP853",
                       "field": J.field_key(),
                       "timing_comparable": False},
            "telemetry": {"n_rhs": int(info["n_rhs"]),
                          "n_accepted_steps": int(info["n_accepted_steps"]),
                          "n_attempted_steps": int(info["n_attempted_steps"]),
                          "n_rejected_trials": int(info["n_rejected_trials"]),
                          "mean_degree_sq": float(rhs.sum_deg_sq / rhs.n_calls),
                          "total_wall_s": float(info["wall_s"])},
            "raw_path": str(cn), "raw_sha256": J.sha256_file(cn)})
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
        for level in NEW_LEVELS:
            for policy in ("reference", "constant", "radial"):
                cj, cn = case_paths(int(row["sobol_index"]), policy, level)
                if cj.exists() and cn.exists():
                    continue
                tasks.append({"row": row, "policy": policy, "level": level})
    log(f"run: {len(tasks)} trajectories to run")
    if not tasks:
        return 0
    # The looser level first: it is cheap, and it exercises every code path
    # before the expensive level commits the machine to it.
    tasks.sort(key=lambda t: NEW_LEVELS.index(t["level"]))
    stop = None
    if deadline:
        from datetime import datetime as _dt
        stop = _dt.fromisoformat(deadline)
    t0 = time.time()
    ok = fail = 0
    failures = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(run_task, t) for t in tasks]
        for fut in as_completed(futs):
            res = fut.result()
            if res["ok"]:
                ok += 1
            else:
                fail += 1
                failures.append({k: res[k] for k in
                                 ("index", "policy", "level")})
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
    J.atomic_json(METRICS / f"rJ3_run_complete{DESIGN_SUFFIX}.json", {
        "schema": "rJ3_run_complete_v1", "created_utc": J.utc_now(),
        "design": DESIGN, "plan_sha256": plan["plan_sha256"],
        "trajectories_run": ok, "failures": failures,
        "complete": fail == 0})
    return 0 if fail == 0 else 1


# ------------------------------------------------------------------- score
LEVEL_PAIRS = (("loose", "tight"), ("tight", "tighter"),
               ("tighter", "tightest"))
MANUSCRIPT_PAIR = ("tight", "tighter")


def _conditional(rows: list[dict]) -> dict:
    """The question behind the question.

    The counts the manuscript reports are counts over *resolved* comparisons,
    so what a reader wants to know is not whether every rho in the sample keeps
    its side of 1 -- the rule already says the low-margin ones will not -- but
    whether the comparisons the rule admits keep theirs. Reported by margin
    stratum as well, because the sample was drawn by stratum and a rate over
    the pooled sample would hide which end of the margin range moved.
    """
    by_role = {}
    for role in ("boundary", "low", "high"):
        sub = [r for r in rows if r["role"] == role and "rho_range" in r]
        if not sub:
            continue
        by_role[role] = {
            "orbits": len(sub),
            "sign_stable_across_all_levels":
                sum(1 for r in sub if r["all_same_side_of_one"]),
            "median_M_res_archived": float(np.median(
                [r["M_res_archived"] for r in sub])),
            "median_rho_spread_dex": float(np.median(
                [r["rho_spread_dex"] for r in sub])),
        }
    admitted = [r for r in rows if r.get("resolved_in_manuscript_pair")]
    agree = [r for r in admitted if r["sign_agrees_across_resolved_pairs"]]
    pair_counts = {}
    for lo, hi in LEVEL_PAIRS:
        key = f"{lo}_{hi}"
        live = [r["pairs"][key] for r in rows if key in r.get("pairs", {})]
        if live:
            pair_counts[key] = {
                "resolved": sum(1 for v in live if v["resolved"]),
                "of": len(live),
                "radial_worse_among_resolved":
                    sum(1 for v in live if v["resolved"] and v["radial_worse"]),
                "median_rho": float(np.median([v["rho"] for v in live])),
            }
    return {
        "question": ("do the comparisons the resolution rule admits keep their "
                     "sign when the rule is evaluated on a different pair of "
                     "tolerance levels?"),
        "by_margin_stratum": by_role,
        "resolution_rule_applied_to_each_adjacent_level_pair": pair_counts,
        "admitted_by_manuscript_pair": len(admitted),
        "of_those_sign_unchanged_in_every_pair_that_resolves": len(agree),
        "reading": ("the sample was drawn to include orbits the rule excludes, "
                    "so instability among them is the rule working, not the "
                    "ranking failing"),
    }


def _pair_verdicts(index: int) -> dict:
    """Re-apply the manuscript's resolution rule on each adjacent level pair.

    The rule needs two levels: the finer one supplies the errors, and the
    difference between the two supplies each policy's numerical envelope. With
    four levels there are three such pairs, so the question "would the reported
    counts change if the criterion had been evaluated on a different pair?" can
    be asked directly instead of inferred.
    """
    out = {}
    for lo, hi in LEVEL_PAIRS:
        Y = {(p, lv): load_states(index, p, lv)
             for p in ("reference", "constant", "radial") for lv in (lo, hi)}
        if any(v is None for v in Y.values()):
            continue
        n = min(v.shape[1] for v in Y.values())
        self_ref = J.self_difference(Y[("reference", lo)][:, :n],
                                     Y[("reference", hi)][:, :n])
        err, env = {}, {}
        for policy in ("constant", "radial"):
            err[policy] = J.trajectory_error(
                Y[(policy, hi)][:, :n], Y[("reference", hi)][:, :n]
            )["J_traj_rms_m"]
            env[policy] = self_ref + J.self_difference(
                Y[(policy, lo)][:, :n], Y[(policy, hi)][:, :n])
        gap = abs(err["constant"] - err["radial"])
        out[f"{lo}_{hi}"] = {
            "levels": [lo, hi],
            "E_constant_m": err["constant"], "E_radial_m": err["radial"],
            "envelope_sum_m": env["constant"] + env["radial"],
            "M_res": gap / (env["constant"] + env["radial"]),
            "resolved": bool(J.resolved(err["constant"], err["radial"],
                                        env["constant"], env["radial"])),
            "rho": err["radial"] / err["constant"],
            "radial_worse": bool(err["radial"] > err["constant"])}
    return out


def command_score() -> int:
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    out = []
    for row in plan["rows"]:
        index = int(row["sobol_index"])
        rec = {"sobol_index": index, "name": row["name"], "role": row["role"],
               "hp_km": row["hp_km"], "incl_deg": row["incl_deg"],
               "M_res_archived": row["M_res"],
               "rho_traj_archived": row["rho_traj_archived"],
               "by_level": {}, "levels_available": []}
        for level in ALL_LEVELS:
            Y = {p: load_states(index, p, level)
                 for p in ("reference", "constant", "radial")}
            if any(v is None for v in Y.values()):
                continue
            n = min(v.shape[1] for v in Y.values())
            e_c = J.trajectory_error(Y["constant"][:, :n],
                                     Y["reference"][:, :n])["J_traj_rms_m"]
            e_r = J.trajectory_error(Y["radial"][:, :n],
                                     Y["reference"][:, :n])["J_traj_rms_m"]
            rec["by_level"][level] = {
                "E_constant_m": e_c, "E_radial_m": e_r,
                "rho": e_r / e_c, "radial_worse": bool(e_r > e_c),
                "rtol": J.LEVELS[level]["rtol"]}
            rec["levels_available"].append(level)
        rhos = [rec["by_level"][lv]["rho"] for lv in rec["levels_available"]]
        if rhos:
            rec["all_same_side_of_one"] = bool(
                all(r > 1.0 for r in rhos) or all(r < 1.0 for r in rhos))
            rec["rho_range"] = [float(min(rhos)), float(max(rhos))]
            rec["rho_spread_dex"] = float(math.log10(max(rhos) / min(rhos)))
        rec["pairs"] = _pair_verdicts(index)
        pairs = rec["pairs"]
        res_pairs = [k for k, v in pairs.items() if v["resolved"]]
        rec["resolved_in_n_pairs"] = len(res_pairs)
        rec["resolved_in_manuscript_pair"] = bool(
            pairs.get("_".join(MANUSCRIPT_PAIR), {}).get("resolved"))
        if res_pairs:
            signs = {pairs[k]["radial_worse"] for k in res_pairs}
            rec["sign_agrees_across_resolved_pairs"] = bool(len(signs) == 1)
        else:
            rec["sign_agrees_across_resolved_pairs"] = None
        out.append(rec)

    complete = [r for r in out if len(r["levels_available"]) == len(ALL_LEVELS)]
    stable = [r for r in complete if r["all_same_side_of_one"]]
    payload = {
        "schema": "rJ3_score_v1", "created_utc": J.utc_now(),
        "plan_sha256": plan["plan_sha256"], "beta": BETA,
        "levels": list(ALL_LEVELS),
        "counts": {"orbits": len(out),
                   "with_all_four_levels": len(complete),
                   "rho_same_side_at_every_level": len(stable),
                   "by_role": {role: sum(1 for r in complete
                                         if r["role"] == role
                                         and r["all_same_side_of_one"])
                               for role in ("boundary", "low", "high")}},
        "aggregates": {
            "max_rho_spread_dex": float(max(
                (r["rho_spread_dex"] for r in complete), default=float("nan"))),
            "median_rho_by_level": {
                lv: float(np.median([r["by_level"][lv]["rho"]
                                     for r in complete]))
                for lv in ALL_LEVELS} if complete else {},
        },
        "verdict": ("A_sign_stable"
                    if complete and len(stable) >= 0.9 * len(complete)
                    else "B_sign_moves"),
        "verdict_note": (
            "the pre-registered threshold is applied to every orbit in the "
            "sample, including the ones the resolution rule itself declares "
            "unresolved and excludes from the reported counts. It is kept "
            "unchanged and reported as found, but it is not the question a "
            "reader is asking: see conditional_on_resolution"),
        "conditional_on_resolution": _conditional(out),
        "rows": out, "provenance": J.provenance(),
    }
    J.atomic_json(SCORE, payload)
    cond = payload["conditional_on_resolution"]
    log(f"score: {payload['verdict']} on the pre-registered threshold "
        f"({len(stable)}/{len(complete)} orbits keep rho on one side of 1 "
        f"across all four levels); conditional on resolution, "
        f"{cond['of_those_sign_unchanged_in_every_pair_that_resolves']}/"
        f"{cond['admitted_by_manuscript_pair']} admitted comparisons keep "
        "their sign; by stratum "
        + ", ".join(f"{k} {v['sign_stable_across_all_levels']}/{v['orbits']}"
                    for k, v in cond["by_margin_stratum"].items()))
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
