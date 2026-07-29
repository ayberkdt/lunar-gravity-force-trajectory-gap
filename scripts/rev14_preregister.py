"""Pre-registration of the O25 fixed-budget Pareto campaign (R14).

Written and frozen BEFORE any aggregate O25 result is inspected. It fixes the
hypotheses, the budget grid, the calibration and comparator rules, the censoring
convention, the rule that governs any later extension of the budget grid, and
the decision logic that maps outcomes onto manuscript claims.

The point of freezing the extension rule is that the trajectory campaign is
staged: only beta = 1 is propagated on both populations up front, and further
budgets are added afterwards. Without a written rule, choosing those budgets
after seeing the beta = 1 aggregate would be result shopping.

Usage:
    python rev14_preregister.py
"""

from __future__ import annotations

import json
from pathlib import Path

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metrics" / "r14_preregistration.json"

BUDGET_GRID = [0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 3.00]

PAYLOAD = {
    "schema": "r14_preregistration_v1",
    "campaign": "R14 -- O25 fixed-budget radial-allocation Pareto study",
    "question": (
        "Given a fixed gravity-computation budget C, does a budget-calibrated "
        "Atallah radial degree allocation produce better accuracy than spending "
        "the same budget on one constant harmonic degree?"),
    "why_the_existing_benchmark_does_not_answer_it": (
        "The archived R12/R13 benchmark evaluates the published rule at the "
        "operating point where its acceleration tolerance is set from the "
        "perilune fidelity of the critical-altitude fixed degree. That point "
        "costs roughly 2.7-2.8 times the critical fixed quadratic work, so the "
        "archived comparisons are A(2.8x) vs F(x) and A(2.8x) vs F(2.8x). The "
        "missing comparison is A(x) vs F(x)."),

    "definitions": {
        "W_crit": "N_crit^2, the per-call quadratic work of the orbit's critical-altitude fixed degree",
        "beta": "W / W_crit, the normalized per-call gravity-work budget",
        "budget_calibrated_atallah": (
            "the published Atallah (2022) radial mapping N_req(r; eps_A), with its "
            "user accuracy parameter eps_A chosen so the resulting degree history "
            "consumes the requested budget. This is not a new algorithm: the online "
            "degree decision remains radial and local, while its tolerance is "
            "calibrated offline to a declared orbit-level compute budget."),
        "W_A_sampled": "<N_A(t)^2> over the archived truth output epochs (primary diagnostic budget)",
        "W_total": "sum over accepted RHS calls of N_j^2 (includes the scheduler's effect on call count)",
        "measured_time": "serial single-core gravity-kernel wall time (architecture dependent)",
    },

    "budget_grid": BUDGET_GRID,
    "additional_non_grid_point": (
        "each orbit's archived Atallah accuracy-target operating point, at its own "
        "measured beta_At_orig = W_At_orig / W_crit; never approximated as 2.8"),

    "calibration_rule": {
        "method": "log-space bisection on eps_A",
        "target": "|W_A_sampled / (beta * N_crit^2) - 1| < 0.01",
        "degree_history": (
            "10-km binned altitude table, identical to the binning convention used "
            "by the archived Atallah campaign and by the paper's other altitude "
            "schedules, so that the calibrated policy is the one actually propagated"),
        "floor": 2,
        "cap": "the orbit's adopted truth degree",
        "on_failure": (
            "if integer-degree discretization prevents the 1% target, take the "
            "nearest attainable point and record the achieved mismatch; never "
            "interpolate degree"),
    },
    "fixed_comparator_rule": {
        "definition": "N_F(beta) = argmin_N |N^2 - beta * N_crit^2| over integers",
        "note": "equivalently round(N_crit*sqrt(beta)) up to the integer tie rule; never interpolated",
        "cap": "the orbit's adopted truth degree",
    },
    "censoring_rule": {
        "statement": (
            "a budget point is censored for an orbit when the comparator or the "
            "calibrated Atallah history would need degrees above the adopted truth "
            "degree, because the defect against that truth is then not a meaningful "
            "measure of truncation error"),
        "handling": (
            "censored points are flagged and counted in every table (cap-bound "
            "orbit count) and are excluded from the population statistics of that "
            "budget, never silently dropped and never joined across in a figure"),
    },

    "hypotheses": {
        "H1": ("At equal critical-fixed quadratic budget (beta = 1), budget-calibrated "
               "Atallah yields a smaller truncation force defect than the fixed degree."),
        "H2": ("At beta = 1, budget-calibrated Atallah yields a smaller seven-day "
               "Cartesian RMS trajectory error on a majority of both independent "
               "populations."),
        "H3": ("There exists a budget range in which radial adaptation lies on a "
               "better cost-accuracy frontier than fixed degree."),
        "H4": ("Force-level dominance does not necessarily imply orbit-by-orbit "
               "trajectory dominance."),
    },
    "hypotheses_frozen": True,

    "staging": {
        "priority_1": "all-128 force-level budget Pareto sweep (archived truth reused)",
        "priority_2": "full design A + design B trajectory test at beta = 1",
        "priority_3": "serial measured-time validation of beta = 1 on the 14-orbit panel",
        "priority_4": "design-A trajectory sweep at beta = 0.75, 1.50, 3.00",
        "priority_5": "design-B replication at beta = 1 always, plus one budget nearest any crossover",
        "priority_6": "forced-variational calibration at beta = 1 on the eight-orbit panel",
        "priority_7": "optional O26 force-allocation oracle",
    },

    "adaptive_extension_rule": {
        "written_before_results": True,
        "rule": (
            "Additional trajectory budgets beyond the pre-declared {0.75, 1.50, 3.00} "
            "may be propagated only when the Phase-A force-level curves or the "
            "already-propagated trajectory budgets place a crossing of the fixed and "
            "Atallah population medians strictly between two adjacent grid budgets "
            "beta_a < beta_b. In that case exactly the two bracketing grid values are "
            "propagated, on design A first and on design B only for the single grid "
            "value nearest the crossing. No budget may be added for any other reason, "
            "and no budget already run may be removed from the reported tables."),
        "crossover_reporting": (
            "a crossover is reported only as an interval between sampled budgets "
            "('the tested Pareto curves cross between beta_a and beta_b'); no "
            "interpolated crossover value is quoted"),
    },

    "resolution_rule": {
        "inherited_from": "R11/R12 truth-inclusive pairwise resolution",
        "criterion": "|E_A - E_B| > (E_self,A + E_self,T) + (E_self,B + E_self,T)",
        "unresolved_handling": "reported as unresolved; never counted as a tie or as a win",
    },
    "numerical_contract": {
        "tight": {"rtol": 1.0e-12, "atol_position_m": 1.0e-5, "atol_velocity_m_s": 1.0e-8},
        "tighter": {"rtol": 1.0e-13, "atol_position_m": 1.0e-6, "atol_velocity_m_s": 1.0e-9},
        "atol_kind": "vector", "max_step_s": 60.0,
        "duration_s": 604800.0, "output_step_s": 120.0,
        "truth_reuse": ("archived R11 truth trajectories are reused only where the "
                        "tolerance, force model, frame, output grid, truth degree, "
                        "initial state, and code identity contract match exactly"),
    },

    "decision_logic": {
        "A_atallah_dominates_at_beta_1": (
            "if force and trajectory population results both favor Atallah at equal "
            "budget, the manuscript may state that radial allocation provides a more "
            "efficient cost-accuracy trade over the tested populations. 'More "
            "efficient' requires actual cost matching. 'Optimal' is not permitted "
            "without the O26 oracle, and not even then at trajectory level."),
        "B_force_favors_atallah_trajectory_mixed": (
            "state that budget-constrained radial allocation suppresses the local "
            "truncation defect more efficiently than constant degree, but that the "
            "force-level advantage does not map monotonically into long-arc state error"),
        "C_curves_cross": (
            "identify the crossover only as a bracketing interval between sampled "
            "budgets; conclude that radial adaptation becomes beneficial once enough "
            "compute is available"),
        "D_fixed_dominates_at_equal_budget": (
            "report it plainly: the accuracy of the published rule at its "
            "accuracy-target operating point is purchased with additional gravity "
            "work, and at held-fixed budget the constant degree is more "
            "trajectory-efficient over the tested population"),
    },

    "prohibited": [
        "adding budget points because they improve the story",
        "removing low-budget failures",
        "selecting only Atallah-favorable orbits",
        "reporting only force-level results when trajectory results disagree",
        "calling an unresolved comparison a tie",
        "converting a force-defect win into a claim of trajectory superiority",
        "describing the O26 oracle as operational or flight-realizable",
        "quoting a crossover budget more precisely than the sampled grid supports",
    ],

    "interpretation_rule_non_negotiable": (
        "Even if the force defect of the budget-calibrated rule is smaller at every "
        "budget and on every orbit, no trajectory-optimality claim follows. Long-arc "
        "displacement is the state-transition-weighted integral of the defect, so "
        "sign, direction, phase and coherence matter: a uniformly smaller defect norm "
        "can still produce a larger state displacement."),

    "stop_condition": (
        "the campaign stops after the full 128-orbit force-level Pareto sweep, the "
        "full 128-orbit beta = 1 trajectory comparison, the numerical-resolution "
        "audit, the serial measured-time control, and at least one additional "
        "trajectory budget on each side of any detected crossover. No third Sobol "
        "population, no other central body, no further integrator or symplectic "
        "campaign, and no further accuracy-target Atallah runs are added unless a "
        "specific result fails an existing qualification criterion."),
}


def main() -> int:
    payload = dict(PAYLOAD)
    payload["created_utc"] = base.utc_now()
    payload["source"] = base.provenance()
    payload["protocol_sha256"] = base.object_hash(
        {k: v for k, v in payload.items() if k not in ("created_utc", "source")})
    base.atomic_json(OUTPUT, payload)
    print(f"[preregistration] {OUTPUT.name}")
    print(f"  protocol_sha256 = {payload['protocol_sha256']}")
    print(f"  budget grid      = {BUDGET_GRID}")
    print(json.dumps(payload["hypotheses"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
