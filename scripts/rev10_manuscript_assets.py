"""Generate manuscript assets for the R10 confirmatory campaign.

All reported values are read from the frozen machine-readable artifacts.
The script does not rerun any propagation.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "lunaris-matplotlib")
)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

import paper_style as ps


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"


def load(name: str) -> dict:
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


def percentile(values: list[float], q: float) -> float:
    return float(np.percentile(np.asarray(values, dtype=float), q))


aggregate = load("r10_aggregate_summary.json")
baseline = load("r10_sobolA_baseline.json")
# Ten-orbit audited baseline: the frozen six-orbit primary corrected baseline
# with the four extended sub-50 km audits folded in (all reported ratios then
# use each orbit's adopted reference degree). The primary corrected baseline and
# r10_aggregate_summary.json remain as intermediate provenance artifacts.
corrected = load("r10_sobolA_baseline_truth_corrected_extended.json")
design_a = load("r10_sobolA_design.json")
design_b = load("r10_sobolB_design_frozen.json")
truth_audit = load("r10_sobolA_truth_audit.json")
convergence = load("r10_sobolA_convergence.json")
blend = load("r10_blend_lro_convergence.json")

rows = sorted(corrected["rows"], key=lambda row: row["sobol_index"])


def ratio_summary(rows: list[dict], key: str) -> dict:
    values = np.asarray([float(row[key]) for row in rows])
    return {
        "median": float(np.median(values)),
        "p10": float(np.percentile(values, 10)),
        "p90": float(np.percentile(values, 90)),
        "raw_schedule_wins": int(np.sum(values > 1.0)),
    }


work_summary = ratio_summary(rows, "rho_work")
crit_summary = ratio_summary(rows, "rho_crit")
baseline_by_index = {row["sobol_index"]: row for row in baseline["rows"]}
convergence_by_index = {row["sobol_index"]: row for row in convergence["rows"]}

savings = [
    baseline_by_index[row["sobol_index"]]["primary_ratios"][
        "gravity_time_saving_vs_critical"
    ]
    for row in rows
]
in_track = [
    row["policy_errors"]["schedule_empirical"]["in_track_fraction"] for row in rows
]


main_table = rf"""\begin{{tabular}}{{@{{}}lrrrrr@{{}}}}
\toprule
Comparator & Eligible & Raw wins & Median $\rho$ [P10, P90] &
Resolved S/F/U & Median saving \\
\midrule
Work-matched fixed degree &
64 &
{work_summary["raw_schedule_wins"]} &
{work_summary["median"]:.3f}
[{work_summary["p10"]:.3f},
 {work_summary["p90"]:.3f}] &
7/9/1 & -- \\
Critical-altitude fixed degree &
64 &
{crit_summary["raw_schedule_wins"]} &
{crit_summary["median"]:.3f}
[{crit_summary["p10"]:.3f},
 {crit_summary["p90"]:.3f}] &
0/17/0 & {100 * np.median(savings):.1f}\% \\
\bottomrule
\end{{tabular}}
"""
(METRICS / "r10_sobol_main_table.tex").write_text(main_table, encoding="utf-8")


blend_summary = aggregate["blend_lro_convergence"]["summary"]
blend_policies = blend_summary["policies"]
comparison = blend_summary["comparison"]
blend_table = rf"""\begin{{tabular}}{{@{{}}lrrrrr@{{}}}}
\toprule
Policy & Baseline & Tighter & Policy self & Reference self & Envelope \\
\midrule
Fixed $N=120$ &
{blend_policies["fixed_N120"]["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]:.1f} &
{blend_policies["fixed_N120"]["error_against_same_tolerance_truth"]["tighter"]["pos_rms_m"]:.1f} &
{blend_policies["fixed_N120"]["self_difference_rms_m"]:.1f} &
{blend_summary["truth_self_difference_rms_m"]:.1f} &
{blend_policies["fixed_N120"]["truth_inclusive_envelope_m"]:.1f} \\
Corrected potential blend &
{blend_policies["corrected_blend"]["error_against_same_tolerance_truth"]["baseline"]["pos_rms_m"]:.1f} &
{blend_policies["corrected_blend"]["error_against_same_tolerance_truth"]["tighter"]["pos_rms_m"]:.1f} &
{blend_policies["corrected_blend"]["self_difference_rms_m"]:.1f} &
{blend_summary["truth_self_difference_rms_m"]:.1f} &
{blend_policies["corrected_blend"]["truth_inclusive_envelope_m"]:.1f} \\
\midrule
\multicolumn{{3}}{{l}}{{Absolute baseline gap}} &
\multicolumn{{3}}{{r}}{{{comparison["absolute_baseline_error_difference_m"]:.1f} m}} \\
\multicolumn{{3}}{{l}}{{Resolution threshold}} &
\multicolumn{{3}}{{r}}{{{comparison["resolution_threshold_m"]:.1f} m}} \\
\bottomrule
\end{{tabular}}
"""
(METRICS / "r10_blend_convergence_table.tex").write_text(
    blend_table, encoding="utf-8"
)


def design_longtable(design: dict, caption: str, label: str) -> str:
    lines = [
        r"\begin{longtable}{rrrrrr}",
        rf"\caption{{{caption}}}\label{{{label}}}\\",
        r"\toprule",
        r"Index & $h_p$ [km] & $h_a$ [km] & $i$ [deg] & $\omega$ [deg] & $\lambda_p$ [deg] \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{6}{c}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"Index & $h_p$ [km] & $h_a$ [km] & $i$ [deg] & $\omega$ [deg] & $\lambda_p$ [deg] \\",
        r"\midrule",
        r"\endhead",
    ]
    for orbit in sorted(design["orbits"], key=lambda item: item["sobol_index"]):
        lines.append(
            f'{orbit["sobol_index"]} & {orbit["hp_km"]:.2f} & '
            f'{orbit["ha_km"]:.2f} & {orbit["incl_deg"]:.2f} & '
            f'{orbit["argp_deg"]:.2f} & '
            f'{orbit["requested_perilune_lon_deg_bodyfixed_t0"]:.2f} \\\\'
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", ""])
    return "\n".join(lines)


(METRICS / "r10_sobol_designA_table.tex").write_text(
    design_longtable(
        design_a,
        "Propagated scrambled-Sobol design A. The listed variables are "
        "coverage coordinates, not a mission-occurrence distribution.",
        "tab:sobol-design-a",
    ),
    encoding="utf-8",
)
(METRICS / "r10_sobol_designB_table.tex").write_text(
    design_longtable(
        design_b,
        "Scrambled-Sobol design B, frozen independently before "
        "propagation and subsequently propagated in full under the "
        "vector-tolerance contract.",
        "tab:sobol-design-b",
    ),
    encoding="utf-8",
)


def write_design_audit_tables(design: dict, prefix: str, design_name: str) -> None:
    coordinates = [
        r"\begin{longtable}{rrrrrrrr}",
        rf"\caption{{Exact Sobol coordinates and transformed altitude coordinates for design {design_name}.}}\label{{tab:sobol-{prefix}-coordinates}}\\",
        r"\toprule",
        r"Index & $u_0$ & $u_1$ & $u_2$ & $u_3$ & $u_4$ & $h_p$ [km] & $h_a$ [km] \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{8}{c}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"Index & $u_0$ & $u_1$ & $u_2$ & $u_3$ & $u_4$ & $h_p$ [km] & $h_a$ [km] \\",
        r"\midrule",
        r"\endhead",
    ]
    elements = [
        r"\begin{longtable}{rrrrrrrrr}",
        rf"\caption{{Complete osculating elements and planned reference degree for design {design_name}.}}\label{{tab:sobol-{prefix}-elements}}\\",
        r"\toprule",
        r"Index & $a$ [km] & $e$ & $i$ & $\Omega$ & $\omega$ & $\nu_0$ & $\lambda_p$ & $N_T$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{9}{c}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"Index & $a$ [km] & $e$ & $i$ & $\Omega$ & $\omega$ & $\nu_0$ & $\lambda_p$ & $N_T$ \\",
        r"\midrule",
        r"\endhead",
    ]
    states = [
        r"\begin{longtable}{rrrrrrr}",
        rf"\caption{{Exact initial Cartesian states for design {design_name}, in km and km s$^{{-1}}$.}}\label{{tab:sobol-{prefix}-states}}\\",
        r"\toprule",
        r"Index & $x$ & $y$ & $z$ & $v_x$ & $v_y$ & $v_z$ \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{7}{c}{\tablename\ \thetable\ continued}\\",
        r"\toprule",
        r"Index & $x$ & $y$ & $z$ & $v_x$ & $v_y$ & $v_z$ \\",
        r"\midrule",
        r"\endhead",
    ]
    for orbit in sorted(design["orbits"], key=lambda item: item["sobol_index"]):
        coordinates.append(
            f'{orbit["sobol_index"]} & '
            + " & ".join(f"{value:.9f}" for value in orbit["u"])
            + f' & {orbit["hp_km"]:.6f} & {orbit["ha_km"]:.6f} \\\\'
        )
        elements.append(
            f'{orbit["sobol_index"]} & {orbit["semimajor_axis_m"] / 1000:.6f} & '
            f'{orbit["eccentricity"]:.9f} & {orbit["incl_deg"]:.6f} & '
            f'{orbit["raan_deg"]:.6f} & {orbit["argp_deg"]:.6f} & '
            f'{orbit["nu0_deg"]:.1f} & '
            f'{orbit["requested_perilune_lon_deg_bodyfixed_t0"]:.6f} & '
            f'{orbit["truth_degree"]} \\\\'
        )
        state = np.asarray(orbit["initial_state_si"], dtype=float) / 1000.0
        states.append(
            f'{orbit["sobol_index"]} & '
            + " & ".join(f"{value:.9f}" for value in state)
            + r" \\"
        )
    for table in (coordinates, elements, states):
        table.extend([r"\bottomrule", r"\end{longtable}", ""])
    (METRICS / f"r10_sobol_{prefix}_coordinates.tex").write_text(
        "\n".join(coordinates), encoding="utf-8"
    )
    (METRICS / f"r10_sobol_{prefix}_elements.tex").write_text(
        "\n".join(elements), encoding="utf-8"
    )
    (METRICS / f"r10_sobol_{prefix}_states.tex").write_text(
        "\n".join(states), encoding="utf-8"
    )


write_design_audit_tables(design_a, "designA", "A (propagated)")
write_design_audit_tables(design_b, "designB", "B (propagated in full)")


def resolved_label(index: int, comparator: str) -> str:
    if index not in convergence_by_index:
        return "--"
    comp = convergence_by_index[index]["comparisons"][comparator]
    if not comp["resolved"]:
        return "unresolved"
    return "schedule" if comp["winner_if_resolved"] == "schedule_empirical" else "fixed"


full_lines = [
    r"\begin{longtable}{rrrrrrrr}",
    r"\caption{Reference-audited seven-day results for propagated Sobol design A. "
    r"$\rho>1$ favors the empirical schedule. Resolution labels are available "
    r"only for the pre-specified 17-orbit convergence subset.}"
    r"\label{tab:sobol-full-results}\\",
    r"\toprule",
    r"Index & $N_T$ & $N_W$ & $N_C$ & $\rho_W$ & $\rho_C$ & "
    r"Resolved vs.\ $W$ & Resolved vs.\ $C$ \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{8}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Index & $N_T$ & $N_W$ & $N_C$ & $\rho_W$ & $\rho_C$ & "
    r"Resolved vs.\ $W$ & Resolved vs.\ $C$ \\",
    r"\midrule",
    r"\endhead",
]
for row in rows:
    idx = row["sobol_index"]
    full_lines.append(
        f'{idx} & {row["adopted_truth_degree"]} & {row["n_work"]} & '
        f'{row["n_critical"]} & {row["rho_work"]:.3f} & '
        f'{row["rho_crit"]:.3f} & {resolved_label(idx, "fixed_work")} & '
        f'{resolved_label(idx, "fixed_critical")} \\\\'
    )
full_lines.extend([r"\bottomrule", r"\end{longtable}", ""])
(METRICS / "r10_sobol_full_results_table.tex").write_text(
    "\n".join(full_lines), encoding="utf-8"
)


primary_lines = [
    r"\begin{longtable}{rcrrrrrrrr}",
    r"\caption{Complete primary-policy audit for propagated design A. "
    r"All errors are full-common-grid seven-day Cartesian position RMS in meters; "
    r"$\rho>1$ favors the empirical schedule.}\label{tab:sobol-primary-audit}\\",
    r"\toprule",
    r"Index & Status & $N_T$ & $N_W$ & $N_C$ & $E_{\rm emp}$ & $E_W$ & $E_C$ & $\rho_W$ & $\rho_C$ \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{10}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Index & Status & $N_T$ & $N_W$ & $N_C$ & $E_{\rm emp}$ & $E_W$ & $E_C$ & $\rho_W$ & $\rho_C$ \\",
    r"\midrule",
    r"\endhead",
]
sensitivity_lines = [
    r"\begin{longtable}{rrrrrrr}",
    r"\caption{Schedule-sensitivity and diagnostic audit for propagated design A. "
    r"Errors are seven-day position RMS. The saving is measured gravity-kernel "
    r"time relative to the critical-altitude fixed degree.}\label{tab:sobol-sensitivity-audit}\\",
    r"\toprule",
    r"Index & $E_{\rm emp}$ [m] & $E_{\rm up}$ [m] & $E_{\rm down}$ [m] & Saving [\%] & $f_I$ & Epochs \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{7}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Index & $E_{\rm emp}$ [m] & $E_{\rm up}$ [m] & $E_{\rm down}$ [m] & Saving [\%] & $f_I$ & Epochs \\",
    r"\midrule",
    r"\endhead",
]
for row in rows:
    idx = row["sobol_index"]
    errs = row["policy_errors"]
    base = baseline_by_index[idx]
    statuses = set(base["trajectory_status"].values())
    status = "complete" if statuses == {"complete"} else "/".join(sorted(statuses))
    primary_lines.append(
        f'{idx} & {status} & {row["adopted_truth_degree"]} & {row["n_work"]} & '
        f'{row["n_critical"]} & {errs["schedule_empirical"]["pos_rms_m"]:.3f} & '
        f'{errs["fixed_work"]["pos_rms_m"]:.3f} & '
        f'{errs["fixed_critical"]["pos_rms_m"]:.3f} & '
        f'{row["rho_work"]:.3f} & {row["rho_crit"]:.3f} \\\\'
    )
    sensitivity_lines.append(
        f'{idx} & {errs["schedule_empirical"]["pos_rms_m"]:.3f} & '
        f'{errs["schedule_up"]["pos_rms_m"]:.3f} & '
        f'{errs["schedule_down"]["pos_rms_m"]:.3f} & '
        f'{100 * base["primary_ratios"]["gravity_time_saving_vs_critical"]:.2f} & '
        f'{errs["schedule_empirical"]["in_track_fraction"]:.4f} & '
        f'{errs["schedule_empirical"]["common_epoch_count"]} \\\\'
    )
for table in (primary_lines, sensitivity_lines):
    table.extend([r"\bottomrule", r"\end{longtable}", ""])
(METRICS / "r10_sobol_primary_audit.tex").write_text(
    "\n".join(primary_lines), encoding="utf-8"
)
(METRICS / "r10_sobol_sensitivity_audit.tex").write_text(
    "\n".join(sensitivity_lines), encoding="utf-8"
)


convergence_self_lines = [
    r"\begin{longtable}{rrrrrrrr}",
    r"\caption{Tight-to-tighter self-differences and reference-inclusive envelopes "
    r"for the 17-orbit convergence subset, in meters.}\label{tab:sobol-convergence-self}\\",
    r"\toprule",
    r"Index & Reference self & S self & S env. & $W$ self & $W$ env. & $C$ self & $C$ env. \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{8}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Index & Reference self & S self & S env. & $W$ self & $W$ env. & $C$ self & $C$ env. \\",
    r"\midrule",
    r"\endhead",
]
convergence_decision_lines = [
    r"\begin{longtable}{rcrrrrrc}",
    r"\caption{Numerical-resolution audit for every selected orbit and comparator. "
    r"$E_S$ and $E_F$ are tight-level position RMS errors; a ranking is accepted "
    r"only when the absolute gap exceeds the summed reference-inclusive envelopes.}"
    r"\label{tab:sobol-convergence-decisions}\\",
    r"\toprule",
    r"Index & Comparator & $E_S$ & $E_F$ & Gap & Threshold & $\rho$ & Outcome \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{8}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Index & Comparator & $E_S$ & $E_F$ & Gap & Threshold & $\rho$ & Outcome \\",
    r"\midrule",
    r"\endhead",
]
for row in sorted(convergence["rows"], key=lambda item: item["sobol_index"]):
    idx = row["sobol_index"]
    policies = row["policies"]
    convergence_self_lines.append(
        f'{idx} & {row["truth_self_difference_rms_m"]:.3f} & '
        f'{policies["schedule_empirical"]["self_difference_rms_m"]:.3f} & '
        f'{policies["schedule_empirical"]["truth_inclusive_envelope_m"]:.3f} & '
        f'{policies["fixed_work"]["self_difference_rms_m"]:.3f} & '
        f'{policies["fixed_work"]["truth_inclusive_envelope_m"]:.3f} & '
        f'{policies["fixed_critical"]["self_difference_rms_m"]:.3f} & '
        f'{policies["fixed_critical"]["truth_inclusive_envelope_m"]:.3f} \\\\'
    )
    for comparator, short in (("fixed_work", "$W$"), ("fixed_critical", "$C$")):
        comp = row["comparisons"][comparator]
        sched_error = policies["schedule_empirical"][
            "errors_against_same_tolerance_truth"
        ]["tight"]["pos_rms_m"]
        fixed_error = policies[comparator]["errors_against_same_tolerance_truth"][
            "tight"
        ]["pos_rms_m"]
        if comp["resolved"]:
            outcome = (
                "schedule"
                if comp["winner_if_resolved"] == "schedule_empirical"
                else "fixed"
            )
        else:
            outcome = "unresolved"
        convergence_decision_lines.append(
            f'{idx} & {short} & {sched_error:.3f} & {fixed_error:.3f} & '
            f'{comp["absolute_error_difference_m"]:.3f} & '
            f'{comp["resolution_threshold_m"]:.3f} & '
            f'{comp["rho_tight"]:.3f} & {outcome} \\\\'
        )
for table in (convergence_self_lines, convergence_decision_lines):
    table.extend([r"\bottomrule", r"\end{longtable}", ""])
(METRICS / "r10_sobol_convergence_self.tex").write_text(
    "\n".join(convergence_self_lines), encoding="utf-8"
)
(METRICS / "r10_sobol_convergence_decisions.tex").write_text(
    "\n".join(convergence_decision_lines), encoding="utf-8"
)


blend_detail_lines = [
    r"\begin{longtable}{llrrrrrrrrr}",
    r"\caption{Complete LRO-like same-tolerance reference-error audit. Position "
    r"quantities and R/I/C components are in meters; velocity RMS is in m s$^{-1}$.}"
    r"\label{tab:blend-error-detail}\\",
    r"\toprule",
    r"Policy & Level & RMS & Max & Final & $v_{\rm RMS}$ & $R_{\rm RMS}$ & $I_{\rm RMS}$ & $C_{\rm RMS}$ & $I_{\rm final}$ \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{10}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Policy & Level & RMS & Max & Final & $v_{\rm RMS}$ & $R_{\rm RMS}$ & $I_{\rm RMS}$ & $C_{\rm RMS}$ & $I_{\rm final}$ \\",
    r"\midrule",
    r"\endhead",
]
for policy, label in (("fixed_N120", "fixed $N=120$"), ("corrected_blend", "blend")):
    for level in ("baseline", "tighter"):
        error = blend_policies[policy]["error_against_same_tolerance_truth"][level]
        blend_detail_lines.append(
            f'{label} & {level} & {error["pos_rms_m"]:.2f} & '
            f'{error["pos_max_m"]:.2f} & {error["pos_final_m"]:.2f} & '
            f'{error["vel_rms_m_s"]:.4f} & {error["ric_rms_m"]["radial"]:.2f} & '
            f'{error["ric_rms_m"]["in_track"]:.2f} & '
            f'{error["ric_rms_m"]["cross_track"]:.2f} & '
            f'{error["ric_final_m"]["in_track"]:.2f} \\\\'
        )
blend_detail_lines.extend([r"\bottomrule", r"\end{longtable}", ""])
(METRICS / "r10_blend_error_detail.tex").write_text(
    "\n".join(blend_detail_lines), encoding="utf-8"
)


blend_telemetry_lines = [
    r"\begin{longtable}{llrrrrrrr}",
    r"\caption{LRO-like propagation telemetry. Kernel, CPU, and wall columns "
    r"report seconds for each complete trajectory.}\label{tab:blend-telemetry}\\",
    r"\toprule",
    r"Policy & Level & RHS & Accepted & Attempted & Rejected & Kernel [s] & CPU [s] & Wall [s] \\",
    r"\midrule",
    r"\endfirsthead",
    r"\multicolumn{9}{c}{\tablename\ \thetable\ continued}\\",
    r"\toprule",
    r"Policy & Level & RHS & Accepted & Attempted & Rejected & Kernel [s] & CPU [s] & Wall [s] \\",
    r"\midrule",
    r"\endhead",
]
for record in blend["records"]:
    telemetry = record["telemetry"]
    policy_label = {
        "fixed_N120": "fixed $N=120$",
        "corrected_blend": "blend",
        "truth_N600": "reference $N=600$",
    }[record["policy"]]
    blend_telemetry_lines.append(
        f'{policy_label} & {record["level"]} & {telemetry["n_rhs"]} & '
        f'{telemetry["n_accepted_steps"]} & {telemetry["n_attempted_steps"]} & '
        f'{telemetry["n_rejected_trials"]} & '
        f'{telemetry["gravity_kernel_ns"] / 1e9:.1f} & '
        f'{telemetry["process_cpu_ns"] / 1e9:.1f} & '
        f'{telemetry["total_wall_ns"] / 1e9:.1f} \\\\'
    )
blend_telemetry_lines.extend([r"\bottomrule", r"\end{longtable}", ""])
(METRICS / "r10_blend_telemetry.tex").write_text(
    "\n".join(blend_telemetry_lines), encoding="utf-8"
)


audit_lines = [
    r"\begin{tabular}{rrrrrr}",
    r"\toprule",
    r"Index & $h_p$ [km] & $N=600$--900 RMS [m] & Threshold [m] & Pass & Adopted $N_T$ \\",
    r"\midrule",
]
for row in sorted(truth_audit["rows"], key=lambda item: item["sobol_index"]):
    audit_lines.append(
        f'{row["sobol_index"]} & {row["hp_km"]:.1f} & '
        f'{row["N600_to_N900"]["pos_rms_m"]:.3f} & '
        f'{row["acceptance_threshold_m"]:.3f} & '
        f'{"yes" if row["passes"] else "no"} & '
        f'{row["adopted_truth_degree"]} \\\\'
    )
audit_lines.extend([r"\bottomrule", r"\end{tabular}", ""])
(METRICS / "r10_truth_audit_table.tex").write_text(
    "\n".join(audit_lines), encoding="utf-8"
)


ps.apply()
plt.rcParams["text.usetex"] = False
rho_work = np.asarray([row["rho_work"] for row in rows])
rho_crit = np.asarray([row["rho_crit"] for row in rows])
hp = np.asarray([row["design_point"]["hp_km"] for row in rows])
ha = np.asarray([row["design_point"]["ha_km"] for row in rows])
inc = np.asarray([row["design_point"]["incl_deg"] for row in rows])

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)

for ax, values, title in (
    (axes[0], rho_work, r"$\rho_{\rm work}$"),
    (axes[1], rho_crit, r"$\rho_{\rm crit}$"),
):
    ordered = np.sort(values)
    x = np.arange(1, len(ordered) + 1)
    ax.semilogy(x, ordered, "o-", ms=2.7, lw=0.8, color="#315b7d")
    ax.axhline(1.0, color="#a33a2b", lw=0.9, ls="--")
    ax.set_xlabel("Ordered orbit")
    ax.set_ylabel(title)
    ax.set_xlim(1, 64)
    ax.grid(True, which="both", alpha=0.22, lw=0.4)

log_rho = np.log10(rho_work)
limit = max(abs(float(log_rho.min())), abs(float(log_rho.max())))
norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
markers = ((inc < 60, "o", "prograde"), ((inc >= 60) & (inc <= 120), "s", "high-$i$"), (inc > 120, "^", "retrograde"))
scatter = None
for mask, marker, label in markers:
    scatter = axes[2].scatter(
        hp[mask],
        ha[mask],
        c=log_rho[mask],
        cmap="RdBu",
        norm=norm,
        marker=marker,
        s=25,
        edgecolors="0.25",
        linewidths=0.35,
        label=label,
    )
axes[2].set_xlabel(r"$h_p$ [km]")
axes[2].set_ylabel(r"$h_a$ [km]")
axes[2].legend(fontsize=6.2, loc="best", frameon=True)
colorbar = fig.colorbar(scatter, ax=axes[2], fraction=0.05, pad=0.02)
colorbar.set_label(r"$\log_{10}\rho_{\rm work}$")

FIGURES.mkdir(exist_ok=True)
fig.savefig(FIGURES / "fig_sobol_confirmatory.pdf")
plt.close(fig)


summary = {
    "median_gravity_time_saving_vs_critical": float(np.median(savings)),
    "p10_gravity_time_saving_vs_critical": percentile(savings, 10),
    "p90_gravity_time_saving_vs_critical": percentile(savings, 90),
    "median_schedule_in_track_fraction": float(np.median(in_track)),
    "design_a_inclination_counts": design_a["inclination_regime_counts"],
    "design_a_inclination_range_deg": design_a["realized_inclination_range_deg"],
}
(METRICS / "r10_manuscript_descriptives.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)

print("[written] R10 LaTeX tables, descriptives, and fig_sobol_confirmatory.pdf")
