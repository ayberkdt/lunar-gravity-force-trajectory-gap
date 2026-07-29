# Claim-to-artifact map

Which script produced which manuscript item, and which campaign manifest covers it. Generated from the LaTeX sources, the compiled label numbering and the campaign manifests; the machine-readable version is [`REPRODUCIBILITY_INDEX.csv`](../REPRODUCIBILITY_INDEX.csv).

Analysis and table passes read the archived records in `metrics/` in place, so most rows can be re-run without propagating anything. Rows marked `(inline)` are tables typed directly in the manuscript from numbers reported in the text; they have no generated file.

## Main text

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **Figure 1** | Per-coefficient RMS spectrum of JGGRX\_1800F with power-law fits over the observationally dominated band (n ) … | `fig_spectrum.pdf` | `make_figures_r1.py` | R1 |
| **Figure 2** | Truncation criteria against the model-relative empirical reference (log--log) | `fig_truncation.pdf` | `make_figures_r1.py` | R1 |
| **Figure 3** | 24-hour RMS position error versus truncation degree for the 100~km circular polar orbit at two tolerance setti… | `fig_orbit_error.pdf` | `make_figures_r1.py` | R1 |
| **Figure 4** | Signed (upper) and absolute-log (lower) relative specific-energy error over the band-crossing orbit | `fig_blend_energy.pdf` | `make_figures_r2.py` | R2 |
| **Figure 5** | Confirmatory scrambled-Sobol design | `fig_sobol_confirmatory.pdf` | `rev10_finalize_manifest.py`<br>`rev10_manuscript_assets.py` | R10 |
| **Figure 6** | Force-defect curves and propagated trajectory operating points for radial and constant-degree truncation, one … | `budget_pareto.pdf` | `make_figures_r14.py`<br>`rev14_deliverables.py`<br>`rev14_finalize_manifest.py` | R14 |
| **Figure 7** | Equal-budget error ratio =E_ fixed | `fig_variational_parity.pdf` | `make_figures_r19.py` | R19 |
| **Figure 8** | The interpolation path at the critical degree's own nominal per-call budget ( =1), plotted against the degree … | `fig_span_curve.pdf` | `make_figures_r20.py` | R20 |
| **Table 1** | Adjacent degree-selection studies and the aspects measured in each | *(inline)* | — | — |
| **Table 2** | What carries what | *(inline)* | — | — |
| **Table 3** | Recommended truncation degree: model-relative empirical reference and closed-form rules | *(inline)* | — | — |
| **Table 4** | Practical static selections | `review_static_selection_table.tex` | `review_postprocess.py` | — |
| **Table 5** | Transfer of the calibration procedure | `r16_transfer_table.tex` | `rev16_finalize_manifest.py`<br>`rev16_tables.py` | R16 |
| **Table 6** | Physics of four degree-transition policies on a common band | *(inline)* | — | — |
| **Table 7** | 7-day eccentric arc (50 300~km polar), tight tolerance: position error against the N=300 adopted reference, qu… | *(inline)* | — | — |
| **Table 8** | Direct Atallah radial-adaptive benchmark on the 64-orbit design-A set, perilune-tolerance matched to the criti… | `r12_atallah_benchmark_table.tex` | `rev12_finalize_manifest.py` | R12 |
| **Table 9** | Fixed-budget force-level comparison; each cell pairs design~A with design~B so the replication reads across | `r14_force_pareto_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **Table 10** | Equal-budget comparison with the budget defined as measured serial gravity-kernel time rather than the quadrat… | `r14_timing_budget_table_compact.tex` | `rev14_timing_budget.py` | R14 |
| **Table 11** | Propagated fixed-budget comparison | `r14_trajectory_pareto_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **Table 12** | The same path in numbers, at =1 | `r18_span_table.tex` | `rev18_finalize_manifest.py`<br>`rev18_tables.py` | R18 |
| **Table 13** | The interior member (k=0.5) against the constant degree under two definitions of an equal budget, at two decla… | `r19_equal_work_table.tex` | `rev19_finalize_manifest.py`<br>`rev19_tables.py` | R19 |
| **Table 14** | The =1 realized-work comparison re-run at a third, ultra-tight tolerance level | `r23_ultra_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **Table 15** | The interpolation path at three declared budgets | `r18_budget_table.tex` | `rev18_budget_table.py` | R18 |
| **Table 16** | Practical recommendations within the scope tested here | *(inline)* | — | — |

## Supplement

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **SFigure S9** | Uniform margin control on the 24-orbit design | `fig_alpha_margin.pdf` | `make_figures_r8.py` | R8 |
| **SFigure S11** | The published radial-adaptive rule against both fixed comparators, design~A (circles) and design~B (triangles) | `fig_atallah_benchmark.pdf` | `make_figures_r13.py`<br>`rev13_finalize_manifest.py` | R13 |
| **SFigure S1** | Dense-altitude calibration objective for the effective tail-budget exponent on JGGRX\_1800F | `fig_pstar_objective.pdf` | `make_figures_supplemental.py` | — |
| **SFigure S2** | Measured single-thread evaluation time of the serial kernel (median and interquartile range over nine timing b… | `fig_timing.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S3** | Two orbital periods of the eccentric arc: altitude, scheduled degree ( =10^ -3 | `fig_schedule_panels.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S4** | Regime map of the stratified 24-orbit design | `fig_doe_regime.pdf` | `make_figures_doe_regime.py` | — |
| **SFigure S5** | Event-aligned accepted-step median and interquartile range (upper panels) and direct rejection probability (lo… | `fig_switch_aggregate.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S6** | Seven-day RMS position error for six degree policies across four orbit/orientation controls | `fig_longarc_matrix.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S7** | Dense fixed-degree sweep at the tight tolerance for the reference phase, apolune-start phase, and 60^ -inclina… | `fig_degree_sweep.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S8** | Position-error envelopes over the tested 7-day arc (tight tolerance) | `fig_longarc_growth.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S10** | One-revolution PEFRL energy control | `fig_symplectic_switch.pdf` | `rev4_robustness_controls.py` | R4 |
| **STable S1** | Field-level acceleration comparison against SHTOOLS, pooled over 30, 100, and 500~km and including near-pole p… | *(inline)* | — | — |
| **STable S2** | Reproducibility contract of the Tudat--Lunaris comparison | *(inline)* | — | — |
| **STable S3** | Three-resolution Tudat RK4 controls (5, 2.5, and 1.25~s) on the common 60-s comparison grid | *(inline)* | — | — |
| **STable S21** | Uniform-margin control, design aggregates over the 24-orbit design | `r8_alpha_margin_table.tex` | `make_figures_r8.py` | R8 |
| **STable S22** | Rank-resolution counts for the aggregate margin analysis | `r8_alpha_margin_resolution_table.tex` | `make_figures_r8.py` | R8 |
| **STable S69** | Twenty-eight-day Cartesian position RMS (m) against the common high-degree reference for the smooth N=30/120 t… | `r9_potential_blend_longarc_table.tex` | `rev9_blend_postprocess.py` | R9 |
| **STable S70** | LRO-like 28-day convergence of fixed N=120 and the corrected potential blend | `r10_blend_convergence_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S63** | Gradient-degree sensitivity of the forced variational solve, over 241 sampled epochs of each panel orbit's arc… | `r21_gradient_sensitivity_table.tex` | `rev21_finalize_manifest.py`<br>`rev21_gradient_sensitivity.py` | R21 |
| **STable S5** | Verification of the per-degree acceleration-RMS identity by exact Gauss--Legendre quadrature (nodes n + | *(inline)* | — | — |
| **STable S6** | The kernel cost curve re-measured under both archived protocols at both archived model degrees, on an idle mac… | `r23_cost_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **STable S7** | Empirical minimum degree N^ emp | `r16_transfer_detail_table.tex` | `rev16_finalize_manifest.py`<br>`rev16_tables.py` | R16 |
| **STable S8** | Normalized band-difference RMS amplitudes (1000 directions per altitude; bootstrap 95\% intervals in brackets) | *(inline)* | — | — |
| **STable S9** | Degree recommendations from the sampled sphere-RMS criterion and directional tail-amplitude percentiles normal… | *(inline)* | — | — |
| **STable S4** | Integration-tolerance regimes and the results quoted from each | *(inline)* | — | — |
| **STable S74** | Sixty-day comparison on the widened geometry set | `r17_longarc60_table.tex` | `rev17_finalize_manifest.py`<br>`rev17_tables.py` | R17 |
| **STable S75** | Sixty-day growth by checkpoint (meters) | `r17_longarc60_growth_table.tex` | `rev17_finalize_manifest.py`<br>`rev17_tables.py` | R17 |
| **STable S76** | Twenty-eight-day Cartesian position RMS E_ RMS | *(inline)* | — | — |
| **STable S10** | Seven-day N=300 minus N=600 reference-degree control (meters), using vector absolute tolerances of 10^ -6 | *(inline)* | — | — |
| **STable S11** | Narrow-gap precision controls | *(inline)* | — | — |
| **STable S12** | Canonical geometries in the 24-orbit design | *(inline)* | — | — |
| **STable S13** | Per-switch acceleration jump: direction-averaged RMS of \\| a | *(inline)* | — | — |
| **STable S14** | Direct DOP853 switching telemetry over 2.2 orbital periods | *(inline)* | — | — |
| **STable S15** | Scalar-tolerance screening convergence; not used for close policy ranking | *(inline)* | — | — |
| **STable S16** | Exploratory scalar-tolerance multi-geometry screening matrix; final narrow-gap rankings use vector tolerances | *(inline)* | — | — |
| **STable S17** | Initial-phase scalar-tolerance screening sweep of the 7-day 50 300~km polar comparison: RMS position error (m)… | *(inline)* | — | — |
| **STable S19** | Weekly checkpoints recomputed from the archived 28-day records for the schedule selected by minimum full-arc p… | `review_stage3_checkpoints.tex` | `review_postprocess.py` | — |
| **STable S20** | Fraction of the seven-day arc spent at each active degree under event-resolved downward scheduling | *(inline)* | — | — |
| **STable S24** | Seven-day robustness controls: Cartesian RMS position error against the same-solver, same-force N=300 referenc… | *(inline)* | — | — |
| **STable S25** | LRO-like seven-day accuracy controls against the same-force N=300 reference | *(inline)* | — | — |
| **STable S26** | LRO-like computational controls | *(inline)* | — | — |
| **STable S27** | Seven-day rotating-frame Jacobi-integral control | *(inline)* | — | — |
| **STable S28** | Independent TudatPy reproduction of dynamic scheduling: seven-day Cartesian RMS position error against each co… | *(inline)* | — | — |
| **STable S29** | LRO-like Tudat step-convergence | *(inline)* | — | — |
| **STable S30** | Confirmatory truth-degree audit, frozen six-orbit primary set | `r10_truth_audit_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S31** | Extended truth-degree audit of the four remaining sub-50~km survivors, using the same acceptance rule | `r10_truth_audit_extended_table.tex` | `rev10_truth_audit_extended_postprocess.py` | R10 |
| **STable S32** | Pre-specified 17-orbit convergence subset and the selection category (or categories) that placed each orbit in… | *(inline)* | — | — |
| **STable S39** | Aggregate confirmatory counts of the initial screening pass | `r10_sobol_main_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S77** | The continuation family member by member, at every declared budget | `r18_by_k_table.tex` | `rev22_supplement_tables.py` | R22 |
| **STable S78** | Design~A: seven-day position RMS (m) at the tighter level for each member of the family, by orbit | `r18_span_detail_table_A.tex` | `rev18_finalize_manifest.py` | R18 |
| **STable S79** | Design~B: seven-day position RMS (m) at the tighter level for each member of the family, by orbit | `r18_span_detail_table_B.tex` | `rev18_finalize_manifest.py` | R18 |
| **STable S80** | The interior member (k=0.5) at =1 against two constant comparators on the 16-orbit oracle panel: the nominated… | `r23_oracle_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **STable S81** | Sixty-day record of the interpolation family on the eight design-A orbits with an archived sixty-day truth, or… | `r20_longarc_detail_table.tex` | `rev22_supplement_tables.py` | R22 |
| **STable S72** | Complete LRO-like same-tolerance truth-error audit | `r10_blend_error_detail.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S73** | LRO-like propagation telemetry | `r10_blend_telemetry.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S34** | Numerical-resolution audit for every selected orbit and comparator | `r10_sobol_convergence_decisions.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S33** | Tight-to-tighter self-differences and truth-inclusive envelopes for the 17-orbit convergence subset, in meters | `r10_sobol_convergence_self.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S42** | Exact Sobol coordinates and transformed altitude coordinates for design A (propagated) | `r10_sobol_designA_coordinates.tex` | — | — |
| **STable S43** | Complete osculating elements and planned truth degree for design A (propagated) | `r10_sobol_designA_elements.tex` | — | — |
| **STable S44** | Exact initial Cartesian states for design A (propagated), in km and km s^ -1 | `r10_sobol_designA_states.tex` | — | — |
| **STable S45** | Exact Sobol coordinates and transformed altitude coordinates for design B (independently frozen before propaga… | `r10_sobol_designB_coordinates.tex` | — | — |
| **STable S46** | Complete osculating elements and planned truth degree for design B (independently frozen before propagation) | `r10_sobol_designB_elements.tex` | — | — |
| **STable S47** | Exact initial Cartesian states for design B (independently frozen before propagation), in km and km s^ -1 | `r10_sobol_designB_states.tex` | — | — |
| **STable S40** | Complete primary-policy audit for propagated design A | `r10_sobol_primary_audit.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S41** | Schedule-sensitivity and diagnostic audit for propagated design A | `r10_sobol_sensitivity_audit.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S71** | 28-day LRO-like (30 216~km polar) fixed-degree versus corrected-blend comparison at position/velocity-split ve… | `r11_blend_vector_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S36** | Independent design-B replication (64-orbit scrambled-Sobol coverage design, seed 20260724) at the same vector … | `r11_designB_convergence_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S38** | Per-orbit design-B vector-tolerance results, columns as in the design-A per-orbit table | `r11_designB_per_orbit_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S35** | Full 64-orbit design-A comparison at position/velocity-split vector tolerance | `r11_full64_convergence_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S37** | Per-orbit design-A vector-tolerance results | `r11_full64_per_orbit_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S18** | Initial-phase dispersion of the seven-day 50 300~km polar scheduling penalty at position/velocity-split vector… | `r11_phase_sweep_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S52** | Bin-resolution control for the Atallah benchmark | `r12_atallah_bincontrol_table.tex` | `rev12_atallah_bincontrol.py`<br>`rev12_finalize_manifest.py` | R12 |
| **STable S50** | Per-orbit matching record for the Atallah benchmark, design~A | `r12_atallah_matching_table_A.tex` | `rev12_finalize_manifest.py` | R12 |
| **STable S51** | Per-orbit matching record for the Atallah benchmark, design~B | `r12_atallah_matching_table_B.tex` | `rev12_finalize_manifest.py` | R12 |
| **STable S48** | Verification of the Atallah selection rule on JGGRX\_1800F (degree 600 verification field), against the accept… | `r12_atallah_verification_table.tex` | `rev12_atallah_verification.py`<br>`rev12_finalize_manifest.py` | R12 |
| **STable S54** | Integration-noise-free comparison at matched work | `r13_force_defect_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_force_defect.py` | R13 |
| **STable S53** | Diagnosis of the unresolved matched-work comparisons | `r13_resolution_diagnosis_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_resolution_diagnosis.py` | R13 |
| **STable S56** | Measured-time-matched comparator | `r13_timing_match_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_timing_match.py` | R13 |
| **STable S55** | Forced-variational calibration of the matched-work comparison | `r13_variational_check_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_variational_check.py` | R13 |
| **STable S67** | Per-orbit equal-budget comparison at = 1, design~A | `r14_beta1_per_orbit_A.tex` | `rev14_finalize_manifest.py` | R14 |
| **STable S68** | Per-orbit equal-budget comparison at = 1, design~B | `r14_beta1_per_orbit_B.tex` | `rev14_finalize_manifest.py` | R14 |
| **STable S60** | Degree-cap audit of the fixed-budget sweep | `r14_cap_audit_table.tex` | `rev14_tables.py` | R14 |
| **STable S61** | Does the per-call budget match survive the integrator? The per-call ratio is the declared budget match, N^2 _A… | `r14_cost_bookkeeping_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **STable S58** | Fixed-budget force-level comparison, full form | `r14_force_pareto_table_full.tex` | `rev14_tables.py` | R14 |
| **STable S66** | O26 Lagrangian-relaxed reference-path allocation benchmark | `r14_oracle_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **STable S59** | Measured-time budget control, full form | `r14_timing_budget_table_full.tex` | `rev14_timing_budget.py` | R14 |
| **STable S62** | Forced-variational mechanism check at equal budget | `r14_variational_budget_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_variational_budget.py` | R14 |
| **STable S49** | Grid convergence of the sampled maximum used to set the published rule's tolerance, and it does not converge | `r15_atallah_grid_table.tex` | `rev15_atallah_grid.py`<br>`rev15_finalize_manifest.py` | R15 |
| **STable S57** | Output-cadence convergence of the force-defect measurement | `r15_cadence_check_table.tex` | `rev15_cadence_check.py`<br>`rev15_finalize_manifest.py` | R15 |
| **STable S65** | Budget calibration without the reference arc | `r15_deployable_table.tex` | `rev15_finalize_manifest.py`<br>`rev15_tables.py` | R15 |
| **STable S64** | Three fixed comparators at = 1 | `r15_fixed_oracle_table.tex` | `rev15_finalize_manifest.py`<br>`rev15_tables.py` | R15 |
| **STable S23** | Per-orbit uniform-margin ladder for the 24-orbit design | `r8_alpha_margin_supplement.tex` | `make_figures_r8.py` | R8 |

## Generated but not used in the manuscript

These artifacts are produced by the pipeline but no manuscript item references them.

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **(not in manuscript)** | Propagated scrambled-Sobol design A | `r10_sobol_designA_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Scrambled-Sobol design B, frozen independently before propagation and subsequently propagated in full under th… | `r10_sobol_designB_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Truth-audited seven-day results for propagated Sobol design A | `r10_sobol_full_results_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Direct Atallah radial-adaptive benchmark on the independent design-B set, perilune-tolerance matched to the cr… | `r12_atallah_benchmark_table_designB.tex` | `rev12_finalize_manifest.py` | R12 |

