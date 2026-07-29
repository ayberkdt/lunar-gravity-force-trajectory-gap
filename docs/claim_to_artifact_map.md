# Claim-to-artifact map

Which script produced which manuscript item, and which campaign manifest covers it. Generated from the LaTeX sources, the compiled label numbering and the campaign manifests; the machine-readable version is [`REPRODUCIBILITY_INDEX.csv`](../REPRODUCIBILITY_INDEX.csv).

Analysis and table passes read the archived records in `metrics/` in place, so most rows can be re-run without propagating anything. Rows marked *(inline)* are tables typed directly in the manuscript from numbers reported in the text; they have no generated file.

## Main text

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **Figure 1** | Per-coefficient RMS spectrum of JGGRX_1800F with power-law fits over the observationally dominated band (nin[10,600])… | `fig_spectrum.pdf` | `make_figures_r1.py` | R1 |
| **Figure 2** | Truncation criteria against the model-relative empirical reference (log--log) | `fig_truncation.pdf` | `make_figures_r1.py` | R1 |
| **Figure 3** | 24-hour RMS position error versus truncation degree for the 100 km circular polar orbit at two tolerance settings, with… | `fig_orbit_error.pdf` | `make_figures_r1.py` | R1 |
| **Figure 4** | Signed (upper) and absolute-log (lower) relative specific-energy error over the band-crossing orbit. The corrected… | `fig_blend_energy.pdf` | `make_figures_r2.py` | R2 |
| **Figure 5** | Confirmatory scrambled-Sobol design. Left and center: ordered truth-degree-audited error ratios against work-matched… | `fig_sobol_confirmatory.pdf` | `rev10_finalize_manifest.py`<br>`rev10_manuscript_assets.py` | R10 |
| **Figure 6** | Force-defect curves and propagated trajectory operating points for radial and constant-degree truncation, one column… | `budget_pareto.pdf` | `make_figures_r14.py`<br>`rev14_deliverables.py`<br>`rev14_finalize_manifest.py` | R14 |
| **Figure 7** | Equal-budget error ratio ρ=E_fixed/E_radial predicted by the forced variational solution against the propagated value,… | `fig_variational_parity.pdf` | `make_figures_r19.py` | R19 |
| **Figure 8** | The interpolation path at the critical degree's own nominal per-call budget (β=1), plotted against the degree span each… | `fig_span_curve.pdf` | `make_figures_r20.py` | R20 |
| **Table 1** | Adjacent degree-selection studies and the aspects measured in each. "Local" = pointwise or per-direction error… | *(inline)* | — | — |
| **Table 2** | What carries what. Each claim of the paper is listed against the single measurement it rests on and the section… | *(inline)* | — | — |
| **Table 3** | Recommended truncation degree: model-relative empirical reference and closed-form rules. The spectrum-tail columns use… | *(inline)* | — | — |
| **Table 4** | Practical static selections. The compact fit summarizes the empirical trend but does not preserve the budget at every… | `review_static_selection_table.tex` | `review_postprocess.py` | — |
| **Table 5** | Transfer of the calibration procedure. N_max is the highest degree present in the distributed file. p_spec is the OLS… | `r16_transfer_table.tex` | `rev16_finalize_manifest.py`<br>`rev16_tables.py` | R16 |
| **Table 6** | Physics of four degree-transition policies on a common band. Position error is against fixed N=120. Per-path… | *(inline)* | — | — |
| **Table 7** | 7-day eccentric arc (50×300 km polar), tight tolerance: position error against the N=300 adopted reference, quadratic… | *(inline)* | — | — |
| **Table 8** | Direct benchmark of the published Atallah radial-adaptive rule on both 64-orbit scrambled-Sobol designs, at the… | `r12_atallah_benchmark_table.tex` | `rev12_finalize_manifest.py` | R12 |
| **Table 9** | Fixed-budget force-level comparison; each cell pairs design A with design B so the replication reads across. At each… | `r14_force_pareto_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **Table 10** | Equal-budget comparison with the budget defined as measured serial gravity-kernel time rather than the quadratic proxy.… | `r14_timing_budget_table_compact.tex` | `rev14_timing_budget.py` | R14 |
| **Table 11** | Propagated fixed-budget comparison. Seven-day Cartesian position RMS against the common adopted truth at the tight… | `r14_trajectory_pareto_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **Table 12** | The same path in numbers, at β=1. Span is the max-to-min degree ratio of the propagated schedule; error is the… | `r18_span_table.tex` | `rev18_finalize_manifest.py`<br>`rev18_tables.py` | R18 |
| **Table 13** | The interior member (k=0.5) against the constant degree under two definitions of an equal budget, at two declared… | `r19_equal_work_table.tex` | `rev19_finalize_manifest.py`<br>`rev19_tables.py` | R19 |
| **Table 14** | The β=1 realized-work comparison re-run at a third, ultra-tight tolerance level. The panel is every comparison that… | `r23_ultra_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **Table 15** | The interpolation path at three declared budgets. "Best k" is the member with the lowest measured error on the most… | `r18_budget_table.tex` | `rev18_budget_table.py` | R18 |
| **Table 16** | Practical recommendations within the scope tested here | *(inline)* | — | — |

## Supplement

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **SFigure S9** | Uniform margin control on the 24-orbit design. Left: median error ratios ρ_crit and re-matched ρ_work^α of the inflated… | `fig_alpha_margin.pdf` | `make_figures_r8.py` | R8 |
| **SFigure S11** | The published radial-adaptive rule against both fixed comparators, design A (circles) and design B (triangles). (a)… | `fig_atallah_benchmark.pdf` | `make_figures_r13.py`<br>`rev13_finalize_manifest.py` | R13 |
| **SFigure S1** | Dense-altitude calibration objective for the effective tail-budget exponent on JGGRX_1800F. The shaded interval is the… | `fig_pstar_objective.pdf` | `make_figures_supplemental.py` | — |
| **SFigure S2** | Measured single-thread evaluation time of the serial kernel (median and interquartile range over nine timing blocks) | `fig_timing.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S3** | Two orbital periods of the eccentric arc: altitude, scheduled degree (ε=10^-3, floor 60), and accepted DOP853 step… | `fig_schedule_panels.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S4** | Regime map of the stratified 24-orbit design. Left and center: the best-of-three-schedule ⟨ N^2⟩-matched ratio ρ_work… | `fig_doe_regime.pdf` | `make_figures_doe_regime.py` | — |
| **SFigure S5** | Event-aligned accepted-step median and interquartile range (upper panels) and direct rejection probability (lower… | `fig_switch_aggregate.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S6** | Seven-day RMS position error for six degree policies across four orbit/orientation controls. The logarithmic scale… | `fig_longarc_matrix.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S7** | Dense fixed-degree sweep at the tight tolerance for the reference phase, apolune-start phase, and 60^circ-inclination… | `fig_degree_sweep.pdf` | `make_figures_r3.py` | R3 |
| **SFigure S8** | Position-error envelopes over the tested 7-day arc (tight tolerance). Scheduled runs grow superlinearly; the tested… | `fig_longarc_growth.pdf` | `make_figures_r1.py` | R1 |
| **SFigure S10** | One-revolution PEFRL energy control. Fixed-degree energy remains at the numerical floor, whereas changing the… | `fig_symplectic_switch.pdf` | `rev4_robustness_controls.py` | R4 |
| **STable S1** | Field-level acceleration comparison against SHTOOLS, pooled over 30, 100, and 500 km and including near-pole points.… | *(inline)* | — | — |
| **STable S2** | Reproducibility contract of the Tudat--Lunaris comparison. The SPICE kernel set is naif0012.tls, pck00011.tpc,… | *(inline)* | — | — |
| **STable S3** | Three-resolution Tudat RK4 controls (5, 2.5, and 1.25 s) on the common 60-s comparison grid. Altitude, latitude, and… | *(inline)* | — | — |
| **STable S21** | Uniform-margin control, design aggregates over the 24-orbit design. ρ_work^α is the median error ratio against the… | `r8_alpha_margin_table.tex` | `make_figures_r8.py` | R8 |
| **STable S22** | Rank-resolution counts for the aggregate margin analysis. Each family contains 24 orbits. Raw schedule wins ignore… | `r8_alpha_margin_resolution_table.tex` | `make_figures_r8.py` | R8 |
| **STable S69** | Twenty-eight-day Cartesian position RMS (m) against the common high-degree reference for the smooth N=30/120 transition… | `r9_potential_blend_longarc_table.tex` | `rev9_blend_postprocess.py` | R9 |
| **STable S70** | LRO-like 28-day convergence of fixed N=120 and the corrected potential blend. All entries are meters. Policy self is… | `r10_blend_convergence_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S63** | Gradient-degree sensitivity of the forced variational solve, over 241 sampled epochs of each panel orbit's archived… | `r21_gradient_sensitivity_table.tex` | `rev21_finalize_manifest.py`<br>`rev21_gradient_sensitivity.py` | R21 |
| **STable S5** | Verification of the per-degree acceleration-RMS identity by exact Gauss--Legendre quadrature (nodes n+10 × 2n+12;… | *(inline)* | — | — |
| **STable S6** | The kernel cost curve re-measured under both archived protocols at both archived model degrees, on an idle machine, in… | `r23_cost_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **STable S7** | Empirical minimum degree N^emp and the p_fit proxy at four points of the dimensionless band, at ε=10^-2. Dashes mark… | `r16_transfer_detail_table.tex` | `rev16_finalize_manifest.py`<br>`rev16_tables.py` | R16 |
| **STable S8** | Normalized band-difference RMS amplitudes (1000 directions per altitude; bootstrap 95% intervals in brackets). Bands… | *(inline)* | — | — |
| **STable S9** | Degree recommendations from the sampled sphere-RMS criterion and directional tail-amplitude percentiles normalized by… | *(inline)* | — | — |
| **STable S4** | Integration-tolerance regimes and the results quoted from each. Absolute tolerances are scalar unless a… | *(inline)* | — | — |
| **STable S74** | Sixty-day comparison on the widened geometry set. N_crit is the critical-altitude empirical degree; "fixed crit." and… | `r17_longarc60_table.tex` | `rev17_finalize_manifest.py`<br>`rev17_tables.py` | R17 |
| **STable S75** | Sixty-day growth by checkpoint (meters). "cum. RMS" is the cumulative Cartesian position RMS over [0,d]; "in-track" is… | `r17_longarc60_growth_table.tex` | `rev17_finalize_manifest.py`<br>`rev17_tables.py` | R17 |
| **STable S76** | Twenty-eight-day Cartesian position RMS E_RMS,28 (m) against a common fixed high-degree reference trajectory (N=600 for… | *(inline)* | — | — |
| **STable S10** | Seven-day N=300 minus N=600 reference-degree control (meters), using vector absolute tolerances of 10^-6 m and 10^-9 m… | *(inline)* | — | — |
| **STable S11** | Narrow-gap precision controls. RMS, maximum, and final R/I/C compare each tighter policy run with the tighter N=300… | *(inline)* | — | — |
| **STable S12** | Canonical geometries in the 24-orbit design. Osculating elements are defined at t=0 with true anomaly ν_0=0^circ. For… | *(inline)* | — | — |
| **STable S13** | Per-switch acceleration jump: direction-averaged RMS of \|a(120)-a(120-q)\| relative to the perturbation RMS (400… | *(inline)* | — | — |
| **STable S14** | Direct DOP853 switching telemetry over 2.2 orbital periods. The controlled-ODE validation satisfies n_rm RHS=2+12n_rm… | *(inline)* | — | — |
| **STable S15** | Scalar-tolerance screening convergence; not used for close policy ranking. Directly measured seven-day differences… | *(inline)* | — | — |
| **STable S16** | Exploratory scalar-tolerance multi-geometry screening matrix; final narrow-gap rankings use vector tolerances. Final… | *(inline)* | — | — |
| **STable S17** | Initial-phase scalar-tolerance screening sweep of the 7-day 50×300 km polar comparison: RMS position error (m) against… | *(inline)* | — | — |
| **STable S19** | Weekly checkpoints recomputed from the archived 28-day records for the schedule selected by minimum full-arc position… | `review_stage3_checkpoints.tex` | `review_postprocess.py` | — |
| **STable S20** | Fraction of the seven-day arc spent at each active degree under event-resolved downward scheduling. The LRO-like orbit… | *(inline)* | — | — |
| **STable S24** | Seven-day robustness controls: Cartesian RMS position error against the same-solver, same-force N=300 reference. Dashes… | *(inline)* | — | — |
| **STable S25** | LRO-like seven-day accuracy controls against the same-force N=300 reference. R/I/C are signed final radial, in-track,… | *(inline)* | — | — |
| **STable S26** | LRO-like computational controls. ⟨ N^2⟩-matched fixed degrees are the nearest integers to sqrt⟨ N^2⟩ of the… | *(inline)* | — | — |
| **STable S27** | Seven-day rotating-frame Jacobi-integral control. For each geometry the fixed-degree relative drift \|Δ H/H\|_max under… | *(inline)* | — | — |
| **STable S28** | Independent TudatPy reproduction of dynamic scheduling: seven-day Cartesian RMS position error against each code's own… | *(inline)* | — | — |
| **STable S29** | LRO-like Tudat step-convergence. Successive maximum-step differences (RMS position, m) for the N=300 reference and each… | *(inline)* | — | — |
| **STable S30** | Confirmatory truth-degree audit, frozen six-orbit primary set. The threshold is the stricter of 5 m and 5% of the… | `r10_truth_audit_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S31** | Extended truth-degree audit of the four remaining sub-50 km survivors, using the same acceptance rule. All four adopt… | `r10_truth_audit_extended_table.tex` | `rev10_truth_audit_extended_postprocess.py` | R10 |
| **STable S32** | Pre-specified 17-orbit convergence subset and the selection category (or categories) that placed each orbit in it: (A)… | *(inline)* | — | — |
| **STable S39** | Aggregate confirmatory counts of the initial screening pass. Eligible and raw counts use all 64 truth-surviving… | `r10_sobol_main_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S77** | The continuation family member by member, at every declared budget. Span, error and the two cost ratios are medians… | `r18_by_k_table.tex` | `rev22_supplement_tables.py` | R22 |
| **STable S78** | Design A: seven-day position RMS (m) at the tighter level for each member of the family, by orbit. k=0 and k=1 are the… | `r18_span_detail_table_A.tex` | `rev18_finalize_manifest.py` | R18 |
| **STable S79** | Design B: seven-day position RMS (m) at the tighter level for each member of the family, by orbit. k=0 and k=1 are the… | `r18_span_detail_table_B.tex` | `rev18_finalize_manifest.py` | R18 |
| **STable S80** | The interior member (k=0.5) at β=1 against two constant comparators on the 16-orbit oracle panel: the nominated… | `r23_oracle_table.tex` | `rev23_finalize_manifest.py`<br>`rev23_tables.py` | R23 |
| **STable S81** | Sixty-day record of the interpolation family on the eight design-A orbits with an archived sixty-day truth, ordered by… | `r20_longarc_detail_table.tex` | `rev22_supplement_tables.py` | R22 |
| **STable S72** | Complete LRO-like same-tolerance truth-error audit | `r10_blend_error_detail.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S73** | LRO-like propagation telemetry | `r10_blend_telemetry.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S34** | Numerical-resolution audit for every selected orbit and comparator | `r10_sobol_convergence_decisions.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S33** | Tight-to-tighter self-differences and truth-inclusive envelopes for the 17-orbit convergence subset, in meters | `r10_sobol_convergence_self.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S42** | Exact Sobol coordinates and transformed altitude coordinates for design A (propagated) | `r10_sobol_designA_coordinates.tex` | — | — |
| **STable S43** | Complete osculating elements and planned truth degree for design A (propagated) | `r10_sobol_designA_elements.tex` | — | — |
| **STable S44** | Exact initial Cartesian states for design A (propagated), in km and km s^ -1 | `r10_sobol_designA_states.tex` | — | — |
| **STable S45** | Exact Sobol coordinates and transformed altitude coordinates for design B (independently frozen before propagation) | `r10_sobol_designB_coordinates.tex` | — | — |
| **STable S46** | Complete osculating elements and planned truth degree for design B (independently frozen before propagation) | `r10_sobol_designB_elements.tex` | — | — |
| **STable S47** | Exact initial Cartesian states for design B (independently frozen before propagation), in km and km s^ -1 | `r10_sobol_designB_states.tex` | — | — |
| **STable S40** | Complete primary-policy audit for propagated design A | `r10_sobol_primary_audit.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S41** | Schedule-sensitivity and diagnostic audit for propagated design A | `r10_sobol_sensitivity_audit.tex` | `rev10_manuscript_assets.py` | R10 |
| **STable S71** | 28-day LRO-like (30×216 km polar) fixed-degree versus corrected-blend comparison at position/velocity-split vector… | `r11_blend_vector_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S36** | Independent design-B replication (64-orbit scrambled-Sobol coverage design, seed 20260724) at the same vector tolerance | `r11_designB_convergence_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S38** | Per-orbit design-B vector-tolerance results, columns as in the design-A per-orbit table | `r11_designB_per_orbit_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S35** | Full 64-orbit design-A comparison at position/velocity-split vector tolerance. All 64 orbits were rerun under the… | `r11_full64_convergence_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S37** | Per-orbit design-A vector-tolerance results | `r11_full64_per_orbit_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S18** | Initial-phase dispersion of the seven-day 50×300 km polar scheduling penalty at position/velocity-split vector… | `r11_phase_sweep_table.tex` | `rev11_finalize_manifest.py`<br>`rev11_manuscript_tables.py` | R11 |
| **STable S52** | Bin-resolution control for the Atallah benchmark. For 10 design-A orbits spanning the perilune range, the published… | `r12_atallah_bincontrol_table.tex` | `rev12_atallah_bincontrol.py`<br>`rev12_finalize_manifest.py` | R12 |
| **STable S50** | Per-orbit matching record for the Atallah benchmark, design A | `r12_atallah_matching_table_A.tex` | `rev12_finalize_manifest.py` | R12 |
| **STable S51** | Per-orbit matching record for the Atallah benchmark, design B | `r12_atallah_matching_table_B.tex` | `rev12_finalize_manifest.py` | R12 |
| **STable S48** | Verification of the Atallah selection rule on JGGRX_1800F (degree 600 verification field), against the acceptance… | `r12_atallah_verification_table.tex` | `rev12_atallah_verification.py`<br>`rev12_finalize_manifest.py` | R12 |
| **STable S54** | Integration-noise-free comparison at matched work. Along each orbit's archived truth trajectory, the truncation… | `r13_force_defect_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_force_defect.py` | R13 |
| **STable S53** | Diagnosis of the unresolved matched-work comparisons. M_res = \|E_At-E_fix\|/(E_num,At+E_num,fix) is the resolution… | `r13_resolution_diagnosis_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_resolution_diagnosis.py` | R13 |
| **STable S56** | Measured-time-matched comparator. N_work is the quadratic-proxy comparator of the campaign and N_time the fixed degree… | `r13_timing_match_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_timing_match.py` | R13 |
| **STable S55** | Forced-variational calibration of the matched-work comparison. For each orbit one augmented integration along the… | `r13_variational_check_table.tex` | `rev13_finalize_manifest.py`<br>`rev13_variational_check.py` | R13 |
| **STable S67** | Per-orbit equal-budget comparison at = 1, design A | `r14_beta1_per_orbit_A.tex` | `rev14_finalize_manifest.py` | R14 |
| **STable S68** | Per-orbit equal-budget comparison at = 1, design B | `r14_beta1_per_orbit_B.tex` | `rev14_finalize_manifest.py` | R14 |
| **STable S60** | Degree-cap audit of the fixed-budget sweep. Wherever the calibrated radial degree reaches the orbit's adopted truth… | `r14_cap_audit_table.tex` | `rev14_tables.py` | R14 |
| **STable S61** | Does the per-call budget match survive the integrator? The per-call ratio is the declared budget match, ⟨ N^2⟩_A /… | `r14_cost_bookkeeping_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **STable S58** | Fixed-budget force-level comparison, full form. The compact version in the main text carries the ordering and the… | `r14_force_pareto_table_full.tex` | `rev14_tables.py` | R14 |
| **STable S66** | O26 Lagrangian-relaxed reference-path allocation benchmark. At each budget the relaxation N_λ(t)=argmin_N[d(t,N)+λ N^2]… | `r14_oracle_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_tables.py` | R14 |
| **STable S59** | Measured-time budget control, full form. The compact version in the main text carries the matched degree and the two… | `r14_timing_budget_table_full.tex` | `rev14_timing_budget.py` | R14 |
| **STable S62** | Forced-variational mechanism check at equal budget. Both policies are carried along one shared reference under one… | `r14_variational_budget_table.tex` | `rev14_finalize_manifest.py`<br>`rev14_variational_budget.py` | R14 |
| **STable S49** | Grid convergence of the sampled maximum used to set the published rule's tolerance, and it does not converge. The… | `r15_atallah_grid_table.tex` | `rev15_atallah_grid.py`<br>`rev15_finalize_manifest.py` | R15 |
| **STable S57** | Output-cadence convergence of the force-defect measurement. Every force-level statistic in this paper is taken on the… | `r15_cadence_check_table.tex` | `rev15_cadence_check.py`<br>`rev15_finalize_manifest.py` | R15 |
| **STable S65** | Budget calibration without the reference arc. The O25 tolerance is found by bisecting on the archived reference… | `r15_deployable_table.tex` | `rev15_finalize_manifest.py`<br>`rev15_tables.py` | R15 |
| **STable S64** | Three fixed comparators at β = 1. N_sat is the budget-saturating degree, the largest whose per-call quadratic work fits… | `r15_fixed_oracle_table.tex` | `rev15_finalize_manifest.py`<br>`rev15_tables.py` | R15 |
| **STable S23** | Per-orbit uniform-margin ladder for the 24-orbit design | `r8_alpha_margin_supplement.tex` | `make_figures_r8.py` | R8 |

## Generated but not used in the manuscript

Produced by the pipeline, but no manuscript item references them.

| Item | Claim | Artifact | Script | Campaign |
|---|---|---|---|---|
| **(not in manuscript)** | Propagated scrambled-Sobol design A | `r10_sobol_designA_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Scrambled-Sobol design B, frozen independently before propagation and subsequently propagated in full under the… | `r10_sobol_designB_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Truth-audited seven-day results for propagated Sobol design A | `r10_sobol_full_results_table.tex` | `rev10_manuscript_assets.py` | R10 |
| **(not in manuscript)** | Direct Atallah radial-adaptive benchmark on the independent design-B set, perilune-tolerance matched to the critical… | `r12_atallah_benchmark_table_designB.tex` | `rev12_finalize_manifest.py` | R12 |

