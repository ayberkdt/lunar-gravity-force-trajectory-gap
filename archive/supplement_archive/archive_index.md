# Archive index

What is in this archive, directory by directory. Generated from the tree itself by `build_archive_docs.py`; if a directory is listed as empty it is empty.

## Directories

| path | files | contents |
|---|---|---|
| `environment/` | 5 | interpreter specifications and dependency locks |
| `source/release_snapshot/` | 1 | pinned source snapshot of the kernel |
| `gravity_models/` | 1 | coefficient-product digests; the products themselves are not redistributed |
| `registrations/original/` | 21 | pre-registrations, one per campaign |
| `registrations/amendments/` | 5 | amendments, each naming its parent |
| `designs/` | 14 | frozen design coordinates and initial states |
| `designs/geometry_strata/` | 5 | frozen sub-box designs, including the four that did not reach a ladder |
| `raw/field_validation/` | 16 | field-level validation records |
| `raw/trajectory_validation/` | 100 | trajectory-level qualification records |
| `raw/calibration/` | 25 | calibration objective and transfer records |
| `raw/atallah_benchmark/` | 24 | published-rule implementation and benchmark |
| `raw/fixed_budget/` | 47 | fixed-budget allocation campaign |
| `raw/interpolation/` | 60 | interpolation-family campaign |
| `raw/operational_elliptical/` | 22 | operational elliptical population |
| `processed/aggregate_tables/` | 89 | aggregate tables as they appear in the PDF |
| `processed/per_orbit_tables/` | 19 | CSV renderings of the per-orbit tables that left the supplement |
| `processed/resolution_envelopes/` | 21 | numerical-resolution envelopes |
| `processed/work_and_timing/` | 37 | realized work and measured kernel time |
| `scripts/reproduce_tables/` | 21 | table generators |
| `scripts/reproduce_figures/` | 13 | figure generators |
| `scripts/validation/` | 11 | integrity and consistency checks |
| `manifests/` | 31 | campaign manifests, campaign map, checksums |
| `supplement_archive/` | 2 | what left the supplement PDF, and why |
| `supplement_archive/full_legacy_tables/` | 19 | the moved tables as rendered LaTeX, byte-identical to what the supplement compiled |

## Campaign manifests

| campaign | registration | scope |
|---|---|---|
| R10 | not recorded | SHA-256 integrity index for the formal R10 confirmatory evidence and the manuscript state that reports it |
| R11 | not recorded | R11 extension: full-population design-A vector-tolerance rerun, independent design-B replication, 24-phase dis |
| R12 | not recorded | R12 extension: faithful implementation and verification of the Atallah (2022) analytical radial-adaptive rule, |
| R13 | not recorded | R13 extension: diagnosis of the unresolved matched-work Atallah comparisons, targeted third-tolerance-level re |
| R14 | not recorded | R14 extension (O25/O26): fixed-budget radial-allocation Pareto study. Pre-registered protocol; 128-orbit integ |
| R15 | not recorded | R15 audit response: output-cadence convergence of the force-defect statistic (R15-D), deployable truth-free bu |
| R16 | not recorded | R16 transfer test: the static-selection calibration repeated on two further lunar solutions (cross-solution re |
| R17 | not recorded | R17 sixty-day long arcs on a widened geometry set, with a two-level vector-tolerance ladder supplying a numeri |
| R18 | not recorded | R18 span sweep: a one-parameter family interpolating geometrically between the equal-budget constant degree an |
| R19 | not recorded | R19: the interior span member (k = 0.5) against a constant degree matched on realized total quadratic work rat |
| R20 | not recorded | R20: the constant-to-radial interpolation family propagated for 60 days at beta = 1 on the eight design-A orbi |
| R21 | not recorded | R21: relative truncation of the reference gravity gradient at the variational solve's degree-120 evaluation, a |
| R22 | not recorded | R22: two supplementary tables derived from archived R18 and R20 records. Propagates nothing; reads the frozen  |
| R23 | not recorded | R23 (O31): the three registered controls on the interior span member -- the realized-work comparison repeated  |
| R24 | not recorded | R24 (O32): two registered controls. O32a re-runs the O31b oracle panel at the third tolerance level, reusing t |
| R25 | not recorded | R25 (O33): the interior-member comparison at the midpoint budget beta = 0.75, run to locate the sign change th |
| R28 | NOT pre-registered | R28 (O34): the post-hoc midpoint beta = 0.62. One Phase-A calibration point added to the frozen grid after the |
| R29 | pre-registered | R29: design C, the third scrambled-Sobol coverage design frozen by R26. Its base, its regenerated accuracy-tar |
| R30 | pre-registered | R30: geometry strata of the frozen factor box. Each stratum is 64 orbits drawn with the pinned orbit map on a  |
| R31 | pre-registered | R31: one 64-orbit population of operational elliptical lunar orbits, perilune 80-120 km with apolune 700-2500  |
| R35 | not pre-registered | R34-R36: derived products over frozen records. R34 assembles the four-rung instrument ladder at beta = 1 -- de |
| R37 | pre-registered | R37: extension of the forced-variational panel of Section 7.2 from the archived eight orbits toward both full  |
| R38 | pre-registered | R38: the R31 operational elliptical population re-propagated with its adopted reference degree raised from 300 |
| R39 | pre-registered | R39: the gradient-degree audit of the enlarged forced-variational panel (O40). The forced solve of (O25) and ( |
| R41 | pre-registered | R41: the reference-degree control (O41). Every truncation error in the paper is measured against an adopted re |

## Tables moved out of the supplement PDF

19 tables. Full record in `supplement_archive/migration_map.csv`.

| table | rows | archived LaTeX | CSV |
|---|---|---|---|
| tab:sobol-convergence-self | 18 | `supplement_archive/full_legacy_tables/r10_sobol_convergence_self.tex` | `processed/per_orbit_tables/r10_sobol_convergence_self.csv` |
| tab:sobol-convergence-decisions | 35 | `supplement_archive/full_legacy_tables/r10_sobol_convergence_decisions.tex` | `processed/per_orbit_tables/r10_sobol_convergence_decisions.csv` |
| tab:full64-per-orbit | 65 | `supplement_archive/full_legacy_tables/r11_full64_per_orbit_table.tex` | `processed/per_orbit_tables/r11_full64_per_orbit_table.csv` |
| tab:designB-per-orbit | 65 | `supplement_archive/full_legacy_tables/r11_designB_per_orbit_table.tex` | `processed/per_orbit_tables/r11_designB_per_orbit_table.csv` |
| tab:sobol-primary-audit | 65 | `supplement_archive/full_legacy_tables/r10_sobol_primary_audit.tex` | `processed/per_orbit_tables/r10_sobol_primary_audit.csv` |
| tab:sobol-sensitivity-audit | 65 | `supplement_archive/full_legacy_tables/r10_sobol_sensitivity_audit.tex` | `processed/per_orbit_tables/r10_sobol_sensitivity_audit.csv` |
| tab:sobol-designA-coordinates | 65 | `supplement_archive/full_legacy_tables/r10_sobol_designA_coordinates.tex` | `processed/per_orbit_tables/r10_sobol_designA_coordinates.csv` |
| tab:sobol-designB-coordinates | 65 | `supplement_archive/full_legacy_tables/r10_sobol_designB_coordinates.tex` | `processed/per_orbit_tables/r10_sobol_designB_coordinates.csv` |
| tab:blend-error-detail | 5 | `supplement_archive/full_legacy_tables/r10_blend_error_detail.tex` | `processed/per_orbit_tables/r10_blend_error_detail.csv` |
| tab:blend-telemetry | 7 | `supplement_archive/full_legacy_tables/r10_blend_telemetry.tex` | `processed/per_orbit_tables/r10_blend_telemetry.csv` |
| tab:alpha-margin-per-orbit | 121 | `supplement_archive/full_legacy_tables/r8_alpha_margin_supplement.tex` | `processed/per_orbit_tables/r8_alpha_margin_supplement.csv` |
| tab:phase-sweep-vector | 5 | `supplement_archive/full_legacy_tables/r11_phase_sweep_table.tex` | `processed/per_orbit_tables/r11_phase_sweep_table.csv` |
| tab:atallah-matching-A | 65 | `supplement_archive/full_legacy_tables/r12_atallah_matching_table_A.tex` | `processed/per_orbit_tables/r12_atallah_matching_table_A.csv` |
| tab:atallah-matching-B | 65 | `supplement_archive/full_legacy_tables/r12_atallah_matching_table_B.tex` | `processed/per_orbit_tables/r12_atallah_matching_table_B.csv` |
| tab:budget-per-orbit-A | 65 | `supplement_archive/full_legacy_tables/r14_beta1_per_orbit_A.tex` | `processed/per_orbit_tables/r14_beta1_per_orbit_A.csv` |
| tab:budget-per-orbit-B | 65 | `supplement_archive/full_legacy_tables/r14_beta1_per_orbit_B.tex` | `processed/per_orbit_tables/r14_beta1_per_orbit_B.csv` |
| tab:span-detail-A | 66 | `supplement_archive/full_legacy_tables/r18_span_detail_table_A.tex` | `processed/per_orbit_tables/r18_span_detail_table_A.csv` |
| tab:span-detail-B | 66 | `supplement_archive/full_legacy_tables/r18_span_detail_table_B.tex` | `processed/per_orbit_tables/r18_span_detail_table_B.csv` |
| tab:span-longarc-detail | 10 | `supplement_archive/full_legacy_tables/r20_longarc_detail_table.tex` | `processed/per_orbit_tables/r20_longarc_detail_table.csv` |

## Not in this archive

Raw trajectory arrays (about 2.1 GB), regenerable from the drivers and hashed in the campaign manifests. Gravity-model coefficient files and SPICE kernels, which are distributed by their own archives and are recorded here by digest only.

