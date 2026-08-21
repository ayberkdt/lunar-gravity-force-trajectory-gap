# Archive index

What is in this archive, directory by directory. Generated from the tree itself by `build_archive_docs.py`; if a directory is listed as empty it is empty.

## Directories

| path | files | contents |
|---|---|---|
| `environment/` | 5 | interpreter specifications and dependency locks |
| `source/release_snapshot/` | 1 | pinned source snapshot of the kernel |
| `gravity_models/` | 1 | coefficient-product digests; the products themselves are not redistributed |
| `registrations/original/` | 43 | pre-registrations, one per campaign |
| `registrations/amendments/` | 6 | amendments, each naming its parent |
| `designs/` | 18 | frozen design coordinates and initial states |
| `designs/geometry_strata/` | 5 | frozen sub-box designs, including the four that did not reach a ladder |
| `raw/field_validation/` | 16 | field-level validation records |
| `raw/trajectory_validation/` | 105 | trajectory-level qualification records |
| `raw/calibration/` | 25 | calibration objective and transfer records |
| `raw/atallah_benchmark/` | 24 | published-rule implementation and benchmark |
| `raw/fixed_budget/` | 88 | fixed-budget allocation campaign |
| `raw/interpolation/` | 138 | interpolation-family campaign |
| `raw/operational_elliptical/` | 22 | operational elliptical population |
| `processed/aggregate_tables/` | 128 | aggregate tables as they appear in the PDF |
| `processed/per_orbit_tables/` | 19 | CSV renderings of the per-orbit tables that left the supplement |
| `processed/resolution_envelopes/` | 26 | numerical-resolution envelopes |
| `processed/work_and_timing/` | 117 | realized work and measured kernel time |
| `scripts/reproduce_tables/` | 36 | table generators |
| `scripts/reproduce_figures/` | 17 | figure generators |
| `scripts/validation/` | 13 | integrity and consistency checks |
| `manifests/` | 52 | campaign manifests, campaign map, checksums |
| `supplement_archive/` | 10 | what left the supplement PDF, and why |
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
| R42 | pre-registered, and inheriting a registration rather than restating on | R42: completion of the R37 forced-variational level chain to its last two levels, 56 and 64 per design (112 an |
| R44 | not recorded | R44 (O42): the interior span member (k = 0.5) against a constant degree matched on realized total quadratic wo |
| R48 | not recorded | R48 (O48): the interior span member (k = 0.5) at beta = 1 against a constant degree matched on measured serial |
| R50 | pre-registered | R50: a paired apolune ladder. Sixteen orbit identities per block, each flown at 300, 600, 1200 and 2400 km wit |
| R51 | pre-registered | R51: block A of the R50 paired radial-span ladder re-propagated with its adopted reference degree raised from  |
| R52 | pre-registered | R52: block B of the R50 paired radial-span ladder re-propagated with its adopted reference degree raised from  |
| R53 | declared before the first propagation in r53_preregistration | R53: the declared post-hoc budget beta = 0.62, already computed on designs A and B and on the two identity blo |
| R56 | not recorded | R56 (O56): the k = 0.5 interior member against a constant degree over sixty days, on the eight Design-A orbits |
| R57 | not recorded | R57: the three unfinished orbits of the (O40) gradient-degree audit, solved at reference degree 900, with R39' |
| R58 | not recorded | R58: the budget-calibrated radial endpoint against a constant degree matched on realized total quadratic work  |
| R59 | not recorded | R59: the budget comparison re-tallied over leading scrambled-Sobol prefixes of each coverage design, and the s |
| R60 | not recorded | R60: the editorial round's derivatives. The claims-ledger entries pinning the realized-work excess, the main-t |
| R61 | not recorded | R61 (O42-ext): the interior span member (k = 0.5) against a constant degree matched on realized total quadrati |
| R62 | not recorded | R62 (O54): the interior span member (k = 0.5) of the paired apolune ladder against a constant degree matched o |
| R63 | not recorded | R63 (O55): the interior span member (k = 0.5) of the ceiling-free apolune ladder against a constant degree mat |
| R64 | not recorded | R64 (O57): the k = 0.5 interior member against a constant degree matched on measured kernel time at the tighte |
| R65 | not recorded | R65 (O58): the sampled interior family, k = 0.25, 0.50 and 0.75, each against its own constant degree matched  |
| R66 | not recorded | R66: the main-text allocation-anatomy figure, drawn from the archived (O14) design-A beta = 1 record, that orb |
| R67 | not recorded | R67: the rank form of the R14 span/switch association and the perilune distribution of the R42 panel's undecid |
| R68 | not recorded | R68: every orbit of coverage designs A and B at beta = 1, each member against a constant degree refined until  |
| RJ | registered before propagation, with one disclosed staging: population  | J1-J3: pre-submission replication block for the JGCD retarget. J1 recalibrates the entire recipe on GSFC GRGM1 |

## Tables moved out of the supplement PDF

31 tables. Full record in `supplement_archive/migration_map.csv`.

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
| tab:screening-block |  | `supplement_archive/screening_block.tex` | `` |
| tab:practical-recommendations |  | `supplement_archive/practical_recommendations_table.tex` | `` |
| tab:related |  | `supplement_archive/adjacent_studies_table.tex` | `` |
| sec:supp-experiment-contract |  | `supplement_archive/experiment_chronology.tex` | `` |
| sec:reproducibility + sec:cost-curve |  | `supplement_archive/provenance_bookkeeping.tex` | `` |
| (O19)-(O30) prose entries |  | `supplement_archive/experiment_registry_prose_primary.tex` | `` |
| (O35)-(O60) prose entries |  | `supplement_archive/experiment_registry_prose.tex` | `` |
| fig:doe-regime |  | `supplement_archive/full_legacy_tables/` | `` |
| tab:canonical-geometries |  | `supplement_archive/full_legacy_tables/` | `` |
| tab:stage3-weekly |  | `supplement_archive/full_legacy_tables/` | `` |
| tab:sobol-main |  | `supplement_archive/full_legacy_tables/` | `` |
| S7.13 discarded-attempt forensics |  | `supplement_archive/timing_family_forensics.tex` | `` |

## Not in this archive

Raw trajectory arrays under `metrics/*_raw/` (about 14 GB), regenerable from the drivers and hashed in the campaign manifests, and the J1-J3 arrays (about 0.7 GB) held outside `metrics/` and hashed in `rJ_final_experiment_manifest.json`; about 15 GB together. Gravity-model coefficient files and SPICE kernels, which are distributed by their own archives and are recorded here by digest only.

