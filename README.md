# lunar-gravity-force-trajectory-gap

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21824029-blue.svg)](https://doi.org/10.5281/zenodo.21824029)

Reproducibility archive for fixed-budget spherical-harmonic degree allocation
and the lunar gravity force–trajectory gap, on the degree-1800 GRAIL
JGGRX\_1800F field.

The study asks how a prescribed nominal per-call spherical-harmonic budget
should be spread along an eccentric lunar arc, and measures where a smaller
trajectory-averaged truncation-force defect stops predicting a smaller
trajectory error.

This archive holds the measurement instrument, the campaign records, the
manifests and the drivers needed to audit every reported result, from the
field-level calibrations through campaign R68 and the J-series cross-field
replication. All 824 recorded digests currently match the files they name,
under the default audit and under `--strict`. The only things fetched from outside are public archive
data products, and `fetch_data.py` retrieves and checksums those for you.

The archived deposit is at <https://doi.org/10.5281/zenodo.21824029>.

## How to read the scientific claims

The archive contains results with different evidential status. They should not
all be read as if they were one confirmatory experiment.

- The two fixed-budget endpoint comparisons are the confirmatory core. Their
  designs, metrics, censoring rules, numerical-resolution rule, and decision
  logic were frozen before their aggregate results were inspected.
- The allocation family between the constant and radial endpoints is
  exploratory. The reported intermediate member was identified after its
  sweep, and the archive labels it accordingly.
- Geometry-localizing populations and the wide-elliptic population are
  registered scope extensions. They test where the endpoint ranking changes;
  they are not probability samples of operational lunar missions.
- Post-hoc controls and scoring amendments remain in the archive with that
  status rather than being promoted to confirmatory evidence.
- The J-series campaigns answer the three standing objections directly: J1
  repeats the comparison on the independent GSFC GRGM1200A solution with its
  own calibration and a fresh population, J2 repeats it under full lunar
  dynamics rather than the isolated gravity-only system, and J3 probes how the
  verdict counts move when the numerical-resolution rule itself is varied.

Every error here is **model-relative**: a truncated evaluation is compared
with a higher-degree evaluation of the same adopted coefficient product. The
field named `truth` in historical records means that adopted numerical
reference, not the unknown true lunar gravity field.

Likewise, “fixed budget” usually means declared nominal per-call work,
`<N^2>`. Realized total quadratic work, kernel time, wall time, right-hand-side
calls, and accepted integrator steps are distinct measurements and are
reported separately. A smaller force-defect norm is not treated as proof of a
smaller trajectory error; that gap is the phenomenon being measured.

The results therefore do **not** establish a universally optimal allocation,
a flight-ready onboard rule, or an absolute mission error bound. They establish
comparisons within the recorded gravity products, policies, geometries,
budgets, force models, integration settings, and arc lengths. The exact status
of each claim is in
[`archive/manifests/campaign_map.csv`](archive/manifests/campaign_map.csv).

## Quick check

Bash:

```bash
pip install -r requirements-smoke-test.txt
python verify_snapshot.py
```

PowerShell:

```powershell
pip install -r requirements-smoke-test.txt
python verify_snapshot.py
```

This needs no download. It loads the Lunar Prospector LP150Q product shipped in
`data/`, evaluates nine spherical-harmonic accelerations across three degrees
and three field points, and compares them against values recorded from the
archived snapshot used for the submitted manuscript. If it passes, those
representative evaluations reproduce within the declared tolerance; it is not
a complete campaign rerun.

The default compares against a declared relative tolerance of `1e-12`. Exact
reproduction to the last bit depends on the LLVM version Numba compiles
through, the CPU's fused-multiply-add behaviour and the platform's libm, so it
is expected on the archived Windows 11 x64 environment and not guaranteed
elsewhere. To assert it there:

```bash
python verify_snapshot.py --strict-bitwise
```

## Maintaining the archive safely

This is an evidence archive, so historical experiment files are not tidied in
place. Existing files under `scripts/`, `metrics/`, `figures/`, `archive/`,
`experiments/`, `environments/`, and `data/` are append-only. The measurement
snapshot under `src/` is fully frozen. Formatting, renaming, or refactoring an
old file can invalidate the digest chain even when its numerical behaviour is
unchanged.

Maintenance work belongs in `tools/`, `tests/`, `.github/workflows/`, and the
explanatory parts of `docs/`. New campaigns should add newly named artifacts
and their own sealed manifest instead of overwriting an earlier campaign.

Before committing maintenance work, run:

```bash
python -m unittest discover -s tests -v
python tools/check_frozen_changes.py --base HEAD --working-tree
python tools/build_source_manifest.py --check
```

The archive guard also runs automatically on pull requests and pushes to
`main`. It rejects edits, removals, and renames of existing evidence while
allowing new campaign artifacts. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for
the contribution rules and [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for
the frozen-evidence/maintenance-layer boundary.

## Where does a given result come from?

[`docs/REPRODUCIBILITY_INDEX.csv`](docs/REPRODUCIBILITY_INDEX.csv) maps every
manuscript figure and table to the artifact file, the script that writes it,
and the campaign manifest that covers it. The readable version is
[`docs/claim_to_artifact_map.md`](docs/claim_to_artifact_map.md), and
[`archive/manifests/campaign_map.csv`](archive/manifests/campaign_map.csv)
gives the same thing one level up: each claim of the paper against the
experiment family, the registration status and the records behind it.

For example, Table 4 of the main text is
`metrics/r68_measured_time_summary_table.tex`, written by
`scripts/rev68_tables.py`, covered by
`metrics/r68_final_experiment_manifest.json`.

## Environments

The campaigns did **not** all run in one environment, and using the wrong one
is the most likely way a reproduction attempt fails for reasons unrelated to
the science.

| Campaign family | Python | Key packages |
|---|---|---|
| Field, timing, budget and long-arc campaigns (R8, R14–R17, and the R18–R68 and J-series campaigns reusing their trajectory machinery) | 3.12.1 | NumPy 2.2.6, SciPy 1.14.1, Numba 0.63.1 |
| Orbit-level verification runs (R11, R12) and the reference side of the external cross-validation | 3.10.20 | NumPy 2.2.6, SciPy 1.15.3 |
| External cross-validation comparator | 3.12.13 (conda-forge) | TudatPy 1.0.0, NumPy 1.26.4 |

Requirement files and the full derivation are in
[`environments/`](environments/README.md); the table there is regenerated from
the archived provenance records by `tools/show_environments.py`.
`requirements-smoke-test.txt` at the root is for the quick check only and is
not the environment of any campaign.

All runs were on Windows 11 x64, an Intel Core i7-9750H at 2.6 GHz. The
numerical kernel is pinned across every campaign to one source snapshot,
release tag `paper-truncation-v1.0`, commit
`27e9ab86ed61d623f78c453ea2054348f1044c23`; that snapshot is what `src/`
contains.

## Layout

| path | contents |
|---|---|
| `src/lunaris/` | the measurement instrument: spherical-harmonic kernel, ephemeris, third-body and solar-radiation models, solid tides, symplectic step, gravity-file loaders |
| `scripts/` | 358 tracked files: 331 Python drivers/analysis passes and 27 shell, PowerShell and cmd launchers |
| `metrics/` | 670 result records, 131 generated LaTeX tables, preregistrations, the claims ledger, and the 46 campaign manifests (`r10`–`r25`, `r28`–`r31`, `r35`, `r37`–`r39`, `r41`–`r42`, `r44`, `r48`, `r50`–`r53`, `r56`–`r68`, `rJ`); per-orbit case configurations in `*_cases/` |
| `archive/` | the reader-facing evidence package: the campaign map, the archive index, the migration record for tables that were moved out of the supplement PDF, and their LaTeX and CSV forms |
| `figures/` | the 34 figure PDFs used in the manuscript |
| `data/` | LP150Q gravity product, the external cross-validation record, and the download target for fetched products |
| `experiments/` | experiment protocol definitions |
| `environments/` | one requirement file per campaign family, and how they were derived |
| `docs/` | the reproducibility index, the claim-to-artifact map, the source manifest and the external-data table |
| `tools/` | generators for the derived files above; see [`tools/README.md`](tools/README.md) |
| `tests/` | standard-library regression tests for archive-maintenance safeguards |

At the root are the three things a reader runs, `verify_snapshot.py`,
`audit_manifest_digests.py` and `fetch_data.py`, plus the audit's input,
`known_stale_digests.json`.

## External data

```bash
python fetch_data.py --list          # what is needed, and where it comes from
python fetch_data.py --group lunar   # the GRAIL fields (365 MB)
python fetch_data.py                 # everything (~430 MB)
```

Each gravity file is verified against the SHA-256 recorded in the campaign
manifest as it lands, so a corrupted or superseded product is caught
immediately instead of turning into a wrong number later. Files that already
verify are skipped, so the command is safe to re-run.

Exact filenames, byte counts, sources and expected paths are in
[`docs/EXTERNAL_DATA.md`](docs/EXTERNAL_DATA.md). Note that for several
products the PDS copy and the Tudat resource copy are not byte-identical, and
only the one listed there verifies against the recorded digest.

## Verifying the archive

```bash
python audit_manifest_digests.py
```

Every driver script is hashed in the manifest of the campaign it belongs to.
This walks all 46 manifests and compares them against the files as they stand.
Digests of everything under `src/` are in
[`docs/SOURCE_MANIFEST.md`](docs/SOURCE_MANIFEST.md).

These checks answer different questions:

| Check | What it establishes | What it does not establish |
|---|---|---|
| `audit_manifest_digests.py --strict` | recorded artifacts still have the bytes sealed by their campaign manifests | that the scientific method or implementation is correct |
| `tools/build_source_manifest.py --check` | the vendored source snapshot has not drifted | that every campaign used an otherwise identical runtime |
| `verify_snapshot.py` | representative harmonic evaluations reproduce the archived snapshot within the declared tolerance | complete trajectory reproduction or physical ground truth |
| JSON parsing and maintenance tests | records are readable and archive safeguards behave as specified | numerical convergence of every historical run |

Together they protect provenance and detect accidental drift. Scientific
support additionally comes from convergence runs, independent SHTOOLS/Tudat
comparisons, registered controls, and the per-claim evidence map; no single
green check substitutes for those layers.

At the current state the audit reports 824 digests matching, no file missing
and no difference of any kind, and `--strict` reports the same. That was not
always so. Thirty-one files once differed from the revisions their manifests
recorded, because the repository had fallen behind the working tree by the
campaigns R24 to R39 while the drivers moved on. Copying the drivers, the
records and the campaign manifests across together closed the gap: a digest and
the file it names describe the same object again. The digests were not
re-hashed to make that true. The resolved entries are kept, with their recorded
and observed values, in
[`known_stale_digests.json`](known_stale_digests.json), and the evidence
gathered while they were open is in
[`docs/DIGEST_STATUS.md`](docs/DIGEST_STATUS.md). The same closure was
repeated when the campaigns R42 to R68 and the J series came across: drivers,
records and manifests moved together, and the strict audit passed on the first
walk afterwards, with no entry added to the stale list.

The audit fails on any mismatch that is not recorded, so new drift is
caught rather than absorbed. `.github/workflows/ci.yml` runs it, the numerical
smoke test, the source-manifest check and a parse of every JSON record; run
those four locally if you want the same assurance without waiting on a hosted
runner. `.github/workflows/archive-guard.yml` separately prevents maintenance
changes from rewriting frozen evidence, even when an artifact and its manifest
would otherwise be changed together.

Because the archive is digest-verified throughout, `.gitattributes` disables
line-ending normalization. Without it every recorded digest would break on
checkout, so please leave it in place.

Drivers are preserved byte-for-byte, including absolute paths from the machine
that produced the results. They are deliberately not tidied, because editing
them would invalidate the digests that make the archive checkable.

## Three levels of reproduction

1. **Audit the deposited evidence.** Run the digest, source-manifest, JSON, and
   numerical smoke checks. This requires no external gravity download.
2. **Regenerate analysis products.** Table and figure builders read the frozen
   records under `metrics/`; these generally do not propagate trajectories.
3. **Re-run a campaign.** Fetch the declared external products, select the
   campaign's recorded Python environment, set `PYTHONPATH=src`, and run the
   preserved driver. Historical absolute paths may need an environment/path
   adapter, but the driver bytes should remain unchanged. Compare a new run
   with the archive using `tools/compare_records.py` instead of overwriting the
   deposited record.

## Reproducing a result

Drivers run with `src/` on the import path.

Bash:

```bash
PYTHONPATH=src python scripts/rev18_span_sweep.py --help
PYTHONPATH=src python scripts/rev19_tables.py
```

PowerShell:

```powershell
$env:PYTHONPATH = "src"
python scripts/rev18_span_sweep.py --help
python scripts/rev19_tables.py
```

Analysis and presentation passes read the archived records in `metrics/` in
place and regenerate tables and figures without propagating anything.

Raw trajectory arrays (about 5.5 GB) are not carried here. They are regenerable
from the drivers, and every file's digest is in the campaign manifests, so an
independent run can be checked against the reported set without shipping it.

## License

MIT, see [`LICENSE`](LICENSE). Citation metadata is in
[`CITATION.cff`](CITATION.cff).
