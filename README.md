# lunar-gravity-force-trajectory-gap

[![CI](https://github.com/ayberkdt/lunar-gravity-force-trajectory-gap/actions/workflows/ci.yml/badge.svg)](https://github.com/ayberkdt/lunar-gravity-force-trajectory-gap/actions/workflows/ci.yml)
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
manifests and the drivers needed to audit every reported result. All 436
recorded digests currently match the files they name, under the default audit
and under `--strict`. The only things fetched from outside are public archive
data products, and `fetch_data.py` retrieves and checksums those for you.

The archived deposit is at <https://doi.org/10.5281/zenodo.21824029>.

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
archived snapshot used for the submitted manuscript. If it passes, the
numerical core is intact.

The default compares against a declared relative tolerance of `1e-12`. Exact
reproduction to the last bit depends on the LLVM version Numba compiles
through, the CPU's fused-multiply-add behaviour and the platform's libm, so it
is expected on the archived Windows 11 x64 environment and not guaranteed
elsewhere. To assert it there:

```bash
python verify_snapshot.py --strict-bitwise
```

## Where does a given result come from?

[`docs/REPRODUCIBILITY_INDEX.csv`](docs/REPRODUCIBILITY_INDEX.csv) maps every
manuscript figure and table to the artifact file, the script that writes it,
and the campaign manifest that covers it. The readable version is
[`docs/claim_to_artifact_map.md`](docs/claim_to_artifact_map.md), and
[`archive/manifests/campaign_map.csv`](archive/manifests/campaign_map.csv)
gives the same thing one level up: each claim of the paper against the
experiment family, the registration status and the records behind it.

For example, Table 13 of the main text is `metrics/r19_equal_work_table.tex`,
written by `scripts/rev19_tables.py`, covered by
`metrics/r19_final_experiment_manifest.json`.

## Environments

The campaigns did **not** all run in one environment, and using the wrong one
is the most likely way a reproduction attempt fails for reasons unrelated to
the science.

| Campaign family | Python | Key packages |
|---|---|---|
| Field, timing, budget and long-arc campaigns (R8, R14–R17, and the R18–R25, R28–R39 and R41 campaigns reusing their trajectories) | 3.12.1 | NumPy 2.2.6, SciPy 1.14.1, Numba 0.63.1 |
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
| `scripts/` | 234 campaign drivers, analysis passes, table builders and figure generators |
| `metrics/` | 325 result records, 90 generated LaTeX tables, preregistrations, and the 25 campaign manifests (`r10`–`r25`, `r28`–`r31`, `r35`, `r37`–`r39`, `r41`); per-orbit case configurations in `*_cases/` |
| `archive/` | the reader-facing evidence package: the campaign map, the archive index, the migration record for tables that were moved out of the supplement PDF, and their LaTeX and CSV forms |
| `figures/` | the 31 figure PDFs used in the manuscript |
| `data/` | LP150Q gravity product, the external cross-validation record, and the download target for fetched products |
| `experiments/` | experiment protocol definitions |
| `environments/` | one requirement file per campaign family, and how they were derived |
| `docs/` | the reproducibility index, the claim-to-artifact map, the source manifest and the external-data table |
| `tools/` | generators for the derived files above; see [`tools/README.md`](tools/README.md) |

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
This walks all 25 manifests and compares them against the files as they stand.
Digests of everything under `src/` are in
[`docs/SOURCE_MANIFEST.md`](docs/SOURCE_MANIFEST.md).

At the current state the audit reports 436 digests matching, no file missing
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
[`docs/DIGEST_STATUS.md`](docs/DIGEST_STATUS.md).

The audit fails on any mismatch that is not recorded, so new drift breaks CI.

Because the archive is digest-verified throughout, `.gitattributes` disables
line-ending normalization. Without it every recorded digest would break on
checkout, so please leave it in place.

Drivers are preserved byte-for-byte, including absolute paths from the machine
that produced the results. They are deliberately not tidied, because editing
them would invalidate the digests that make the archive checkable.

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
