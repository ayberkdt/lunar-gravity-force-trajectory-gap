# lunar-gravity-force-trajectory-gap

[![CI](https://github.com/ayberkdt/lunar-gravity-force-trajectory-gap/actions/workflows/ci.yml/badge.svg)](https://github.com/ayberkdt/lunar-gravity-force-trajectory-gap/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Reproducibility archive for fixed-budget spherical-harmonic degree allocation
and the lunar gravity force–trajectory gap, on the degree-1800 GRAIL
JGGRX\_1800F field.

The study asks how a prescribed nominal per-call spherical-harmonic budget
should be spread along an eccentric lunar arc, and measures where a smaller
trajectory-averaged truncation-force defect stops predicting a smaller
trajectory error.

The archive documents the measurement instrument, the campaign records, the
manifests and the available drivers needed to audit the reported results.
Thirty-one files currently differ from the revisions their campaign manifests
record — 16 drivers, 13 tables and 2 result records. Nineteen of those have
been shown to reproduce the archived output byte-for-byte; the rest are listed
explicitly in [`docs/DIGEST_STATUS.md`](docs/DIGEST_STATUS.md) and must be
resolved before an archival release.

The only things fetched from outside are public archive data products, and
`fetch_data.py` retrieves and checksums those for you.

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

[`docs/REPRODUCIBILITY_INDEX.csv`](docs/REPRODUCIBILITY_INDEX.csv) maps every manuscript
figure and table to the artifact file, the script that writes it, and the
campaign manifest that covers it. The readable version is
[`docs/claim_to_artifact_map.md`](docs/claim_to_artifact_map.md).

For example, Table 13 of the main text is `metrics/r19_equal_work_table.tex`,
written by `scripts/rev19_tables.py`, covered by
`metrics/r19_final_experiment_manifest.json`.

## Environments

The campaigns did **not** all run in one environment, and using the wrong one
is the most likely way a reproduction attempt fails for reasons unrelated to
the science.

| Campaign family | Python | Key packages |
|---|---|---|
| Field, timing, budget and long-arc campaigns (R8, R14–R17, and the R18–R23 campaigns reusing their trajectories) | 3.12.1 | NumPy 2.2.6, SciPy 1.14.1, Numba 0.63.1 |
| Orbit-level verification runs (R11, R12) and the reference side of the external cross-validation | 3.10.20 | NumPy 2.2.6, SciPy 1.15.3 |
| External cross-validation comparator | 3.12.13 (conda-forge) | TudatPy 1.0.0, NumPy 1.26.4 |

Requirement files and the full derivation are in
[`environments/`](environments/README.md); the table there is regenerated from
the archived provenance records by `tools/show_environments.py`.
`requirements-smoke-test.txt` at the root is for the quick check only and is
not the environment of any campaign.

All runs were on Windows 11 x64, an Intel Core i7-9750H at 2.6 GHz. The
numerical kernel is pinned across R10–R23 to one source snapshot, release tag
`paper-truncation-v1.0`, commit `27e9ab86ed61d623f78c453ea2054348f1044c23`;
that snapshot is what `src/` contains.

## Layout

| path | contents |
|---|---|
| `src/lunaris/` | the measurement instrument: spherical-harmonic kernel, ephemeris, third-body and solar-radiation models, solid tides, symplectic step, gravity-file loaders |
| `scripts/` | 136 campaign drivers, analysis passes, table builders and figure generators |
| `metrics/` | 245 result records, generated LaTeX tables, preregistrations, and the fourteen campaign manifests `r10`–`r23`; per-orbit case configurations in `*_cases/` |
| `figures/` | the 29 figure PDFs used in the manuscript |
| `data/` | LP150Q gravity product, the external cross-validation record, and the download target for fetched products |
| `experiments/` | experiment protocol definitions |
| `environments/` | one requirement file per campaign family, and how they were derived |
| `docs/` | the reproducibility index, the claim-to-artifact map, the source manifest and the external-data table |
| `tools/` | generators for the derived files above; see [`tools/README.md`](tools/README.md) |

At the root are the three things a reader runs — `verify_snapshot.py`,
`audit_manifest_digests.py` and `fetch_data.py` — plus the audit's input,
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
[`docs/EXTERNAL_DATA.md`](docs/EXTERNAL_DATA.md). Note that for several products the PDS
copy and the Tudat resource copy are not byte-identical, and only the one
listed there verifies against the recorded digest.

## Verifying the archive

```bash
python audit_manifest_digests.py
```

Every driver script is hashed in the manifest of the campaign it belongs to.
This walks all fourteen manifests and compares them against the scripts as they
stand. Digests of everything under `src/` are in
[`docs/SOURCE_MANIFEST.md`](docs/SOURCE_MANIFEST.md).

**Known open issue.** 143 digests match; **31 do not** — 16 driver scripts, 13
generated tables and 2 result records, across R11–R19. The originals are
unrecoverable: no copy matching any recorded digest survives anywhere on the
authoring machine, across 53,720 Python files, the editor's local history and
every available git history.

Re-running settles 19 of the 31. Four presentation passes and the full R16
cross-body calibration driver regenerate their committed artifacts
byte-for-byte, so those edits provably did not change the numbers; the R12
table differs only by caption line-wrapping. Eleven drivers cannot be checked
that way, because re-running them means re-running the propagation campaign.
The evidence, the search, and what is left open are set out in
[`docs/DIGEST_STATUS.md`](docs/DIGEST_STATUS.md); the entries themselves are in
[`known_stale_digests.json`](known_stale_digests.json).

The digests are deliberately **not** refreshed: re-hashing would record that the
current file is the one the manifest describes, which is a different claim from
the evidence above and is unsupported for the eleven. They are resolved by
re-running the affected campaign with the current driver and refreshing the
manifest from that run. The audit fails on any mismatch that is *not* recorded,
so new drift breaks CI; `--strict` fails on the known ones too and is the state
to reach before an archival release.

Because the archive is digest-verified throughout, `.gitattributes` disables
line-ending normalization. Without it every recorded digest would break on
checkout, so please leave it in place.

Drivers whose digests match their manifests are preserved byte-for-byte,
including absolute paths from the machine that produced the results; they are
deliberately not tidied, because editing them would invalidate the digests that
make the archive checkable. The known exceptions are documented in
[`known_stale_digests.json`](known_stale_digests.json).

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

Raw trajectory arrays (~2.1 GB) are not carried here. They are regenerable from
the drivers, and every file's digest is in the campaign manifests, so an
independent run can be checked against the reported set without shipping it.

## License

MIT — see [`LICENSE`](LICENSE). Citation metadata is in
[`CITATION.cff`](CITATION.cff).
