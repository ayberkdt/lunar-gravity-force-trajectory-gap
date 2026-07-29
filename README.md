# lunar-gravity-force-trajectory-gap

Source, experiment configurations, manifests, analysis scripts and
reproducibility records for a study of the force–trajectory gap in high-degree
lunar gravity modelling, on the degree-1800 GRAIL JGGRX\_1800F field.

The study asks how a prescribed nominal per-call spherical-harmonic budget
should be spread along an eccentric lunar arc, and measures where a smaller
trajectory-averaged truncation-force defect stops predicting a smaller
trajectory error.

This repository is self-contained: the measurement instrument, the drivers that
produced every published number, and the records they wrote are all here. The
only things fetched from outside are public archive data products, and
`fetch_data.py` retrieves and checksums those for you.

## Getting started

```bash
pip install -r requirements.txt
python verify_snapshot.py
```

`verify_snapshot.py` needs no download. It loads the Lunar Prospector LP150Q
product shipped in `data/`, evaluates nine spherical-harmonic accelerations
across three degrees and three field points, and asserts they are bitwise
identical to values recorded from the source tree that produced the published
results. If that passes, the numerical core is intact.

To run the experiments themselves you need the archive data products:

```bash
python fetch_data.py --list          # what is needed, and where it comes from
python fetch_data.py --group lunar   # the GRAIL fields (365 MB)
python fetch_data.py                 # everything (~430 MB)
```

Each gravity file is verified against the SHA-256 recorded in the campaign
manifest as it lands, so a corrupted or superseded product is caught
immediately instead of turning into a wrong number later. Files that already
verify are skipped, so the command is safe to re-run.

## Layout

| path | contents |
|---|---|
| `src/lunaris/` | the measurement instrument: spherical-harmonic kernel, ephemeris, third-body and solar-radiation models, solid tides, symplectic step, gravity-file loaders |
| `scripts/` | 136 campaign drivers, analysis passes, table builders and figure generators |
| `metrics/` | 245 result records, generated LaTeX tables, preregistrations, and the fourteen campaign manifests `r10`–`r23`; per-orbit case configurations in `*_cases/` |
| `figures/` | the 29 figure PDFs used in the paper |
| `data/` | LP150Q gravity product, the Tudat cross-validation record, and the download target for fetched products |
| `experiments/` | experiment protocol definitions |

Campaign manifests (`metrics/r*_final_experiment_manifest.json`) record each
campaign's driver scripts, result records, generated tables and figures, its
input products with digests, and per-trajectory sidecar hashes with rolled-up
digests of the raw state arrays.

## Verifying the archive

```bash
python audit_manifest_digests.py
```

Every driver script is hashed in the manifest of the campaign it belongs to.
This walks all fourteen manifests and reports any script whose current bytes
differ from what was recorded — the signal that a driver was edited after its
manifest was finalized and that the manifest needs refreshing with the
campaign's own `revNN_finalize_manifest.py`.

Digests of everything under `src/` are listed in
[`SOURCE_MANIFEST.md`](SOURCE_MANIFEST.md).

Because the archive is digest-verified throughout, `.gitattributes` disables
line-ending normalization. Without it every recorded digest would break on
checkout, so please leave it in place.

Driver scripts are kept exactly as they ran, including absolute paths from the
machine that produced the results. They are deliberately not tidied: editing
them would invalidate the digests that make the archive checkable.

## Reproducing a result

Drivers are grouped by revision prefix and run with `src/` on the import path:

```bash
PYTHONPATH=src python scripts/rev18_span_sweep.py --help
```

Analysis and presentation passes read the archived records in `metrics/` in
place and regenerate tables and figures without propagating anything:

```bash
PYTHONPATH=src python scripts/rev19_tables.py
PYTHONPATH=src python scripts/make_figures_r1.py
```

Raw trajectory arrays (~2.1 GB) are not carried here. They are regenerable from
the drivers, and every file's digest is in the campaign manifests, so an
independent run can be checked against the published set without shipping it.

## Environment

Python 3.12.1 on Windows 11 x64; pinned versions in `requirements.txt`. Numba
compiles the spherical-harmonic kernel on first call, so the first evaluation
in a process is slow.

## Data products

| product | body | role | source |
|---|---|---|---|
| JGGRX\_1800F | Moon | primary field | PDS Geosciences Node |
| GRGM1200A, GGGRX\_1200L | Moon | solution transfer | PDS Geosciences Node |
| LP150Q | Moon | shipped in `data/` | PDS Geosciences Node |
| GOCO05c, EGM96 | Earth | cross-body transfer | Tudat resource archive |
| JGMRO120D | Mars | cross-body transfer | Tudat resource archive |
| JGMESS\_160A | Mercury | cross-body transfer | Tudat resource archive |
| SHGJ180U | Venus | cross-body transfer | Tudat resource archive |
| DE440 kernels | — | ephemeris and lunar orientation | NAIF |

## License

MIT — see [`LICENSE`](LICENSE).

## Citation

Ayberk Demirkanat, Department of Astronautical Engineering, Istanbul Technical
University. Citation details will be added once the associated paper is
published.
