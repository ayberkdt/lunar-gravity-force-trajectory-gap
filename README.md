# lunar-gravity-force-trajectory-gap

Paper-specific source snapshot, experiment configurations, manifests, analysis
scripts, and reproducibility records for the lunar-gravity force–trajectory gap
study on the degree-1800 GRAIL JGGRX\_1800F field.

The study asks how a prescribed nominal per-call spherical-harmonic budget
should be spread along an eccentric lunar arc, and measures where a smaller
trajectory-averaged truncation-force defect stops predicting a smaller
trajectory error.

## Relationship to Lunaris

The spherical-harmonic kernel used as the measurement instrument comes from
**Lunaris**, a general lunar orbit propagation and surrogate-gravity toolkit
that is considerably broader than this paper. Rather than pointing at that
project, this repository carries a **verbatim vendored subset** of it: the
transitive closure of exactly the modules the published experiments import —
22 files, about 9.7k lines, roughly 6% of the Lunaris source tree — under
`src/lunaris/`.

Every vendored module is byte-for-byte the file that produced the published
results; SHA-256 digests are listed in [`VENDOR_MANIFEST.md`](VENDOR_MANIFEST.md).
The only two exceptions are `src/lunaris/core/__init__.py` and
`src/lunaris/core/propagation/__init__.py`, which upstream eagerly import the
full propagation stack. They are replaced by trimmed initializers that expose
nothing; both are marked as such in the manifest and in the file header. No
computational module was modified.

`verify_snapshot.py` checks this claim numerically — see below.

## Layout

| path | contents |
|---|---|
| `src/lunaris/` | vendored measurement instrument (SH kernel, ephemeris, third-body, SRP, solid tides, symplectic step, gravity-file loaders) |
| `scripts/` | 155 campaign drivers, analysis passes, table builders and figure generators (`rev2`–`rev23`, `make_figures*`, `check_*`) |
| `metrics/` | 253 result records, generated LaTeX tables, preregistrations, and the fourteen campaign manifests `r10`–`r23`; plus per-orbit case configurations in `*_cases/` |
| `figures/` | 29 generated figure PDFs as used in the manuscript |
| `data/` | Lunar Prospector LP150Q product and the Tudat cross-validation record |
| `experiments/` | experiment protocol definitions |

Campaign manifests (`metrics/r*_final_experiment_manifest.json`) record each
campaign's driver scripts, result records, generated tables and figures, and
per-trajectory sidecar hashes with rolled-up digests of the raw state arrays.

## What is deliberately not here

- **Raw trajectory arrays** (`metrics/*_raw/`, ~2.1 GB). Regenerable from the
  drivers; every file's digest is recorded in the campaign manifests, so an
  independent run can be checked against the published set without shipping it.
- **The primary gravity field.** GRAIL JGGRX\_1800F is a 189 MB public PDS
  product, above GitHub's per-file limit. Fetch it yourself (below).
- **SPICE kernels.** Public NAIF products, fetched separately.
- **The manuscript sources.** This repository is the code and record archive;
  the paper is distributed through its own channel.

## External data you need to fetch

The drivers resolve the gravity file through
`lunaris.common.lunar_data.resolve_lunar_gravity_path`, which searches a data
root. Place the products where that resolver can find them, or set the data
root explicitly.

| product | source | note |
|---|---|---|
| GRAIL JGGRX\_1800F SHADR | NASA PDS Geosciences Node | primary field; SHA-256 begins `d2a55206`, ends `6738fec` |
| GSFC GGGRX\_1200L, GRGM1200A | NASA PDS | solution-transfer tests |
| Earth GOCO05c / EGM96, Mars JGMRO120D, Mercury JGMESS\_160A, Venus SHGJ180U | respective archives | cross-body transfer tests |
| `de440s.bsp`, `gm_de440.tpc`, `moon_de440_250416.tf`, `moon_pa_de440_200625.bpc`, `naif0012.tls` | NAIF | ephemeris and lunar orientation |

Lunar Prospector LP150Q (`data/jgl150q1.sha`) is small enough to ship and is
included, so the verification below needs no download.

## Quick check

```bash
pip install -r requirements.txt
python verify_snapshot.py
```

This loads the bundled LP150Q product through the vendored loader, evaluates
nine spherical-harmonic accelerations across three degrees and three points,
and asserts they are bitwise identical to values recorded from the Lunaris
working tree that produced the published results. It requires no external
download and no GPU.

## Reproducing a result

Campaign drivers are grouped by revision prefix and are run from `scripts/`
with `src/` on the import path:

```bash
PYTHONPATH=src python scripts/rev18_span_sweep.py --help
```

Analysis and presentation passes read the archived records in `metrics/` in
place and regenerate the tables and figures without propagating anything:

```bash
PYTHONPATH=src python scripts/rev19_tables.py
```

Note that `scripts/make_figures.py` is a superseded generator kept for
provenance; the current static figures come from `make_figures_r1.py`.

## Environment

Python 3.12.1 on Windows 11 x64; pinned package versions in
`requirements.txt`. Numba compiles the spherical-harmonic kernel on first call,
so the first evaluation in a process is slow.

## License

MIT — see [`LICENSE`](LICENSE). The vendored Lunaris subset is covered by the
same license and copyright.

## Citation

Ayberk Demirkanat, Department of Astronautical Engineering, Istanbul Technical
University. Citation details will be added once the associated paper is
published.
