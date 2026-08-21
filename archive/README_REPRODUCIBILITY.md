# Reproducing the results

This archive is the complete machine-readable record behind the manuscript.
The Supplementary Information carries the definitions, the contracts and the
aggregate results; everything per-orbit, per-case and per-file is here.

Cite it as <https://doi.org/10.5281/zenodo.21824029>. Paths quoted in the manuscript and the
supplement are relative to the root of this deposit.

## What to read first

| you want to | read |
|---|---|
| see which record carries which claim | `manifests/campaign_map.csv` |
| find a directory | `supplement_archive/archive_index.md` |
| find a table that left the supplement PDF | `supplement_archive/migration_map.csv` |
| check that nothing has been altered | `manifests/` and the procedure below |

## Verifying integrity

Every campaign owns one manifest, `metrics/rNN_final_experiment_manifest.json`,
which records a SHA-256 for each file it claims and a seal over itself. The
check that verifies the manifests against the files is
`python_codes/check_manifest_integrity.py`. It recomputes every recorded digest,
recomputes each manifest seal, confirms that no trajectory record is claimed by
two manifests, and confirms that every table and figure the documents pull in is
hashed somewhere. It exits non-zero on any failure.

    python python_codes/check_manifest_integrity.py

The trajectory partition is a property of the archive, not an assertion: a
record appears under exactly one manifest.

## Regenerating a table that left the supplement PDF

Twelve per-orbit and per-case tables were moved out of the supplement because no
text pointed at them and they cannot be read on a page. Each is archived three
ways: the rendered LaTeX the supplement compiled, a CSV rendering, and the
source record plus the generating script, both named and hashed in
`supplement_archive/migration_map.csv`.

To regenerate one, run its generating script from `python_codes/`. The scripts
read frozen records and propagate nothing.

The CSV is a convenience rendering of the typeset table. Where a column heading
was set in mathematics, the CSV carries a plain-text approximation of it; the
archived LaTeX and the source record are the authoritative forms.

## Re-running a campaign

Campaign drivers live in `python_codes/` and are named `revNN_*.py`. Each reads
its registration from `metrics/rNN_preregistration.json` where one exists,
refuses to run if the registration is missing, and writes only its own records.
Two environments are used and both ship as a specification and a full lock:
`environment-py310.yml` with `requirements-py310.lock` for R10-R13, and
`environment-py312.yml` with `requirements-py312.lock` for everything later.

Propagation campaigns are expensive. The gravity-model files are not
redistributed here; obtain each from its archive and verify the digest recorded
in the R16 manifest and in Supplementary Section S19 before use.

## What is deliberately not here

Raw trajectory arrays under `metrics/*_raw/` (about 14 GB) are regenerable from
the drivers and are excluded; their per-file digests are recorded in the
campaign manifests. The J1-J3 arrays (about 0.7 GB) are held outside `metrics/`
at the location recorded in `rJ_final_experiment_manifest.json` and are excluded
on the same terms, which together are the about 15 GB the manuscript reports as
not redistributed. Large external coefficient products and SPICE kernels are not
redistributed.
