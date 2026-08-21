# Contributing without changing the evidence

This repository is a reproducibility archive. Its historical campaign files
are evidence, not a conventional application source tree. A cleanup that
changes their bytes can invalidate the SHA-256 records used by the paper.

## Safe changes

The following changes are normally safe and welcome:

- new documentation under `docs/`;
- new non-mutating maintenance utilities under `tools/`;
- tests for those utilities under `tests/`;
- CI checks that read the archive without rewriting it; and
- new campaign artifacts with new filenames and their own manifest.

## Frozen and append-only areas

Treat these paths as append-only evidence collections:

- `scripts/`
- `metrics/`
- `figures/`
- `archive/`
- `experiments/`
- `environments/`
- `data/`

Existing files in those directories must not be reformatted, renamed, moved,
or regenerated in place. New campaign files may be added with new names. The
one exception is `environments/README.md`, which is documentation regenerated
by `tools/show_environments.py` rather than sealed evidence.

The measurement snapshot under `src/` is stricter: do not add, edit, remove,
or rename anything there. It is pinned by `docs/SOURCE_MANIFEST.md` and is the
source exercised by `verify_snapshot.py`.

The root files `.gitattributes`, `audit_manifest_digests.py`, `fetch_data.py`,
`known_stale_digests.json`, `requirements-smoke-test.txt`, and
`verify_snapshot.py` are also part of the verification contract. Changes to
them require an explicit evidence revision, not a maintenance cleanup.

## Before committing

Run the checks from the repository root:

```powershell
python -m unittest discover -s tests -v
python tools/check_frozen_changes.py --base HEAD --working-tree
python tools/build_source_manifest.py --check
python audit_manifest_digests.py --strict
python verify_snapshot.py
```

The archive digest and numerical checks can take longer than the maintenance
tests. Do not update `known_stale_digests.json` merely to make a check pass.

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the boundary between the
frozen evidence and the maintenance layer.
