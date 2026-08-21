# Archive architecture

The repository has two deliberately different layers.

## 1. Frozen evidence layer

This layer contains the objects used to support the manuscript:

| Path | Role |
|---|---|
| `src/` | commit-pinned Lunaris measurement snapshot |
| `scripts/` | historical campaign drivers and analysis passes |
| `metrics/` | registrations, result records, tables, and campaign manifests |
| `figures/` | rendered figures used by the manuscript |
| `archive/` | reader-facing evidence maps and migrated tables |
| `experiments/` | frozen experiment protocols |
| `environments/` | recorded runtime specifications |
| `data/` | shipped small inputs and targets for verified external products |

Historical files are preserved byte-for-byte. Absolute paths, old formatting,
and duplicated helper functions are therefore provenance, even when they are
not ideal software engineering. Refactoring them in place would make the code
look cleaner while weakening the evidence chain.

New campaigns extend this layer by adding newly named files and a new sealed
manifest. They do not overwrite an earlier campaign.

## 2. Maintenance layer

The maintenance layer makes the archive easier to inspect without changing
the evidence:

| Path | Role |
|---|---|
| `tools/` | non-mutating validators and documentation generators |
| `tests/` | regression tests for maintenance utilities and archive rules |
| `.github/workflows/` | automated verification |
| `docs/` | navigation and explanatory documentation |

Maintenance utilities should be standard-library-only where practical, read
files in place, and write only when their command explicitly says it is a
generator. A validator must never repair or absorb drift automatically.

## Change rule

For ordinary maintenance work:

1. Existing evidence files may not change path or bytes.
2. Additions under append-only evidence directories are allowed.
3. Any change under `src/` is rejected.
4. Documentation, tests, and new read-only tools may evolve normally.
5. A genuine scientific correction uses a separately reviewed evidence
   revision and explains which prior record it supersedes.

`tools/check_frozen_changes.py` enforces the first three rules against a Git
base revision. The existing digest and numerical checks then verify content,
source identity, JSON readability, and kernel behaviour.
