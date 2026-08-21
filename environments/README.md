# Environments

The campaigns in this archive did **not** all run in one environment. Three
distinct interpreters appear in the provenance records the drivers wrote, and
mixing them up would make a reproduction attempt fail for reasons that have
nothing to do with the science.

Every version below is read from the archived provenance records themselves
(`metrics/*.json`), not from a recollection of how the machine was set up.
`../tools/show_environments.py` regenerates this table from those records.

| Campaign family | Python | Key packages | Records |
|---|---|---|---|
| Field, timing, budget and long-arc campaigns (R8, R14–R17, and the R18–R68 and J-series campaigns that reuse their trajectory machinery) | 3.12.1 | NumPy 2.2.6, SciPy 1.14.1, Numba 0.63.1 | 162 |
| Orbit-level verification runs (R11, R12) and the reference side of the external cross-validation | 3.10.20 | NumPy 2.2.6, SciPy 1.15.3 | 3 |
| External cross-validation comparator, and the LRO convergence robustness test | 3.12.13 (conda-forge) | TudatPy 1.0.0, NumPy 1.26.4 | 3 |

All runs were on Windows 11 x64, an Intel Core i7-9750H at 2.6 GHz.

The numerical kernel is pinned across every campaign to one source snapshot,
release tag `paper-truncation-v1.0`, commit
`27e9ab86ed61d623f78c453ea2054348f1044c23`. That snapshot is what `src/`
contains.

## Files

| file | use |
|---|---|
| `field-and-timing-py312.txt` | the bulk of the campaigns; install under Python 3.12.1 |
| `orbit-campaign-py310.txt` | the orbit-level verification runs; install under Python 3.10.20 |
| `tudat-validation-py312.txt` | the independent comparator; install under Python 3.12.13 |

For the quick check only, use `../requirements-smoke-test.txt`, which is
smaller than any of these and is not tied to a campaign.

## A note on what these files can and cannot promise

These pin the packages whose versions the drivers recorded. They are not
lockfiles: transitive dependencies were not captured at run time, so a fresh
install may resolve a different BLAS or LLVM build. That is the usual reason a
long propagation reproduces to a slightly different last bit rather than
exactly, and it is why the trajectory comparisons in the paper are decided
against a measured numerical envelope instead of against exact equality.
