# Digest status

The campaign manifests record a SHA-256 for each driver script and for the
tables, figures and result records those drivers wrote. Most of them still
match. This page states exactly which do not, what was done to find out
whether it matters, and what remains open.

Run `python audit_manifest_digests.py` for the live version of these counts.

## Where the archive stands

| | count |
|---|---|
| digests matching their manifest | 143 |
| differences, reproduction verified | 19 |
| differences, not verified | 12 |
| files missing | 0 |

The 31 differences span 16 driver scripts, 13 generated tables and 2 result
records, across campaigns R11–R19. They are enumerated with recorded and
observed values in [`../known_stale_digests.json`](../known_stale_digests.json).

## The originals could not be recovered

The first thing to try is restoring the exact file each manifest describes. A
digest-matching copy would close the question outright. The following were
searched, hashing every candidate against the 16 recorded script digests:

- the working tree and every sibling project directory on both drives,
  53,720 Python files in total;
- the editor's local-history store, 1,932 revisions — none of these drivers
  was ever tracked there;
- the git histories available on the machine — none has ever contained a file
  matching these driver names.

**No copy matching any recorded digest exists.** Only the current revision of
each file survives, so restoration is not available and the question has to be
settled by reproduction instead.

## What reproduction shows

Where a driver could be re-run against the archived inputs, the test is
whether it regenerates the committed artifact byte-for-byte. That does not
prove the current file is the one the manifest names, but it does establish
the thing that matters scientifically: the edit did not change the numbers.

| driver | campaign | result |
|---|---|---|
| `rev14_tables.py` | R14 | regenerates its tables byte-for-byte |
| `rev15_tables.py` | R15 | regenerates its tables byte-for-byte |
| `rev19_tables.py` | R19 | regenerates its table byte-for-byte |
| `rev16_multibody_calibration.py` | R16 | regenerates `r16_multibody_calibration.json` byte-for-byte |
| `rev12_atallah_tables.py` | R12 | every value reproduces; the committed table differs only by caption line-wrapping, so it was hand-rewrapped after generation |

R16 is a full campaign driver rather than a presentation pass, so its
reproduction covers the cross-body calibration end to end.

The reading for these five is therefore narrow and firm: the manifest digests
are stale relative to a later, internally coherent re-run, not evidence that
the reported numbers came from something unavailable.

## What remains open

Eleven drivers could not be checked this way, because re-running them means
re-running the propagation campaign itself, which needs the raw trajectory
trees that are not carried in this archive:

`rev11_full_convergence.py`, `rev11_geometric_verification.py`,
`rev11_phase_sweep.py`, `rev13_variational_check.py`, `rev13_timing_match.py`,
`rev14_timing_budget.py`, `rev14_variational_budget.py`,
`rev15_fixed_oracle.py`, `rev17_longarc60.py`, `rev18_span_sweep.py`,
`rev19_equal_total_work.py`.

`rev14_timing_budget.py` is a special case: it measures wall-clock cost, so
re-running it cannot reproduce an earlier record byte-for-byte on any machine.

Each is resolved by re-running its campaign with the current driver and
refreshing the manifest from that run, which is the only route that leaves the
chain closed rather than merely re-hashed.

## Why the digests are not simply refreshed

Re-hashing would record that the current file is what the manifest describes.
That is a different claim from the one the evidence supports, and for the
eleven unverified drivers there is no evidence at all. `record_known_stale.py`
exists to document differences, not to absorb them, and it prints whatever it
newly takes in so the decision shows up in the diff.
