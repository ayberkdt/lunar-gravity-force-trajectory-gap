# Tools

Generators for the derived files in this archive. None of them touch the
campaign records; they only rebuild documentation from what is already there.

| tool | regenerates | needs the manuscript? |
|---|---|---|
| `build_source_manifest.py` | `docs/SOURCE_MANIFEST.md` | no |
| `build_claim_map.py` | `docs/claim_to_artifact_map.md` | yes, for the captions |
| `show_environments.py` | the table in `environments/README.md` | no |
| `record_known_stale.py` | `known_stale_digests.json` | no |

Run them from the repository root:

```bash
python tools/build_source_manifest.py
python tools/build_claim_map.py --manuscript-root ../codebase
python tools/show_environments.py
```

```powershell
python tools\build_source_manifest.py
python tools\build_claim_map.py --manuscript-root ..\codebase
python tools\show_environments.py
```

`build_source_manifest.py --check` verifies without rewriting and is what CI
runs, so a change under `src/` that is not reflected in the manifest fails the
build.

## Two of these deserve care

`docs/REPRODUCIBILITY_INDEX.csv` has no generator here. It was built from the
compiled manuscript — the `.aux` numbering, the float definitions and the
campaign manifests — and is checked in as a frozen artifact. `build_claim_map.py`
consumes it. Rebuilding it requires the manuscript sources, which are not part
of this archive.

`record_known_stale.py` records every driver whose bytes differ from its
manifest digest, and running it absorbs whatever mismatches exist at that
moment. That is the opposite of what the audit is for, so run it only after
deciding, case by case, that a difference is understood. It prints which
entries are newly absorbed so the decision is visible in the diff.
