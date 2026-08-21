"""Gate the evidence package: do the campaign manifests still describe the repo?

The manifests are written once per campaign and then go stale silently, because
nothing re-reads them. Two failures have actually happened here: a digest that
no longer matches the file it names (any later edit to a script or a generated
table does this), and a record written as ``{"missing": true}`` because the
generator guessed a filename that was never produced under that name.

A third failure is coverage: a table or figure the manuscript \\input{}s or
\\includegraphics{}es but no manifest ever hashed, so the evidence package
would ship the claim without its provenance.

A fourth is the partition. The supplement states that a trajectory appears
under exactly one manifest, which is what lets the fourteen be read as a
partition rather than as overlapping inventories. A campaign that reuses an
earlier campaign's driver writes its sidecars into the earlier campaign's
directory tree, and a generator that sweeps that tree blindly then claims them
twice. Only trajectory records are held to this; a script or a figure may
legitimately be hashed under more than one manifest.

This checks all four, and then asks the archive-side questions the forward
scan cannot: reverse ownership (every non-exempt artifact in metrics/ is
owned by some manifest), twin agreement (a same-name copy under archive/
hashes like its metrics/ original), and derivative freshness (a tracked
derivative is newer than every input it is produced from). Exits non-zero
on any failure, so it can gate a submission build. It writes nothing.

Usage:  python check_manifest_integrity.py [--quiet] [--skip-mtime]

Usage:  python check_manifest_integrity.py [--quiet]
"""

from __future__ import annotations

import glob
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
FIGURES = ROOT / "figures"
CH = ROOT / "chapters"

HEX = re.compile(r"^[0-9a-f]{64}$")

# Digest fields that seal a manifest or name an external object rather than a
# file in this repository. Hashing them against a path would be meaningless.
NON_FILE_KEYS = {"manifest_sha256", "protocol_sha256", "raw_rollup_sha256",
                 "lunaris_commit", "sha256"}

# Where a section's bare filenames live. Sections not listed fall back to the
# chain below, which is why an unrecognised section still resolves rather than
# being reported as absent.
SECTION_BASE = {
    "scripts": CODE,
    "generated_figures": FIGURES,
}
FALLBACK_BASES = (METRICS, CODE, FIGURES, ROOT)

FILE_SUFFIX = re.compile(
    r"\.(tex|py|json|npz|npy|pdf|csv|bib|md|txt|png|ps1|log)$")

# A trajectory record: the per-case sidecar and the raw state array. These are
# the records the supplement's partition claim is about. The replication
# campaigns' sidecars live under rJ-prefixed trees and are held to the same
# partition.
TRAJECTORY = re.compile(r"^metrics/r(?:\d+|J\d*)_(cases|raw)/")

_digests: dict[Path, str] = {}


def sha(path: Path) -> str:
    if path not in _digests:
        d = hashlib.sha256()
        with path.open("rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                d.update(blk)
        _digests[path] = d.hexdigest()
    return _digests[path]


def rel(path: Path) -> str:
    """Repo-relative when possible. Some manifests hash an input that lives
    outside the repository (a gravity coefficient product), and those are
    verified in place rather than skipped."""
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def resolve(name: str, section: str) -> Path | None:
    n = name.replace("\\", "/")
    if "/" in n:
        # repo-relative (R10 entries), metrics-relative (sidecar keys), or a
        # path reaching outside the repo (external coefficient products)
        for cand in (ROOT / n, METRICS / n, Path(n)):
            if cand.is_file():
                return cand
        return None
    bases = ([SECTION_BASE[section]] if section in SECTION_BASE else [])
    for base in bases + list(FALLBACK_BASES):
        cand = base / n
        if cand.is_file():
            return cand
    return None


def collect(node, section: str, out: list) -> None:
    """Pull (name, digest, flagged_missing, section) out of any manifest shape.

    Three record shapes are in use across the fourteen manifests:
        "name.ext": "<sha>"
        "name.ext": {"sha256": "<sha>", "bytes": n}
        {"path": "repo/relative", "sha256": "<sha>", ...}
    """
    if isinstance(node, dict):
        digest, name = node.get("sha256"), (
            node.get("path") or node.get("file") or node.get("relpath"))
        if isinstance(digest, str) and HEX.match(digest) and isinstance(name, str):
            out.append((name, digest, bool(node.get("missing")), section))
        for key, val in node.items():
            sub = key if key not in ("entries",) else section
            if key in NON_FILE_KEYS and not isinstance(val, dict):
                continue
            if isinstance(val, str):
                if HEX.match(val) and FILE_SUFFIX.search(key):
                    out.append((key, val, False, section))
            elif isinstance(val, dict):
                inner = val.get("sha256")
                if isinstance(inner, str) and HEX.match(inner) and \
                        FILE_SUFFIX.search(key):
                    out.append((key, inner, bool(val.get("missing")), section))
                elif val.get("missing") is True and FILE_SUFFIX.search(key):
                    out.append((key, None, True, section))
                collect(val, sub, out)
            elif isinstance(val, list):
                collect(val, sub, out)
    elif isinstance(node, list):
        for val in node:
            collect(val, section, out)


def seal_ok(payload: dict) -> bool | None:
    """Recompute the self-seal. None when the manifest carries no seal field."""
    sealed = payload.get("manifest_sha256")
    if not sealed:
        return None
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    return hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest() == sealed


def manuscript_inventory() -> dict[str, set[str]]:
    """What the manuscript actually pulls in, as repo-relative paths.

    preamble.tex is scanned too: it \\input{}s generated files (the campaign
    number macros), and leaving it out once let every macro the abstract
    quotes sit outside the audit without even an 'uncovered' report.
    """
    tex = sorted(glob.glob(str(CH / "*.tex"))) + [
        str(ROOT / "main.tex"), str(ROOT / "supplement.tex"),
        str(ROOT / "preamble.tex")]
    text = "\n".join(Path(p).read_text(encoding="utf-8") for p in tex)

    inputs = set()
    for target in re.findall(r"\\input\{([^}]+)\}", text):
        t = target.strip()
        for cand in (t, t + ".tex"):
            if (ROOT / cand).is_file():
                inputs.add(cand)
                break

    figures = set()
    for target in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}",
                             text):
        t = target.strip()
        for cand in (t, t + ".pdf", t + ".png",
                     f"figures/{t}", f"figures/{t}.pdf", f"figures/{t}.png"):
            if (ROOT / cand).is_file():
                figures.add(cand)
                break

    sources = {p for p in ("main.tex", "supplement.tex", "preamble.tex",
                           "references.bib") if (ROOT / p).is_file()}
    sources |= {f"chapters/{Path(p).name}" for p in glob.glob(str(CH / "*.tex"))}
    return {"input": inputs, "figure": figures, "source": sources}


def main() -> int:
    quiet = "--quiet" in sys.argv
    problems: list[str] = []
    covered: set[str] = set()
    external: list[tuple[str, str, str]] = []
    claimed_by: dict[str, list[str]] = {}
    totals = {"match": 0, "mismatch": 0, "absent": 0, "flagged": 0}

    paths = sorted(glob.glob(
        str(METRICS / "r[0-9]*_final_experiment_manifest.json"))
        + glob.glob(str(METRICS / "rJ_final_experiment_manifest.json")))
    if not paths:
        print("no campaign manifests found", file=sys.stderr)
        return 1

    for mp in paths:
        name = Path(mp).name
        payload = json.loads(Path(mp).read_text(encoding="utf-8"))

        seal = seal_ok(payload)
        if seal is False:
            problems.append(f"{name}: manifest_sha256 does not recompute "
                            f"(hand-edited after it was sealed)")

        records, seen = [], set()
        collect(payload, "", records)
        per = {"match": 0, "mismatch": 0, "absent": 0, "flagged": 0}
        detail: list[str] = []

        # Inputs distributed by an external archive are not in the package, so
        # their digest is the whole of the provenance and there is nothing here
        # to compare it against. They are counted and reported, never failed on
        # absence -- failing would make the gate pass only on the one machine
        # that happens to hold the downloads.
        for key, rec in (payload.get("input_products") or {}).items():
            if isinstance(rec, dict) and rec.get("shipped_in_package") is False:
                external.append((name, key, rec.get("sha256", "")))

        for rec_name, digest, flagged, section in records:
            if section == "input_products":
                continue
            key = (rec_name, digest)
            if key in seen:
                continue
            seen.add(key)
            if flagged or digest is None:
                per["flagged"] += 1
                detail.append(f"    flagged-missing  {rec_name}")
                problems.append(f"{name}: '{rec_name}' recorded as missing")
                continue
            found = resolve(rec_name, section)
            if found is None:
                per["absent"] += 1
                detail.append(f"    not-on-disk      {rec_name}")
                problems.append(f"{name}: '{rec_name}' is not on disk")
                continue
            covered.add(rel(found))
            if TRAJECTORY.match(rel(found)):
                owners = claimed_by.setdefault(rel(found), [])
                if name not in owners:
                    owners.append(name)
            if sha(found) == digest:
                per["match"] += 1
            else:
                per["mismatch"] += 1
                detail.append(f"    stale-digest     {rel(found)}")
                problems.append(f"{name}: digest for '{rel(found)}' is stale")

        for k in totals:
            totals[k] += per[k]
        if not quiet:
            seal_txt = {True: "sealed", False: "SEAL BROKEN",
                        None: "unsealed"}[seal]
            print(f"{name}: {per['match']} ok, {per['mismatch']} stale, "
                  f"{per['absent']} absent, {per['flagged']} flagged "
                  f"({len(seen)} records, {seal_txt})")
            for line in detail[:20]:
                print(line)
            if len(detail) > 20:
                print(f"    ... +{len(detail) - 20} more")

    if external:
        bad = [f"{m}:{k}" for m, k, d in external if not HEX.match(d or "")]
        if bad:
            problems.append("external inputs recorded without a usable digest: "
                            + ", ".join(bad))
        if not quiet:
            print(f"\n[external] {len(external)} distributed inputs recorded by "
                  f"digest, not shipped in the package "
                  f"({len(external) - len(bad)} carry a valid SHA-256)")

    shared = {p: owners for p, owners in claimed_by.items() if len(owners) > 1}
    if not quiet:
        print()
    if shared:
        pairs: dict[tuple, int] = {}
        for owners in shared.values():
            pairs[tuple(sorted(owners))] = pairs.get(tuple(sorted(owners)), 0) + 1
        for owners, n in sorted(pairs.items(), key=lambda kv: -kv[1]):
            problems.append(
                f"{n} trajectory records claimed by more than one manifest: "
                + ", ".join(o.replace("_final_experiment_manifest.json", "")
                            for o in owners))
        if not quiet:
            print(f"[overlap] {len(shared)} trajectory records under two or "
                  f"more manifests:")
            for p in sorted(shared)[:10]:
                print(f"    {p}  <- {', '.join(shared[p])}")
            if len(shared) > 10:
                print(f"    ... +{len(shared) - 10} more")
    elif not quiet:
        print(f"[ok] trajectory records form a partition "
              f"({len(claimed_by)} records, no record under two manifests)")

    # Reverse ownership: every result record in metrics/ must be indexed by
    # some manifest. The forward scan alone (manuscript -> covered) let two
    # full campaign records go unsealed while feeding a printed table,
    # because nothing asked the archive-side question. Legacy-package records
    # are owned by the legacy manifest, which is parsed here for its names
    # rather than trusted as a campaign manifest.
    legacy_owned: set[str] = set()
    legacy = METRICS / "r3_experiment_manifest.json"
    if legacy.exists():
        legacy_owned = {Path(m).name for m in re.findall(
            r"[\w./-]+\.(?:json|npz|npy|tex|csv)",
            legacy.read_text(encoding="utf-8"))}
    orphan_exempt = {
        # scheduling/progress state, not evidence
        "rJ_supervisor_progress.json",
        # The claims ledger is an instrument, not a result: every value in it
        # is derived from records that are themselves sealed. Owning it in a
        # manifest would make the two chase each other, since re-sealing a
        # record would stale the ledger and re-pinning the ledger would stale
        # the manifest that owned it. claims_ledger.py checks it instead.
        "claims_ledger.json",
        "claims_ledger_macros.tex",
        # The re-seal journal is this checker's counterpart, and owning it has
        # the same circularity: it records the digests a manifest was given,
        # so sealing it inside one would mean every entry it writes stales the
        # manifest that holds it. It is evidence about the instruments rather
        # than about the Moon, and it is readable in the package.
        "reseal_log.json",
    }
    # classes the R10 manifest policy excludes by name: smoke tests,
    # active/checkpoint files, and campaign progress state
    exempt_suffix = re.compile(r"_(?:smoke|active|progress)\.(?:json|tex|csv)$")
    orphans = []
    # scope is defined by exclusion, not by a name-prefix allowlist: an
    # earlier r*/rJ* prefix filter silently passed fifteen records (the
    # e-series and the robustness_* qualification records), one of which was
    # the sole source of a printed paragraph. It is likewise defined over
    # every artifact class the package ships (.json records, generated .tex
    # tables, .csv lookups) rather than an extension allowlist: a .json-only
    # sweep let two leftover .tex tables sit unowned. Everything in metrics/
    # that is not a manifest and not explicitly exempt must be owned.
    for p in sorted(q for pat in ("*.json", "*.tex", "*.csv")
                    for q in METRICS.glob(pat)):
        n = p.name
        if n.endswith("manifest.json"):
            continue
        if n in orphan_exempt or n in legacy_owned or exempt_suffix.search(n):
            continue
        if rel(p) not in covered and n not in covered:
            orphans.append(n)
    if orphans:
        problems.append("result records owned by no manifest: "
                        + ", ".join(orphans))
        if not quiet:
            print(f"[orphan] {len(orphans)} result records under no manifest:")
            for n in orphans[:20]:
                print(f"    {n}")
            if len(orphans) > 20:
                print(f"    ... +{len(orphans) - 20} more")
    elif not quiet:
        print("[ok] every result record in metrics/ is owned by a manifest")

    # Same-basename twin copies. The curated archive/ tree mirrors metrics/
    # artifacts into a structured layout, and every stale-derivative incident
    # to date has had the same signature: a fresh file in metrics/ beside an
    # old same-name copy that kept circulating. metrics/ is canonical; a twin
    # under archive/ that hashes differently is a blocker until synced.
    archive_root = ROOT / "archive"
    twin_bad: list[str] = []
    if archive_root.exists():
        met_files = {p.name: p for p in METRICS.iterdir()
                     if p.is_file() and p.name not in orphan_exempt
                     and not exempt_suffix.search(p.name)}
        for p in sorted(archive_root.rglob("*")):
            if not p.is_file() or exempt_suffix.search(p.name):
                continue
            twin = met_files.get(p.name)
            if twin is not None and sha(p) != sha(twin):
                twin_bad.append(str(p.relative_to(ROOT)).replace("\\", "/"))
    if twin_bad:
        problems.append(f"{len(twin_bad)} archive/ twin copies disagree "
                        "with their metrics/ original")
        if not quiet:
            print(f"[twin] {len(twin_bad)} stale archive/ copies "
                  f"(metrics/ is canonical):")
            for n in twin_bad[:20]:
                print(f"    {n}")
            if len(twin_bad) > 20:
                print(f"    ... +{len(twin_bad) - 20} more")
    elif not quiet and archive_root.exists():
        print("[ok] every archive/ twin matches its metrics/ original")

    # The test above compares the twins that exist; it cannot see one that was
    # never made. R65 shipped for a round with no archive copy of either its
    # manifest or its registration, and the check printed [ok] throughout,
    # because a missing file has nothing to disagree with. Every sealed
    # manifest and every registration must have a twin.
    twin_missing: list[str] = []
    if archive_root.exists():
        arch_names = {p.name for p in archive_root.rglob("*") if p.is_file()}
        for p in sorted(METRICS.glob("r*_final_experiment_manifest.json")):
            if p.name not in arch_names:
                twin_missing.append(f"archive/manifests/{p.name}")
        for p in sorted(METRICS.glob("r*_preregistration.json")):
            if p.name not in arch_names:
                twin_missing.append(f"archive/registrations/original/{p.name}")
    if twin_missing:
        problems.append(f"{len(twin_missing)} sealed records have no archive "
                        "twin")
        if not quiet:
            print(f"[twin] {len(twin_missing)} sealed records not twinned:")
            for n in twin_missing[:20]:
                print(f"    {n}")
    elif not quiet and archive_root.exists():
        print("[ok] every sealed manifest and registration has an archive twin")

    # Derivative freshness. Three times a generator was fixed or an input
    # completed and the derivative was not re-produced -- and once a figure
    # was regenerated before its input run had finished, so the seal blessed
    # an output impossible to produce from its sealed input. mtime(output) >=
    # mtime(newest input) is a cheap necessary condition for "produced after
    # its inputs". It guards the working tree the package is built from; a
    # fresh checkout that flattens mtimes can pass --skip-mtime.
    if "--skip-mtime" not in sys.argv:
        derived = [
            # figure/table and record pairs written by one script are
            # siblings, not stages: both trace to the propagation records
            (METRICS / "r36_regime_map.json",
             ("r19_equal_total_work_*.json", "r14_trajectory_*_beta_*.json")),
            (FIGURES / "fig_regime_map.pdf",
             ("r19_equal_total_work_*.json", "r14_trajectory_*_beta_*.json")),
            # The main text's rendering of the same map, watched because it is
            # the one that went stale: it is drawn by a second invocation of
            # the same script and therefore against a second state of the
            # disk, and on 14 August it spent four hours showing a cell the
            # supplement's rendering correctly left blank. Watching only the
            # supplement's copy watched the copy that was right.
            (FIGURES / "fig_regime_map_main.pdf",
             ("r19_equal_total_work_*.json", "r14_trajectory_*_beta_*.json")),
            # Added after budget_pareto.pdf was found three days stale: it had
            # been drawn before design B's 2.00 and 3.00 cells landed, and the
            # rule below did not watch it, so the audit passed a figure the
            # records no longer supported. The patterns are anchored to designs
            # A and B, because these two figures are about the coverage designs
            # and a wider glob would make every later population mark them
            # stale without cause.
            (FIGURES / "budget_pareto.pdf",
             ("r14_trajectory_[AB]_beta_*.json", "r14_budget_pareto.json")),
            # The parity figure is drawn from the completed R42 panel (the
            # make_figures_r42.py view), not the R37 files it once read;
            # watching the R37 pair kept passing a figure the R42 records
            # could have outrun.
            (FIGURES / "fig_variational_parity.pdf",
             ("r42_variational_completion.json", "r42_panel_verdict.json")),
            # The allocation-interior figure (R60) reads the frozen (O28)
            # family records at beta = 1 on the two coverage designs.
            (FIGURES / "fig_allocation_interior.pdf",
             ("r18_span_sweep_[AB]_beta_1.00.json",)),
            # Anchored to designs A and B for the same reason budget_pareto
            # is: these derivatives are about the coverage designs, and a bare
            # wildcard makes every later population's records look like an
            # input to them. The paired ladder's ceiling-free run tripped
            # exactly that, marking seven unrelated tables stale.
            # The cut-sensitivity pair now covers the post-hoc budget's
            # populations too, so its inputs are no longer only A and B.
            # the ladder-geometry pair carries four rows of numbers into a
            # sentence, so it is watched like the tables that do the same
            (METRICS / "r50_ladder_geometry.json",
             ("r50_span_ladder_[ab]_design_frozen.json",)),
            (METRICS / "r50_ladder_geometry_table.tex",
             ("r50_span_ladder_[ab]_design_frozen.json",)),
            (METRICS / "r27_threshold_sensitivity.json",
             ("r19_equal_total_work_[AB]_beta_*.json",
              "r19_equal_total_work_C_beta_*.json",
              "r19_equal_total_work_SP_beta_*.json",
              "r19_equal_total_work_OEU_beta_*.json")),
            (METRICS / "r27_threshold_sensitivity_table.tex",
             ("r19_equal_total_work_[AB]_beta_*.json",
              "r19_equal_total_work_C_beta_*.json",
              "r19_equal_total_work_SP_beta_*.json",
              "r19_equal_total_work_OEU_beta_*.json")),
            (METRICS / "r28_equal_work_table.tex",
             ("r19_equal_total_work_[AB]_beta_*.json",)),
            (METRICS / "r28_midpoint_table.tex",
             ("r19_equal_total_work_[AB]_beta_*.json",)),
            (METRICS / "r19_equal_work_table.tex",
             ("r19_equal_total_work_[AB]_beta_*.json",)),
            (METRICS / "r18_span_table.tex",
             ("r18_span_sweep_[AB]_beta_1.00.json",)),
            (METRICS / "r18_by_k_table.tex",
             ("r18_span_sweep_[AB]_beta_*.json",)),
            (METRICS / "r44_equal_work_table.tex",
             ("r44_equal_work_tighter_*.json",
              "r18_span_sweep_[AB]_beta_*.json")),
            (METRICS / "r44_manuscript_descriptives.json",
             ("r44_equal_work_tighter_*.json",)),
            (METRICS / "r48_interior_timing_table.tex",
             ("r48_interior_timing.json",)),
            # The ladder's generated tables and the macro file the prose
            # counts come from. The macro file is watched too: it is not a
            # table but it carries numbers into sentences, which is the same
            # exposure with none of the visibility.
            (METRICS / "r50_span_ladder_primary_table.tex",
             ("r14_trajectory_RS[12]_beta_*.json",)),
            (METRICS / "r50_span_ladder_secondary_table.tex",
             ("r19_equal_total_work_RS[12]_beta_*.json",)),
            (METRICS / "r50_span_ladder_macros.tex",
             ("r14_trajectory_RS[12]_beta_*.json",
              "r19_equal_total_work_RS[12]_beta_*.json")),
            (METRICS / "r51_capped_vs_uncapped_table.tex",
             ("r51_verdict.json",)),
            (METRICS / "r51_capped_vs_uncapped_macros.tex",
             ("r51_verdict.json",)),
        ]
        stale_out: list[str] = []
        for out, in_pats in derived:
            if not out.exists():
                continue
            newest: tuple[float, str] | None = None
            for pat in in_pats:
                for q in METRICS.glob(pat):
                    m = q.stat().st_mtime
                    if newest is None or m > newest[0]:
                        newest = (m, q.name)
            if newest and out.stat().st_mtime < newest[0]:
                stale_out.append(f"{out.name} older than input {newest[1]}")
        if stale_out:
            problems.append("derivative older than its input: "
                            + "; ".join(stale_out))
            if not quiet:
                print(f"[stale] {len(stale_out)} derivatives older than "
                      f"their inputs:")
                for s in stale_out:
                    print(f"    {s}")
        elif not quiet:
            print("[ok] tracked derivatives are newer than their inputs")

    inventory = manuscript_inventory()
    for kind, wanted in inventory.items():
        gap = sorted(w for w in wanted if w not in covered)
        if gap:
            problems.append(f"manuscript {kind} not covered by any manifest: "
                            + ", ".join(gap))
            if not quiet:
                print(f"[gap] {len(gap)}/{len(wanted)} {kind} targets "
                      f"uncovered:")
                for g in gap[:20]:
                    print(f"    {g}")
                if len(gap) > 20:
                    print(f"    ... +{len(gap) - 20} more")
        elif not quiet:
            print(f"[ok] all {len(wanted)} manuscript {kind} targets covered")

    print(f"\n{len(paths)} manifests, {totals['match']} digests verified, "
          f"{totals['mismatch']} stale, {totals['absent']} absent, "
          f"{totals['flagged']} flagged-missing")
    if problems:
        print(f"blocking manifest problems: {len(problems)}")
        for p in problems[:40]:
            print(f"  - {p}")
        if len(problems) > 40:
            print(f"  ... +{len(problems) - 40} more")
        return 1
    print("blocking manifest problems: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
