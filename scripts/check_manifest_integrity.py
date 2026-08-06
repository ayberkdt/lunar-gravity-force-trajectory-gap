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

This checks all four and exits non-zero on any of them, so it can gate a
submission build. It writes nothing.

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
# the records the supplement's partition claim is about.
TRAJECTORY = re.compile(r"^metrics/r\d+_(cases|raw)/")

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
    """What the manuscript actually pulls in, as repo-relative paths."""
    tex = sorted(glob.glob(str(CH / "*.tex"))) + [
        str(ROOT / "main.tex"), str(ROOT / "supplement.tex")]
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
        str(METRICS / "r[0-9]*_final_experiment_manifest.json")))
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
