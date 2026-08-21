"""Refresh the digests a manifest already indexes, and nothing else.

A manifest goes stale whenever a file it hashes is legitimately edited: a
generated table regenerated after a symbol was renamed, a figure script
amended, a chapter revised. The obvious repair, re-running the campaign's own
finalizer, is not available for every manifest here. Several finalizers
discover their contents by globbing the disk, and the disk has since acquired
records belonging to later campaigns: re-running R14's finalizer today would
sweep R53's beta = 0.62 records and case trees into R14's inventory and break
the partition the integrity check enforces. The finalizer is the right tool the
day a campaign ends and the wrong one afterwards.

So this tool does the narrow thing instead. It walks a manifest's existing
entries, recomputes the digest of each file that is already indexed, and writes
the manifest back with its self-seal recomputed. It never adds an entry, never
removes one and never touches ownership: the set of files a manifest claims
before and after is identical, and the only thing that can change is what those
files hash to. Anything that would change ownership is a finalizer's job and is
left to fail loudly at the gate.

Every refreshed entry is printed with both digests, because a re-seal that says
nothing about what moved is indistinguishable from one that hid something.

Usage:
    python reseal_stale_digests.py                 # report, change nothing
    python reseal_stale_digests.py --write r14 r35 # refresh those manifests
    python reseal_stale_digests.py --write --all
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import check_manifest_integrity as chk

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

HEX = re.compile(r"^[0-9a-f]{64}$")

# What this tool must never touch. A propagation output is the evidence; if its
# digest no longer matches, the answer is that something rewrote a result, and
# the gate is supposed to say so rather than be talked out of it. Derived
# summaries, generated tables, figures and scripts are a different class: they
# are recomputed from the evidence, and a re-seal after a symbol rename or an
# amended comment is bookkeeping. The line between the two is drawn here, by
# name, and refusing is not a warning: the run stops.
PROTECTED = (
    re.compile(r"(^|/)r(?:\d+|J\d*)_(?:cases|raw)/"),
    re.compile(r"(^|/)r(?:\d+|J\d*)_trajectory_[^/]*\.json$"),
    re.compile(r"(^|/)r(?:\d+|J\d*)_span_sweep_[^/]*\.json$"),
    re.compile(r"(^|/)r(?:\d+|J\d*)_equal_total_work_[^/]*\.json$"),
    re.compile(r"(^|/)r(?:\d+|J\d*)_(?:budget_pareto|preregistration|verdict)"
               r"[^/]*\.json$"),
    # A departure record states what a campaign decided when its registration
    # did not cover the outcome. That is evidence, and it was missing from
    # this list until a hand edit to R51's record passed through the tool as
    # bookkeeping. Records are rewritten by their own generator or not at all.
    re.compile(r"(^|/)r(?:\d+|J\d*)_registration_departure\.json$"),
)


def protected(name: str) -> bool:
    n = name.replace("\\", "/")
    return any(p.search(n) for p in PROTECTED)


def digest(path: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest(), path.stat().st_size


def refresh(node, section: str, changes: list, refused: list) -> None:
    """Mirror of check_manifest_integrity.collect, writing instead of reading."""
    if isinstance(node, list):
        for val in node:
            refresh(val, section, changes, refused)
        return
    if not isinstance(node, dict):
        return

    name = node.get("path") or node.get("file") or node.get("relpath")
    if isinstance(node.get("sha256"), str) and HEX.match(node["sha256"]) \
            and isinstance(name, str):
        found = chk.resolve(name, section)
        if found is not None:
            new, size = digest(found)
            if new != node["sha256"]:
                changes.append((name, node["sha256"], new))
                if protected(name):
                    refused.append(name)
                else:
                    node["sha256"] = new
                    if "bytes" in node:
                        node["bytes"] = size

    for key, val in list(node.items()):
        sub = key if key not in ("entries",) else section
        if key in chk.NON_FILE_KEYS and not isinstance(val, dict):
            continue
        if isinstance(val, str):
            if HEX.match(val) and chk.FILE_SUFFIX.search(key):
                found = chk.resolve(key, section)
                if found is not None:
                    new, _ = digest(found)
                    if new != val:
                        changes.append((key, val, new))
                        if protected(key):
                            refused.append(key)
                        else:
                            node[key] = new
        elif isinstance(val, dict):
            inner = val.get("sha256")
            if isinstance(inner, str) and HEX.match(inner) and \
                    chk.FILE_SUFFIX.search(key):
                found = chk.resolve(key, section)
                if found is not None:
                    new, size = digest(found)
                    if new != inner:
                        changes.append((key, inner, new))
                        if protected(key):
                            refused.append(key)
                        else:
                            val["sha256"] = new
                            if "bytes" in val:
                                val["bytes"] = size
            refresh(val, sub, changes, refused)
        elif isinstance(val, list):
            refresh(val, sub, changes, refused)


def indent_of(text: str) -> int:
    """The file's own indentation, so a re-seal is not also a reformat."""
    for line in text.split("\n")[1:]:
        if line.strip():
            return len(line) - len(line.lstrip(" "))
    return 1


def process(path: Path, write: bool, journal: list) -> tuple[int, int]:
    text = path.read_text(encoding="utf-8")
    payload = json.loads(text)
    changes: list = []
    refused: list = []
    refresh(payload, "", changes, refused)
    if not changes:
        print(f"[ok] {path.name}: nothing stale")
        return 0, 0
    print(f"[{'reseal' if write else 'would reseal'}] {path.name}: "
          f"{len(changes)} entries")
    for name, old, new in changes:
        mark = "  REFUSED (result record)" if name in refused else ""
        print(f"    {name}{mark}\n      {old[:16]} -> {new[:16]}")
    if refused:
        return len(changes), len(refused)
    if not write:
        return len(changes), 0
    if "manifest_sha256" in payload:
        body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
        payload["manifest_sha256"] = hashlib.sha256(
            json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
    path.write_text(json.dumps(payload, indent=indent_of(text)),
                    encoding="utf-8")
    journal.append({
        "manifest": path.name,
        "manifest_sha256_after": payload.get("manifest_sha256"),
        "entries": [{"name": n, "from": o, "to": w} for n, o, w in changes]})
    return len(changes), 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("registries", nargs="*")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--stamp", default="",
                    help="timestamp recorded in the log; the tool does not "
                         "read a clock, so a dated log needs one passed in")
    a = ap.parse_args()

    if a.all or not a.registries:
        paths = sorted(METRICS.glob("*_final_experiment_manifest.json"))
    else:
        paths = [METRICS / f"{r}_final_experiment_manifest.json"
                 for r in a.registries]
    total, refused, journal = 0, 0, []
    for p in paths:
        if not p.exists():
            print(f"[absent] {p.name}")
            return 2
        n, r = process(p, a.write, journal)
        total += n
        refused += r
    if refused:
        print(f"[abort] {refused} stale entries are result records. A "
              f"propagation output whose digest moved is not a bookkeeping "
              f"problem; nothing was written.")
        return 1
    if a.write and journal:
        log = METRICS / "reseal_log.json"
        prior = []
        if log.exists():
            prior = json.loads(log.read_text(encoding="utf-8")).get("runs", [])
        prior.append({"stamp": a.stamp, "manifests": journal})
        log.write_text(json.dumps(
            {"schema": "reseal_log_v1",
             "what_this_is": (
                 "every digest this tool refreshed, with both values. It "
                 "refreshes only entries a manifest already indexes and "
                 "refuses result records, so this log is the whole of what a "
                 "re-seal can do."),
             "runs": prior}, indent=1), encoding="utf-8")
        print(f"[written] {log.name}")
    print(f"{total} entries {'refreshed' if a.write else 'stale'} "
          f"across {len(paths)} manifests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
