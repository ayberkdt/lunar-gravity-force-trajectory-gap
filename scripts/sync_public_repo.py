"""Bring the public repository up to date with this working tree.

The repository's copies of the shared drivers are not a fork: they are older
versions of the same files, so the sync overwrites them rather than merging.
What makes that safe is that the campaign manifests travel with them, so the
digests the manifests record and the files they name move in one step and the
chain closes again on the other side. After an apply, run the repository's
``audit_manifest_digests.py --strict`` to confirm the chain did close.

  copied      every driver and every metrics record the working tree has and
              the repository does not, plus every shared file whose contents
              differ. Figures too.

  excluded    what the repository's own .gitignore already excludes: raw
              trajectory arrays under metrics/*_raw/, run logs, lock
              directories, and the large external coefficient products.

  added       the evidence-archive deliverables that did not exist before:
              the migration map, the archive index, the campaign map, the
              reproducibility README, the archive manifest and checksums, and
              the tables that left the supplement PDF. These go under
              archive/ and are not duplicated from metrics/.

Nothing is deleted from the repository. A file the repository has and the
working tree does not is left alone and reported, because that is a question
for a person rather than a script.

Usage:  python sync_public_repo.py --dry-run
        python sync_public_repo.py --apply
"""

from __future__ import annotations

import argparse
import filecmp
import re
import shutil
from pathlib import Path

CB = Path(__file__).resolve().parents[1]
RP = CB.parent / "lunar-gravity-force-trajectory-gap"

PAIRS = [
    (CB / "python_codes", RP / "scripts"),
    (CB / "metrics", RP / "metrics"),
    (CB / "figures", RP / "figures"),
]
ARCHIVE_ONLY = [
    "supplement_archive/migration_map.csv",
    "supplement_archive/archive_index.md",
    "manifests/campaign_map.csv",
    "manifests/MANIFEST.json",
    "manifests/CHECKSUMS.sha256",
    "README_REPRODUCIBILITY.md",
]


# Quarantined partial records: a supervisor renames a half-written JSON to
# NAME.invalid.<UTC timestamp> instead of deleting it. Operational debris, not
# manifest-sealed evidence, so it stays out of the public archive.
QUARANTINED = re.compile(r"\.invalid\.\d{8}T\d{6}Z$")


def excluded(rel: str) -> bool:
    return ("_raw/" in rel or rel.endswith((".log", ".err", ".out"))
            or "__pycache__" in rel
            or rel.endswith(".pyc") or "/.lock" in rel
            or QUARANTINED.search(rel) is not None
            or Path(rel).name.startswith("morning_"))


def plan(src: Path, dst: Path):
    new, changed, only_repo = [], [], []
    for p in sorted(src.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(src).as_posix()
        if excluded(rel):
            continue
        q = dst / rel
        if not q.exists():
            new.append(rel)
        elif p.stat().st_size != q.stat().st_size or not filecmp.cmp(p, q, shallow=False):
            changed.append(rel)
    for q in sorted(dst.rglob("*")):
        if not q.is_file():
            continue
        rel = q.relative_to(dst).as_posix()
        if excluded(rel):
            continue
        if not (src / rel).exists():
            only_repo.append(rel)
    return new, changed, only_repo


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    if not RP.exists():
        print(f"[abort] repository not found at {RP}")
        return 2

    total_new = total_changed = 0
    for src, dst in PAIRS:
        new, changed, only_repo = plan(src, dst)
        total_new += len(new)
        total_changed += len(changed)
        print(f"\n{dst.name}/  new {len(new)}, changed {len(changed)}, "
              f"only-in-repo {len(only_repo)}")
        for rel in only_repo[:10]:
            print(f"   [only in repo, left alone] {rel}")
        if a.apply:
            for rel in new + changed:
                d = dst / rel
                d.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src / rel, d)

    arch_dst = RP / "archive"
    n_arch = 0
    for rel in ARCHIVE_ONLY:
        s = CB / "archive" / rel
        if s.exists():
            n_arch += 1
            if a.apply:
                d = arch_dst / rel
                d.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(s, d)
    legacy = sorted((CB / "archive" / "supplement_archive"
                     / "full_legacy_tables").glob("*.tex"))
    csvs = sorted((CB / "archive" / "processed" / "per_orbit_tables").glob("*.csv"))
    for group, sub in ((legacy, "supplement_archive/full_legacy_tables"),
                       (csvs, "processed/per_orbit_tables")):
        for s in group:
            n_arch += 1
            if a.apply:
                d = arch_dst / sub / s.name
                d.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(s, d)
    print(f"\narchive/  {n_arch} deliverables "
          f"({len(legacy)} moved tables, {len(csvs)} CSV renderings)")

    print(f"\ntotal: {total_new} new, {total_changed} updated, {n_arch} archive")
    print("[dry run, nothing written]" if not a.apply else "[applied]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
