"""Compare re-run result records against the archived ones.

    python tools/compare_records.py ARCHIVE_DIR RERUN_DIR
    python tools/compare_records.py ARCHIVE_DIR RERUN_DIR --rtol 1e-12
    python tools/compare_records.py ARCHIVE_DIR RERUN_DIR --per-orbit

Answers the only question that matters when a driver is re-run: did any
reported number change? Wall-clock telemetry, timestamps, interpreter and
platform strings and absolute script paths are expected to differ on any
re-run and are ignored; everything else is compared value by value.

``--per-orbit`` is for re-runs that are still in progress or that stopped at a
deadline. A campaign covering a subset of the orbits has legitimately different
aggregates, which would swamp the report, so this mode compares the per-orbit
rows on the orbits both sides have and reports the coverage separately. That is
the comparison that establishes whether the driver still behaves identically;
the aggregates follow once coverage is complete.

Exit status is 0 when no scientific difference was found.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

# Keys whose values legitimately change on any re-run.
IGNORE_KEYS = {
    "created_utc", "recorded_utc", "generated_utc", "timestamp",
    "started_utc", "ended_utc",
    "repo_commit_sha", "repo_working_tree_clean",
    "source", "provenance", "environment", "telemetry",
    "wall_time_s", "elapsed_s", "runtime_s", "duration_s",
    "experiment_script", "platform", "python", "host", "machine",
    # The driver's own digest is the thing under investigation, not evidence
    # about the numbers, so comparing it here would only restate the question.
    "script_sha256",
}

# Suffixes that mark a timing measurement.
IGNORE_SUFFIXES = ("_ns", "_wall_s", "_cpu_s", "_walltime", "_seconds")


def ignored(key: str) -> bool:
    return key in IGNORE_KEYS or key.endswith(IGNORE_SUFFIXES)


class Diff:
    def __init__(self) -> None:
        self.numeric: list[tuple[str, float, float, float]] = []
        self.structural: list[str] = []

    def __bool__(self) -> bool:
        return bool(self.numeric or self.structural)


def compare(a, b, path: str, rtol: float, out: Diff) -> None:
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if ignored(key):
                continue
            if key not in a:
                out.structural.append(f"{path}/{key}: only in re-run")
            elif key not in b:
                out.structural.append(f"{path}/{key}: only in archive")
            else:
                compare(a[key], b[key], f"{path}/{key}", rtol, out)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.structural.append(f"{path}: length {len(a)} vs {len(b)}")
            return
        for i, (u, v) in enumerate(zip(a, b)):
            compare(u, v, f"{path}[{i}]", rtol, out)
        return
    if isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            out.structural.append(f"{path}: {a} vs {b}")
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a == b:
            return
        if math.isnan(a) and math.isnan(b):
            return
        scale = max(abs(a), abs(b), 1e-300)
        rel = abs(a - b) / scale
        if rel > rtol:
            out.numeric.append((path, float(a), float(b), rel))
        return
    if a != b:
        out.structural.append(f"{path}: {str(a)[:70]!r} vs {str(b)[:70]!r}")


ROW_KEYS = ("sobol_index", "orbit_index", "index", "name", "orbit")


def row_index(record: dict) -> tuple[str, dict] | None:
    """Return (key name, {key value: row}) if the record has per-orbit rows."""
    rows = record.get("rows")
    if not isinstance(rows, list) or not rows:
        return None
    if not all(isinstance(r, dict) for r in rows):
        return None
    for key in ROW_KEYS:
        if all(key in r for r in rows):
            return key, {r[key]: r for r in rows}
    return None


def compare_per_orbit(a: dict, b: dict, rtol: float, diff: Diff) -> tuple[int, int, int] | None:
    """Compare only the rows both sides carry. Returns coverage counts."""
    left, right = row_index(a), row_index(b)
    if left is None or right is None or left[0] != right[0]:
        return None
    key, archive_rows = left
    _, rerun_rows = right
    common = sorted(set(archive_rows) & set(rerun_rows), key=str)
    for value in common:
        compare(archive_rows[value], rerun_rows[value], f"{key}={value}",
                rtol, diff)
    return len(archive_rows), len(rerun_rows), len(common)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("rerun", type=Path)
    parser.add_argument("--rtol", type=float, default=0.0,
                        help="relative tolerance for floats (default: exact)")
    parser.add_argument("--quiet", action="store_true",
                        help="only report files that differ")
    parser.add_argument("--per-orbit", action="store_true",
                        help="compare per-orbit rows on the shared orbits and "
                             "report coverage, instead of whole records")
    args = parser.parse_args()

    names = sorted(p.name for p in args.archive.glob("*.json"))
    identical, differing, skipped = [], [], []

    for name in names:
        rerun_path = args.rerun / name
        if not rerun_path.exists():
            skipped.append(name)
            continue
        try:
            a = json.loads((args.archive / name).read_text(encoding="utf-8"))
            b = json.loads(rerun_path.read_text(encoding="utf-8"))
        except Exception as exc:
            differing.append((name, Diff()))
            print(f"[error]  {name}: {exc}")
            continue
        diff = Diff()
        coverage = None
        if args.per_orbit:
            coverage = compare_per_orbit(a, b, args.rtol, diff)
        if coverage is None:
            compare(a, b, "", args.rtol, diff)
        if diff:
            differing.append((name, diff))
        else:
            identical.append(name)
        if coverage and not args.quiet:
            n_archive, n_rerun, n_common = coverage
            partial = "" if n_common == n_archive else "  (partial re-run)"
            print(f"[rows]    {name}: archive {n_archive}, re-run {n_rerun}, "
                  f"compared {n_common}{partial}")
        elif coverage:
            n_archive, n_rerun, n_common = coverage
            if n_common != n_archive:
                print(f"[rows]    {name}: compared {n_common} of "
                      f"{n_archive} orbit(s)")

    for name, diff in differing:
        print(f"[differs] {name}")
        for path, x, y, rel in diff.numeric[:8]:
            print(f"    {path}: archive {x!r} vs re-run {y!r}  rel={rel:.3e}")
        if len(diff.numeric) > 8:
            print(f"    ... {len(diff.numeric) - 8} more numeric difference(s)")
        for line in diff.structural[:8]:
            print(f"    {line}")
        if len(diff.structural) > 8:
            print(f"    ... {len(diff.structural) - 8} more structural difference(s)")

    if not args.quiet:
        for name in identical:
            print(f"[same]    {name}")

    print(f"\ncompared {len(identical) + len(differing)} record(s): "
          f"{len(identical)} identical, {len(differing)} differing, "
          f"{len(skipped)} not re-run")
    return 1 if differing else 0


if __name__ == "__main__":
    raise SystemExit(main())
