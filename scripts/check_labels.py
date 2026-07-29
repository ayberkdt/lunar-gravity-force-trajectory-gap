"""Audit LaTeX labels and references across both documents.

Catches what a compile warning does not: labels defined twice, labels nothing
points at, cross-document references whose target no longer exists, and
references to a label defined in the wrong document.

The two documents are linked by `xr`: main.tex reads supplement.aux with the
prefix `supp-`, so a reference `\\ref*{supp-X}` in main resolves against a label
`X` defined in the supplement. The supplement does not read main's labels at all.

Usage:  python check_labels.py
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = [ROOT / "main.tex"]
SUPP = [ROOT / "supplement.tex"]
CHAPTERS = ROOT / "chapters"
METRICS = ROOT / "metrics"

RE_LABEL = re.compile(r"\\label\{([^}]+)\}")
RE_REF = re.compile(r"\\(?:ref|eqref|autoref|cref|Cref)\*?\{([^}]+)\}")
RE_INPUT = re.compile(r"\\input\{([^}]+?)(?:\.tex)?\}")


def expand(entry: Path, seen=None) -> list[Path]:
    """Follow \\input from a root file, so each document sees only its own tree."""
    seen = seen if seen is not None else set()
    if entry in seen or not entry.exists():
        return []
    seen.add(entry)
    out = [entry]
    for m in RE_INPUT.findall(entry.read_text(encoding="utf-8")):
        cand = ROOT / (m + ".tex")
        if cand.exists():
            out += expand(cand, seen)
    return out


def collect(files):
    labels, refs = Counter(), Counter()
    where = {}
    for f in files:
        s = f.read_text(encoding="utf-8")
        for m in RE_LABEL.findall(s):
            labels[m] += 1
            where.setdefault(m, []).append(f.name)
        for m in RE_REF.findall(s):
            refs[m] += 1
    return labels, refs, where


def main() -> int:
    main_files = expand(ROOT / "main.tex")
    supp_files = expand(ROOT / "supplement.tex")
    ml, mr, mw = collect(main_files)
    sl, sr, sw = collect(supp_files)

    problems = 0

    def report(title, items, blocking=True):
        nonlocal problems
        print(f"\n{title}: {len(items)}")
        for i in items:
            print("   ", i)
        if blocking:
            problems += len(items)

    report("duplicate labels in main",
           [f"{k} (x{v}) in {', '.join(mw[k])}" for k, v in ml.items() if v > 1])
    report("duplicate labels in supplement",
           [f"{k} (x{v}) in {', '.join(sw[k])}" for k, v in sl.items() if v > 1])

    # references inside main that are neither local nor a supp- prefixed target
    dangling_main = []
    for r in mr:
        if r in ml:
            continue
        if r.startswith("supp-") and r[len("supp-"):] in sl:
            continue
        dangling_main.append(r)
    report("references in main resolving nowhere", sorted(dangling_main))

    dangling_supp = [r for r in sr if r not in sl]
    report("references in supplement resolving nowhere", sorted(dangling_supp))

    # a supplement label that main points at with the wrong prefix
    mis = [r for r in mr if r not in ml and not r.startswith("supp-") and r in sl]
    report("main references a supplement label without the supp- prefix",
           sorted(mis))

    unused_main = sorted(k for k in ml if mr[k] == 0)
    unused_supp = sorted(k for k in sl
                         if sr[k] == 0 and mr.get("supp-" + k, 0) == 0)
    report("labels in main nothing points at", unused_main, blocking=False)
    report("labels in supplement nothing points at", unused_supp, blocking=False)

    print(f"\nblocking label problems: {problems}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
