"""Find passages the main text and the supplement both carry.

Some overlap is correct: the main text states a result and the supplement
documents how it was obtained. Overlap becomes redundancy when the same
explanation is written out twice at the same level of detail, which is what long
shared word sequences indicate.

Usage:  python find_main_supp_overlap.py [--n 12]
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ["01_introduction", "02_related", "03_kernel", "04_truncation",
        "05_setup", "06_qualification_main", "07_results", "08_discussion",
        "09_conclusion"]


def strip_tex(s: str) -> str:
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\(begin|end)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", s)
    s = re.sub(r"[${}~\\&_^]", " ", s)
    return re.sub(r"\s+", " ", s).lower()


def index(files, n):
    g = defaultdict(set)
    for f in files:
        if not f.exists():
            continue
        w = strip_tex(f.read_text(encoding="utf-8")).split()
        for i in range(len(w) - n):
            g[" ".join(w[i:i + n])].add(f.stem)
    return g


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    a = ap.parse_args()
    ch = ROOT / "chapters"
    mg = index([ch / f"{m}.tex" for m in MAIN], a.n)
    sg = index(sorted(ch.glob("supp_*.tex")) + [ch / "06_qualification.tex"], a.n)
    shared = set(mg) & set(sg)
    print(f"passages of {a.n}+ words present in both documents: {len(shared)}\n")
    shown = []
    for g in sorted(shared, key=len, reverse=True):
        if any(g in prev for prev in shown):
            continue
        shown.append(g)
        print(f"  main[{','.join(sorted(mg[g]))}] supp[{','.join(sorted(sg[g]))}]")
        print(f"    {g[:130]}\n")
        if len(shown) >= 12:
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
