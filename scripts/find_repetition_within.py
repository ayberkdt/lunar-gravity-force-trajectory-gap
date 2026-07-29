"""Find repeated phrasing *within* one section.

Cross-section repetition is caught by find_repetition.py. This catches the other
kind: a long section edited incrementally, where the same qualification gets
restated in three different paragraphs because each edit was made in isolation.

Usage:  python find_repetition_within.py chapters/07_results.tex [--n 8]
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def strip_tex(s: str) -> str:
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\(begin|end)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", s)
    s = re.sub(r"[${}~\\&_^]", " ", s)
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--n", type=int, default=8)
    a = ap.parse_args()
    text = strip_tex((ROOT / a.path).read_text(encoding="utf-8"))
    paras = [p for p in re.split(r"\n\s*\n", text) if len(p.split()) > 25]

    grams = defaultdict(set)
    for i, p in enumerate(paras):
        w = re.sub(r"\s+", " ", p).lower().split()
        for j in range(len(w) - a.n):
            grams[" ".join(w[j:j + a.n])].add(i)

    hits = {g: s for g, s in grams.items() if len(s) > 1}
    shown = []
    print(f"{len(paras)} paragraphs; {len(hits)} repeated {a.n}-word sequences\n")
    for g, s in sorted(hits.items(), key=lambda kv: (-len(kv[1]), -len(kv[0]))):
        if any(g in prev for prev in shown):
            continue
        shown.append(g)
        print(f"  paragraphs {sorted(s)}:  {g[:110]}")
        if len(shown) >= 20:
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
