"""Find repeated content across the main manuscript.

Two passes. The first reports long shared word sequences between different
sections, which is where copy-adjacent restatement shows up. The second reports
how often a set of load-bearing claims is restated, so a claim that has to
appear three times can be told apart from one that drifted into six.

Usage:  python find_repetition.py [--n 12] [--supplement]
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

# claims that legitimately recur, to count rather than to flag blindly
CLAIMS = {
    "2.7-2.8x operating point": r"2\.7.{0,12}2\.8",
    "nominal per-call budget": r"nominal per-call",
    "force vs trajectory disagreement": r"(force|defect).{0,60}(not|does not).{0,40}traject",
    "unresolved is not a tie": r"unresolved|never counted as (a )?tie",
    "29% total work excess": r"29\\?%",
    "median 4.6 and 5.0": r"4\.59|4\.6.{0,10}5\.0|\$4\.6\$",
    "crossing bracket 0.50-0.75": r"0\.50.{0,30}0\.75|bracket",
    "cap audit / truth degree": r"touch(es)? the cap|cap-free|adopted truth degree",
    "resolution rule statement": r"truth-inclusive|resolution rule",
    "1.5 decades offset": r"decade",
}


def strip_tex(s: str) -> str:
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\(begin|end)\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", s)
    s = re.sub(r"[${}~\\&_^]", " ", s)
    return re.sub(r"\s+", " ", s)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--supplement", action="store_true")
    a = ap.parse_args()
    names = (sorted(p.stem for p in (ROOT / "chapters").glob("supp_*.tex"))
             if a.supplement else MAIN)
    texts = {}
    for n in names:
        p = ROOT / "chapters" / f"{n}.tex"
        if p.exists():
            texts[n] = strip_tex(p.read_text(encoding="utf-8"))

    grams = defaultdict(set)
    for name, t in texts.items():
        w = t.lower().split()
        for i in range(len(w) - a.n):
            grams[" ".join(w[i:i + a.n])].add(name)

    cross = {g: s for g, s in grams.items() if len(s) > 1}
    print(f"shared {a.n}-word sequences across sections: {len(cross)}")
    seen = []
    for g, s in sorted(cross.items(), key=lambda kv: -len(kv[0])):
        if any(g in prev for prev in seen):
            continue
        seen.append(g)
        print(f"\n  [{', '.join(sorted(s))}]")
        print(f"    {g[:150]}")
        if len(seen) >= 15:
            break

    print("\n\nrestatement counts of load-bearing claims:")
    print(f"{'claim':38s} " + " ".join(f"{n[:6]:>7s}" for n in texts))
    for label, pat in CLAIMS.items():
        counts = {n: len(re.findall(pat, t, re.I)) for n, t in texts.items()}
        tot = sum(counts.values())
        if tot == 0:
            continue
        print(f"{label:38s} " + " ".join(f"{counts[n]:7d}" for n in texts)
              + f"   total {tot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
