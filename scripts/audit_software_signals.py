"""Locate software-paper signals in the main manuscript.

Reports every occurrence of the branding and provenance vocabulary, with its
file, line and context, so each can be judged individually: needed for the
science, needed to distinguish independent implementations, or provenance-only
and better placed in the supplement.

Usage:  python audit_software_signals.py [--supplement]
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TERMS = ["Lunaris", "stack", "framework", "repository", "commit", "package",
         "release", "manifest", "SHA-256", "SHA", "module", "production"]

MAIN_CHAPTERS = ["01_introduction", "02_related", "03_kernel", "04_truncation",
                 "05_setup", "06_qualification_main", "07_results",
                 "08_discussion", "09_conclusion"]


def files(supplement: bool):
    if supplement:
        return sorted((ROOT / "chapters").glob("supp_*.tex")) + \
               [ROOT / "chapters" / "06_qualification.tex", ROOT / "supplement.tex"]
    return [ROOT / "main.tex"] + [ROOT / "chapters" / f"{c}.tex"
                                  for c in MAIN_CHAPTERS]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--supplement", action="store_true")
    a = ap.parse_args()
    counts = Counter()
    per_file = Counter()
    print(f"{'term':12s} {'file':26s} line  context")
    print("-" * 100)
    for f in files(a.supplement):
        if not f.exists():
            continue
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            for t in TERMS:
                if re.search(rf"\b{re.escape(t)}", line):
                    counts[t] += 1
                    per_file[(t, f.stem)] += 1
                    ctx = line.strip()
                    print(f"{t:12s} {f.stem:26s} {i:5d}  {ctx[:64]}")
                    break
    print("\ntotals by term:")
    for t, n in counts.most_common():
        print(f"  {t:12s} {n}")
    print("\nLunaris by file:")
    for (t, stem), n in sorted(per_file.items()):
        if t == "Lunaris":
            print(f"  {stem:26s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
