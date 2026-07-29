"""Cross-check generated tables and figures against what the sources actually use.

Catches three things before submission: generated tables that no longer appear in
either document, `\\input` lines pointing at files that no longer exist, and
figure files that are referenced but missing (or shipped but orphaned).

Usage:  python check_assets.py
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = list((ROOT / "chapters").glob("*.tex")) + [ROOT / "main.tex",
                                                 ROOT / "supplement.tex"]

RE_INPUT = re.compile(r"\\input\{metrics/([^}]+?)(?:\.tex)?\}")
RE_FIG = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}")


def main() -> int:
    inputs, figs = set(), set()
    for p in SRC:
        s = p.read_text(encoding="utf-8")
        inputs |= {m + ".tex" for m in RE_INPUT.findall(s)}
        figs |= set(RE_FIG.findall(s))

    have_tab = {p.name for p in (ROOT / "metrics").glob("*.tex")}
    have_fig = {p.name for p in (ROOT / "figures").glob("*.pdf")}
    figs_norm = {f if f.endswith(".pdf") else f + ".pdf" for f in figs}

    problems = 0
    for title, items in (
            ("generated tables not used by either document", sorted(have_tab - inputs)),
            ("inputs with no file on disk", sorted(inputs - have_tab)),
            ("figures referenced but missing", sorted(figs_norm - have_fig)),
            ("figure files never referenced", sorted(have_fig - figs_norm))):
        print(f"\n{title}: {len(items)}")
        for i in items:
            print("   ", i)
        if "missing" in title or "no file" in title:
            problems += len(items)
    print(f"\nblocking problems: {problems}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
