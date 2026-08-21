#!/usr/bin/env python3
"""Report per-page layout facts that the .log file cannot see.

Two questions the LaTeX log never answers. First, does a float interrupt a
sentence: a top float on page N+1 is only a defect when page N stopped
mid-sentence, because then the reader meets half a page of figure between a
subject and its verb. Second, where are the pages that are mostly white.

Prose here means running text: not figure labels, not table cells, not
captions. All three are set at the text margin or centred and would otherwise
be read as the sentence a page ends on.

Usage:
  python inspect_layout.py main.pdf                 # float interruptions
  python inspect_layout.py supplement.pdf --gaps    # white space
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parent.parent
# Page text carries Greek and dashes the Windows console codepage cannot map.
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# A line ending a sentence, a heading or a display equation is a clean break;
# anything else carries the reader into the next page mid-thought. The trailing
# \s+\d+ is the running page number, which extraction sometimes appends to the
# last line; the space is required. Written \s*\d+ it also matched the tail of a
# decimal, so a page ending "...design B's 0.75" read as a finished sentence and
# the float that then opened the next page went unreported. That is the fourth
# defect found in this file and the second that hid an interruption rather than
# inventing one, so a change here is measured against a page known to break
# mid-sentence, not only against the total.
CLEAN_END = re.compile(r"(?:[.!?:](?:\s+\d+)?|\)|\})\s*$")
# A caption is a label pattern set in bold. The pattern alone is not enough:
# running text prints "Table 1. What the displacement is made of ..." and a
# caption prints "Table S1. Field-level acceleration ...", which are the same
# shape. Both document classes set the caption label bold and neither sets
# prose bold, so the typeface decides. Getting this wrong is expensive in both
# directions: a prose line read as a caption invents float interruptions, and a
# caption read as prose hides them and leaves caption text in the prose stream.
CAPTION_LABEL = re.compile(r"^\s*(?:Fig\.|Figure|Table)\s*S?\d+\b")
# Running text in this double-spaced single column reaches roughly 65
# characters; anything much shorter is a label, a tag or a caption tail.
BODY_CHARS = 55


class Line:
    __slots__ = ("y0", "y1", "x0", "text", "block", "size", "bold")

    def __init__(self, y0, y1, x0, text, block, size, bold=False):
        self.y0, self.y1, self.x0 = y0, y1, x0
        self.text, self.block, self.size = text, block, size
        self.bold = bold


def is_caption(line: "Line") -> bool:
    """A bold float label. See the note at CAPTION_LABEL."""
    return line.bold and bool(CAPTION_LABEL.match(line.text))


def page_lines(page: fitz.Page) -> list[Line]:
    """Every text line, with the folio dropped."""
    foot = page.rect.height - 72
    out = []
    for bi, block in enumerate(page.get_text("dict")["blocks"]):
        if block["type"] != 0:
            continue
        for line in block["lines"]:
            text = "".join(s["text"] for s in line["spans"]).strip()
            if not text:
                continue
            if line["bbox"][1] > foot and text.isdigit():
                continue
            size = max(s["size"] for s in line["spans"])
            first = line["spans"][0]
            bold = bool(first["flags"] & 16) or any(
                tag in first["font"] for tag in ("Bold", "BX", "bx"))
            out.append(Line(line["bbox"][1], line["bbox"][3],
                            line["bbox"][0], text, bi, size, bold))
    out.sort(key=lambda l: (l.y0, l.x0))
    return out


def caption_blocks(lines: list[Line]) -> set[int]:
    return {l.block for l in lines if is_caption(l)}


def prose_lines(lines: list[Line]) -> list[Line]:
    """Running text: long, starting at the modal left edge, outside a caption.

    A wide label inside a figure is long too, but it is centred, so the modal
    left edge separates prose from graphics; captions share that edge, so they
    are removed by block instead.
    """
    caps = caption_blocks(lines)
    long_ = [l for l in lines if len(l.text) >= BODY_CHARS]
    if not long_:
        return []
    edges = [round(l.x0) for l in long_]
    margin = max(set(edges), key=edges.count)
    # A paragraph's first line is one indent to the right of the margin, and
    # dropping it made every page that ends on a paragraph opening look as if
    # it ended on the line before -- which reported clean breaks that were not.
    indents = sorted({e for e in edges if 5 <= e - margin <= 30},
                     key=edges.count, reverse=True)
    body = max(set(round(l.size) for l in long_),
               key=[round(l.size) for l in long_].count)
    out = []
    for l in lines:
        if l.block in caps or round(l.size) != body or len(l.text) < 4:
            continue
        if abs(l.x0 - margin) < 3:
            # No length floor at the margin: the line that matters most is
            # often a paragraph's short last line.
            out.append(l)
        elif indents and abs(l.x0 - indents[0]) < 3 and len(l.text) >= BODY_CHARS:
            # At the indent, require full width. A paragraph opening always is;
            # a figure's axis label is not, and reading one as prose invents
            # interruptions.
            out.append(l)
    return out


def page_bands(page: fitz.Page, lines: list[Line]) -> list[tuple[float, float]]:
    bands = [(l.y0, l.y1) for l in lines]
    for block in page.get_text("dict")["blocks"]:
        if block["type"] == 1:
            bands.append((block["bbox"][1], block["bbox"][3]))
    # Vector art (matplotlib output, table rules) is not an image block.
    for d in page.get_drawings():
        r = d["rect"]
        if r.height > 2 and r.width > 20:
            bands.append((r.y0, r.y1))
    return bands


def biggest_gap(bands, top: float, bottom: float):
    """Largest vertical band carrying neither text nor graphics."""
    if not bands:
        return bottom - top, top, bottom
    bands = sorted(bands)
    merged = [list(bands[0])]
    for y0, y1 in bands[1:]:
        if y0 <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], y1)
        else:
            merged.append([y0, y1])
    best = (0.0, top, top)
    for (_, a1), (b0, _) in zip(merged, merged[1:]):
        if b0 - a1 > best[0]:
            best = (b0 - a1, a1, b0)
    if bottom - merged[-1][1] > best[0]:
        best = (bottom - merged[-1][1], merged[-1][1], bottom)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--pages", default=None, help="1-based range, e.g. 14-22")
    ap.add_argument("--gaps", action="store_true",
                    help="list pages whose largest empty band exceeds --min-gap")
    ap.add_argument("--min-gap", type=float, default=150.0)
    ap.add_argument("--all", action="store_true", help="print every page")
    args = ap.parse_args()

    path = Path(args.pdf)
    if not path.is_absolute():
        path = ROOT / path
    doc = fitz.open(path)

    if args.pages:
        lo, _, hi = args.pages.partition("-")
        rng = range(int(lo) - 1, int(hi or lo))
    else:
        rng = range(doc.page_count)

    prev_last = None
    prev_no = None
    splits = 0
    for i in rng:
        page = doc[i]
        lines = page_lines(page)
        if not lines:
            continue
        prose = prose_lines(lines)
        bands = page_bands(page, lines)
        top = min(b[0] for b in bands)
        gap, g0, g1 = biggest_gap(bands, top, page.rect.height - 78)

        cap_y = next((l.y0 for l in lines if is_caption(l)), None)
        opens_with_float = (
            cap_y is not None
            and not any(l.y0 < cap_y - 6 for l in prose)
            and any(l.y0 > cap_y + 6 for l in prose))
        split = (prev_last is not None and opens_with_float
                 and not CLEAN_END.search(prev_last))
        splits += split

        if args.gaps:
            if gap >= args.min_gap:
                print(f"p{i+1:>3}  largest empty band {gap:6.1f}pt "
                      f"(y {g0:.0f}-{g1:.0f})   opens: {lines[0].text[:48]!r}")
        elif split or args.all:
            print(f"p{i+1:>3}  opens with a float; p{prev_no} ended "
                  f"{prev_last[-52:]!r}" if split else
                  f"p{i+1:>3}  {lines[0].text[:60]!r}")
        if prose:
            prev_last, prev_no = prose[-1].text, i + 1
    if not args.gaps:
        print(f"float interruptions: {splits}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
