"""Emit a design group as a rotated spanning label instead of a banner row.

Several tables in this paper are two blocks of rows, one per coverage design.
The original form put a full-width ``Design A'' row above each block. That
costs a line, interrupts the column rules twice, and still leaves the reader to
work out where the group ends. A rotated label in a leading column spans the
block and shows the boundary directly.

The transformation is the same in every table that needs it, so it lives here
once rather than five times. Callers pass the rows of each block already
formatted, minus the leading ampersand; this returns the body lines, and
``prefix_colspec`` gives the ``c`` column to add at the front of the tabular
spec and the empty leading cell for the header row.

Requires ``multirow`` and ``graphicx`` in the preamble.
"""

from __future__ import annotations


def prefix_colspec(spec: str) -> str:
    """Insert the label column at the front of a tabular column spec."""
    if spec.startswith("@{}"):
        return "@{}c " + spec[3:]
    return "c " + spec


def header_row(header: str) -> str:
    """Push an existing header row right by one empty cell."""
    return "& " + header


def blocks(groups, rule: str = "\\midrule") -> list[str]:
    """Body lines for a table grouped by design.

    ``groups`` is a sequence of ``(label, rows)``. Each row is the LaTeX for
    one line \\emph{without} the leading label column and without a leading
    ampersand; the row is expected to end in ``\\\\``. Empty groups are
    skipped, so a table whose second design has not been propagated yet
    degrades to a single block rather than emitting an empty span.
    """
    out: list[str] = []
    real = [(label, rows) for label, rows in groups if rows]
    for i, (label, rows) in enumerate(real):
        if i:
            out.append(rule)
        span = (f"\\multirow{{{len(rows)}}}{{*}}"
                f"{{\\rotatebox[origin=c]{{90}}{{\\emph{{{label}}}}}}}")
        out.append(f"{span} & {rows[0]}")
        out += [f" & {r}" for r in rows[1:]]
    return out
