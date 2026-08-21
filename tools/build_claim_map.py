"""Render the human-readable claim-to-artifact map.

    python tools/build_claim_map.py [--manuscript-root PATH]

Reads the machine-readable index in ``docs/REPRODUCIBILITY_INDEX.csv`` for the
item numbering, artifact, script and campaign of every manuscript float, and
takes the caption text from the manuscript sources so the maths survives as
readable Unicode rather than being flattened.

The manuscript sources are not part of this archive, so this needs a checkout
of them; ``--manuscript-root`` points at it and defaults to ``../codebase``.
Without it the captions already stored in the CSV are used, which are shorter
and have had their maths stripped.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "docs" / "REPRODUCIBILITY_INDEX.csv"
OUT = ROOT / "docs" / "claim_to_artifact_map.md"

GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "Delta": "Δ",
    "epsilon": "ε", "varepsilon": "ε", "rho": "ρ", "sigma": "σ", "Sigma": "Σ",
    "mu": "μ", "nu": "ν", "omega": "ω", "Omega": "Ω", "phi": "φ", "Phi": "Φ",
    "theta": "θ", "lambda": "λ", "tau": "τ", "chi": "χ", "psi": "ψ",
    "leq": "≤", "geq": "≥", "approx": "≈", "times": "×", "cdot": "·",
    "pm": "±", "to": "→", "ldots": "…", "infty": "∞", "propto": "∝",
    "langle": "⟨", "rangle": "⟩", "nabla": "∇", "partial": "∂", "in": "∈",
}

# purely presentational macros that should vanish rather than leave their name
DROP = {"hat", "bar", "tilde", "vec", "left", "right", "big", "Big", "bigl",
        "bigr", "quad", "qquad", "noindent", "centering", "small",
        "footnotesize", "textwidth", "linewidth", "par", "hfill", "medskip",
        "smallskip", "bigskip", "protect", "relax", "displaystyle"}

SUBSCRIPT = re.compile(r"_\{?\\mathrm\{([^}]*)\}\}?|_\{([^}]*)\}|_([A-Za-z0-9])")


def demath(text: str) -> str:
    """Turn a LaTeX fragment into readable plain text."""
    text = re.sub(r"\\(?:cite[a-z]*|ref\*?|label|nolinkurl)\{[^}]*\}", "", text)
    for _ in range(3):
        text = re.sub(
            r"\\(?:emph|textbf|textit|code|text|mathrm|mathbf)\s*\{([^{}]*)\}",
            r"\1", text)
    text = re.sub(r"\\frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}", r"\1/\2", text)
    text = text.replace("``", '"').replace("''", '"')
    # an unknown alphabetic macro keeps its own name, which is what
    # operator-style ones such as \max, \min and \deg want
    text = re.sub(r"\\([A-Za-z]+)",
                  lambda m: ("" if m.group(1) in DROP
                             else GREEK.get(m.group(1), m.group(1))), text)
    text = SUBSCRIPT.sub(
        lambda m: "_" + (m.group(1) or m.group(2) or m.group(3) or ""), text)
    text = re.sub(r"\^\{?([^{}\s]*)\}?", r"^\1", text)
    text = text.replace("$", "").replace("{", "").replace("}", "")
    text = text.replace("~", " ").replace("\\", "")
    text = re.sub(r"\s+([,.;:])", r"\1", text)
    text = re.sub(r" =(?=\d)", "=", text)
    return re.sub(r"\s+", " ", text).strip(" .,")


def raw_captions(manuscript_root: Path) -> dict[str, str]:
    """label -> raw caption, straight from the LaTeX sources."""
    out: dict[str, str] = {}
    sources: list[Path] = []
    chapters = manuscript_root / "chapters"
    if chapters.is_dir():
        sources += sorted(chapters.glob("*.tex"))
    sources += sorted((ROOT / "metrics").glob("*.tex"))
    for path in sources:
        text = path.read_text(encoding="utf-8", errors="replace")
        for env in (r"figure\*?", r"table\*?"):
            for block in re.findall(rf"\\begin{{{env}}}(.*?)\\end{{{env}}}",
                                    text, re.S):
                label = re.search(r"\\label\{([^}]+)\}", block)
                caption = re.search(r"\\caption\{(.*)", block, re.S)
                if not (label and caption):
                    continue
                depth, collected = 1, []
                for ch in caption.group(1):
                    if ch == "{":
                        depth += 1
                    elif ch == "}":
                        depth -= 1
                        if depth == 0:
                            break
                    collected.append(ch)
                out[label.group(1)] = "".join(collected)
    return out


def cell(text: str, limit: int = 120) -> str:
    text = (text or "").replace("|", "\\|").strip()
    if len(text) > limit:
        return text[:limit].rsplit(" ", 1)[0] + "…"
    return text or "—"


def code_list(text: str) -> str:
    if not text:
        return "—"
    return "<br>".join(f"`{p.strip()}`" for p in text.split(";") if p.strip())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-root", type=Path,
                        default=ROOT.parent / "codebase",
                        help="checkout of the manuscript sources "
                             "(default: ../codebase)")
    args = parser.parse_args()

    rows = list(csv.DictReader(CSV_PATH.open(encoding="utf-8")))
    captions = raw_captions(args.manuscript_root)
    if not captions:
        print(f"[warn] no captions found under {args.manuscript_root}; "
              f"falling back to the stripped captions in the CSV")

    lines = [
        "# Claim-to-artifact map\n\n",
        "Which script produced which manuscript item, and which campaign ",
        "manifest covers it. Generated by ",
        "[`tools/build_claim_map.py`](../tools/build_claim_map.py) from the ",
        "manuscript sources, the compiled label numbering and the campaign ",
        "manifests; the machine-readable version is ",
        "[`REPRODUCIBILITY_INDEX.csv`](REPRODUCIBILITY_INDEX.csv).\n\n",
        "Analysis and table passes read the archived records in `metrics/` in ",
        "place, so most rows can be re-run without propagating anything. Rows ",
        "marked *(inline)* are tables typed directly in the manuscript from ",
        "numbers reported in the text; they have no generated file.\n\n",
    ]

    for document, title in (("main", "Main text"),
                            ("supplement", "Supplement"),
                            ("none", "Generated but not used in the manuscript")):
        subset = [r for r in rows if r["document"] == document]
        if not subset:
            continue
        lines.append(f"## {title}\n\n")
        if document == "none":
            lines.append("Produced by the pipeline, but no manuscript item "
                         "references them.\n\n")
        lines.append("| Item | Claim | Artifact | Script | Campaign |\n")
        lines.append("|---|---|---|---|---|\n")
        for r in subset:
            claim = demath(captions.get(r["label"], r["claim"]))
            artifact = (f"`{r['artifact']}`" if r["artifact"] != "(inline)"
                        else "*(inline)*")
            lines.append(
                f"| **{r['item']}** | {cell(claim)} | {artifact} | "
                f"{code_list(r['script'])} | {r['campaign'] or '—'} |\n")
        lines.append("\n")

    OUT.write_text("".join(lines), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
