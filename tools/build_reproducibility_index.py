"""Regenerate docs/REPRODUCIBILITY_INDEX.csv from the manuscript and manifests.

    python tools/build_reproducibility_index.py [--manuscript-root PATH]
    python tools/build_reproducibility_index.py --check

Walks every figure and table environment in the manuscript sources, takes the
item numbering and the caption from the compiled ``.aux`` files, and attributes
each generated artifact to the campaign manifests that seal it. The scripts
column lists the manifest-registered scripts whose source mentions the
artifact by name; where no manifest seals an artifact (the early field-level
figures), the attribution recorded in the previous index is carried forward.

Artifacts that a manifest seals but no manuscript item references are listed
at the end with document ``none``, so nothing sealed can silently disappear
from the map.

This needs a checkout of the manuscript sources with fresh ``main.aux`` and
``supplement.aux``; ``--manuscript-root`` points at it and defaults to
``../codebase``. ``--check`` rewrites nothing and exits non-zero if the
committed index no longer matches what would be generated.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "REPRODUCIBILITY_INDEX.csv"

FIELDS = ["item", "document", "label", "claim", "artifact", "script",
          "campaign", "manifest", "defined_in"]

NEWLABEL = re.compile(r"\\newlabel\{([^}]+)\}\{")
INPUT_CHAPTER = re.compile(r"\\input\{(chapters/[^}]+?)(?:\.tex)?\}")
GRAPHICS = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
TABLE_INPUT = re.compile(r"\\input\{((?:\.\./)?metrics/[^}]+?)(?:\.tex)?\}")


def balanced_groups(text: str, start: int, count: int) -> list[str]:
    """Read ``count`` consecutive {...} groups with brace matching."""
    groups: list[str] = []
    i = start
    for _ in range(count):
        while i < len(text) and text[i] != "{":
            i += 1
        depth, begin = 0, i
        while i < len(text):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        groups.append(text[begin + 1:i])
        i += 1
    return groups


def aux_entries(aux_path: Path) -> dict[str, dict[str, str]]:
    """label -> {number, caption, kind} from a compiled .aux file."""
    entries: dict[str, dict[str, str]] = {}
    text = aux_path.read_text(encoding="utf-8", errors="replace")
    for match in NEWLABEL.finditer(text):
        label = match.group(1)
        if label.endswith("@cref") or "." in label.split(":")[0]:
            continue
        wrapper = balanced_groups(text, match.end() - 1, 1)[0]
        number, _page, caption, anchor, _ = balanced_groups(wrapper, 0, 5)
        kind = ("figure" if anchor.startswith("figure")
                else "table" if anchor.startswith("table") else "")
        if kind and label not in entries:
            entries[label] = {"number": number, "caption": caption,
                              "kind": kind}
    return entries


def chapter_files(manuscript_root: Path, top: str) -> list[Path]:
    """The top-level file plus every chapter it inputs, in input order."""
    top_path = manuscript_root / top
    files = [top_path]
    text = top_path.read_text(encoding="utf-8", errors="replace")
    for rel in INPUT_CHAPTER.findall(text):
        path = manuscript_root / f"{rel}.tex"
        if path.exists():
            files.append(path)
    return files


def float_blocks(path: Path) -> list[dict[str, str]]:
    """Every labeled figure/table environment referenced from one source file.

    A table environment may live either in the chapter itself (inline tables,
    and tables that ``\\input`` their tabular body) or entirely inside a
    generated file under ``metrics/`` that the chapter inputs at top level.
    Both shapes end up as one block with the chapter as the referencing file.
    """
    blocks: list[dict[str, str]] = []
    text = path.read_text(encoding="utf-8", errors="replace")
    inside_env: set[str] = set()
    for env in ("figure", "table"):
        pattern = rf"\\begin{{{env}\*?}}(.*?)\\end{{{env}\*?}}"
        for body in re.findall(pattern, text, re.S):
            label = re.search(r"\\label\{([^}]+)\}", body)
            if not label:
                continue
            if env == "figure":
                artifacts = [Path(g).name for g in GRAPHICS.findall(body)]
            else:
                artifacts = [f"{Path(g).name}.tex"
                             for g in TABLE_INPUT.findall(body)]
                inside_env.update(artifacts)
            blocks.append({"label": label.group(1), "kind": env,
                           "artifacts": "; ".join(artifacts) or "(inline)"})
    for name in (f"{Path(g).name}.tex" for g in TABLE_INPUT.findall(text)):
        if name in inside_env:
            continue
        generated = ROOT / "metrics" / name
        if not generated.exists():
            continue
        body = generated.read_text(encoding="utf-8", errors="replace")
        label = re.search(r"\\label\{([^}]+)\}", body)
        if label:
            blocks.append({"label": label.group(1), "kind": "table",
                           "artifacts": name})
    return blocks


def demath(text: str) -> str:
    from build_claim_map import demath as claim_demath
    return claim_demath(text)


def load_manifest_seals(metrics: Path) -> dict[str, list[str]]:
    """artifact basename -> manifest filenames that seal it."""
    seals: dict[str, list[str]] = {}
    for manifest_path in sorted(metrics.glob("r*_final_experiment_manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for section in ("generated_tables", "generated_figures"):
            for name in (manifest.get(section) or {}):
                seals.setdefault(Path(name).name, []).append(manifest_path.name)
    return seals


def campaign_of(manifest_name: str) -> str:
    stem = manifest_name.split("_", 1)[0]
    body = stem[1:] if stem.startswith("r") else stem
    return f"R{body}" if body.isdigit() else body.upper()


class ScriptAttributor:
    """Which registered scripts mention an artifact by name."""

    def __init__(self, scripts_dir: Path, metrics: Path) -> None:
        self.sources: dict[str, str] = {}
        for path in sorted(scripts_dir.glob("*.py")):
            self.sources[path.name] = path.read_text(encoding="utf-8",
                                                     errors="replace")
        self.registered: dict[str, list[str]] = {}
        for manifest_path in sorted(
                metrics.glob("r*_final_experiment_manifest.json")):
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.registered[manifest_path.name] = [
                Path(name).name for name in (manifest.get("scripts") or {})]

    VALIDATOR_PREFIXES = ("audit_", "campaign_ownership", "check_",
                          "compare_", "record_known", "submission_",
                          "test_", "verify_")

    def attribute(self, artifact: str, manifests: list[str]) -> str:
        candidates: list[str] = []
        for manifest_name in manifests:
            candidates += self.registered.get(manifest_name, [])
        if not manifests:
            candidates = [name for name in self.sources
                          if not name.startswith(self.VALIDATOR_PREFIXES)]
        candidates = list(dict.fromkeys(candidates))
        hits = [name for name in candidates
                if artifact in self.sources.get(name, "")]
        if not hits:
            stem = Path(artifact).stem
            hits = [name for name in candidates
                    if stem in self.sources.get(name, "")]
        if not hits and manifests:
            # The writer may predate the sealing campaign (a re-sealed
            # artifact), so widen to every non-validator script.
            return self.attribute(artifact, [])
        return "; ".join(sorted(hits))


def build_rows(manuscript_root: Path) -> list[dict[str, str]]:
    previous = {row["artifact"]: row
                for row in csv.DictReader(OUT.open(encoding="utf-8"))
                if row["artifact"] not in ("", "(inline)")}
    seals = load_manifest_seals(ROOT / "metrics")
    attributor = ScriptAttributor(ROOT / "scripts", ROOT / "metrics")

    rows: list[dict[str, str]] = []
    referenced: set[str] = set()
    for document, top, aux_name in (("main", "main.tex", "main.aux"),
                                    ("supplement", "supplement.tex",
                                     "supplement.aux")):
        aux = aux_entries(manuscript_root / aux_name)
        doc_rows: list[tuple[tuple[int, int], dict[str, str]]] = []
        seen: set[str] = set()
        for path in chapter_files(manuscript_root, top):
            for block in float_blocks(path):
                label = block["label"]
                entry = aux.get(label)
                if entry is None or label in seen:
                    if entry is None:
                        print(f"[warn] {label} in {path.name} has no "
                              f"{aux_name} entry; skipped")
                    continue
                seen.add(label)
                number = entry["number"]
                kind_title = block["kind"].capitalize()
                item = (f"{kind_title} {number}" if document == "main"
                        else f"S{kind_title} {number}")
                artifact = block["artifacts"]
                manifests: list[str] = []
                script = campaign = manifest_cell = ""
                if artifact != "(inline)":
                    for name in artifact.split("; "):
                        referenced.add(name)
                        manifests += [m for m in seals.get(name, [])
                                      if m not in manifests]
                    if manifests:
                        manifest_cell = "; ".join(
                            f"metrics/{m}" for m in manifests)
                        campaign = ", ".join(
                            dict.fromkeys(campaign_of(m) for m in manifests))
                        script = attributor.attribute(artifact.split("; ")[0],
                                                      manifests)
                    else:
                        old = previous.get(artifact.split("; ")[0], {})
                        script = old.get("script", "")
                        campaign = old.get("campaign", "")
                        manifest_cell = old.get("manifest", "")
                        if not script:
                            script = attributor.attribute(
                                artifact.split("; ")[0], [])
                sort_key = (0 if block["kind"] == "figure" else 1,
                            int(re.sub(r"\D", "", number) or 0))
                doc_rows.append((sort_key, {
                    "item": item, "document": document, "label": label,
                    "claim": demath(entry["caption"]), "artifact": artifact,
                    "script": script, "campaign": campaign,
                    "manifest": manifest_cell,
                    "defined_in": path.relative_to(manuscript_root).as_posix(),
                }))
        rows += [row for _key, row in sorted(doc_rows, key=lambda r: r[0])]

    for name in sorted(seals):
        if name in referenced or not name.endswith((".tex", ".pdf")):
            continue
        manifests = seals[name]
        rows.append({
            "item": "(not in manuscript)", "document": "none", "label": "",
            "claim": previous.get(name, {}).get("claim", ""),
            "artifact": name,
            "script": attributor.attribute(name, manifests),
            "campaign": ", ".join(
                dict.fromkeys(campaign_of(m) for m in manifests)),
            "manifest": "; ".join(f"metrics/{m}" for m in manifests),
            "defined_in": "",
        })
    return rows


def render(rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-root", type=Path,
                        default=ROOT.parent / "codebase",
                        help="checkout of the manuscript sources "
                             "(default: ../codebase)")
    parser.add_argument("--check", action="store_true",
                        help="verify without rewriting; non-zero on drift")
    args = parser.parse_args()

    for aux_name in ("main.aux", "supplement.aux"):
        if not (args.manuscript_root / aux_name).exists():
            print(f"[abort] {aux_name} not found under "
                  f"{args.manuscript_root}; compile the manuscript first")
            return 2

    rendered = render(build_rows(args.manuscript_root))
    if args.check:
        if OUT.read_text(encoding="utf-8") != rendered:
            print(f"[stale] {OUT.relative_to(ROOT)} no longer matches the "
                  f"manuscript; rerun without --check")
            return 1
        print(f"[ok] {OUT.relative_to(ROOT)} matches the manuscript")
        return 0
    OUT.write_text(rendered, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {rendered.count(chr(10)) - 1} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
