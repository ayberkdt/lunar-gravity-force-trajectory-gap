"""Build the final R10 evidence-package integrity manifest.

This script hashes existing evidence and manuscript files. It does not run
propagations or alter scientific outputs.

R10 carries the manuscript-state half of the package: the sources, the two
PDFs, and the assets the document pulls in that belong to no later campaign.
That last group is derived from the manuscript rather than listed by hand,
because listing it by hand is what left the figures and the surviving pre-R10
tables hashed by nothing at all while the supplement claimed the manifests were
jointly exhaustive. Generated tables whose names carry a later campaign's
prefix stay with that campaign and are not duplicated here. Figures and
figure-generation code are indexed here even where a campaign also lists its
own copy: the non-duplication convention the supplement states is about
trajectories, and a file hashed under two manifests has to satisfy both.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metrics" / "r10_final_experiment_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_formal_r10(path: Path) -> bool:
    lowered = path.as_posix().lower()
    return (
        "smoke" not in lowered
        and "active" not in lowered
        and not lowered.endswith(".invalid")
        and path != OUTPUT
    )


selected: set[Path] = set()

for pattern in (
    "metrics/r10*.json",
    "metrics/r10*.tex",
    "metrics/r10_cases/**/*.json",
    "metrics/r10_raw/**/*.npz",
    "python_codes/rev10*.py",
    "python_codes/run_r10*.ps1",
):
    selected.update(path for path in ROOT.glob(pattern) if is_formal_r10(path))

for fixed in (
    "experiments/protocols/sobol_confirmatory_protocol.json",
    "figures/fig_sobol_confirmatory.pdf",
    "metrics/r3_experiment_manifest.json",
    "main.tex",
    "supplement.tex",
    "preamble.tex",
    "references.bib",
    "main.pdf",
    "supplement.pdf",
    "revision/CLAIM_EVIDENCE_MATRIX.md",
    "revision/CHANGE_LIST.md",
    "revision/CODE_CORRECTNESS_AUDIT_R10.md",
    "revision/FINAL_SUBMISSION_REPORT_R10.md",
):
    path = ROOT / fixed
    if path.exists():
        selected.add(path)

selected.update(ROOT.glob("chapters/*.tex"))

# Manuscript assets no later campaign owns. The \input targets are filtered by
# prefix so an r14_* or r23_* table stays under its own manifest; the figures
# and the figure-generation code have no campaign home at all, and the
# supplement's source-snapshot statement promises both.
OWNED_BY_LATER_CAMPAIGN = re.compile(r"^r(1[0-9]|2[0-9])_")

_tex = [ROOT / "main.tex", ROOT / "supplement.tex"] + \
    sorted(ROOT.glob("chapters/*.tex"))
_body = "\n".join(p.read_text(encoding="utf-8") for p in _tex if p.is_file())

for target in re.findall(r"\\input\{(metrics/[^}]+)\}", _body):
    stem = Path(target.strip()).name
    if OWNED_BY_LATER_CAMPAIGN.match(stem):
        continue
    for candidate in (ROOT / target.strip(), ROOT / (target.strip() + ".tex")):
        if candidate.is_file():
            selected.add(candidate)
            break

for target in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", _body):
    t = target.strip()
    for candidate in (ROOT / t, ROOT / f"{t}.pdf", ROOT / f"{t}.png",
                      ROOT / "figures" / t, ROOT / "figures" / f"{t}.pdf",
                      ROOT / "figures" / f"{t}.png"):
        if candidate.is_file():
            selected.add(candidate)
            break

selected.update(ROOT.glob("python_codes/make_figures*.py"))


def category(path: Path) -> str:
    relative = path.relative_to(ROOT).as_posix()
    if relative.startswith("metrics/r10_raw/"):
        return "raw_trajectory"
    if relative.startswith("metrics/r10_cases/"):
        return "case_sidecar"
    if relative.startswith("metrics/"):
        return "metric_or_generated_table"
    if relative.startswith("python_codes/"):
        return "script"
    if relative.startswith("figures/"):
        return "figure"
    if relative.startswith("revision/"):
        return "report"
    if relative.startswith("experiments/"):
        return "protocol"
    return "manuscript"


entries = []
for path in sorted(selected, key=lambda item: item.relative_to(ROOT).as_posix()):
    if not path.is_file():
        continue
    entries.append(
        {
            "path": path.relative_to(ROOT).as_posix(),
            "category": category(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    )

category_counts: dict[str, int] = {}
for entry in entries:
    category_counts[entry["category"]] = category_counts.get(entry["category"], 0) + 1

legacy = ROOT / "metrics" / "r3_experiment_manifest.json"
payload = {
    "schema": "r10_final_experiment_manifest_v1",
    "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "purpose": (
        "SHA-256 integrity index for the formal R10 confirmatory evidence and "
        "the manuscript state that reports it"
    ),
    "protocol_sha256": (
        "ebed6ac042e85fd1399b931bc849fb78a9f167cfc343f2165dde2e254d32e691"
    ),
    "lunaris_release": {
        "tag": "paper-truncation-v1.0",
        "commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        "url": "https://github.com/ayberkdt/lunaris/releases/tag/paper-truncation-v1.0",
        "is_doi_archive": False,
    },
    "legacy_manifest": {
        "path": "metrics/r3_experiment_manifest.json",
        "sha256": sha256(legacy) if legacy.exists() else None,
    },
    "selection": {
        "included": (
            "formal R10 summaries, case sidecars, raw trajectories, generated "
            "tables, scripts, protocol, confirmatory figure, manuscript "
            "sources and PDFs, and final reporting documents"
        ),
        "excluded": (
            "smoke-test, active/checkpoint, and invalid artifacts; these do "
            "not contribute to manuscript claims"
        ),
    },
    "file_count": len(entries),
    "category_counts": category_counts,
    "total_bytes": sum(entry["bytes"] for entry in entries),
    "entries": entries,
}

OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(
    f"[written] {OUTPUT.relative_to(ROOT)}: "
    f"{len(entries)} files, {payload['total_bytes']} bytes"
)
