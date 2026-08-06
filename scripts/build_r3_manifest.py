"""Build a hash manifest for the manuscript evidence package."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metrics" / "r3_experiment_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path, role: str) -> dict[str, object]:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "role": role,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def collect() -> list[dict[str, object]]:
    groups = [
        (sorted((ROOT / "metrics").glob("r[1234]_*.json")), "raw metric"),
        (sorted((ROOT / "metrics").glob("r8_*.json")), "raw metric"),
        (sorted((ROOT / "metrics").glob("r8_*.tex")), "generated table fragment"),
        (sorted((ROOT / "metrics").glob("external_*.json")), "external-validation evidence"),
        (sorted((ROOT / "metrics").glob("supplemental_*.json")), "supplemental control metric"),
        (sorted((ROOT / "python_codes").glob("rev[234]_*.py")), "experiment script"),
        (sorted((ROOT / "python_codes").glob("rev8_*.py")), "experiment script"),
        (sorted((ROOT / "python_codes").glob("make_figures_r[1238].py")), "figure script"),
        (sorted((ROOT / "python_codes").glob("supplemental_*.py")), "supplemental control script"),
        ([ROOT / "python_codes" / "make_figures_supplemental.py"], "figure script"),
        ([ROOT / "python_codes" / "paper_style.py"], "figure style"),
        (sorted((ROOT / "figures").glob("*.pdf")), "figure"),
        (sorted((ROOT / "chapters").glob("*.tex")), "manuscript source"),
        ([ROOT / "main.tex", ROOT / "preamble.tex", ROOT / "references.bib"], "manuscript source"),
        ([ROOT / "main.pdf"], "compiled manuscript"),
    ]
    records: list[dict[str, object]] = []
    seen: set[Path] = set()
    for paths, role in groups:
        for path in paths:
            path = path.resolve()
            if path == OUTPUT.resolve() or path in seen or not path.is_file():
                continue
            seen.add(path)
            records.append(artifact(path, role))
    return sorted(records, key=lambda item: str(item["path"]))


def metric_commits() -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    paths = list((ROOT / "metrics").glob("r[34]_*.json"))
    paths += list((ROOT / "metrics").glob("r8_*.json"))
    paths += list((ROOT / "metrics").glob("supplemental_*.json"))
    for path in sorted(paths):
        if path.resolve() == OUTPUT.resolve():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        records[path.name] = {
            "repo_commit_sha": payload.get("repo_commit_sha"),
            "repo_working_tree_clean": payload.get("repo_working_tree_clean"),
        }
    return records


def main() -> None:
    artifacts = collect()
    manifest = {
        "schema": "lunaris.manuscript-evidence-manifest.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Metrics, scripts, figures, and manuscript sources present in this package.",
        "source_repository": {
            "path_at_generation": "D:/Masaustu/LUNAR_SIMULATION",
            "archived_tag": "paper-truncation-v1.0",
            "archived_tag_commit_sha": "48dd98d9fae4aa2c7cbfa26262099444a44d86c7",
            "computational_evidence_commit_sha": "c63de18580bbcc108f8825c4c271cbbfeae10123",
            "release_status": "A new immutable release tag/DOI for the complete evidence set is still required.",
        },
        "metric_provenance": metric_commits(),
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
    }
    OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
