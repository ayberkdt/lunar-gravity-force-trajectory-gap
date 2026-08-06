"""Write the three reader-facing archive documents.

  campaign_map.csv          every claim -> experiment family, registration
                            status, source record, aggregate table, archive path
  archive_index.md          what is in the archive, directory by directory
  README_REPRODUCIBILITY.md how to verify a digest, regenerate a table, and
                            re-run a campaign

Everything that can be read from the records is read from them: the manifest
scopes, the registration statuses, the digests and the directory tree are not
retyped here. The claim rows are the final experiment matrix of Supplementary
Section S1.1, which is the manuscript's own statement of what carries what.

Usage:  python build_archive_docs.py
"""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
ARCHIVE = ROOT / "archive"

# The final experiment matrix of S1.1: claim, the family that carries it in this
# version, the evidential status it is claimed at, and the campaign that owns
# the records. Campaign attributions are the ones the reproducibility section
# states; none is inferred here.
CLAIMS = [
    # The calibration objective predates the R10 campaign series and is owned by
    # the legacy evidence package, indexed in r3_experiment_manifest.json.
    ("The exponent a compact tail rule needs is not the spectral slope",
     "F2", "confirmatory, field level", "R3", "spectrum_pfit|local_tail", "7.1"),
    ("The procedure transfers across products; the value does not",
     "F2", "descriptive, model relative", "R16", "multibody|transfer", "7.1"),
    ("Acceleration-level degree blending is non-conservative",
     "F5, F6", "confirmatory, field level", "R10, R11", "blend", "7.1"),
    ("Cheap empirical altitude-only schedules lose to fixed degrees over long arcs",
     "O19, O20, O27", "confirmatory", "R11, R17", "convergence|long", "7.1"),
    ("Given an accuracy target the published radial rule is accurate and expensive",
     "O23 resolved by O24", "confirmatory", "R12, R13",
     "atallah_campaign|atallah_verification", "7.1"),
    ("At a declared budget, force and trajectory disagree",
     "O25 against O26", "confirmatory", "R14", "pareto|variational", "7.2"),
    ("The disagreement is dynamical filtering, not numerical noise",
     "O38, gradient degree audited by O40",
     "mechanism; resolved-only scoring declared post hoc", "R37, R39",
     "variational|gradient", "7.2"),
    ("Within the tested family, an intermediate span beats both endpoints",
     "O28, O29, O35", "exploratory; budget dependent below beta 0.75",
     "R18, R19, R29", "span_sweep|equal_total_work", "7.3"),
    ("Where inside the factor box that allocation result holds",
     "O36", "locating, one stratum of five declared", "R30",
     "high_apolune|verdict", "S17"),
    ("Outside the factor box the endpoint comparison reverses",
     "O37 with ceiling control O39", "registered scope extension", "R31, R38",
     "uncapped|operational_elliptical", "S18"),
]

README = """# Reproducing the results

This archive is the complete machine-readable record behind the manuscript.
The Supplementary Information carries the definitions, the contracts and the
aggregate results; everything per-orbit, per-case and per-file is here.

## What to read first

| you want to | read |
|---|---|
| see which record carries which claim | `manifests/campaign_map.csv` |
| find a directory | `supplement_archive/archive_index.md` |
| find a table that left the supplement PDF | `supplement_archive/migration_map.csv` |
| check that nothing has been altered | `manifests/` and the procedure below |

## Verifying integrity

Every campaign owns one manifest, `metrics/rNN_final_experiment_manifest.json`,
which records a SHA-256 for each file it claims and a seal over itself. The
check that verifies the manifests against the files is
`python_codes/check_manifest_integrity.py`. It recomputes every recorded digest,
recomputes each manifest seal, confirms that no trajectory record is claimed by
two manifests, and confirms that every table and figure the documents pull in is
hashed somewhere. It exits non-zero on any failure.

    python python_codes/check_manifest_integrity.py

The trajectory partition is a property of the archive, not an assertion: a
record appears under exactly one manifest.

## Regenerating a table that left the supplement PDF

Twelve per-orbit and per-case tables were moved out of the supplement because no
text pointed at them and they cannot be read on a page. Each is archived three
ways: the rendered LaTeX the supplement compiled, a CSV rendering, and the
source record plus the generating script, both named and hashed in
`supplement_archive/migration_map.csv`.

To regenerate one, run its generating script from `python_codes/`. The scripts
read frozen records and propagate nothing.

The CSV is a convenience rendering of the typeset table. Where a column heading
was set in mathematics, the CSV carries a plain-text approximation of it; the
archived LaTeX and the source record are the authoritative forms.

## Re-running a campaign

Campaign drivers live in `python_codes/` and are named `revNN_*.py`. Each reads
its registration from `metrics/rNN_preregistration.json` where one exists,
refuses to run if the registration is missing, and writes only its own records.
Two environments are used and both ship as a specification and a full lock:
`environment-py310.yml` with `requirements-py310.lock` for R10-R13, and
`environment-py312.yml` with `requirements-py312.lock` for everything later.

Propagation campaigns are expensive. The gravity-model files are not
redistributed here; obtain each from its archive and verify the digest recorded
in the R16 manifest and in Supplementary Section S19 before use.

## What is deliberately not here

Raw trajectory arrays (about 2.1 GB) are regenerable from the drivers and are
excluded; their per-file digests are recorded in the campaign manifests. Large
external coefficient products and SPICE kernels are not redistributed.
"""


def manifest_records(d: dict) -> list[str]:
    """Every result record a manifest claims, whatever shape it stores them in.
    R10 keeps a flat `entries` inventory; later manifests use `result_json`."""
    names = list((d.get("result_json") or {}).keys())
    # the legacy package indexes its records under metric_provenance/artifacts
    names += list((d.get("metric_provenance") or {}).keys())
    for e in (d.get("entries") or []) + (d.get("artifacts") or []):
        n = e.get("path") or e.get("name") or ""
        # keep the repo-relative path: the basename alone does not resolve for
        # records that live in a campaign's case tree
        if n.endswith(".json") and "manifest" not in n:
            names.append(n)
    return names


def read_manifests() -> dict:
    out = {}
    for p in sorted(list(METRICS.glob("r*_final_experiment_manifest.json"))
                    + [METRICS / "r3_experiment_manifest.json"]):
        if not p.exists():
            continue
        name = p.name.split("_")[0].upper()
        d = json.loads(p.read_text(encoding="utf-8"))
        out[name] = {
            "scope": re.sub(r"\s+", " ", d.get("scope", d.get("purpose", ""))),
            "registration_status": re.sub(
                r"\s+", " ", d.get("registration_status", "not recorded")),
            "manifest_sha256": d.get("manifest_sha256", ""),
            "records": manifest_records(d),
            "file": p.name,
        }
    return out


def resolve_records(man: dict, campaigns: str, hint: str) -> str:
    """Name the records the claim rests on by looking them up in the campaigns'
    own manifests. Nothing is invented: if the hint matches nothing, the row
    says so and points at the manifest instead."""
    pat = re.compile(hint, re.I)
    hits = []
    for c in [x.strip() for x in campaigns.split(",")]:
        for rec in man.get(c, {}).get("records", []):
            if pat.search(rec) and rec not in hits:
                hits.append(rec)
    if not hits:
        return "see campaign manifest"
    return "; ".join(hits[:4]) + (f" (+{len(hits) - 4} more)"
                                  if len(hits) > 4 else "")


def write_campaign_map(man: dict) -> Path:
    out = ARCHIVE / "manifests" / "campaign_map.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["claim", "experiment_family", "evidential_status", "campaign",
              "registration_status", "source_record", "aggregate_table_section",
              "manifest_file", "manifest_sha256", "archive_path"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for claim, fam, status, camps, hint, section in CLAIMS:
            first = camps.split(",")[0].strip()
            m = man.get(first, {})
            w.writerow({
                "claim": claim,
                "experiment_family": fam,
                "evidential_status": status,
                "campaign": camps,
                "registration_status": m.get("registration_status",
                                             "see campaign manifest"),
                "source_record": resolve_records(man, camps, hint),
                "aggregate_table_section": section,
                "manifest_file": m.get("file", ""),
                "manifest_sha256": m.get("manifest_sha256", ""),
                "archive_path": f"manifests/{m.get('file', '')}",
            })
    return out


def write_index(man: dict) -> Path:
    lines = ["# Archive index", "",
             "What is in this archive, directory by directory. Generated from "
             "the tree itself by `build_archive_docs.py`; if a directory is "
             "listed as empty it is empty.", ""]
    lines += ["## Directories", "",
              "| path | files | contents |", "|---|---|---|"]
    described = {
        "environment": "interpreter specifications and dependency locks",
        "source/release_snapshot": "pinned source snapshot of the kernel",
        "gravity_models": "coefficient-product digests; the products themselves "
                          "are not redistributed",
        "registrations/original": "pre-registrations, one per campaign",
        "registrations/amendments": "amendments, each naming its parent",
        "designs": "frozen design coordinates and initial states",
        "designs/geometry_strata": "frozen sub-box designs, including the four "
                                   "that did not reach a ladder",
        "raw/field_validation": "field-level validation records",
        "raw/trajectory_validation": "trajectory-level qualification records",
        "raw/calibration": "calibration objective and transfer records",
        "raw/atallah_benchmark": "published-rule implementation and benchmark",
        "raw/fixed_budget": "fixed-budget allocation campaign",
        "raw/interpolation": "interpolation-family campaign",
        "raw/operational_elliptical": "operational elliptical population",
        "processed/aggregate_tables": "aggregate tables as they appear in the PDF",
        "processed/per_orbit_tables": "CSV renderings of the per-orbit tables "
                                      "that left the supplement",
        "processed/resolution_envelopes": "numerical-resolution envelopes",
        "processed/work_and_timing": "realized work and measured kernel time",
        "scripts/reproduce_tables": "table generators",
        "scripts/reproduce_figures": "figure generators",
        "scripts/validation": "integrity and consistency checks",
        "manifests": "campaign manifests, campaign map, checksums",
        "supplement_archive": "what left the supplement PDF, and why",
        "supplement_archive/full_legacy_tables": "the moved tables as rendered "
                                                 "LaTeX, byte-identical to what "
                                                 "the supplement compiled",
    }
    for rel, desc in described.items():
        d = ARCHIVE / rel
        n = len([p for p in d.glob("*") if p.is_file()]) if d.exists() else 0
        lines.append(f"| `{rel}/` | {n if n else 'empty'} | {desc} |")

    lines += ["", "## Campaign manifests", "",
              "| campaign | registration | scope |", "|---|---|---|"]
    for name, m in man.items():
        reg = m["registration_status"].split(".")[0][:70]
        lines.append(f"| {name} | {reg} | {m['scope'][:110]} |")

    mig = ARCHIVE / "supplement_archive" / "migration_map.csv"
    if mig.exists():
        rows = list(csv.DictReader(mig.open(encoding="utf-8")))
        lines += ["", "## Tables moved out of the supplement PDF", "",
                  f"{len(rows)} tables. Full record in "
                  "`supplement_archive/migration_map.csv`.", "",
                  "| table | rows | archived LaTeX | CSV |", "|---|---|---|---|"]
        for r in rows:
            lines.append(f"| {r['old_label']} | {r['rows_in_csv']} | "
                         f"`{r['archive_path_tex']}` | `{r['archive_path_csv']}` |")

    lines += ["", "## Not in this archive", "",
              "Raw trajectory arrays (about 2.1 GB), regenerable from the "
              "drivers and hashed in the campaign manifests. Gravity-model "
              "coefficient files and SPICE kernels, which are distributed by "
              "their own archives and are recorded here by digest only.", ""]
    out = ARCHIVE / "supplement_archive" / "archive_index.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> int:
    man = read_manifests()
    cm = write_campaign_map(man)
    ix = write_index(man)
    rd = ARCHIVE / "README_REPRODUCIBILITY.md"
    rd.write_text(README, encoding="utf-8")
    (ARCHIVE / "manifests" / "archive_docs_build.json").write_text(json.dumps({
        "schema": "archive_docs_build_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "claims": len(CLAIMS), "campaigns_indexed": len(man),
    }, indent=2), encoding="utf-8")
    for p in (cm, ix, rd):
        print(f"[written] {p.relative_to(ROOT)}  ({p.stat().st_size} bytes)")
    print(f"[campaign_map] {len(CLAIMS)} claims, {len(man)} campaigns indexed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
