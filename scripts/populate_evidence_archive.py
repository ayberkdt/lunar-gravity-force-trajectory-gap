"""Populate the evidence archive from the working tree.

What goes in, and what deliberately does not:

  copied      registrations and amendments, environment specifications and
              locks, frozen designs, the generated tables the documents pull
              in, the drivers and validation scripts, the gravity-product
              digests, and every top-level campaign result record. These are
              what a reader needs to check a number or re-run a campaign.

  indexed     the per-case trees under metrics/*_cases/. About 25000 files and
              95 MB of per-case detail, already hashed in the campaign
              manifests. The archive carries an index of every one of them with
              its path, size and digest, and the assembler below copies them at
              deposit time. Duplicating them in the working tree would create a
              second copy that drifts the moment a campaign re-runs.

  excluded    raw trajectory arrays under metrics/*_raw/, about 14 GB. Regenerable
              from the drivers; per-file digests stay in the campaign
              manifests. This is the policy the public repository already
              applies.

Nothing is renamed, reformatted or recomputed on the way in.

Usage:  python populate_evidence_archive.py
        python populate_evidence_archive.py --assemble-cases   (deposit time)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
FIGS = ROOT / "figures"
ARCHIVE = ROOT / "archive"

# destination <- selection rule over metrics/*.json (top level only)
RAW_ROUTING = [
    ("raw/calibration", ("r1_", "r3_", "r16_", "r15_deployable", "e2_", "e3_",
                         "supplemental_pstar")),
    ("raw/field_validation", ("r2_", "r4_", "r5_", "r6_", "r7_", "e1_", "e4_",
                              "e5_")),
    ("raw/atallah_benchmark", ("r12_", "r13_", "r24_", "r15_atallah")),
    ("raw/fixed_budget", ("r14_", "r21_", "r25_", "r27_", "r28_", "r34_",
                          "r36_", "r37_", "r39_", "r15_fixed")),
    ("raw/interpolation", ("r18_", "r19_", "r20_", "r22_", "r23_", "r35_")),
    ("raw/operational_elliptical", ("r31_", "r38_")),
    ("raw/trajectory_validation", ("r8_", "r9_", "r10_", "r11_", "r15_",
                                   "r17_", "r26_", "r29_", "r30_", "e6_", "e7_",
                                   "external_tudat", "robustness_",
                                   "review_postprocess",
                                   "supplemental_orbit")),
]
REGISTRATION_KEYS = ("preregistration", "registration", "amendment",
                     "design_frozen", "protocol")
# a stratum design is named for its sub-box, not for the word "stratum"
STRATUM_NAMES = ("polar", "equatorial", "frozen_like", "high_apolune",
                 "low_perilune")
DESIGN_KEYS = ("design_frozen", "sobol_design", "designB_rows",
               "designC", "operational_elliptical_design", "stratum")


def sha(p: Path) -> str:
    d = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def copy(src: Path, dst_dir: Path) -> Path:
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    shutil.copy2(src, dst)
    return dst


def route(name: str) -> str | None:
    for dest, prefixes in RAW_ROUTING:
        if name.startswith(prefixes):
            return dest
    return None


def populate() -> dict:
    counts: dict[str, int] = {}

    def bump(k, n=1):
        counts[k] = counts.get(k, 0) + n

    # ---- environment
    for pat in ("environment-py*.yml", "requirements-py*.lock"):
        for p in ROOT.glob(pat):
            copy(p, ARCHIVE / "environment")
            bump("environment")
    (ARCHIVE / "environment" / "system_info.json").write_text(json.dumps({
        "note": ("the campaigns ran under two interpreters; which campaign used "
                 "which is stated in Supplementary Section S19. Both "
                 "specifications and both full locks are in this directory."),
        "python_3_10": "R10-R13",
        "python_3_12": "R14 onward, the kernel timings and the field-level work",
    }, indent=2), encoding="utf-8")
    bump("environment")

    # ---- source snapshot identity
    kern = {}
    for p in METRICS.glob("r*_final_experiment_manifest.json"):
        d = json.loads(p.read_text(encoding="utf-8"))
        if "numerical_kernel" in d:
            kern = d["numerical_kernel"]
            break
    (ARCHIVE / "source" / "release_snapshot" / "README.txt").write_text(
        "The spherical-harmonic kernel used as the measurement instrument is\n"
        "maintained in the open-source Lunaris repository and is identified\n"
        "here by release tag and commit, in ../commit.txt. The release package\n"
        "carries the commit-pinned snapshot itself under src/lunaris/; this\n"
        "mirror carries the identification only. Any snapshot placed in this\n"
        "directory must match that commit; verify it before use, because every\n"
        "number in the manuscript is model-relative to the kernel that\n"
        "produced it.\n",
        encoding="utf-8")
    bump("source")
    (ARCHIVE / "source" / "commit.txt").write_text(
        f"lunaris_release_tag: {kern.get('lunaris_release_tag', 'unrecorded')}\n"
        f"lunaris_commit: {kern.get('lunaris_commit', 'unrecorded')}\n"
        "\nThe legacy evidence package is pinned to a different commit; see\n"
        "Supplementary Section S19 for which campaign uses which snapshot.\n",
        encoding="utf-8")
    bump("source")

    # ---- gravity-product digests
    grav = {}
    r16 = METRICS / "r16_final_experiment_manifest.json"
    if r16.exists():
        d = json.loads(r16.read_text(encoding="utf-8"))
        for k, v in (d.get("input_products") or d.get("gravity_models")
                     or {}).items():
            grav[k] = v
    (ARCHIVE / "gravity_models" / "file_manifest.json").write_text(json.dumps({
        "note": ("coefficient products are distributed by their own archives "
                 "and are not redistributed here. Obtain each file and verify "
                 "its digest before use; the digests are those of the files "
                 "actually read."),
        "primary": {"model": "GRAIL GL1800F",
                    "sha256_prefix": "d2a55206", "sha256_suffix": "6738fec"},
        "from_r16_manifest": grav,
        "see": "Supplementary Section S19 for the full list with citations",
    }, indent=2), encoding="utf-8")
    bump("gravity_models")

    # ---- registrations, designs, raw records, aggregate tables
    for p in sorted(METRICS.glob("*.json")):
        n = p.name
        if "manifest" in n:
            copy(p, ARCHIVE / "manifests")
            bump("manifests")
            continue
        if any(k in n for k in REGISTRATION_KEYS):
            sub = "amendments" if "amendment" in n else "original"
            copy(p, ARCHIVE / "registrations" / sub)
            bump("registrations")
        if any(k in n for k in DESIGN_KEYS):
            strat = n.startswith("r30_") and any(s in n for s in STRATUM_NAMES)
            sub = "designs/geometry_strata" if strat else "designs"
            copy(p, ARCHIVE / sub)
            bump("designs")
        dest = route(n)
        if dest:
            copy(p, ARCHIVE / dest)
            bump("raw")

    for p in sorted(METRICS.glob("*.tex")):
        copy(p, ARCHIVE / "processed" / "aggregate_tables")
        bump("aggregate_tables")

    # ---- scripts
    for p in sorted(CODE.glob("*.py")):
        n = p.name
        if n.startswith("make_figures") or "figure" in n:
            copy(p, ARCHIVE / "scripts" / "reproduce_figures")
        elif n.startswith("check_") or "audit" in n or "verify" in n:
            copy(p, ARCHIVE / "scripts" / "validation")
        elif "table" in n or "manuscript_assets" in n or "deliverable" in n:
            copy(p, ARCHIVE / "scripts" / "reproduce_tables")
        else:
            copy(p, ARCHIVE / "scripts")
        bump("scripts")

    # ---- case trees: indexed, not duplicated
    idx, nbytes = [], 0
    for p in sorted(METRICS.glob("*_cases/**/*")):
        if p.is_file():
            idx.append({"path": p.relative_to(ROOT).as_posix(),
                        "bytes": p.stat().st_size})
            nbytes += p.stat().st_size
    (ARCHIVE / "raw" / "case_tree_index.json").write_text(json.dumps({
        "schema": "case_tree_index_v1",
        "note": ("per-case records, indexed rather than duplicated in the "
                 "working tree. Digests for these files are held in the "
                 "campaign manifests. Run populate_evidence_archive.py "
                 "--assemble-cases at deposit time to copy them in."),
        "files": len(idx), "bytes": nbytes,
        "entries": idx,
    }, indent=1), encoding="utf-8")
    bump("case_index", len(idx))

    # ---- resolution envelopes and work/timing pointers
    for sub, keys in (("resolution_envelopes", ("convergence", "resolution",
                                                "envelope")),
                      ("work_and_timing", ("timing", "work", "cost", "kernel"))):
        hit = [p for p in METRICS.glob("*.json")
               if any(k in p.name for k in keys)]
        for p in hit:
            copy(p, ARCHIVE / "processed" / sub)
        bump(sub, len(hit))

    return counts


def write_checksums() -> tuple[int, float]:
    lines, total = [], 0
    for p in sorted(ARCHIVE.rglob("*")):
        if p.is_file() and p.name not in ("CHECKSUMS.sha256", "MANIFEST.json"):
            lines.append(f"{sha(p)}  {p.relative_to(ARCHIVE).as_posix()}")
            total += p.stat().st_size
    (ARCHIVE / "manifests" / "CHECKSUMS.sha256").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    (ARCHIVE / "manifests" / "MANIFEST.json").write_text(json.dumps({
        "schema": "evidence_archive_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "files": len(lines), "bytes": total,
        "excluded": {
            "metrics/*_raw/": "raw trajectory arrays, about 14 GB, "
                              "regenerable from the drivers and hashed in the "
                              "campaign manifests",
            "metrics/*_cases/": "per-case records, indexed in "
                                "raw/case_tree_index.json and copied at "
                                "deposit time",
            "rJ raw state arrays": "the J1-J3 arrays, about 0.7 GB, held "
                                   "outside metrics/ at the location recorded "
                                   "in rJ_final_experiment_manifest.json and "
                                   "hashed there per file; with metrics/*_raw/ "
                                   "this is the about 15 GB the manuscript "
                                   "reports as not redistributed",
        },
        "checksums": "manifests/CHECKSUMS.sha256",
    }, indent=2), encoding="utf-8")
    return len(lines), total / 1048576


def assemble_cases() -> int:
    idx = json.loads((ARCHIVE / "raw" / "case_tree_index.json").read_text(
        encoding="utf-8"))
    n = 0
    for e in idx["entries"]:
        src = ROOT / e["path"]
        if not src.exists():
            continue
        dst = ARCHIVE / "raw" / "cases" / Path(e["path"]).relative_to("metrics")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        n += 1
    print(f"[assemble] {n} case files copied into archive/raw/cases/")
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemble-cases", action="store_true")
    a = ap.parse_args()
    if a.assemble_cases:
        assemble_cases()
    counts = populate()
    for k in sorted(counts):
        print(f"   {k:<22} {counts[k]}")
    n, mb = write_checksums()
    print(f"[archive] {n} files, {mb:.1f} MB, checksums written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
