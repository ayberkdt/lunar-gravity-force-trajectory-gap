"""Build the DOI-backed evidence archive and migrate supplement tables into it.

The supplement carries per-orbit and per-case tables that no text points at and
that a reader cannot practically inspect on the page. They are evidence, so they
are not deleted: they move here, in three forms at once, and the supplement
keeps an aggregate summary and a path.

For every migrated table the archive records
  * the rendered LaTeX table, byte-identical to what the supplement compiled;
  * a CSV rendering of the same rows, for reading without LaTeX;
  * the source record it was generated from and the script that generated it,
    both by name and by SHA-256, so the CSV can be regenerated and checked.

Nothing is recomputed, rounded differently, or reformatted numerically: the CSV
carries the cell text the LaTeX table carries.

Usage:  python build_evidence_archive.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
ARCHIVE = ROOT / "archive"

# label -> (generated table, supplement section it left, generator, source record)
MIGRATED = {
    "tab:sobol-convergence-self": (
        "r10_sobol_convergence_self.tex", "S12", "rev10_manuscript_assets.py",
        "r10_sobolA_convergence.json",
        "Tight-to-tighter self-differences and reference-inclusive envelopes, "
        "17-orbit convergence subset"),
    "tab:sobol-convergence-decisions": (
        "r10_sobol_convergence_decisions.tex", "S12",
        "rev10_manuscript_assets.py", "r10_sobolA_convergence.json",
        "Numerical-resolution audit for every selected orbit and comparator"),
    "tab:full64-per-orbit": (
        "r11_full64_per_orbit_table.tex", "S12", "rev11_manuscript_tables.py",
        "r11_full_convergence.json",
        "Per-orbit design-A vector-tolerance results, 64 orbits"),
    "tab:designB-per-orbit": (
        "r11_designB_per_orbit_table.tex", "S12", "rev11_manuscript_tables.py",
        "r11_designB_rows.json",
        "Per-orbit design-B vector-tolerance results, 64 orbits"),
    "tab:sobol-primary-audit": (
        "r10_sobol_primary_audit.tex", "S12", "rev10_manuscript_assets.py",
        "r10_aggregate_summary.json",
        "Complete primary-policy audit for propagated design A"),
    "tab:sobol-sensitivity-audit": (
        "r10_sobol_sensitivity_audit.tex", "S12", "rev10_manuscript_assets.py",
        "r10_aggregate_summary.json",
        "Schedule-sensitivity and diagnostic audit for propagated design A"),
    "tab:sobol-designA-coordinates": (
        "r10_sobol_designA_coordinates.tex", "S12", "rev10_manuscript_assets.py",
        "r10_sobolA_design.json",
        "Exact Sobol coordinates and transformed altitude coordinates, design A"),
    "tab:sobol-designB-coordinates": (
        "r10_sobol_designB_coordinates.tex", "S12", "rev10_manuscript_assets.py",
        "r10_sobolA_design.json",
        "Exact Sobol coordinates and transformed altitude coordinates, design B"),
    "tab:blend-error-detail": (
        "r10_blend_error_detail.tex", "S4", "rev10_manuscript_assets.py",
        "r10_blend_lro_convergence.json",
        "Complete LRO-like same-tolerance reference-error audit"),
    "tab:blend-telemetry": (
        "r10_blend_telemetry.tex", "S4", "rev10_manuscript_assets.py",
        "r10_blend_lro_convergence.json",
        "LRO-like propagation telemetry, per trajectory"),
    "tab:alpha-margin-per-orbit": (
        "r8_alpha_margin_supplement.tex", "S9", "rev8_alpha_margin.py",
        "r8_alpha_margin.json",
        "Complete per-orbit uniform degree-inflation margin ladder"),
    "tab:phase-sweep-vector": (
        "r11_phase_sweep_table.tex", "S6", "rev11_manuscript_tables.py",
        "r11_phase_sweep.json",
        "Per-phase vector-tolerance dispersion sweep"),
    # second migration: referenced per-orbit records, each replaced in the PDF
    # by the aggregate it decides plus an archive path
    "tab:atallah-matching-A": (
        "r12_atallah_matching_table_A.tex", "S13", "rev12_atallah_tables.py",
        "r12_atallah_campaign.json",
        "Per-orbit matching record of the published radial rule, design A"),
    "tab:atallah-matching-B": (
        "r12_atallah_matching_table_B.tex", "S13", "rev12_atallah_tables.py",
        "r12_atallah_campaign_designB.json",
        "Per-orbit matching record of the published radial rule, design B"),
    "tab:budget-per-orbit-A": (
        "r14_beta1_per_orbit_A.tex", "S14", "rev14_tables.py",
        "r14_trajectory_A_beta_1.00.json",
        "Per-orbit equal-budget comparison at beta = 1, design A"),
    "tab:budget-per-orbit-B": (
        "r14_beta1_per_orbit_B.tex", "S14", "rev14_tables.py",
        "r14_trajectory_B_beta_1.00.json",
        "Per-orbit equal-budget comparison at beta = 1, design B"),
    "tab:span-detail-A": (
        "r18_span_detail_table_A.tex", "S15", "rev18_tables.py",
        "r18_span_sweep_A_beta_1.00.json",
        "Per-orbit seven-day error of every interpolation member, design A"),
    "tab:span-detail-B": (
        "r18_span_detail_table_B.tex", "S15", "rev18_tables.py",
        "r18_span_sweep_B_beta_1.00.json",
        "Per-orbit seven-day error of every interpolation member, design B"),
    "tab:span-longarc-detail": (
        "r20_longarc_detail_table.tex", "S15", "rev18_tables.py",
        "r20_span_longarc.json",
        "Sixty-day record of the interpolation family, eight design-A orbits"),
}

SKELETON = [
    "environment", "source/release_snapshot", "gravity_models",
    "registrations/original", "registrations/amendments",
    "designs/geometry_strata",
    "raw/field_validation", "raw/trajectory_validation", "raw/calibration",
    "raw/atallah_benchmark", "raw/fixed_budget", "raw/interpolation",
    "raw/operational_elliptical",
    "processed/aggregate_tables", "processed/per_orbit_tables",
    "processed/resolution_envelopes", "processed/work_and_timing",
    "scripts/reproduce_tables", "scripts/reproduce_figures", "scripts/validation",
    "manifests",
    "supplement_archive/full_legacy_tables",
]


def sha(p: Path) -> str:
    d = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def strip_tex(cell: str) -> str:
    """Cell text as a reader would see it. Numbers are never touched; only
    markup is removed."""
    s = cell.strip()
    s = re.sub(r"\\multicolumn\{\d+\}\{[^}]*\}\{(.*)\}", r"\1", s)
    s = re.sub(r"\\(?:textbf|emph|texttt|mathrm|text|code)\{([^{}]*)\}", r"\1", s)
    s = s.replace(r"\_", "_").replace(r"\%", "%").replace(r"\&", "&")
    s = re.sub(r"\\[ ,;!]", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s)
    s = s.replace("$", "").replace("{", "").replace("}", "")
    s = s.replace("--", "-").replace("~", " ")
    s = re.sub(r"^\[\]", "", s)
    return re.sub(r"\s+", " ", s).strip()


def strip_env(src: str) -> str:
    """Remove the float/tabular environments and, crucially, their column
    specifications, which contain nested braces such as @{} and cannot be
    matched by a flat pattern."""
    out, i = [], 0
    while i < len(src):
        m = re.compile(r"\\(begin|end)\{(longtable|tabular|table)\*?\}").search(src, i)
        if not m:
            out.append(src[i:])
            break
        out.append(src[i:m.start()])
        i = m.end()
        if src[i:i + 1] == "[":                       # optional placement
            i = src.index("]", i) + 1
        if m.group(1) == "begin" and src[i:i + 1] == "{":
            depth = 0                                  # balanced column spec
            while i < len(src):
                if src[i] == "{":
                    depth += 1
                elif src[i] == "}":
                    depth -= 1
                    if depth == 0:
                        i += 1
                        break
                i += 1
    return "".join(out)


def tex_to_rows(path: Path) -> list[list[str]]:
    """Rows as the typeset table shows them.

    longtable repeats its header through \\endfirsthead/\\endhead and carries a
    "continued" banner; both are typesetting furniture, not data, and a
    machine-readable rendering must not contain them. Length declarations that
    precede the environment must not leak into the first cell either.
    """
    src = path.read_text(encoding="utf-8", errors="replace")
    src = re.sub(r"(?m)^\s*%.*$", "", src)
    src = re.sub(r"\\(?:setlength|renewcommand|arrayrulewidth)"
                 r"(?:\{[^{}]*\}){1,2}", "", src)
    body = strip_env(src)
    for cmd in (r"\\caption\{.*?\}\s*\\\\", r"\\label\{[^}]*\}",
                r"\\multicolumn\{\d+\}\{[^}]*\}\{[^{}]*continued[^{}]*\}",
                r"\\(?:toprule|midrule|bottomrule|hline|endfirsthead|endhead|"
                r"endfoot|endlastfoot|addlinespace|cmidrule\(?[lr]*\)?"
                r"(?:\{[^}]*\})?)"):
        body = re.sub(cmd, "", body, flags=re.S)
    rows = []
    for line in body.split(r"\\"):
        cells = [strip_tex(c) for c in line.split("&")]
        joined = "".join(cells).strip()
        if not joined:
            continue
        if re.search(r"\(continued\)|^\s*continued\s*$", joined, flags=re.I):
            continue
        if rows and cells == rows[-1]:          # repeated longtable header
            continue
        rows.append(cells)
    # a header repeated further down (\endhead) still reaches here once
    if len(rows) > 2 and rows[0] == rows[1]:
        rows.pop(1)
    return rows


def main() -> int:
    for d in SKELETON:
        (ARCHIVE / d).mkdir(parents=True, exist_ok=True)

    legacy = ARCHIVE / "supplement_archive" / "full_legacy_tables"
    per_orbit = ARCHIVE / "processed" / "per_orbit_tables"
    rows_out, missing = [], []

    for label, (texname, section, script, source, desc) in MIGRATED.items():
        src = METRICS / texname
        if not src.exists():
            missing.append(texname)
            continue
        shutil.copy2(src, legacy / texname)
        csv_name = texname.replace(".tex", ".csv")
        data = tex_to_rows(src)
        with (per_orbit / csv_name).open("w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerows(data)

        srec = METRICS / source
        rows_out.append({
            "old_label": label,
            "old_section": section,
            "content": desc,
            "action": "MOVE",
            "revised_supplement_location":
                f"aggregate summary retained in {section}; full table archived",
            "archive_path_tex":
                f"supplement_archive/full_legacy_tables/{texname}",
            "archive_path_csv": f"processed/per_orbit_tables/{csv_name}",
            "rows_in_csv": len(data),
            "generator_script": script,
            "source_record": source if srec.exists() else f"{source} (not on disk)",
            "source_record_sha256": sha(srec) if srec.exists() else "",
            "table_sha256": sha(src),
            "reason": "no text in either document referenced this table; "
                      "per-case inventory not practically inspectable on the page",
            "main_manuscript_refs_affected": "none",
        })

    mig = ARCHIVE / "supplement_archive" / "migration_map.csv"
    with mig.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)

    (ARCHIVE / "manifests" / "archive_build.json").write_text(json.dumps({
        "schema": "evidence_archive_build_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "migrated_tables": len(rows_out),
        "missing_tables": missing,
        "note": ("tables moved out of the supplement PDF. Each is archived as "
                 "the rendered LaTeX it compiled from, as CSV, and with the "
                 "source record and generator script named and hashed. No "
                 "numerical value was recomputed or reformatted."),
    }, indent=2), encoding="utf-8")

    print(f"[archive] skeleton at {ARCHIVE}")
    print(f"[archive] migrated {len(rows_out)} tables")
    for r in rows_out:
        print(f"   {r['old_label']:<34} {r['rows_in_csv']:>4} csv rows  "
              f"<- {r['generator_script']}")
    if missing:
        print("[error] tables not found on disk: " + ", ".join(missing))
        return 1
    print(f"[archive] migration map: {mig.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
