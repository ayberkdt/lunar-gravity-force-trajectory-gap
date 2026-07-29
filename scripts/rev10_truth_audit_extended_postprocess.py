"""Fold the extended low-perilune N=900 audit into the truth-corrected aggregate.

The frozen primary audit (``r10_sobolA_truth_audit.json``, six cases) is the
only audit read by ``rev10_sobol_convergence`` and ``rev10_queue_postprocess``.
The reviewer-requested extension audits the four remaining sub-50 km survivors
(``r10_sobolA_truth_audit_extended.json``). This script merges both audits and
recomputes, with the *same* frozen functions used for the primary aggregate:

  * the 64-orbit truth-corrected baseline (adopted degree and N=900 ratios
    now taken from whichever audit covered each orbit);
  * the rho_work / rho_crit population medians, 10/90 percentiles, and raw-win
    counts, reported both primary-only and after the extension so the change is
    explicit;
  * the pre-specified selection bands (raw wins and 0.75--1.33 close
    comparisons) for the four extended orbits, to detect any new contested case
    that would require an additional convergence rerun.

Nothing is propagated here; only the already-computed N=900 policy differences
are re-aggregated. Frozen artifacts are not overwritten -- outputs carry the
``_extended`` suffix.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev10_sobol_convergence as conv
import rev10_queue_postprocess as qp


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
BASELINE = METRICS / "r10_sobolA_baseline.json"
PRIMARY_AUDIT = METRICS / "r10_sobolA_truth_audit.json"
EXTENDED_AUDIT = METRICS / "r10_sobolA_truth_audit_extended.json"
CORRECTED_OUT = METRICS / "r10_sobolA_baseline_truth_corrected_extended.json"
REPORT_OUT = METRICS / "r10_truth_audit_extended_aggregate.json"
TABLE_OUT = METRICS / "r10_truth_audit_extended_table.tex"

EXTENDED_INDICES = (11, 35, 43, 52)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def summarize(rows: list[dict], key: str) -> dict:
    return qp.ratio_summary(rows, key)


def main() -> int:
    baseline = load(BASELINE)
    primary = load(PRIMARY_AUDIT)
    extended = load(EXTENDED_AUDIT)
    if not extended.get("complete") or len(extended.get("rows", [])) != 4:
        raise RuntimeError("extended audit is incomplete")
    if not primary.get("complete") or len(primary.get("rows", [])) != 6:
        raise RuntimeError("primary audit is incomplete")

    merged_audit = {"rows": list(primary["rows"]) + list(extended["rows"])}

    # Same frozen builder used for the primary manuscript aggregate.
    corrected_primary = conv.corrected_rows(baseline, primary)
    corrected_merged = conv.corrected_rows(baseline, merged_audit)

    # Persist the extended corrected baseline (frozen file left untouched).
    corrected_payload = {
        "schema": "r10_truth_corrected_baseline_extended_v1",
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "baseline_sha256": base.file_hash(BASELINE),
        "primary_audit_sha256": base.file_hash(PRIMARY_AUDIT),
        "extended_audit_sha256": base.file_hash(EXTENDED_AUDIT),
        "rows": corrected_merged,
    }
    corrected_payload["artifact_sha256"] = base.object_hash(corrected_payload)
    base.atomic_json(CORRECTED_OUT, corrected_payload)

    # Selection bands (frozen thresholds) before and after the extension.
    parts_primary = conv.selection_parts(corrected_primary)
    parts_merged = conv.selection_parts(corrected_merged)
    by_index = {r["sobol_index"]: r for r in corrected_merged}
    primary_by_index = {r["sobol_index"]: r for r in corrected_primary}

    extended_detail = {}
    for idx in EXTENDED_INDICES:
        row = by_index[idx]
        prev = primary_by_index[idx]
        extended_detail[idx] = {
            "adopted_truth_degree": row["adopted_truth_degree"],
            "truth_source": row["truth_source"],
            "rho_work_N600": prev["rho_work"],
            "rho_crit_N600": prev["rho_crit"],
            "rho_work_audited": row["rho_work"],
            "rho_crit_audited": row["rho_crit"],
            "is_raw_win": idx in set(parts_merged["raw_win"]),
            "in_close_band": idx in set(parts_merged["close"]),
            "was_selected_primary": idx in set(parts_primary["selected"]),
            "is_selected_merged": idx in set(parts_merged["selected"]),
        }

    new_contested = sorted(
        idx for idx in EXTENDED_INDICES
        if (extended_detail[idx]["is_raw_win"]
            or extended_detail[idx]["in_close_band"])
        and not extended_detail[idx]["was_selected_primary"]
    )

    report = {
        "schema": "r10_truth_audit_extended_aggregate_v1",
        "created_utc": base.utc_now(),
        "protocol_sha256": base.protocol_payload()["protocol_sha256"],
        "baseline_sha256": base.file_hash(BASELINE),
        "primary_audit_sha256": base.file_hash(PRIMARY_AUDIT),
        "extended_audit_sha256": base.file_hash(EXTENDED_AUDIT),
        "corrected_extended_sha256": corrected_payload["artifact_sha256"],
        "extended_adopted_N900_indices": [
            idx for idx in EXTENDED_INDICES
            if by_index[idx]["adopted_truth_degree"] == 900
        ],
        "aggregate_primary_only": {
            "rho_work": summarize(corrected_primary, "rho_work"),
            "rho_crit": summarize(corrected_primary, "rho_crit"),
        },
        "aggregate_with_extension": {
            "rho_work": summarize(corrected_merged, "rho_work"),
            "rho_crit": summarize(corrected_merged, "rho_crit"),
        },
        "extended_case_detail": {str(k): v for k, v in extended_detail.items()},
        "new_contested_indices": new_contested,
        "requires_additional_convergence": bool(new_contested),
    }
    report["report_sha256"] = base.object_hash(report)
    base.atomic_json(REPORT_OUT, report)

    # Extended audit table, same column format as the primary audit table.
    table_lines = [
        r"\begin{tabular}{rrrrrr}",
        r"\toprule",
        r"Index & $h_p$ [km] & $N=600$--900 RMS [m] & Threshold [m] & "
        r"Pass & Adopted $N_T$ \\",
        r"\midrule",
    ]
    for row in sorted(extended["rows"], key=lambda item: item["sobol_index"]):
        table_lines.append(
            f'{row["sobol_index"]} & {row["hp_km"]:.1f} & '
            f'{row["N600_to_N900"]["pos_rms_m"]:.3f} & '
            f'{row["acceptance_threshold_m"]:.3f} & '
            f'{"yes" if row["passes"] else "no"} & '
            f'{row["adopted_truth_degree"]} \\\\'
        )
    table_lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    TABLE_OUT.write_text("\n".join(table_lines), encoding="utf-8")

    aw = report["aggregate_with_extension"]
    ap = report["aggregate_primary_only"]
    print(json.dumps({
        "extended_adopted_N900": report["extended_adopted_N900_indices"],
        "rho_work_median_primary": round(ap["rho_work"]["median"], 4),
        "rho_work_median_extended": round(aw["rho_work"]["median"], 4),
        "rho_crit_median_primary": round(ap["rho_crit"]["median"], 4),
        "rho_crit_median_extended": round(aw["rho_crit"]["median"], 4),
        "rho_work_raw_wins_primary": ap["rho_work"]["raw_schedule_wins"],
        "rho_work_raw_wins_extended": aw["rho_work"]["raw_schedule_wins"],
        "rho_crit_raw_wins_primary": ap["rho_crit"]["raw_schedule_wins"],
        "rho_crit_raw_wins_extended": aw["rho_crit"]["raw_schedule_wins"],
        "new_contested_indices": new_contested,
        "requires_additional_convergence": report["requires_additional_convergence"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
