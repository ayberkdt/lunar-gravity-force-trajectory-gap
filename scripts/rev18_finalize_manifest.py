"""SHA-256 integrity manifest for the R18 span sweep.

R18 propagates only the interior of the path between the two archived
endpoints, so the manifest records which inputs it reuses unchanged (the R14
tolerances and comparator degrees, the R11 truths) alongside what it produced.

Usage:  python rev18_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import campaign_ownership as own

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r18_final_experiment_manifest.json"

SCRIPTS = ["rev18_span_sweep.py", "rev18_tables.py",
           "rev18_finalize_manifest.py"]
# One record per (design, budget) pair. The sweep has never written an
# undecorated "r18_span_sweep_A.json"; that name was recorded as missing rather
# than failing, which hid all five real records from the manifest.
RESULT_JSON = ["r18_span_sweep_A_beta_0.50.json",
               "r18_span_sweep_A_beta_1.00.json",
               "r18_span_sweep_A_beta_1.50.json",
               "r18_span_sweep_B_beta_0.50.json",
               "r18_span_sweep_B_beta_1.00.json",
               "r18_manuscript_descriptives.json",
               "r18_budget_descriptives.json"]
TABLES = ["r18_span_table.tex", "r18_budget_table.tex",
          "r18_span_detail_table_A.tex", "r18_span_detail_table_B.tex"]
REUSED = ["r14_trajectory_A_beta_1.00.json",
          "r14_trajectory_B_beta_1.00.json", "r14_budget_pareto.json",
          "r10_sobolA_baseline_truth_corrected.json",
          "r11_designB_rows.json"]


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path) -> dict:
    out = {}
    for n in names:
        p = base / n
        out[n] = ({"sha256": sha(p), "bytes": p.stat().st_size}
                  if p.exists() else {"missing": True})
    return out


def index_tree() -> dict:
    """Index only the budgets this campaign propagated.

    Later campaigns reuse this driver with a budget argument, so their
    trajectories land under the same prefix; indexing them here would put the
    same records under two manifests. Ownership is declared in
    ``campaign_ownership``.
    """
    sidecars = {}
    for p in sorted((METRICS / "r18_cases").rglob("*.json")):
        if not own.owned_by_r18(p):
            continue
        sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    for p in sorted((METRICS / "r18_raw").rglob("*.npz")):
        if not own.owned_by_r18(p):
            continue
        roll.update(sha(p).encode())
        n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    payload = {
        "schema": "r18_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R18 span sweep: a one-parameter family interpolating "
                  "geometrically between the equal-budget constant degree and "
                  "the budget-calibrated radial rule, every member rescaled to "
                  "the same per-call budget, propagated at both vector "
                  "tolerance levels on both populations."),
        "relationship_to_r14": (
            "additive. The two endpoints are not re-propagated: k=0 and k=1 are "
            "read from the frozen R14 beta=1 record, and the per-orbit Atallah "
            "tolerance, comparator degree and truth trajectories are reused "
            "unchanged. R18 adds only the interior of the path."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R17"},
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": index_tree(),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    t = payload["trajectory_tree"]
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  {t['n_sidecars']} sidecars, {t['n_raw_arrays']} raw arrays")
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "reused_inputs")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
