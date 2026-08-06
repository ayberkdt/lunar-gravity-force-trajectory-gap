"""SHA-256 integrity manifest for the R24 registered controls (O32).

R24 is two controls, each closing something the previous round left open:

  O32a  the (O31b) oracle panel re-run at the third tolerance level, because
        eight of its sixteen comparisons were undecided and the record showed
        the interior member's own envelope, not its error, was what withheld
        the verdict
  O32b  the bin-resolution control of (O23) applied at the budget-calibrated
        point, where the constructive claim lives and where the response to
        review had narrowed the control's declared scope rather than run it

Neither campaign recomputes an error that is already archived. O32a reuses the
(O31b) tighter-level errors verbatim and rebuilds only the envelopes, so no
quantity already reported can move; O32b holds every parameter of the interior
member at its archived value so that quantization is the only difference.

Usage:  python rev24_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r24_final_experiment_manifest.json"

SCRIPTS = ["rev24_oracle_ultra.py", "rev24_bin_control.py",
           "rev24_tables.py", "rev24_finalize_manifest.py"]

REGISTRATION = ["r24_preregistration.json",
                "r24_bin_control_preregistration.json"]

RESULT_JSON = ["r24_oracle_ultra.json", "r24_bin_control.json",
               "r24_manuscript_descriptives.json"]

TABLES = ["r24_oracle_ultra_table.tex", "r24_bin_control_table.tex"]

# Inputs read but not produced here. Hashed so that a later change to any of
# them is visible as a change to what R24 was computed from.
REUSED = ["r23_oracle_vs_interior.json",
          "r23_ultratight_span.json",
          "r19_equal_total_work_A.json",
          "r19_equal_total_work_B.json",
          "r18_span_sweep_A_beta_1.00.json",
          "r18_span_sweep_B_beta_1.00.json",
          "r15_fixed_oracle.json"]

TREES = {
    "O32a_oracle_panel_ultra": [
        "r24_cases/oracle_ultra", "r24_raw/oracle_ultra"],
    "O32b_bin_control_budget_calibrated": [
        "r24_cases/bin_control", "r24_raw/bin_control"],
}


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path) -> dict:
    """A named file that is not on disk is an error, not a footnote.

    The earlier generators wrote {"missing": true} and exited zero, which is
    how a manifest came to certify records that were not there. Here the caller
    is told and the exit code carries it.
    """
    out, absent = {}, []
    for n in names:
        p = base / n
        if p.exists():
            out[n] = {"sha256": sha(p), "bytes": p.stat().st_size}
        else:
            out[n] = {"missing": True}
            absent.append(n)
    return out, absent


def index_trees() -> dict:
    out = {}
    for name, rels in TREES.items():
        sidecars, roll, n_raw = {}, hashlib.sha256(), 0
        for rel in rels:
            base = METRICS / rel
            if not base.exists():
                continue
            for p in sorted(base.rglob("*.json")):
                key = str(p.relative_to(METRICS)).replace("\\", "/")
                sidecars[key] = sha(p)
            for p in sorted(base.rglob("*.npz")):
                roll.update(sha(p).encode())
                n_raw += 1
        out[name] = {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
                     "sidecar_sha256": sidecars,
                     "raw_rollup_sha256": roll.hexdigest()}
    return out


def completeness() -> dict:
    """Read from the records rather than restated, so it cannot drift."""
    out = {}
    for key, fname in (("O32a", "r24_oracle_ultra.json"),
                       ("O32b", "r24_bin_control.json")):
        p = METRICS / fname
        if not p.exists():
            out[key] = {"missing": True}
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        out[key] = d.get("panel_completeness", {})
    return out


def main() -> int:
    scripts, a1 = index_files(SCRIPTS, CODE)
    registration, a2 = index_files(REGISTRATION, METRICS)
    results, a3 = index_files(RESULT_JSON, METRICS)
    tables, a4 = index_files(TABLES, METRICS)
    reused, a5 = index_files(REUSED, METRICS)
    absent = a1 + a2 + a3 + a4 + a5

    payload = {
        "schema": "r24_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R24 (O32): two registered controls. O32a re-runs the O31b "
                  "oracle panel at the third tolerance level, reusing the "
                  "archived tighter-level errors and rebuilding only the "
                  "envelopes. O32b applies the O23 bin-resolution control at "
                  "the budget-calibrated point, holding every parameter of the "
                  "interior member at its archived value so that only "
                  "quantization differs."),
        "why": ("O31b left eight of sixteen comparisons undecided, and in six "
                "of those the interior member held the smaller raw error while "
                "its own envelope withheld the verdict; O23's bin control had "
                "been run only at the accuracy-targeted point, and the "
                "response to review narrowed its declared scope instead of "
                "running it where the constructive claim lives."),
        "relationship_to_r23": (
            "additive. No archived error is recomputed and no R23 record is "
            "superseded: O32a reuses the R23 tighter-level errors verbatim, so "
            "the median ratios R23 reports are unchanged by construction and "
            "only the resolution of each comparison can move."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "result_json": results,
        "generated_tables": tables,
        "panel_completeness": completeness(),
        "trajectory_tree": index_trees(),
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")
                   ).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in payload["trajectory_tree"].values())
    n_raw = sum(t["n_raw_arrays"] for t in payload["trajectory_tree"].values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}  "
          f"{n_side} sidecars, {n_raw} raw arrays")
    if absent:
        print("[error] named files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
