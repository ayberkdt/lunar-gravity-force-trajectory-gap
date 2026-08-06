"""SHA-256 integrity manifest for the R23 registered controls (O31).

R23 collects the three controls applied to the constructive interior-member
claim after that claim was known, each registered before it was run:

  O31a  the R19 realized-work comparison repeated at beta = 0.5 on both designs
  O31b  the interior member scored against the constant-degree ladder of the
        16-orbit oracle panel at beta = 1
  O31c  the surviving beta = 1 comparison re-run at a third, ultra-tight
        tolerance level on the frozen 83-orbit panel

and one measurement that propagates nothing: the kernel cost curve measured
under both archived timing protocols at both archived model degrees in a single
session.

One naming exception is recorded rather than renamed away. O31a reuses the R19
driver with a budget argument added, so its outputs keep the r19 prefix and
carry a beta_0.50 suffix. Those files belong to R23 and are indexed here only;
the R19 manifest was sealed before they existed and its sidecar inventory
covers the beta = 1 record alone.

Usage:  python rev23_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r23_final_experiment_manifest.json"

SCRIPTS = ["rev23_preregister.py", "rev23_preregister_amendment.py",
           "rev19_equal_total_work.py", "rev23_oracle_vs_interior.py",
           "rev23_ultratight_span.py", "rev23_cost_curve_unified.py",
           "rev23_run_when_idle.py", "rev23_tables.py",
           "rev23_audit_claims.py", "rev23_finalize_manifest.py"]

REGISTRATION = ["r23_preregistration.json",
                "r23_preregistration_amendment.json"]

RESULT_JSON = ["r19_equal_total_work_A_beta_0.50.json",
               "r19_equal_total_work_B_beta_0.50.json",
               "r23_oracle_vs_interior.json",
               "r23_ultratight_span.json",
               "r23_ultra_panel.json",
               "r23_cost_curve_unified.json",
               "r23_cost_curve_unified_pilot.json",
               "r23_cost_curve_reproducibility.json",
               "r23_manuscript_descriptives.json"]

TABLES = ["r23_ultra_table.tex", "r23_oracle_table.tex",
          "r23_oracle_per_orbit_table.tex", "r23_cost_table.tex"]

REUSED = ["r18_span_sweep_A_beta_0.50.json", "r18_span_sweep_B_beta_0.50.json",
          "r18_span_sweep_A_beta_1.00.json", "r18_span_sweep_B_beta_1.00.json",
          "r14_trajectory_A_beta_0.50.json", "r14_trajectory_B_beta_0.50.json",
          "r19_equal_total_work_A.json", "r19_equal_total_work_B.json",
          "r15_fixed_oracle.json"]

# Sub-campaign -> the case/raw subtrees it owns, relative to metrics/.
TREES = {
    "O31a_realized_work_beta_0.50": [
        "r19_cases/A_workmatched_beta_0.50",
        "r19_cases/B_workmatched_beta_0.50",
        "r19_raw/A_workmatched_beta_0.50",
        "r19_raw/B_workmatched_beta_0.50"],
    "O31b_oracle_panel": [
        "r23_cases/oracle_vs_interior", "r23_raw/oracle_vs_interior"],
    "O31c_ultra_panel": [
        "r23_cases/ultra_span", "r23_raw/ultra_span"],
}


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


def index_trees() -> dict:
    """Per-sub-campaign sidecar hashes with a rolled-up digest of raw arrays.

    Sidecars are hashed individually because a reader may want to check one
    orbit; the raw state arrays are rolled up because their per-file hashes
    would dwarf the manifest without being separately actionable.
    """
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


def censoring() -> dict:
    """Read from the results rather than restated, so it cannot drift."""
    out = {}
    for design in ("A", "B"):
        p = METRICS / f"r19_equal_total_work_{design}_beta_0.50.json"
        if not p.exists():
            out[design] = {"missing": True}
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        rows = d.get("rows", [])
        out[design] = [r.get("sobol_index") for r in rows
                       if r.get("status") == "censored"]
    return out


def main() -> int:
    payload = {
        "schema": "r23_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R23 (O31): the three registered controls on the interior "
                  "span member -- the realized-work comparison repeated at "
                  "beta = 0.5 (O31a), the interior member against the "
                  "constant-degree ladder of the 16-orbit oracle panel "
                  "(O31b), and the beta = 1 realized-work comparison re-run "
                  "at a third ultra-tight tolerance level on the frozen "
                  "83-orbit panel (O31c) -- together with a single-session "
                  "re-measurement of the kernel cost curve under both "
                  "archived timing protocols, which propagates nothing."),
        "why": ("The interpolation family was not pre-registered and the "
                "reported member was chosen from the design medians after the "
                "sweep had run. These controls were registered before they "
                "were run so that the mapping from outcome to manuscript "
                "wording was fixed before the numbers existed."),
        "relationship_to_r18_r19": (
            "additive. O31a reuses the R19 protocol and driver at a different "
            "budget and propagates one comparator trajectory per orbit per "
            "tolerance level. O31b propagates the interior member and the "
            "declared ladder on the 16-orbit panel. O31c re-propagates the "
            "already-scored objects at one further tolerance level, reading "
            "their degree tables out of the archived R18 and R19 sidecars "
            "rather than rebuilding them. No R18 or R19 trajectory is "
            "re-propagated or re-indexed."),
        "naming_exception": (
            "O31a outputs carry an r19 prefix and a beta_0.50 suffix because "
            "the R19 driver was reused with a budget argument added, leaving "
            "the beta = 1 configurations and their hashes untouched. Those "
            "files belong to R23 and are indexed here only; the R19 manifest "
            "was sealed before they existed."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R22"},
        "censored_orbits_beta_0.50": censoring(),
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": index_trees(),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    for name, t in payload["trajectory_tree"].items():
        print(f"  {name}: {t['n_sidecars']} sidecars, "
              f"{t['n_raw_arrays']} raw arrays")
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "reused_inputs", "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
