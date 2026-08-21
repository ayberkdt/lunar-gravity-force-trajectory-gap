"""SHA-256 integrity manifest for the R44 (O42) tighter-level equal-work
rematch.

R44 changes only the level the realized-work match is made at. It reuses the
R18 span members and their tighter-level telemetry unchanged and propagates
one new constant-degree trajectory per orbit per tolerance level, matched on
realized total quadratic work at the tighter level, in eight registered
design--budget cells.

Usage:  python rev44_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r44_final_experiment_manifest.json"

CELLS = [("A", "0.50"), ("A", "0.75"), ("A", "1.00"), ("A", "1.25"),
         ("B", "0.50"), ("B", "0.75"), ("B", "1.00"), ("B", "1.50")]

SCRIPTS = ["rev44_equal_work_tighter.py", "rev44_tables.py",
           "rev44_finalize_manifest.py", "run_r44_overnight.ps1",
           "run_overnight_20260809.ps1"]
RESULT_JSON = ([f"r44_equal_work_tighter_{d}_beta_{b}.json" for d, b in CELLS]
               + ["r44_manuscript_descriptives.json"])
REGISTRATION = ["r44_preregistration.json"]
TABLES = ["r44_equal_work_table.tex"]
REUSED = [f"r18_span_sweep_{d}_beta_{b}.json" for d, b in CELLS]


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
    """R44 owns its whole tree: every case directory carries the
    workmatched_tighter tag, so no suffix partition is needed."""
    sidecars = {}
    root = METRICS / "r44_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r44_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    payload = {
        "schema": "r44_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R44 (O42): the interior span member (k = 0.5) against a "
                  "constant degree matched on realized total quadratic work "
                  "at the tighter tolerance, the level every error in the "
                  "comparison is scored at, in eight registered design-"
                  "budget cells."),
        "why": ("Council review flagged the R19 match as level-inconsistent: "
                "cost matched at the tight level, errors read at the tighter "
                "level, where the R19 comparator underspends the member by "
                "about 10% at beta = 1. R44 re-establishes the match at the "
                "scoring level; the R19 record is unchanged and remains "
                "archived under its own manifest."),
        "relationship_to_r18_r19": (
            "additive. The member trajectories are not re-propagated; their "
            "archived tighter-level telemetry supplies the work target. R44 "
            "adds one constant-degree trajectory per orbit per tolerance "
            "level in each cell."),
        "comparator_rule": (
            "N* = round(N_0 * sqrt(W_k/W_0)) with both works at the tighter "
            "level; where the beta-specific fixed_budget sidecar is not "
            "addressable through r14.reuse_paths, n_RHS for W_0 comes from "
            "the orbit's fixed_critical run and the source policy is "
            "recorded in the case config; the achieved ratio is measured "
            "from the propagated runs at both levels. Comparators at or "
            "above the adopted reference degree are censored, not clamped; "
            "no orbit required this in any cell."),
        "declared_outcome_returned": "B (see r44_preregistration.json)",
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R19"},
        "registration": index_files(REGISTRATION, METRICS),
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
                             "reused_inputs", "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
