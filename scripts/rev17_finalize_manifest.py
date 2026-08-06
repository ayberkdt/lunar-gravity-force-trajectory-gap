"""SHA-256 integrity manifest for the R17 sixty-day campaign.

R17 extends the month-long stage in three ways at once -- 60-day arcs, a
widened geometry set drawn by a rule fixed in advance, and a two-level
tolerance ladder that gives every comparison its own resolution envelope. The
manifest indexes the driver, the record, the generated tables, and the
per-trajectory sidecars with a rolled-up digest of the raw state arrays.

Usage:  python rev17_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r17_final_experiment_manifest.json"

SCRIPTS = ["rev17_longarc60.py", "rev17_tables.py",
           "rev17_finalize_manifest.py"]
RESULT_JSON = ["r17_longarc60.json"]
TABLES = ["r17_longarc60_table.tex", "r17_longarc60_growth_table.tex"]


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
    case_root = METRICS / "r17_cases" / "longarc60"
    raw_root = METRICS / "r17_raw" / "longarc60"
    sidecars = {}
    if case_root.exists():
        for p in sorted(case_root.rglob("run_*/*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n_raw = 0
    if raw_root.exists():
        for p in sorted(raw_root.rglob("run_*/*.npz")):
            roll.update(sha(p).encode())
            n_raw += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def campaign_scope() -> dict:
    src = METRICS / "r17_longarc60.json"
    if not src.exists():
        return {}
    d = json.loads(src.read_text(encoding="utf-8"))
    s = d.get("summary", {})
    return {
        "geometry_rule": d["scenario"]["geometry_rule"],
        "duration_days": d["scenario"]["duration_days"],
        "levels": d["scenario"]["levels"],
        "orbits_attempted": s.get("orbits_attempted"),
        "orbits_reaching_full_arc": s.get("orbits_reaching_60_days"),
        "orbits_terminated_early": s.get("orbits_terminated_early"),
        "orbits_with_work_matched_comparator": s.get(
            "orbits_with_work_matched_comparator"),
        "reconstructed_from_disk": d.get("reconstructed_from_disk", False),
    }


def main() -> int:
    payload = {
        "schema": "r17_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R17 sixty-day long arcs on a widened geometry set, with a "
                  "two-level vector-tolerance ladder supplying a numerical "
                  "envelope for every policy and a resolution rule for every "
                  "comparison."),
        "relationship_to_r7_stage3": (
            "extends rather than supersedes. The 28-day stage-3 result stands; "
            "R17 lengthens the arc, widens the geometry set by a rule fixed "
            "before propagation, and adds the resolution envelope that the "
            "single-tolerance 28-day stage lacked. Geometries with perilune "
            "below 50 km need an N=600 truth and are not extended, so the "
            "LRO-like case of stage 3 is absent here by cost."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R15"},
        "campaign": campaign_scope(),
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
    missing = [k for sec in ("scripts", "result_json", "generated_tables")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
