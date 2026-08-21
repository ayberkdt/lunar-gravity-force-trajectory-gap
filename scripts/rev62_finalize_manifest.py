"""SHA-256 integrity manifest for R62 (O54): the controlled apolune ladder's
interior panel, re-matched at the scoring tolerance.

R62 changes neither the method of (O42)/(O53) nor the ladder of (O49). It
applies the former to the latter's interior member, on both identity blocks,
at the two apolune levels where the reference-degree ceiling does not bind.

Partition. The ladder's own records (R50/R51) and the earlier rematches
(R44/R61) are untouched and no file of theirs is written: the driver
redirects the case root, the raw root and the record path to r62_*, and the
subset filter moves R44's hard-coded censored-list name into the r62_ prefix
as well. The R19 ladder records this campaign pairs against are read-only
inputs, indexed here and owned by their own manifests.

Usage:  python rev62_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r62_final_experiment_manifest.json"

KEYS = ("RS1", "RS2")
BETAS = ("0.50", "0.75", "1.00")
CELLS = [(k, b) for b in BETAS for k in KEYS]

SCRIPTS = ["rev62_ladder_interior_rematch.py", "rev62_preregister.py",
           "rev62_campaign.py", "rev62_tables.py",
           "rev62_finalize_manifest.py", "rev44_equal_work_tighter.py",
           "population_registry.py"]
RESULT_JSON = ([f"r62_ladder_interior_{k}_beta_{b}.json" for k, b in CELLS]
               + ["r62_manuscript_descriptives.json",
                  "r62_campaign_progress.json"])
REGISTRATION = ["r62_preregistration.json"]
TABLES = ["r62_ladder_interior_table.tex"]
REUSED = ([f"r18_span_sweep_{k}_beta_{b}.json" for k, b in CELLS]
          + [f"r19_equal_total_work_{k}_beta_{b}.json" for k, b in CELLS]
          + ["r50_span_ladder_a_rows.json", "r50_span_ladder_b_rows.json"])


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
    sidecars = {}
    root = METRICS / "r62_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r62_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome() -> dict:
    p = METRICS / "r62_manuscript_descriptives.json"
    if not p.exists():
        return {"missing": True}
    return json.loads(p.read_text(encoding="utf-8"))["outcome"]


def main() -> int:
    payload = {
        "schema": "r62_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R62 (O54): the interior span member (k = 0.5) of the "
                  "paired apolune ladder against a constant degree matched "
                  "on realized total quadratic work at the tighter "
                  "tolerance, on both identity blocks, at the 300 and 600 km "
                  "apolune levels, in six registered block-budget cells."),
        "why": ("(O53) established that the tolerance the realized-work "
                "match is made at changes which populations share the "
                "crossing bracket. That left the geometry strata carrying "
                "the level-consistent match while the controlled test of the "
                "radial-span direction still carried the older convention on "
                "its interior panel, holding the uncontrolled evidence to a "
                "stricter accounting than the controlled evidence."),
        "scope_bound": (
            "the 1200 and 2400 km levels are deliberately not run: (O49) "
            "reports the reference-degree ceiling beginning to bind at 1200 "
            "km and binding on every orbit at 2400 km, so a rematch there "
            "would confound the accounting change with the ceiling. The plan "
            "subcommand measured 0 censored comparators at 300 and 600 km in "
            "all six cells."),
        "relationship_to_r49_r50": (
            "additive and read-only. The ladder members are not "
            "re-propagated; their archived tighter-level telemetry supplies "
            "the work target, and only the matched constant-degree "
            "comparator is propagated, at both tolerance levels."),
        "statistics_rule": (
            "every tally is computed per apolune level and never pooled "
            "across levels: the campaign's question is the difference "
            "between the levels, which pooling would average away."),
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R61"},
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
