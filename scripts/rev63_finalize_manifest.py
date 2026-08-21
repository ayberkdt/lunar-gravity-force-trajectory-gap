"""SHA-256 integrity manifest for R63 (O55): the ceiling-free apolune ladder,
interior panel, re-matched at the scoring tolerance.

R63 is (O54)'s construction on the ceiling-free blocks of (O50) and (O51), at
all four apolune levels. On the two narrow levels it repeats a measurement
(O54) already made, independently propagated; on the two wide levels it makes
one (O54) could not, because on the capped blocks the reference-degree
ceiling binds there.

Partition: r63_* records and trees. R44, R61, R62 and the ladder campaigns
are untouched, and R44's hard-coded censored-list name is moved into the r63_
prefix by the driver rather than left in the sealed campaign's file family.

Usage:  python rev63_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r63_final_experiment_manifest.json"

KEYS = ("RS1U", "RS2U")
BETAS = ("1.00", "0.75", "0.50")
CELLS = [(k, b) for b in BETAS for k in KEYS]

SCRIPTS = ["rev63_ladder_uncapped_rematch.py", "rev63_preregister.py",
           "rev63_campaign.py", "rev63_tables.py",
           "rev63_finalize_manifest.py", "rev44_equal_work_tighter.py",
           "population_registry.py"]
RESULT_JSON = ([f"r63_ladder_uncapped_{k}_beta_{b}.json" for k, b in CELLS]
               + ["r63_manuscript_descriptives.json",
                  "r63_campaign_progress.json"])
REGISTRATION = ["r63_preregistration.json"]
TABLES = ["r63_ladder_uncapped_table.tex"]
REUSED = ([f"r18_span_sweep_{k}_beta_{b}.json" for k, b in CELLS]
          + [f"r19_equal_total_work_{k}_beta_{b}.json" for k, b in CELLS]
          + ["r51_span_ladder_a_uncapped_rows.json",
             "r52_span_ladder_b_uncapped_rows.json"])


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
    root = METRICS / "r63_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r63_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome() -> dict:
    p = METRICS / "r63_manuscript_descriptives.json"
    if not p.exists():
        return {"missing": True}
    return json.loads(p.read_text(encoding="utf-8"))["outcome"]


def main() -> int:
    payload = {
        "schema": "r63_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R63 (O55): the interior span member (k = 0.5) of the "
                  "ceiling-free apolune ladder against a constant degree "
                  "matched on realized total quadratic work at the tighter "
                  "tolerance, on both ceiling-free identity blocks, at all "
                  "four apolune levels, in six registered block-budget "
                  "cells."),
        "why": ("(O54) carried the level-consistent match to the controlled "
                "ladder but could only run the two narrow levels, because on "
                "the capped blocks the reference-degree ceiling binds above "
                "them. The ceiling-free blocks remove that confound, so this "
                "campaign both repeats (O54) where the two overlap and "
                "extends it to the two widest levels."),
        "relationship_to_r62": (
            "the 300 and 600 km levels of this campaign and of R62 are the "
            "same measurement made twice. Their adopted reference degrees are "
            "identical at those levels, so the ceiling never distinguished "
            "them, and the two campaigns' interior errors agree to about four "
            "significant figures; they differ only in that each propagated "
            "its arcs independently. One block-B comparison at 300 km and "
            "beta = 0.50 sits close enough to its resolution threshold to "
            "fall on different sides of it in the two runs, which is a "
            "measurement of how reproducible a borderline comparison is "
            "rather than a disagreement about a verdict."),
        "relationship_to_r50_r51_r52": (
            "additive and read-only. The ladder members are not "
            "re-propagated; their archived tighter-level telemetry supplies "
            "the work target, and only the matched constant-degree "
            "comparator is propagated, at both tolerance levels."),
        "statistics_rule": (
            "every tally is per block, per budget and per apolune level, "
            "never pooled across levels or blocks: the campaign's question is "
            "the shape of the margin across levels, which pooling destroys."),
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R62"},
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
