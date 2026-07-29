"""SHA-256 integrity manifest for the R19 equal-realized-work comparison.

R19 changes only what is held equal. It reuses the R18 span members and the R14
tolerances, comparator degrees and truths unchanged, and propagates one new
trajectory per orbit: the constant degree matched to the realized total
quadratic work of the k = 0.5 member.

Usage:  python rev19_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r19_final_experiment_manifest.json"

SCRIPTS = ["rev19_equal_total_work.py", "rev19_tables.py",
           "rev19_finalize_manifest.py"]
RESULT_JSON = ["r19_equal_total_work_A.json", "r19_equal_total_work_B.json",
               "r19_manuscript_descriptives.json"]
TABLES = ["r19_equal_work_table.tex"]
REUSED = ["r18_span_sweep_A_beta_1.00.json", "r18_span_sweep_B_beta_1.00.json",
          "r14_trajectory_A_beta_1.00.json", "r14_trajectory_B_beta_1.00.json"]


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
    root = METRICS / "r19_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r19_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def censoring() -> dict:
    out = {}
    for design in ("A", "B"):
        p = METRICS / f"r19_censored_{design}.json"
        out[design] = (json.loads(p.read_text(encoding="utf-8"))
                       if p.exists() else [])
    return out


def main() -> int:
    payload = {
        "schema": "r19_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R19: the interior span member (k = 0.5) against a constant "
                  "degree matched on realized total quadratic work rather than "
                  "nominal per-call work, on both coverage designs."),
        "why": ("The per-call budget the span sweep holds equal is the same "
                "quantity the budget campaign shows the integrator does not "
                "preserve. R19 removes that objection by changing what is held "
                "equal and measuring, rather than assuming, the match achieved."),
        "relationship_to_r18": (
            "additive. The k = 0.5 trajectories are not re-propagated; their "
            "archived telemetry supplies the work target. R19 adds one constant"
            "-degree trajectory per orbit per tolerance level."),
        "comparator_rule": (
            "N* = round(N_0 * sqrt(W_k / W_0)), then propagated and its "
            "realized work ratio measured. Orbits whose N* reaches the adopted "
            "truth degree are censored, not clamped."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R18"},
        "censored_orbits": censoring(),
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
        print("[note] not produced: " + ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
