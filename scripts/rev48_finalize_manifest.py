"""SHA-256 integrity manifest for the R48 (O48) interior measured-time panel.

R48 applies the R13 measured-time protocol to the interior span member at
beta = 1: a serial contention-free member baseline at the tight level, a
comparator matched on measured total kernel time, and errors scored at the
tighter level against the archived truths.

Usage:  python rev48_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r48_final_experiment_manifest.json"

SCRIPTS = ["rev48_interior_timing.py", "rev48_tables.py",
           "rev48_finalize_manifest.py", "run_overnight_20260809.ps1"]
RESULT_JSON = ["r48_interior_timing.json",
               "r48_interior_timing_selection.json"]
REGISTRATION = ["r48_preregistration.json"]
TABLES = ["r48_interior_timing_table.tex"]
REUSED = ["r12_kernel_cost_curve.json",
          "r18_span_sweep_A_beta_1.00.json",
          "r18_span_sweep_B_beta_1.00.json"]


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
    root = METRICS / "r48_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r48_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    payload = {
        "schema": "r48_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R48 (O48): the interior span member (k = 0.5) at beta = 1 "
                  "against a constant degree matched on measured serial "
                  "total kernel time, on a 14-orbit registered panel, seven "
                  "per design spread over perilune with the extremes "
                  "retained."),
        "why": ("The measured-time control existed only for the radial "
                "endpoint; the continuation family was compared under "
                "operation-count proxies whose low-degree flattening the "
                "cost record shows biased in the member's favor. R48 "
                "measures that residual instead of bounding it by "
                "argument."),
        "protocol": (
            "R13 protocol applied to the member: serial member baseline at "
            "the tight level on an otherwise idle machine (every timed "
            "stage refuses to start while other python processes are "
            "alive); first-pass comparator degree from the call-weighted "
            "mean per-call measured cost of the member's archived degree "
            "histogram; one refinement pass on the measured total-kernel-"
            "time ratio with c(N) ~ N^2, capped at the adopted truth "
            "degree; refined comparator propagated at both tolerance "
            "levels; errors scored at the tighter level against the "
            "archived truths under the campaign envelope rule."),
        "declared_outcome_returned": "C (see r48_preregistration.json)",
        "label_note": ("O48 is the label the registration was filed under; "
                       "O43-O47 are unassigned."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R44"},
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
