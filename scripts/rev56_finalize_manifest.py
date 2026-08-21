"""SHA-256 integrity manifest for R56 (O56): the interior member at sixty
days, recalibrated on the sixty-day arc and matched at the scoring tolerance.

R56 adds no population and no reference. It recalibrates one policy on the
arc it is scored over, sizes a constant comparator from that policy's own
tighter-level telemetry, and propagates both over the archived (O27)
references. The (O30) records it is read against are inputs, not outputs.

Usage:  python rev56_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r56_final_experiment_manifest.json"

SCRIPTS = ["rev56_longarc_interior.py", "rev56_preregister.py",
           "rev56_campaign.py", "rev56_finalize_manifest.py",
           "rev20_span_longarc.py", "rev18_span_sweep.py",
           "rev17_longarc60.py"]
RESULT_JSON = ["r56_longarc_interior.json", "r56_campaign_progress.json"]
REGISTRATION = ["r56_preregistration.json"]
REUSED = ["r17_longarc60.json", "r20_span_longarc.json",
          "r18_span_sweep_A_beta_1.00.json"]


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
    root = METRICS / "r56_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r56_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome() -> dict:
    p = METRICS / "r56_longarc_interior.json"
    if not p.exists():
        return {"missing": True}
    return json.loads(p.read_text(encoding="utf-8"))["summary"]


def main() -> int:
    extra = []
    if (METRICS / "r56_censored.json").exists():
        extra.append("r56_censored.json")
    payload = {
        "schema": "r56_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R56 (O56): the k = 0.5 interior member against a constant "
                  "degree over sixty days, on the eight Design-A orbits "
                  "carrying a sixty-day reference, at beta = 1. The member's "
                  "degree table is recalibrated so that its mean squared "
                  "degree meets the budget over the sixty-day reference "
                  "epochs, and the comparator is matched on realized total "
                  "quadratic work read at the tighter tolerance."),
        "why": ("(O30) reused the frozen seven-day table over the sixty-day "
                "arc, where the member overspends the per-call budget by a "
                "median factor of 1.22, and scored the comparison on the "
                "nominal per-call accounting. Its negative result therefore "
                "confounds a horizon effect with a calibration mismatch and "
                "an accounting convention. This campaign removes both."),
        "relationship_to_r17_r20": (
            "additive and read-only. The sixty-day references come from (O27) "
            "unchanged at both tolerance levels; the archived (O30) constant "
            "endpoint supplies W_0 at the tighter level. Neither record is "
            "written."),
        "locked_choices": (
            "panel, budget and interior index were fixed in "
            "r56_preregistration.json before any propagation and are not "
            "revisited: eight orbits, beta = 1, k = 0.5."),
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R63"},
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON + extra, METRICS),
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
    missing = [k for sec in ("scripts", "result_json", "reused_inputs",
                             "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
