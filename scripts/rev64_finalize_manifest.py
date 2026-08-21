"""SHA-256 integrity manifest for R64 (O57): the measured-kernel-time
comparator matched at the tolerance the errors are scored at.

R64 adds no orbit and no reference. It repeats (O48)'s protocol on (O48)'s own
panel with every timed stage and the first-pass call histogram moved to the
tighter level. (O48)'s records are inputs and are not written.

The campaign carries the number 64 because 57 was already held by the (O40)
gradient-degree completion; the registration records that renumbering and its
timestamp still precedes every propagation here.

Usage:  python rev64_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r64_final_experiment_manifest.json"

SCRIPTS = ["rev64_interior_timing_tighter.py", "rev64_preregister.py",
           "rev64_finalize_manifest.py", "rev48_interior_timing.py",
           "rev13_timing_match.py"]
RESULT_JSON = ["r64_interior_timing_tighter.json",
               "r64_interior_timing_state.json"]
REGISTRATION = ["r64_preregistration.json"]
REUSED = ["r48_interior_timing_selection.json", "r48_interior_timing.json",
          "r12_kernel_cost_curve.json", "r18_span_sweep_A_beta_1.00.json",
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
    root = METRICS / "r64_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r64_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome() -> dict:
    p = METRICS / "r64_interior_timing_tighter.json"
    if not p.exists():
        return {"missing": True}
    return json.loads(p.read_text(encoding="utf-8"))["summary"]


def main() -> int:
    payload = {
        "schema": "r64_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R64 (O57): the k = 0.5 interior member against a constant "
                  "degree matched on measured kernel time at the tighter "
                  "tolerance, on (O48)'s fourteen-orbit panel at beta = 1."),
        "why": ("(O48) matched on kernel time measured at the tight level "
                "while scoring every error at the tighter level, the same "
                "level inconsistency (O42) removed for realized work. This "
                "campaign moves the timing to the scoring level and changes "
                "nothing else."),
        "numbering_note": (
            "registered as r57 and renumbered to r64 before the manifest was "
            "written: r57 is held by the (O40) gradient-degree completion. "
            "The registration's timestamp is the original and precedes every "
            "propagation indexed here."),
        "relationship_to_r48": (
            "same panel, same member, same budget, same refinement rule. The "
            "member's contention-free re-run, the first-pass comparator and "
            "the refined comparator are timed at the tighter level instead of "
            "the tight one, and the first-pass degree is inverted on the "
            "member's tighter call histogram. The refined comparator is still "
            "propagated at both levels because the envelope needs the pair; "
            "the tight run is untimed. R48's records are read, never "
            "written."),
        "timing_contract": (
            "every timed stage refuses to start while other python processes "
            "are alive, so the kernel times are contention-free. Achieved "
            "time ratios are measured from the propagated runs and reported "
            "per orbit, including the two that fall outside the 0.90-1.10 "
            "band the protocol aims at."),
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R63"},
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
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
