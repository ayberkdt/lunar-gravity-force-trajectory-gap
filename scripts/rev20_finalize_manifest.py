"""SHA-256 manifest for the R20 sixty-day span check."""
from __future__ import annotations
import hashlib, json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS, CODE = ROOT / "metrics", ROOT / "python_codes"
OUT = METRICS / "r20_final_experiment_manifest.json"
SCRIPTS = ["rev20_span_longarc.py", "rev20_finalize_manifest.py"]
RESULT = ["r20_span_longarc.json"]
REUSED = ["r18_span_sweep_A_beta_1.00.json", "r17_longarc60.json",
          "r14_trajectory_A_beta_1.00.json"]


def sha(p: Path) -> str:
    d = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def index(names, base):
    return {n: ({"sha256": sha(base / n), "bytes": (base / n).stat().st_size}
                if (base / n).exists() else {"missing": True}) for n in names}


def tree():
    sc = {str(p.relative_to(METRICS)).replace("\\", "/"): sha(p)
          for p in sorted((METRICS / "r20_cases").rglob("*.json"))}
    roll, n = hashlib.sha256(), 0
    for p in sorted((METRICS / "r20_raw").rglob("*.npz")):
        roll.update(sha(p).encode()); n += 1
    return {"n_sidecars": len(sc), "n_raw_arrays": n,
            "sidecar_sha256": sc, "raw_rollup_sha256": roll.hexdigest()}


payload = {
    "schema": "r20_final_experiment_manifest_v1",
    "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "scope": ("R20: the constant-to-radial interpolation family propagated for "
              "60 days at beta = 1 on the eight design-A orbits carrying an "
              "archived 60-day truth."),
    "relationship_to_r18": ("the degree tables are read from the frozen R18 "
                            "seven-day sidecars and reused verbatim; only the "
                            "arc length differs. Both endpoints are propagated "
                            "here because R14 is a seven-day campaign."),
    "numerical_kernel": {"lunaris_release_tag": "paper-truncation-v1.0",
                         "lunaris_commit":
                         "27e9ab86ed61d623f78c453ea2054348f1044c23"},
    "reused_inputs": index(REUSED, METRICS),
    "scripts": index(SCRIPTS, CODE),
    "result_json": index(RESULT, METRICS),
    "trajectory_tree": tree(),
}
body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
payload["manifest_sha256"] = hashlib.sha256(
    json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
t = payload["trajectory_tree"]
print(f"[written] {OUT.name}  {t['n_sidecars']} sidecars, "
      f"{t['n_raw_arrays']} raw arrays")
missing = [k for sec in ("scripts", "result_json", "reused_inputs")
           for k, v in payload[sec].items() if v.get("missing")]
if missing:
    print("[error] recorded as missing: " + ", ".join(missing))
    raise SystemExit(1)
