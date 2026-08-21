"""SHA-256 integrity manifests for R58 and R59.

Two small campaigns are sealed here because neither owns a trajectory tree the
other could be confused with, and both were added in the same revision.

R58 is the post-hoc realized-work control on the wide-elliptic endpoint: it
propagates one constant-degree comparator per orbit per tolerance level, at the
degree that matches the radial endpoint's realized total quadratic work, and
reuses the R18 span records and the R11 references unchanged.

R59 propagates nothing at all. It re-tallies archived R14 rows over leading
Sobol prefixes, so its inputs are records other manifests already own and its
only outputs are one record and one table. It is sealed anyway: the manifest
policy is that everything in metrics/ is owned by some manifest, and an
analysis that prints numbers into the manuscript is exactly the kind of file
that should not be exempt from that.

Usage:  python rev58_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"

R58_OUT = METRICS / "r58_final_experiment_manifest.json"
R59_OUT = METRICS / "r59_final_experiment_manifest.json"

BETAS = ("0.50", "0.62", "0.75", "1.00")
POPS = ("OE", "OEU")

R58_SCRIPTS = ["rev58_endpoint_equal_work.py", "rev58_chain.py",
               "rev58_tables.py", "rev58_claims.py",
               "rev58_finalize_manifest.py"]
R58_RESULTS = [f"r58_endpoint_equal_work_{p}_beta_{b}.json"
               for p in POPS for b in BETAS]
R58_TABLES = ["r58_endpoint_equal_work_table.tex"]
R58_REUSED = ([f"r18_span_sweep_{p}_beta_{b}.json" for p in POPS for b in BETAS]
              + ["r31_operational_elliptical_rows.json",
                 "r38_operational_elliptical_uncapped_rows.json"])

R59_SCRIPTS = ["rev59_design_size.py"]
R59_RESULTS = ["r59_design_size_convergence.json"]
R59_TABLES = ["r59_design_size_table.tex"]
R59_REUSED = [f"r14_trajectory_{d}_beta_1.00.json" for d in ("A", "B", "C")]


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
    """R58 owns its whole case and raw trees; no other campaign writes there."""
    sidecars = {}
    root = METRICS / "r58_cases"
    if root.exists():
        for p in sorted(root.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    raw = METRICS / "r58_raw"
    if raw.exists():
        for p in sorted(raw.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def censoring() -> dict:
    out = {}
    for p in POPS:
        for b in BETAS:
            f = METRICS / f"r58_censored_{p}_beta_{b}.json"
            if f.exists():
                out[f"{p}_beta_{b}"] = json.loads(
                    f.read_text(encoding="utf-8"))
    return out


def stamp(payload: dict, out: Path) -> dict:
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {out.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    return payload


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def main() -> int:
    r58 = stamp({
        "schema": "r58_final_experiment_manifest_v1",
        "created_utc": now(),
        "scope": ("R58: the budget-calibrated radial endpoint against a "
                  "constant degree matched on realized total quadratic work "
                  "rather than nominal per-call work, on the wide-elliptic "
                  "population under both reference degrees, at four budgets."),
        "why": ("At the declared budget the endpoint wins that population "
                "while spending about a third more realized work than its "
                "comparator, so the win could be read as bought rather than "
                "earned. This control pays the comparator the same work and "
                "re-scores."),
        "status": ("post hoc. It was added after the registered outcomes of "
                   "(O37) and (O39) were scored, changes no registered "
                   "verdict, and is reported as its own record."),
        "relationship_to_r18": (
            "additive. The endpoint trajectories are not re-propagated; their "
            "archived telemetry supplies the work target. R58 adds one "
            "constant-degree trajectory per orbit per tolerance level."),
        "comparator_rule": (
            "N* = round(N_0 * sqrt(W_rad / W_0)), then propagated and its "
            "realized work ratio measured. Orbits whose N* reaches the "
            "adopted reference degree are censored, not clamped."),
        "design_registry_note": (
            "OE and OEU are not in rev14's static registry; rev58 registers "
            "both unconditionally at import so that spawned pool workers "
            "resolve the same rows as the parent."),
        "censored_orbits": censoring(),
        "reused_inputs": index_files(R58_REUSED, METRICS),
        "scripts": index_files(R58_SCRIPTS, CODE),
        "result_json": index_files(R58_RESULTS, METRICS),
        "generated_tables": index_files(R58_TABLES, METRICS),
        "trajectory_tree": index_tree(),
    }, R58_OUT)

    stamp({
        "schema": "r59_final_experiment_manifest_v1",
        "created_utc": now(),
        "scope": ("R59: the budget comparison re-tallied over leading "
                  "scrambled-Sobol prefixes of each coverage design, and the "
                  "spread of the three independent scrambles at full size."),
        "why": ("Design size and scramble-to-scramble spread were asserted "
                "rather than measured. Both are answerable from archived rows."),
        "propagation": ("none. No orbit is propagated and none is re-scored: "
                        "the resolution flag counted is the one the original "
                        "campaign wrote."),
        "reused_inputs": index_files(R59_REUSED, METRICS),
        "scripts": index_files(R59_SCRIPTS, CODE),
        "result_json": index_files(R59_RESULTS, METRICS),
        "generated_tables": index_files(R59_TABLES, METRICS),
    }, R59_OUT)

    t = r58["trajectory_tree"]
    print(f"  R58: {t['n_sidecars']} sidecars, {t['n_raw_arrays']} raw arrays")
    missing = []
    for out in (R58_OUT, R59_OUT):
        p = json.loads(out.read_text(encoding="utf-8"))
        for sec in ("scripts", "result_json", "generated_tables",
                    "reused_inputs"):
            missing += [f"{out.name}:{k}" for k, v in p.get(sec, {}).items()
                        if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
