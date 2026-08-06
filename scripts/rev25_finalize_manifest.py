"""SHA-256 integrity manifest for the R25 crossover campaign (O33).

R25 measures the budget at which the interior-member result changes sign. The
claim holds at beta = 1 and reverses at beta = 0.5, and the interval between
them was never measured, so the midpoint was run with the archived drivers and
no new machinery:

  design A   the span sweep and the realized-work comparison at beta = 0.75
  design B   the budget campaign first, since design B had no archived R14
             record at that budget, then the same two stages on top of it

Naming exception, recorded rather than renamed away. Every stage reuses an
archived driver called with a budget argument it already accepts, so the
outputs keep their original prefixes and carry a beta_0.75 suffix. Files named
r14_*, r18_* and r19_* at beta = 0.75 belong to R25 and are indexed here only;
the R14, R18 and R19 manifests were sealed before these existed and their
inventories cover the budgets they were sealed on. The same convention was used
for the beta = 0.5 outputs held by R23.

Whatever a truncated stage left on disk is indexed as it stands. A stage that
did not run contributes nothing and is reported as absent rather than faked.

Usage:  python rev25_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r25_final_experiment_manifest.json"

BETA_TAG = "beta_0.75"

SCRIPTS = ["rev25_preregister.py", "rev25_amendment.py",
           "rev25_crossover_queue.py", "rev25_designB_queue.py",
           "rev25_bisect_queue.py", "rev25_supervisor.py",
           "rev25_finalize_manifest.py", "rev27_threshold_sensitivity.py",
           "rev14_budget_trajectory.py", "rev18_span_sweep.py",
           "rev19_equal_total_work.py", "rev19_tables.py"]

REGISTRATION = ["r25_preregistration.json",
                "r25_preregistration_amendment.json"]

# Records this campaign produces. Design B's three only exist if the second
# queue got its window.
RESULT_JSON = [
    f"r18_span_sweep_A_{BETA_TAG}.json",
    f"r19_equal_total_work_A_{BETA_TAG}.json",
    f"r14_trajectory_B_{BETA_TAG}.json",
    f"r18_span_sweep_B_{BETA_TAG}.json",
    f"r19_equal_total_work_B_{BETA_TAG}.json",
]

TABLES = ["r19_equal_work_table.tex", "r19_manuscript_descriptives.json",
          # the resolution-cut sensitivity table is derived from the same
          # equal-work records this campaign extends, propagates nothing, and
          # has no campaign of its own to belong to
          "r27_threshold_sensitivity_table.tex",
          "r27_threshold_sensitivity.json"]

# Inputs read but not produced here.
REUSED = ["r14_trajectory_A_beta_0.75.json",
          "r19_equal_total_work_A.json",
          "r19_equal_total_work_B.json",
          "r18_span_sweep_A_beta_1.00.json",
          "r18_span_sweep_B_beta_1.00.json"]

# Budgets this campaign propagated, and whether it had to build the budget
# record first. Design A already had one at 0.75; design B did not, and neither
# design had one at the budget above the anchor it reached.
PRODUCED = [
    ("A", "0.75", False),
    ("B", "0.75", True),
    ("A", "1.25", True),
    ("B", "1.50", True),
]
K_INTERIOR = ("0.25", "0.50", "0.75")


def _tree_for(design: str, beta: str, built_budget: bool) -> list[str]:
    tag = f"beta_{beta}"
    rels: list[str] = []
    if built_budget:
        rels += [f"r14_cases/{design}_{tag}", f"r14_raw/{design}_{tag}"]
    for k in K_INTERIOR:
        rels += [f"r18_cases/{design}_{tag}_k_{k}",
                 f"r18_raw/{design}_{tag}_k_{k}"]
    rels += [f"r19_cases/{design}_workmatched_{tag}",
             f"r19_raw/{design}_workmatched_{tag}"]
    return rels


TREES = {
    f"O33_design_{design}_beta_{beta}": _tree_for(design, beta, built)
    for design, beta, built in PRODUCED
}


def sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            d.update(blk)
    return d.hexdigest()


def index_files(names, base: Path, required: bool):
    """Absent files are an error only where the campaign promises them."""
    out, absent = {}, []
    for n in names:
        p = base / n
        if p.exists():
            out[n] = {"sha256": sha(p), "bytes": p.stat().st_size}
        else:
            out[n] = {"absent": True}
            if required:
                absent.append(n)
    return out, absent


def index_trees() -> dict:
    out = {}
    for name, rels in TREES.items():
        sidecars, roll, n_raw = {}, hashlib.sha256(), 0
        for rel in rels:
            base = METRICS / rel
            if not base.exists():
                continue
            for p in sorted(base.rglob("*.json")):
                sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
            for p in sorted(base.rglob("*.npz")):
                roll.update(sha(p).encode())
                n_raw += 1
        out[name] = {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
                     "sidecar_sha256": sidecars,
                     "raw_rollup_sha256": roll.hexdigest()}
    return out


def completeness() -> dict:
    """Read from the records, so a truncated stage cannot be quoted as full."""
    out = {}
    for design in ("A", "B"):
        entry = {}
        p = METRICS / f"r14_trajectory_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["budget_campaign"] = {
                "complete": d.get("complete"),
                "stopped_for_deadline": d.get("stopped_for_deadline")}
        p = METRICS / f"r18_span_sweep_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["span_sweep_orbits"] = len(d.get("rows", []))
        p = METRICS / f"r19_equal_total_work_{design}_{BETA_TAG}.json"
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            entry["equal_work"] = d.get("summary")
        if entry:
            out[design] = entry
        else:
            out[design] = {"not_run": True}
    return out


def main() -> int:
    scripts, a1 = index_files(SCRIPTS, CODE, required=True)
    registration, a2 = index_files(REGISTRATION, METRICS, required=True)
    results, _ = index_files(RESULT_JSON, METRICS, required=False)
    tables, a4 = index_files(TABLES, METRICS, required=True)
    reused, a5 = index_files(REUSED, METRICS, required=True)
    absent = a1 + a2 + a4 + a5

    payload = {
        "schema": "r25_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R25 (O33): the interior-member comparison at the midpoint "
                  "budget beta = 0.75, run to locate the sign change that the "
                  "archived beta = 1 and beta = 0.5 results could only "
                  "bracket. Design A reuses its archived budget record; "
                  "design B required that record to be built first."),
        "why": ("the constructive claim held at one budget and reversed at "
                "half of it, so it was stated on an interval whose width had "
                "never been measured"),
        "naming_exception": (
            "every stage reuses an archived driver with a budget argument it "
            "already accepts, so outputs keep the r14_, r18_ and r19_ prefixes "
            "and carry a beta_0.75 suffix. Those files belong to R25 and are "
            "indexed here only; the R14, R18 and R19 manifests were sealed "
            "before they existed. This is the convention R23 used for the "
            "beta = 0.5 outputs."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
        },
        "registration": registration,
        "reused_inputs": reused,
        "scripts": scripts,
        "result_json": results,
        "generated_tables": tables,
        "panel_completeness": completeness(),
        "trajectory_tree": index_trees(),
    }
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")
                   ).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    n_side = sum(t["n_sidecars"] for t in payload["trajectory_tree"].values())
    n_raw = sum(t["n_raw_arrays"] for t in payload["trajectory_tree"].values())
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}  "
          f"{n_side} sidecars, {n_raw} raw arrays")
    for design, info in payload["panel_completeness"].items():
        print(f"  design {design}: {json.dumps(info)[:110]}")
    if absent:
        print("[error] required files not on disk: " + ", ".join(absent))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
