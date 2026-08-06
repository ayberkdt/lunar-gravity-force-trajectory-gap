"""SHA-256 integrity manifest for the R14 fixed-budget Pareto campaign.

R14 asks one question the R12/R13 benchmark could not: at a gravity-work budget
declared in advance, is a constant degree or a radial degree history the better
way to spend it? It contains the frozen pre-registration, the 128-orbit
force-level budget sweep, the propagated fixed-budget trajectory comparisons,
the measured-time control, the forced-variational check at equal budget, and the
O26 allocation bound.

The numerical kernel is unchanged: every R14 trajectory sidecar records Lunaris
tag ``paper-truncation-v1.0`` at commit ``27e9ab86...`` with the same kernel and
gravity-file hashes as R10--R13. Truth and critical-degree trajectories are
reused from the R11 trees and remain indexed there; at beta = 1 the fixed
comparator IS the critical degree, so that row adds no new comparator
trajectories.

Usage:  python rev14_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import campaign_ownership as own

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r14_final_experiment_manifest.json"

SCRIPTS = [
    "rev14_preregister.py",
    "rev14_budget_pareto.py",
    "rev14_budget_trajectory.py",
    "rev14_oracle.py",
    "rev14_timing_budget.py",
    "rev14_variational_budget.py",
    "rev14_tables.py",
    "make_figures_r14.py",
    "run_r14_trajectories.sh",
    "rev14_finalize_manifest.py",
]

RESULT_JSON = [
    "r14_preregistration.json",
    "r14_budget_pareto.json",
    "r14_oracle.json",
    "r14_timing_selection.json",
    "r14_timing_budget.json",
    "r14_variational_budget.json",
    "r14_descriptives.json",
]

# The timing table and the force Pareto table are each emitted twice, as the
# main-text extract and the full supplementary listing; an earlier version of
# this list named a single "r14_timing_budget_table.tex" that rev14_tables.py
# has never written, and recorded it as missing instead of failing.
TABLES = [
    "r14_force_pareto_table.tex",
    "r14_force_pareto_table_full.tex",
    "r14_trajectory_pareto_table.tex",
    "r14_cost_bookkeeping_table.tex",
    "r14_cap_audit_table.tex",
    "r14_oracle_table.tex",
    "r14_timing_budget_table_compact.tex",
    "r14_timing_budget_table_full.tex",
    "r14_variational_budget_table.tex",
    "r14_beta1_per_orbit_A.tex",
    "r14_beta1_per_orbit_B.tex",
]

FIGURES = ["budget_pareto.pdf"]


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


def index_tree(case_dir: Path, raw_dir: Path) -> dict:
    """Index only the budgets this campaign propagated.

    Later campaigns reuse this driver with a budget argument, so their
    trajectories land under the same prefix; indexing them here would put the
    same records under two manifests. Ownership is declared in
    ``campaign_ownership``.
    """
    sidecars = {}
    for p in sorted(case_dir.rglob("*.json")):
        if not own.owned_by_r14(p):
            continue
        sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n_raw = 0
    for p in sorted(raw_dir.rglob("*.npz")):
        if not own.owned_by_r14(p):
            continue
        roll.update(sha(p).encode())
        n_raw += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n_raw,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def main() -> int:
    traj = sorted(METRICS.glob("r14_trajectory_*.json"))
    payload = {
        "schema": "r14_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "scope": ("R14 extension (O25/O26): fixed-budget radial-allocation Pareto "
                  "study. Pre-registered protocol; 128-orbit integration-noise-free "
                  "force-defect sweep over the budget grid; propagated fixed-budget "
                  "trajectory comparisons; measured-serial-time budget control; "
                  "forced-variational check at equal budget; and a "
                  "trajectory-informed force-defect allocation bound."),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R13"},
        "reused_evidence": (
            "archived R11 truth trajectories are the reference states for the "
            "force-level sweep and the truth for every R14 comparison; at beta = 1 "
            "the equal-budget comparator degree equals the critical-altitude degree "
            "on every orbit, so the archived R11 fixed_critical trajectories are "
            "reused unchanged. Both are indexed in the R11 manifest."),
        "preregistration": {
            "path": "metrics/r14_preregistration.json",
            "note": ("hypotheses, budget grid, calibration/comparator/censoring "
                     "rules, staging order, adaptive grid-extension rule and "
                     "decision logic, all frozen before any aggregate result was "
                     "inspected")},
        "timing_note": ("trajectory-campaign kernel times were recorded under five "
                        "concurrent workers and are not comparable; the serial runs "
                        "of r14_timing_budget.json are the only timing reference"),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "trajectory_result_json": index_files([p.name for p in traj], METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "generated_figures": index_files(FIGURES, ROOT / "figures"),
        "trajectory_trees": {},
    }
    case_root, raw_root = METRICS / "r14_cases", METRICS / "r14_raw"
    if case_root.exists():
        for case_dir in sorted(case_root.iterdir()):
            if not case_dir.is_dir():
                continue
            raw_dir = raw_root / case_dir.name
            payload["trajectory_trees"][case_dir.name] = index_tree(case_dir, raw_dir)
            t = payload["trajectory_trees"][case_dir.name]
            print(f"[tree] {case_dir.name}: {t['n_sidecars']} sidecars, "
                  f"{t['n_raw_arrays']} raw arrays")
    prereg = METRICS / "r14_preregistration.json"
    if prereg.exists():
        payload["preregistration"]["protocol_sha256"] = json.loads(
            prereg.read_text(encoding="utf-8"))["protocol_sha256"]
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {OUT.name}  manifest_sha256={payload['manifest_sha256'][:16]}")
    missing = [k for sec in ("scripts", "result_json", "generated_tables",
                             "generated_figures")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
