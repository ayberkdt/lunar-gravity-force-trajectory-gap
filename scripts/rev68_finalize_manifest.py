"""SHA-256 integrity manifest for R68 (O59/O60): the full-design comparison at
equal measured kernel time.

R68 adds no orbit, no reference and no policy. It takes the measured-time
construction the panel campaigns established -- (O13) for the radial endpoint,
(O48)/(O57)/(O58) for the interior member -- and runs it on the whole of both
coverage designs, so that the compute-allocation framing is matched on compute
as the machine spends it rather than on an operation-count proxy.

Two things this manifest indexes that its neighbours do not. Both arms ran in
one exclusive 11.56-hour session, so the contention audit is a real question
rather than a formality and its record is indexed here. And because the
member's error is reused from the archived campaign while its time is measured
fresh, the reproduction check that licenses that reuse is indexed as evidence
rather than asserted in prose.

Usage:  python rev68_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r68_final_experiment_manifest.json"

SCRIPTS = ["rev68_timing_full.py", "rev68_chain.py", "rev68_preregister.py",
           "rev68_band_audit.py", "rev68_member_reproduction.py",
           "rev68_tables.py", "rev68_quiet_audit.py",
           "rev68_finalize_manifest.py", "probe_kernel_quiet.py",
           "rev65_quiet_audit.py", "rev65_timing_family.py",
           "rev48_interior_timing.py", "rev13_timing_match.py",
           "rev14_budget_trajectory.py"]
TABLES = ["r68_measured_time_summary_table.tex", "r68_band_table.tex",
          "r68_apolune_posthoc_table.tex"]
RESULT_JSON = ["r68_timing_full_endpoint.json",
               "r68_timing_full_interior.json",
               "r68_timing_full_endpoint_state.json",
               "r68_timing_full_interior_state.json",
               "r68_band_admissibility.json",
               "r68_member_reproduction.json",
               "r68_quiet_audit.json"]
REGISTRATION = ["r68_preregistration.json"]
# r65_quiet_probe_log.json is R65's file: probe_kernel_quiet.py appends to it,
# so this campaign's pre-start probe is recorded there and the file is indexed
# as a reused input rather than claimed.
REUSED = ["r12_kernel_cost_curve.json",
          "r14_trajectory_A_beta_1.00.json", "r14_trajectory_B_beta_1.00.json",
          "r18_span_sweep_A_beta_1.00.json", "r18_span_sweep_B_beta_1.00.json",
          "r48_interior_timing.json", "r64_interior_timing_tighter.json",
          "r65_timing_family.json", "r65_quiet_probe_log.json"]

OWNED_PREFIX = "r68_"


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


def index_tree(cases: str, raw: str) -> dict:
    if not cases.startswith(OWNED_PREFIX) or not raw.startswith(OWNED_PREFIX):
        raise ValueError(f"refusing to claim {cases}/{raw}: not {OWNED_PREFIX}")
    sidecars = {}
    croot = METRICS / cases
    if croot.exists():
        for p in sorted(croot.rglob("*.json")):
            sidecars[str(p.relative_to(METRICS)).replace("\\", "/")] = sha(p)
    roll = hashlib.sha256()
    n = 0
    rroot = METRICS / raw
    if rroot.exists():
        for p in sorted(rroot.rglob("*.npz")):
            roll.update(sha(p).encode())
            n += 1
    return {"n_sidecars": len(sidecars), "n_raw_arrays": n,
            "sidecar_sha256": sidecars, "raw_rollup_sha256": roll.hexdigest()}


def outcome() -> dict:
    out = {}
    for arm in ("endpoint", "interior"):
        p = METRICS / f"r68_timing_full_{arm}.json"
        if not p.exists():
            out[arm] = {"missing": True}
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        out[arm] = {"member_k": d["member_k"], "timing_band": d["timing_band"],
                    "by_design": d["by_design"]}
    return out


def main() -> int:
    payload = {
        "schema": "r68_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R68: every orbit of coverage designs A and B at beta = 1, "
                  "each member against a constant degree refined until its "
                  "measured gravity-kernel time falls within 0.95-1.05 of the "
                  "member's. (O59) is the budget-calibrated radial endpoint "
                  "k = 1.00; (O60) is the sampled interior member k = 0.50. "
                  "128 cells per arm, none selected and none dropped."),
        "why": ("the principal population comparison is matched on a per-call "
                "quadratic proxy calibrated on the reference arc, and "
                "measured machine time had been matched only on declared "
                "fourteen-orbit panels that are stated not to be population "
                "tallies."),
        "evidential_status": {
            "O59": ("prospectively registered; it repeats an existing "
                    "population comparison under a different cost match and "
                    "returns the registered outcome A"),
            "O60": ("prospectively registered full-design measured-time "
                    "validation of an exploratory candidate. The member "
                    "k = 0.5 was adopted after the nominal family sweep had "
                    "been seen, and this campaign does not make it "
                    "confirmatory: what is registered in advance is this "
                    "test, not the choice of member. The interior family "
                    "keeps its exploratory status."),
        },
        "matching_convention": (
            "the match and the scoring share the tighter level, which is "
            "R44's level-consistency requirement applied to time rather than "
            "to work. The member is re-propagated serially for a "
            "contention-free kernel time; its error is reused from the "
            "archived campaign and metrics/r68_member_reproduction.json shows "
            "the fresh trajectory reproduces the archived error exactly on "
            "all 256 cells."),
        "admissibility": (
            "the registration prescribes the band and states that a cell "
            "still outside it keeps its nearest integer match and is flagged. "
            "Misses are therefore not an exclusion criterion and the primary "
            "counts include them; metrics/r68_band_admissibility.json carries "
            "the without-misses tally as a sensitivity and the side each miss "
            "fell on. Ten cells of 256 are out of band, and the only resolved "
            "out-of-band cell that favours a variable member is B040 at "
            "1.051, where the comparator had more time than the member, not "
            "less."),
        "post_hoc": (
            "the apolune stratification in r68_apolune_posthoc_table.tex was "
            "not registered. It is a diagnostic for why a perilune-stratified "
            "fourteen-orbit panel and the full designs read differently, and "
            "carries no causal or preregistered geometry claim."),
        "machine": (
            "both arms ran in one exclusive session of 11.56 h with no other "
            "python process on the machine. probe_kernel_quiet.py reproduced "
            "the archived idle-machine cost curve to 1.1 per cent before the "
            "start, and metrics/r68_quiet_audit.json gives degree-normalized "
            "within-orbit throughput spreads of 1.19 and 1.23 against the "
            "accepted R65 panel's 1.27, with no case above 1.25 times the "
            "opening baseline."),
        "numbering_note": ("r66 is held by the allocation-anatomy figure and "
                           "r67 by the review-measurement record, so this "
                           "campaign takes r68; the observation labels are "
                           "(O59) and (O60)."),
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R67"},
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": index_tree("r68_cases", "r68_raw"),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    t = payload["trajectory_tree"]
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  trajectories: {t['n_sidecars']} sidecars, "
          f"{t['n_raw_arrays']} raw arrays")
    missing = [k for sec in ("scripts", "result_json", "reused_inputs",
                             "generated_tables", "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("  !! missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
