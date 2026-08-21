"""SHA-256 integrity manifest for R65 (O58): the sampled interior family under
level-consistent measured time.

R65 adds no orbit and no reference. It takes (O57)'s protocol and asks it of
the family rather than of one member: k = 0.25 and k = 0.75 are propagated
here, each against its own constant degree matched on measured kernel time at
the tighter tolerance, and the k = 0.50 column is reused from R64 unchanged.
R48's and R64's records are inputs and are not written.

Two things this manifest indexes that its neighbours do not. The campaign
carries the number 65 because 58 is held by the equal-work endpoint control.
And a first attempt ran to completion on a machine that did not stay idle and
was discarded against the same registration; the attempt's artifacts, its own
outcome and the measurements that condemn them are indexed here so the discard
is auditable rather than merely asserted.

Usage:  python rev65_finalize_manifest.py
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CODE = ROOT / "python_codes"
OUT = METRICS / "r65_final_experiment_manifest.json"

SCRIPTS = ["rev65_timing_family.py", "rev65_preregister.py",
           "rev65_record_discarded.py", "rev65_quiet_audit.py",
           "rev65_finalize_manifest.py", "make_accounting_ladder_table.py",
           "probe_kernel_quiet.py", "rev64_interior_timing_tighter.py",
           "rev48_interior_timing.py", "rev13_timing_match.py"]
# The accounting ladder is emitted from this campaign's own timing record and
# from R18's and R44's, so it is this campaign's product to own.
TABLES = ["r65_accounting_ladder_table.tex"]
RESULT_JSON = ["r65_timing_family.json", "r65_timing_family_state.json",
               "r65_discarded_attempt1.json", "r65_quiet_audit.json",
               "r65_quiet_probe_log.json"]
REGISTRATION = ["r65_preregistration.json"]
REUSED = ["r48_interior_timing_selection.json",
          "r64_interior_timing_tighter.json", "r12_kernel_cost_curve.json",
          "r18_span_sweep_A_beta_1.00.json", "r18_span_sweep_B_beta_1.00.json"]

# Only records under these prefixes are claimed. Re-running a finalizer that
# globs a shared tree has taken ownership of a later campaign's records once
# before; the prefix test is what stops it.
OWNED_PREFIX = "r65_"


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
    p = METRICS / "r65_timing_family.json"
    if not p.exists():
        return {"missing": True}
    d = json.loads(p.read_text(encoding="utf-8"))
    return {"by_k": d["by_k"], "protocol": d.get("protocol"),
            "timing_band": d.get("timing_band")}


def main() -> int:
    payload = {
        "schema": "r65_final_experiment_manifest_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "scope": ("R65 (O58): the sampled interior family, k = 0.25, 0.50 and "
                  "0.75, each against its own constant degree matched on "
                  "measured kernel time at the tighter tolerance, on (O48)'s "
                  "fourteen-orbit panel at beta = 1."),
        "why": ("(O48) and (O57) both time k = 0.5, which was adopted as the "
                "most frequent sampled minimum of the nominal sweep and not "
                "because the family has an optimum there. A negative "
                "measured-time result for that member is not one for interior "
                "allocation, and this campaign asks the family-level "
                "question."),
        "numbering_note": (
            "the records carry the prefix r65 because r58 is held by the "
            "equal-work endpoint control; the observation label is (O58)."),
        "relationship_to_r64": (
            "same panel, same budget, same scoring rule, same timing band. "
            "k = 0.25 and k = 0.75 are propagated here; the k = 0.50 column "
            "is reused from r64 unchanged and is marked as reused in the "
            "result record, so its two out-of-band cells are r64's and are "
            "reported rather than re-run. R48's and R64's records are read, "
            "never written."),
        "refinement_contract": (
            "each comparator degree is refined over integers until the "
            "measured time ratio enters 0.90-1.10 or the integer step cannot "
            "improve it, where r64 took a single step. Any cell still outside "
            "the band is flagged in the record rather than absorbed."),
        "no_oracle": (
            "all three members are reported. No per-orbit minimum over k is "
            "taken and no deployable k is claimed from this campaign; the "
            "reporting rule was fixed in the registration before propagation."),
        "discarded_attempt": {
            "status": "ran to completion, discarded, superseded by the "
                      "campaign indexed here; its outcome is reported, not "
                      "withheld",
            "outcome_it_returned": ("4-5, 1-5 and 0-9, which is the same "
                                    "registered class B the accepted run "
                                    "returned"),
            "reason": ("other work joined the machine mid-run, and the "
                       "campaign matches its comparator on measured kernel "
                       "time. Of the 28 freshly propagated cells the "
                       "comparator degree differs from the accepted run's on "
                       "19, and four cells change verdict: the selected "
                       "degrees are not the ones an idle machine selected, "
                       "and that degree is the quantity the campaign is built "
                       "on. Every cell still reported an in-band ratio, "
                       "because the ratio divides the two measurements by "
                       "each other and a load common to both cancels."),
            "criterion_was_set_after_the_run": (
                "the registration fixes a mechanism, an idle check at the "
                "start of the timed pipeline, not a continuing condition and "
                "not a numeric threshold; the thresholds were set after the "
                "attempt's data was on disk. The cell-level divergence above "
                "is what carries the discard and needs no threshold."),
            "evidence_record": "metrics/r65_discarded_attempt1.json",
            "archive_dir": "metrics/r65_discarded_attempt1",
            "post_hoc_quiet_audit": (
                "metrics/r65_quiet_audit.json applies the same diagnostic to "
                "the accepted run from telemetry the case records already "
                "carry. Degree-normalized, the accepted run has a median "
                "within-orbit throughput spread of 1.03 and no case above "
                "1.25x its opening baseline, against 1.36 and thirteen for "
                "the discarded attempt; the degree model leaves 3.7 per cent "
                "residual on the accepted run and 32.5 per cent on the "
                "discarded one. The raw figures in the discard record are "
                "larger because the raw measure is not flat in the degree."),
            "guard_added": ("python_codes/probe_kernel_quiet.py, run before "
                            "the campaign indexed here and reproducing the "
                            "archived R12 idle-machine cost curve to within "
                            "1 per cent. That run left no artifact; the probe "
                            "now appends to metrics/r65_quiet_probe_log.json, "
                            "so the precondition is auditable from the next "
                            "campaign onward and not for this one."),
        },
        "declared_outcome_returned": outcome(),
        "numerical_kernel": {
            "lunaris_release_tag": "paper-truncation-v1.0",
            "lunaris_commit": "27e9ab86ed61d623f78c453ea2054348f1044c23",
            "note": "unchanged from R10-R64"},
        "registration": index_files(REGISTRATION, METRICS),
        "reused_inputs": index_files(REUSED, METRICS),
        "scripts": index_files(SCRIPTS, CODE),
        "result_json": index_files(RESULT_JSON, METRICS),
        "generated_tables": index_files(TABLES, METRICS),
        "trajectory_tree": index_tree("r65_cases", "r65_raw"),
        "discarded_tree": index_tree("r65_discarded_attempt1/cases",
                                     "r65_discarded_attempt1/raw"),
    }
    body = {k: v for k, v in payload.items() if k != "manifest_sha256"}
    payload["manifest_sha256"] = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    t, dt = payload["trajectory_tree"], payload["discarded_tree"]
    print(f"[written] {OUT.name}  manifest_sha256="
          f"{payload['manifest_sha256'][:16]}")
    print(f"  live:      {t['n_sidecars']} sidecars, {t['n_raw_arrays']} raw")
    print(f"  discarded: {dt['n_sidecars']} sidecars, {dt['n_raw_arrays']} raw")
    missing = [k for sec in ("scripts", "result_json", "reused_inputs",
                             "registration")
               for k, v in payload[sec].items() if v.get("missing")]
    if missing:
        print("[error] recorded as missing: " + ", ".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
