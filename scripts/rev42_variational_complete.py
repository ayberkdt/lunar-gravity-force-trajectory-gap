"""R42: finish the R37 level chain from where its wall clock stopped it.

R37 left 102 orbits in the record and a completed panel of 80. The 26 orbits of
levels 56 and 64 that it never submitted are the whole of this run. Nothing here
is a new method and almost nothing here is new code: the selection rule, the
level construction, the per-orbit post-processing and the scoring are imported
from rev37_variational_extend, which imports the solve itself from
rev14_variational_budget. An orbit solved here and an orbit solved by R37 are
the same computation called from a different __main__.

What this file adds is the one thing R37's driver cannot do: continue. That
driver rebuilds its row set from R14's archived eight every time it starts, so
re-running it to reach the last 26 orbits would recompute the 94 already solved
-- about a hundred core-hours to reproduce a sealed file. This one carries the
sealed record forward instead, and pays for the right to do so with two checks
made before the first new orbit is submitted:

  digest      metrics/r37_variational_extension.json must still match the
              sha256 sealed in the R37 manifest, byte for byte. It is
              re-checked at the end, so a concurrent writer would be caught
              rather than silently mixed in.

  recompute   R37's own admissibility orbit B005 is solved again and its
              predicted ratio compared with the archived value under the same
              1e-3 threshold. The digest says the carried file is the sealed
              one; this says the current source still computes what produced
              it.

The carried rows are never re-scored and never re-selected. They are written out
with the provenance flag they arrived with, plus carried_from_r37, so the record
this run produces says exactly which of its rows it computed.

The sealed R37 record is read and never written. This run writes only
metrics/r42_variational_completion.json.

Usage:
    python rev42_variational_complete.py --workers 4 --deadline 2026-08-07T08:00
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import rev10_sobol_confirmatory as base
import rev14_variational_budget as vb
import rev37_variational_extend as ext
from rev13_variational_check import GRADIENT_DEGREE

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PREREG = METRICS / "r42_preregistration.json"
R37_PREREG = METRICS / "r37_preregistration.json"
R37_MANIFEST = METRICS / "r37_final_experiment_manifest.json"
R37_RECORD = METRICS / "r37_variational_extension.json"
OUTPUT = METRICS / "r42_variational_completion.json"


def sealed_digest() -> str:
    man = json.loads(R37_MANIFEST.read_text(encoding="utf-8"))
    return man["result_json"]["r37_variational_extension.json"]["sha256"]


def level_membership(pareto, levels):
    """R37's level bookkeeping, rebuilt from R37's own selection function."""
    counts = {"members": {}}
    for n_per in levels:
        members = set()
        for design in ("A", "B"):
            for p in ext.selection(pareto, design, n_per):
                members.add((design, int(p["sobol_index"])))
        counts["members"][n_per] = members
        counts[n_per] = len(members)
    return counts


def checkpoint(rows, levels, level_counts, meta):
    payload = {"schema": "r42_variational_completion_v1",
               "created_utc": base.utc_now(),
               "beta": vb.BETA, "gradient_degree": GRADIENT_DEGREE,
               **meta,
               "rows": sorted(rows, key=lambda r: (r["design"],
                                                   r["sobol_index"])),
               "summary": ext.summarize(rows, levels, level_counts)}
    base.atomic_json(OUTPUT, payload)
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--deadline", required=True,
                    help="local ISO time after which no new orbit is submitted")
    a = ap.parse_args()

    prereg = json.loads(PREREG.read_text(encoding="utf-8"))
    levels = json.loads(R37_PREREG.read_text(encoding="utf-8")
                        )["selection_rule"]["levels_per_design"]
    deadline = datetime.fromisoformat(a.deadline).timestamp()

    # ---- admissibility, part one: the carried record is the sealed record
    sealed = sealed_digest()
    on_disk = base.file_hash(R37_RECORD)
    if on_disk != sealed:
        print(f"[r42] ABORT: {R37_RECORD.name} is not the sealed record "
              f"({on_disk[:16]} vs {sealed[:16]}); nothing written.")
        return 1
    record = json.loads(R37_RECORD.read_text(encoding="utf-8"))
    carried = [dict(r, carried_from_r37=True) for r in record["rows"]]
    have = {(r["design"], r["sobol_index"]) for r in carried}
    print(f"[r42] carrying {len(carried)} rows from {R37_RECORD.name} "
          f"(digest {sealed[:16]} verified)", flush=True)

    # ---- resume: orbits an earlier process of this same campaign solved.
    # Only its own output qualifies, and only if it carried the same sealed
    # record; otherwise the two processes are not the same campaign.
    previous_runs, resumed = [], []
    if OUTPUT.exists():
        prev = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if prev.get("carried", {}).get("sha256") == sealed:
            previous_runs = prev.get("runs", [])
            resumed = [dict(r, carried_from_earlier_r42_process=True)
                       for r in prev["rows"]
                       if (r["design"], r["sobol_index"]) not in have]
            carried += resumed
            have |= {(r["design"], r["sobol_index"]) for r in resumed}
            print(f"[r42] resuming: {len(resumed)} orbits already solved by an "
                  f"earlier process of this campaign", flush=True)
        else:
            print(f"[r42] ABORT: {OUTPUT.name} exists but carries a different "
                  f"sealed record; refusing to mix two campaigns.")
            return 1

    pareto = json.loads(vb.PARETO.read_text(encoding="utf-8"))
    rows_by_design = {
        d: {int(r["sobol_index"]): r
            for r in json.loads(vb.ROWS[d].read_text())["rows"]}
        for d in ("A", "B")}

    entries = ext.build_levels(pareto, levels)
    level_counts = level_membership(pareto, levels)
    pending = [e for e in entries
               if (e["design"], e["sobol_index"]) not in have]
    print(f"[r42] {len(entries)} orbits declared over levels {levels}; "
          f"{len(pending)} to solve, deadline {a.deadline}, "
          f"{a.workers} workers", flush=True)
    if not pending:
        print("[r42] nothing to do; the chain is already complete.")
        return 0

    # ---- admissibility, part two: the source still computes what it did
    check_name = prereg["admissibility_self_check"]["recomputed_orbit"]
    cd, ci = check_name[0], int(check_name[1:])
    old = float(prereg["admissibility_self_check"]["archived_predicted_ratio"])
    spec = next(e["spec"] for e in entries
                if e["design"] == cd and e["sobol_index"] == ci)
    print(f"[r42] self-check: recomputing {check_name}", flush=True)
    t0 = time.time()
    got = vb.worker(ext.make_task({"design": cd, "sobol_index": ci,
                                   "spec": spec}, rows_by_design))
    if got["status"] != "complete":
        print(f"[r42] ABORT: self-check orbit failed: {got.get('message')}")
        return 1
    new = got["predicted_ratio_fixed_over_atallah"]
    rel = abs(new - old) / abs(old)
    thr = float(prereg["admissibility_self_check"]["abort_threshold_rel"])
    print(f"[r42] self-check {check_name}: archived {old:.10g}, recomputed "
          f"{new:.10g}, rel {rel:.2e} (threshold {thr:g}), "
          f"{time.time()-t0:.0f}s", flush=True)
    if rel > thr:
        print("[r42] ABORT: the carried record and the current source no "
              "longer describe the same computation; nothing written.")
        return 1

    meta = {"registration_sha256": prereg.get("preregistration_sha256"),
            "parent_registration_sha256": json.loads(
                R37_PREREG.read_text(encoding="utf-8")
            )["preregistration_sha256"],
            "deadline_local": a.deadline,
            "workers": a.workers,
            "levels_per_design": levels,
            "carried": {"file": R37_RECORD.name, "sha256": sealed,
                        "rows": len(carried),
                        "from_earlier_r42_process": len(resumed)},
            "runs": previous_runs + [
                {"started_local": datetime.now().isoformat(timespec="seconds"),
                 "workers": a.workers, "deadline_local": a.deadline,
                 "carried_from_earlier_process": len(resumed)}],
            "restart_note": (
                "a process of this campaign may be stopped and started again; "
                "when it is, it carries the orbits the earlier process solved "
                "rather than resolving them, and every process is listed under "
                "runs with the worker count it ran at. Which orbits are "
                "attempted, and in what order, does not depend on any of it."),
            "self_check": {"orbit": check_name, "archived": old,
                           "recomputed": new, "rel": rel, "threshold": thr,
                           "passed": True},
            "provenance_note": (
                "rows carrying carried_from_r37 were computed by R37 and are "
                "reproduced here unchanged; rows without it were solved by "
                "this run with the same imported worker. The eight rows that "
                "also carry reused_from_r14 came to R37 from R14."),
            "completed_levels_note": (
                "highest_complete_level_orbits is the panel; any orbits beyond "
                "it belong to a level the wall clock truncated and are "
                "reported as such.")}

    rows = list(carried)
    checkpoint(rows, levels, level_counts, meta)

    t0 = time.time()
    submitted, completed, skipped = 0, 0, 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {}
        it = iter(pending)
        for _ in range(min(a.workers, len(pending))):
            e = next(it, None)
            if e is None:
                break
            futs[pool.submit(vb.worker, ext.make_task(e, rows_by_design))] = e
            submitted += 1
        while futs:
            for fut in as_completed(list(futs)):
                e = futs.pop(fut)
                rec = fut.result()
                if rec["status"] == "complete":
                    rec = ext.finish(rec)
                    rec["level"] = e["level"]
                    rec["reused_from_r14"] = False
                    rec["carried_from_r37"] = False
                    rows.append(rec)
                    completed += 1
                    agr = ext.sign_agrees(rec)
                    print(f"  [{completed}/{len(pending)}] L{e['level']:02d} "
                          f"{e['design']}{e['sobol_index']:03d} "
                          f"hp={e['hp_km']:6.1f} "
                          f"pred={rec['predicted_ratio_fixed_over_atallah']:.4g} "
                          f"sign={'ok' if agr else 'DISAGREE' if agr is False else 'n/a'} "
                          f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
                else:
                    print(f"  !! {e['design']}{e['sobol_index']:03d} "
                          f"{rec.get('message')}", flush=True)
                checkpoint(rows, levels, level_counts, meta)
                if time.time() < deadline:
                    nxt = next(it, None)
                    if nxt is not None:
                        futs[pool.submit(vb.worker,
                                         ext.make_task(nxt, rows_by_design))] = nxt
                        submitted += 1
                else:
                    skipped = len(pending) - submitted
                break

    # ---- admissibility, part one again: nobody rewrote the carried record
    after = base.file_hash(R37_RECORD)
    meta["carried"]["sha256_after_run"] = after
    if after != sealed:
        print(f"[r42] WARNING: {R37_RECORD.name} changed during the run "
              f"({after[:16]} vs {sealed[:16]}); the carried rows are the ones "
              f"read at the start and the record says so.")
    payload = checkpoint(rows, levels, level_counts, meta)

    s = payload["summary"]
    print(f"\n[r42] {s['orbits']} orbits in the record "
          f"({completed} solved by this process, "
          f"{len(carried) - len(resumed)} carried from R37, "
          f"{len(resumed)} from an earlier process of this campaign, "
          f"{skipped} never submitted)")
    print(f"[r42] highest complete level: "
          f"{s['highest_complete_level_orbits']} orbits "
          f"(level {s['highest_complete_level_per_design']} per design)")
    print(f"[r42] sign agreement {s['sign_agreement']}/{s['sign_scored']}; "
          f"predicted favors radial on {s['predicted_favors_radial']}")
    if s["sign_disagreements"]:
        print("[r42] sign disagreements:")
        for d in s["sign_disagreements"]:
            print(f"    {d['design']}{d['sobol_index']:03d} hp={d['hp_km']:.1f} "
                  f"pred={d['predicted']:.4g} meas={d['measured']:.4g}")
    c = s["calibration_ratio_fixed"]
    if c:
        print(f"[r42] calibration channel median {c['median']:.4f} "
              f"[{c['min']:.3f}, {c['max']:.3f}], "
              f"{len(c['outside_archived_band'])} outside the 0.94-1.04 band")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
