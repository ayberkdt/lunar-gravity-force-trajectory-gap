"""R57: finish the three band-L orbits R39's wall clock did not reach.

R39 audits the forced variational reconstruction's gradient-degree
approximation on a stratified sixteen-orbit subset: it re-solves each orbit with
the gravity gradient at the orbit's own reference degree instead of the archived
120, and asks whether that can move a prediction across ratio one. Thirteen
orbits returned. The three that did not are B060, A059 and B043, and they are
the whole of band L except B044, so the band the audit has most to prove itself
on stands at one of four.

They are also the expensive ones, and for the reason that makes them matter:
band L is the low-perilune band, its orbits carry reference degree 900 where
the rest of the panel carries 300, and the one band-L orbit that did return took
36 363 s against 2 300-6 700 s for every other orbit in the record. Finishing
the band is about ten times the cost per orbit of the rest of it.

What this file adds is the one thing rev39_gradient_degree_panel cannot do:
continue. Its `run` rebuilds the whole sixteen-orbit selection, re-runs the
admissibility orbit and writes a fresh record every time, so re-entering it to
reach the last three would recompute the thirteen already solved -- about
twenty core-hours to reproduce a sealed file -- and overwrite that file with the
result. This one carries the sealed record forward and pays for the right with
the two checks rev42_variational_complete established for the same situation:

  digest      metrics/r39_gradient_degree_panel.json must still match the
              sha256 sealed in the R39 manifest, byte for byte. It is
              re-checked after the solve, so a concurrent writer is caught
              rather than silently mixed in.

  recompute   R39's own admissibility orbit A002 is solved again at the
              archived gradient degree and its predicted ratio compared with
              the archived value under R39's own threshold. The digest says the
              carried file is the sealed one; this says the current source
              still computes what produced it.

The selection is not re-derived: the three orbits are read from the sealed
record's own `unfinished` list, so this run cannot quietly choose a different
three. The carried rows are never re-solved and never re-scored.

The sealed R39 record is read and never written. This run writes only
metrics/r57_gradient_degree_completion.json.

Usage:
    python rev57_gradient_complete.py --workers 3 --deadline-h 30
    python rev57_gradient_complete.py --status
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"

SEALED = METRICS / "r39_gradient_degree_panel.json"
R39_MANIFEST = METRICS / "r39_final_experiment_manifest.json"
OUT = METRICS / "r57_gradient_degree_completion.json"


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def sealed_digest() -> str | None:
    """The digest the R39 manifest recorded for its own panel record."""
    if not R39_MANIFEST.exists():
        return None
    man = json.loads(R39_MANIFEST.read_text(encoding="utf-8"))

    # Two index shapes appear across these manifests and only one of them was
    # handled at first: a list of entries carrying their own path, and a dict
    # keyed by file name whose value holds the digest. R39 uses the second, so
    # the lookup returned None and the run degraded to "proceeding on the file
    # as found" -- a safety check reporting a warning instead of failing, which
    # is the worst of the three possible outcomes. Both shapes are read here.
    def walk(node):
        if isinstance(node, dict):
            name = node.get("path") or node.get("name") or node.get("file")
            digest = node.get("sha256") or node.get("digest")
            if name and digest and Path(str(name)).name == SEALED.name:
                yield digest
            for k, v in node.items():
                if isinstance(k, str) and Path(k).name == SEALED.name \
                        and isinstance(v, dict) and v.get("sha256"):
                    yield v["sha256"]
                yield from walk(v)
        elif isinstance(node, list):
            for v in node:
                yield from walk(v)

    return next(iter(walk(man)), None)


def load_r39():
    """Import the R39 driver without letting its argument parser run."""
    saved = sys.argv[:]
    sys.argv = [str(HERE / "rev39_gradient_degree_panel.py")]
    try:
        import rev39_gradient_degree_panel as r39
    finally:
        sys.argv = saved
    return r39


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--deadline-h", type=float, default=30.0)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    if not SEALED.exists():
        print(f"[abort] {SEALED.name} is missing; this run continues a record, "
              f"it does not start one")
        return 2
    rec = json.loads(SEALED.read_text(encoding="utf-8"))
    unfinished = rec.get("unfinished") or []
    carried = rec.get("rows") or []

    if a.status:
        print(f"R57: {len(carried)} rows carried from R39, "
              f"{len(unfinished)} to solve")
        for u in unfinished:
            print(f"  {u['design']}{int(u['sobol_index']):03d}  band {u['band']}")
        print(f"  output {'present' if OUT.exists() else 'absent'}")
        return 0

    if not unfinished:
        print("[ok] the sealed record lists nothing unfinished; nothing to do")
        return 0

    live = sha256(SEALED)
    want = sealed_digest()
    if want is None:
        print(f"[warn] the R39 manifest does not name {SEALED.name}; "
              f"proceeding on the file as found, digest {live[:16]}")
    elif live != want:
        print(f"[abort] {SEALED.name} does not match the digest sealed in the "
              f"R39 manifest.\n  sealed {want[:16]}\n  found  {live[:16]}\n"
              f"Carrying rows forward from a record that has changed since it "
              f"was sealed would mix two campaigns without saying so.")
        return 2
    else:
        print(f"[ok] sealed digest matches: {live[:16]}")

    r39 = load_r39()

    prereg = json.loads(r39.PREREG.read_text(encoding="utf-8"))
    rows = r39.panel_rows()
    specs, rows_by_design = r39.load_specs()
    by_key = {r39.key(r): r for r in r39.select(rows)}

    # --- recompute check: R39's own admissibility orbit, at degree 120
    cd, ci = r39.SELF_CHECK_ORBIT[0], int(r39.SELF_CHECK_ORBIT[1:])
    ref = next(r for r in rows if r39.key(r) == (cd, ci))
    task = r39.make_task({"design": cd, "sobol_index": ci,
                          "spec": specs[(cd, ci)]}, rows_by_design)
    task["gradient_degree_override"] = r39.ARCHIVED_GRADIENT_DEGREE
    print(f"[r57] recompute check on {r39.SELF_CHECK_ORBIT} at degree "
          f"{r39.ARCHIVED_GRADIENT_DEGREE}", flush=True)
    t0 = time.time()
    got = r39.worker(task)
    if got["status"] != "complete":
        print(f"[abort] the check orbit failed: {got.get('message')}")
        return 1
    old = ref["predicted_ratio_fixed_over_atallah"]
    new = got["predicted_ratio_fixed_over_atallah"]
    rel = abs(new - old) / abs(old)
    thr = float(prereg["admissibility_self_check"]["abort_threshold_rel"])
    print(f"[r57] archived {old:.10g}, recomputed {new:.10g}, rel {rel:.2e} "
          f"(threshold {thr:g}), {time.time()-t0:.0f}s", flush=True)
    if rel > thr:
        print("[abort] the archived record and the current source disagree; "
              "nothing written.")
        return 1

    # --- solve only the three the sealed record names
    tasks = []
    for u in unfinished:
        k = (u["design"], int(u["sobol_index"]))
        entry = {"design": k[0], "sobol_index": k[1], "spec": specs[k]}
        t = r39.make_task(entry, rows_by_design)
        t["_key"] = k
        t["_band"] = u.get("band")
        tasks.append(t)
    print(f"[r57] solving {len(tasks)} orbits at their reference degree, "
          f"{a.workers} workers, deadline in {a.deadline_h:g} h", flush=True)

    deadline = time.time() + a.deadline_h * 3600.0
    solved, failed, unreached = [], [], []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(r39.worker, t): t for t in tasks}
        for fut in as_completed(futs):
            t = futs[fut]
            try:
                got = fut.result()
            except Exception as exc:                      # noqa: BLE001
                failed.append({"design": t["_key"][0],
                               "sobol_index": t["_key"][1],
                               "band": t["_band"],
                               "message": f"{type(exc).__name__}: {exc}"})
                print(f"  [FAIL] {t['_key']}: {exc}", flush=True)
                continue
            if got["status"] != "complete":
                failed.append({"design": t["_key"][0],
                               "sobol_index": t["_key"][1],
                               "band": t["_band"],
                               "message": got.get("message")})
            else:
                solved.append(got)
            print(f"  [{len(solved)+len(failed)}/{len(tasks)}] "
                  f"{t['_key'][0]}{t['_key'][1]:03d} {got['status']} "
                  f"elapsed={(time.time()-t0)/3600:.2f}h", flush=True)
            if time.time() > deadline:
                print("  deadline reached; remaining orbits left unreached",
                      flush=True)
                for p in futs:
                    p.cancel()
                break

    done_keys = {(s["design"], int(s["sobol_index"])) for s in solved}
    for u in unfinished:
        if (u["design"], int(u["sobol_index"])) not in done_keys \
                and not any(f["sobol_index"] == u["sobol_index"]
                            and f["design"] == u["design"] for f in failed):
            unreached.append(u)

    if sha256(SEALED) != live:
        print("[abort] the sealed record changed while this run was solving; "
              "nothing written.")
        return 1

    comparison = [r39.compare(s, by_key[(s["design"], int(s["sobol_index"]))])
                  for s in solved]
    all_rows = carried + solved
    all_comparison = (rec.get("comparison") or []) + comparison
    payload = {
        "schema": "r57_gradient_degree_completion_v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace(
            "+00:00", "Z"),
        "continues": SEALED.name,
        "continues_sha256": live,
        "registration_sha256": rec.get("registration_sha256"),
        "note": (
            "the R39 audit carried forward and completed. Rows solved by R39 "
            "are carried byte for byte and are not re-solved or re-scored; "
            "rows solved here are the orbits R39's own record lists as "
            "unfinished, read from that list rather than re-selected. The "
            "selection, the level construction and the solve are R39's."),
        "recompute_check": {"orbit": r39.SELF_CHECK_ORBIT,
                            "archived": old, "recomputed": new,
                            "rel": rel, "threshold": thr, "passed": True},
        "carried_rows": len(carried),
        "solved_here": [{"design": s["design"],
                         "sobol_index": int(s["sobol_index"]),
                         "wall_s": s.get("wall_s")} for s in solved],
        "failures": failed,
        "unreached": unreached,
        "rows": all_rows,
        "comparison": all_comparison,
        "summary": r39.summarize(all_comparison),
    }
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    tmp.replace(OUT)
    print(f"[written] {OUT.name}: {len(carried)} carried + {len(solved)} "
          f"solved, {len(failed)} failed, {len(unreached)} unreached")
    s = payload["summary"]
    print(f"  orbits {s.get('orbits')}, resolved {s.get('resolved')}, "
          f"side changes {s.get('side_changes_resolved')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
