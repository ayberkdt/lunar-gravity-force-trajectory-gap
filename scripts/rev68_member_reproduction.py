"""Does R68's serially re-propagated member reproduce the archived member?

The campaign measures the member's kernel time on a fresh serial run but takes
the member's *error* from the archived campaign, because the archived run was
made with concurrent workers and only its timing is unusable. That reuse is
legitimate only if the fresh run is the same trajectory. The dynamics are
deterministic, so it should be; this file checks it, because the one defect
that would not show up anywhere else -- a mis-read degree table, the endpoint
arm's table living under a different key from the interior arm's -- shows up
here immediately.

Every scored cell is re-scored: the fresh member trajectory against the same
reference, compared with the error the archived record carries.

Usage:  python rev68_member_reproduction.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import rev10_sobol_confirmatory as base      # noqa: E402
import rev14_budget_trajectory as r14        # noqa: E402
import rev68_timing_full as r68              # noqa: E402

ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"
OUT = METRICS / "r68_member_reproduction.json"


def load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def main() -> int:
    payload = {
        "schema": "r68_member_reproduction_v1",
        "created_utc": base.utc_now(),
        "question": ("the member's kernel time is measured on a fresh serial "
                     "run and its error is reused from the archived campaign; "
                     "this checks that the two are the same trajectory"),
        "arms": {},
    }
    for arm in sorted(r68.ARMS):
        p = r68.out_path(arm)
        if not p.exists():
            continue
        rec = json.loads(p.read_text(encoding="utf-8"))
        rows, missing = [], 0
        for r in rec["rows"]:
            design, idx = r["design"], r["sobol_index"]
            raw = r68.raw_dir(arm, design, idx) / "member_tighter.npz"
            _, tr = r14.reuse_paths(design, idx, "truth", "tighter")
            if not (raw.exists() and tr.exists()):
                missing += 1
                continue
            gt, gy = load(raw)
            tt, ty = load(tr)
            err = base.common_error(gt, gy, tt, ty)["pos_rms_m"]
            arch = r["member_error_m"]
            rows.append({"design": design, "sobol_index": idx,
                         "archived_error_m": arch, "fresh_error_m": err,
                         "relative_difference":
                             abs(err - arch) / arch if arch else None})
        rel = sorted(x["relative_difference"] for x in rows
                     if x["relative_difference"] is not None)
        payload["arms"][arm] = {
            "member_k": rec["member_k"],
            "cells_checked": len(rows),
            "cells_without_a_raw_pair": missing,
            "median_relative_difference":
                float(np.median(rel)) if rel else None,
            "worst_relative_difference": rel[-1] if rel else None,
            "exact_matches": sum(1 for x in rel if x == 0.0),
            "rows": rows}
        print(f"{arm}: {len(rows)} cells, exact on "
              f"{payload['arms'][arm]['exact_matches']}, worst relative "
              f"{rel[-1]:.2e}" if rel else f"{arm}: nothing to check")
    base.atomic_json(OUT, payload)
    print(f"[written] {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
