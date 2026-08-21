"""R51: measure the degree demand of the span ladder before paying for a control.

R38 established the rule this follows: decide the reference degree of a
cap-lifted control by measuring the demand, not by preference, and measure it
before any propagation. The same question arises for R50, and more sharply,
because the widest levels of the span ladder produce the largest error ratios in
the campaign and the ceiling is the first thing a referee will suspect of
producing them.

The probe re-runs the calibration bisection of rev14_budget_pareto with the cap
raised from 300 to 600, on the archived R50 altitude histories. Nothing is
propagated and no trajectory is computed. The only file it writes is its own
record, r51_ceiling_probe.json, which the freeze script requires: a control
whose reference degree is not backed by a measurement on disk should not be
possible to register.

Two numbers decide whether R51 is worth eight hours:

  fraction_at_cap at 300   how much of the calibrated radial schedule is
                           clamped today. If this is zero at the wide levels,
                           the ceiling cannot be what produces their ratios and
                           the control is unnecessary.
  max requested degree     what the schedule asks for once the cap is 600. If
                           it clears 600, then 600 is the right reference for
                           the control, as it was for R38.

Both are reported per apolune level, because the whole point of R50 is that the
levels differ.

Usage:
    python rev51_ceiling_probe.py --workers 2
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

import rev10_sobol_confirmatory as base                       # noqa: E402
import rev14_budget_pareto as pareto                          # noqa: E402

# The ladder has two identity blocks and the ceiling has to be located on the
# one being controlled: the identities differ between them, so block A's probe
# is not evidence about block B. --block selects which parent is read; nothing
# under metrics/r50_* is ever written by this script either way.
BLOCKS = {"a": ("span_ladder_a", "RS1"), "b": ("span_ladder_b", "RS2")}
BLOCK, PARENT_KEY = BLOCKS["a"]
PROBE_KEY = PARENT_KEY     # the same trees are read; nothing is written
BETAS = [1.00, 0.75, 0.50]
RAISED = 600
PARENT = 300


BLOCK_LETTER = "a"


def bind_block(letter: str) -> None:
    global BLOCK, PARENT_KEY, PROBE_KEY, BLOCK_LETTER
    BLOCK_LETTER = letter
    BLOCK, PARENT_KEY = BLOCKS[letter]
    PROBE_KEY = PARENT_KEY
    pareto.DESIGNS[PARENT_KEY] = {
        "rows": METRICS / f"r50_{BLOCK}_rows.json",
        "r12_case": METRICS / "r50_cases" / f"atallah_{BLOCK}",
        "r11_raw": METRICS / "r11_raw" / f"stratum_{BLOCK}_convergence",
    }


bind_block("a")


def probe_worker(task: dict) -> dict:
    # ProcessPoolExecutor spawns on Windows: the child re-imports this module
    # and gets the import-time binding, not the one main() made. Block A's
    # first run happened to be the import-time default and so passed; block B
    # came back with every record incomplete and no cells at all. The block
    # travels with the task so the child binds it before doing any work.
    bind_block(task.pop("block", BLOCK_LETTER))
    return pareto.worker(task)


def levels() -> dict:
    d = json.loads((METRICS / f"r50_{BLOCK}_design_frozen.json")
                   .read_text(encoding="utf-8"))
    return {o["sobol_index"]: o["apolune_level_km"] for o in d["orbits"]}


def summarize(records: list[dict], level_of: dict, label: str) -> dict:
    cells: dict[float, dict] = {}
    for rec in records:
        if rec["status"] != "complete":
            continue
        lvl = level_of[rec["sobol_index"]]
        c = cells.setdefault(lvl, {"orbits": 0, "max_degree": [],
                                   "fraction_at_cap": [], "at_cap": 0})
        c["orbits"] += 1
        worst_frac = 0.0
        worst_deg = 0
        for key, entry in rec["budgets"].items():
            if key == "original":
                continue
            alloc = entry["atallah"]["allocation"]
            worst_frac = max(worst_frac, float(alloc["fraction_at_cap"]))
            worst_deg = max(worst_deg, int(alloc["max_degree"]))
        c["max_degree"].append(worst_deg)
        c["fraction_at_cap"].append(worst_frac)
        if worst_frac > 0.0:
            c["at_cap"] += 1
    print(f"\n{label}")
    print("  ha_km   orbits  max requested  median frac at cap  orbits touching")
    out = {}
    for lvl in sorted(cells):
        c = cells[lvl]
        out[lvl] = {"orbits": c["orbits"], "max_degree": max(c["max_degree"]),
                    "median_fraction_at_cap": median(c["fraction_at_cap"]),
                    "orbits_touching_cap": c["at_cap"]}
        print(f"  {lvl:6.0f}   {c['orbits']:4d}    {max(c['max_degree']):9d}"
              f"      {median(c['fraction_at_cap']):14.4f}"
              f"      {c['at_cap']:3d}/{c['orbits']}")
    return out


def run(cap: int, workers: int, level_of: dict,
        betas: list[float]) -> list[dict]:
    rows = json.loads((METRICS / f"r50_{BLOCK}_rows.json")
                      .read_text(encoding="utf-8"))["rows"]
    tasks = []
    for r in rows:
        q = dict(r)
        q["adopted_truth_degree"] = cap
        tasks.append({"design": PARENT_KEY, "row": q, "betas": betas,
                      "block": BLOCK_LETTER})
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(probe_worker, tasks))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    # The probe that gated the freeze covered the three budgets the control was
    # going to be scored at. The ladder itself runs six, so the demand question
    # can be asked at all six for the price of a bisection; that wider sweep is
    # written to its own file rather than over the record the freeze read.
    ap.add_argument("--betas", type=float, nargs="+", default=BETAS)
    ap.add_argument("--block", choices=sorted(BLOCKS), default="a")
    ap.add_argument("--out", default="r51_ceiling_probe.json")
    a = ap.parse_args()
    bind_block(a.block)
    betas = sorted(a.betas, reverse=True)

    level_of = levels()
    # Hash the two files this probe reads, not every r50_* file. The first run
    # of this probe globbed r50_* by mtime and reported the parent records as
    # touched, because the campaign was writing block B's rows and operating
    # point in another process at the time. A provenance check has to name what
    # it is checking, or it reports the neighbours' work as its own side effect.
    watched = [METRICS / f"r50_{BLOCK}_design_frozen.json",
               METRICS / f"r50_{BLOCK}_rows.json"]
    before = {p.name: base.file_hash(p) for p in watched}

    at_parent = run(PARENT, a.workers, level_of, betas)
    parent = summarize(at_parent, level_of,
                       f"cap {PARENT} (what the campaign ran)")
    at_raised = run(RAISED, a.workers, level_of, betas)
    raised = summarize(at_raised, level_of, f"cap {RAISED} (the control)")

    after = {p.name: base.file_hash(p) for p in watched}
    untouched = before == after
    print(f"\nparent inputs unchanged ({', '.join(before)}): {untouched}")

    worst = max(v["max_degree"] for v in raised.values())
    binding = {lvl: v["orbits_touching_cap"] for lvl, v in parent.items()}
    print(f"\nmax requested degree at cap {RAISED}: {worst}")
    print(f"orbits clamped at cap {PARENT}, by level: {binding}")
    if worst >= RAISED:
        verdict = (f"{RAISED} does not clear the demand; a control at "
                   f"{RAISED} would still be capped")
    elif not any(binding.values()):
        verdict = (f"the cap never binds at {PARENT}: the ratios are not a "
                   f"ceiling effect and the control buys nothing")
    else:
        verdict = (f"{RAISED} clears the demand and the cap does bind today: "
                   f"the control is both meaningful and correctly referenced")
    print(f"  -> {verdict}")

    record = {
        "schema": "r51_ceiling_probe_v1",
        "created_utc": base.utc_now(),
        "population": BLOCK,
        "design_key": PARENT_KEY,
        "method": ("the calibration bisection of rev14_budget_pareto re-run at "
                   f"cap {RAISED} on the archived R50 altitude histories; no "
                   "propagation, no new trajectory, and nothing in "
                   "metrics/r50_* modified"),
        "budgets_probed": betas,
        "parent_reference_degree": PARENT,
        "raised_reference_degree": RAISED,
        "by_level_at_parent_cap": {f"{k:.0f}": v for k, v in parent.items()},
        "by_level_at_raised_cap": {f"{k:.0f}": v for k, v in raised.items()},
        "max_requested_degree_at_raised_cap": worst,
        "orbits_clamped_at_parent_cap_by_level": {f"{k:.0f}": v
                                                  for k, v in binding.items()},
        "parent_inputs_watched": sorted(before),
        "parent_inputs_sha256": after,
        "parent_inputs_unchanged": untouched,
        "conclusion": verdict,
    }
    out = METRICS / a.out
    out.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(f"[written] {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
