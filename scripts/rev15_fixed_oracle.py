"""Budget-saturating and best-under-budget fixed comparators (R15-A).

O25 compares the radial rule with the integer degree nearest the budget in N^2.
That is a legitimate, deployable comparator, but it is not the fixed-degree
Pareto envelope: this paper's own dense sweep shows the seven-day error is not
monotone in N, with order-of-magnitude swings between neighboring degrees from
cancellation that does not transfer across geometry or phase. "The constant
degree at this budget" and "the best constant degree available at this budget"
are therefore different comparators, and only the second can support a claim
about the fixed-degree frontier.

This measures all three at beta = 1, where the budget is B = N_crit^2:

  F_near   argmin_N |N^2 - B|          -- the O25 comparator
  F_sat    max{N : N^2 <= B}           -- the budget-saturating degree, which
                                          never overspends
  F_oracle argmin_{N^2 <= B} E_traj(N) -- post-hoc, over a fine ladder below the
                                          budget; a lower envelope, not a policy

F_oracle can only improve the fixed side, so it can only strengthen a conclusion
that already favors the constant degree. It is reported as a post-hoc oracle and
never as something a user could select in advance.

Usage:
    python rev15_fixed_oracle.py run --orbits 16 --workers 5
"""

from __future__ import annotations

import argparse
import json
import math
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
from rev14_budget_pareto import DESIGNS, _model
from rev14_budget_trajectory import LEVELS, MAX_STEP, DURATION, OUTPUT_STEP

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RAW_ROOT = METRICS / "r15_raw" / "fixed_oracle"
BETA = 1.00
# offsets below the budget-saturating degree; fine near the top because the
# non-monotonicity the sweep found lives at the few-degree scale
OFFSETS = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24]


def worker(task: dict) -> dict:
    design, row, prov = task["design"], task["row"], task["provenance"]
    index = int(row["sobol_index"])
    try:
        adopted = int(row["adopted_truth_degree"])
        n_crit = int(row["n_critical"])
        y0 = np.asarray(row["design_point"]["initial_state_si"], dtype=float)
        model, args = _model(adopted)
        grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
        B = float(n_crit ** 2)

        n_near = min([int(math.floor(math.sqrt(B))), int(math.ceil(math.sqrt(B)))],
                     key=lambda k: abs(k ** 2 - B))
        n_sat = int(math.floor(math.sqrt(B)))
        cands = sorted({max(2, n_sat - d) for d in OFFSETS})

        t_ref, y_ref = base.load_raw(DESIGNS[design]["r11_raw"]
                                     / f"sobolA_{index:03d}" / "truth_tight.npz")
        ladder = {}
        for n in cands:
            t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, lambda _t, _h, k=n: k, args,
                LEVELS["tight"]["rtol"], LEVELS["tight"]["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "where": f"fixed_{n}", "message": fail}
            raw = RAW_ROOT / design / f"sobolA_{index:03d}" / f"fixed_{n}_tight.npz"
            base.atomic_npz(raw, t_s=t, state_si=y)
            e = base.common_error(t, y, t_ref, y_ref)["pos_rms_m"]
            ladder[n] = {"error_tight_m": e,
                         "quadratic_work_per_call": float(n ** 2),
                         "beta": (n ** 2) / B,
                         "n_rhs": int(tel["n_rhs"]),
                         "total_quadratic_work": float(n ** 2) * int(tel["n_rhs"])}
        n_oracle = min(ladder, key=lambda k: ladder[k]["error_tight_m"])
        return {"index": index, "design": design, "status": "complete",
                "n_critical": n_crit, "adopted_truth_degree": adopted,
                "hp_km": float(row["design_point"]["hp_km"]),
                "n_near": n_near, "n_sat": n_sat, "n_oracle": n_oracle,
                "ladder": ladder, "provenance_sha": prov.get("kernel_sha256")}
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


def stat(v):
    a = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)), "p90": float(np.percentile(a, 90)),
            "min": float(a.min()), "max": float(a.max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run",))
    ap.add_argument("--orbits", type=int, default=16)
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    prov = base.provenance()
    tasks = []
    for design in ("A", "B"):
        rows = json.loads(DESIGNS[design]["rows"].read_text())["rows"]
        rows.sort(key=lambda r: r["design_point"]["hp_km"])
        picks = [rows[int(i)] for i in
                 np.linspace(0, len(rows) - 1, a.orbits // 2).round()]
        tasks += [{"design": design, "row": r, "provenance": prov} for r in picks]
    print(f"[fixed-oracle] {len(tasks)} orbits x {len(OFFSETS)} degrees at "
          f"beta={BETA}", flush=True)
    t0 = time.time()
    done, fails = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(worker, t): t for t in tasks}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                fails.append(rec)
                print(f"  !! {rec['index']:03d} {rec.get('message')}", flush=True)
                continue
            done.append(rec)
            l = rec["ladder"]
            print(f"  [{n}/{len(tasks)}] {rec['design']}{rec['index']:03d} "
                  f"N_sat={rec['n_sat']} N_oracle={rec['n_oracle']} "
                  f"gain={l[rec['n_sat']]['error_tight_m'] / l[rec['n_oracle']]['error_tight_m']:.2f}x "
                  f"elapsed={(time.time()-t0)/60:.1f}min", flush=True)
    done.sort(key=lambda r: (r["design"], r["index"]))

    # how much does the best-under-budget degree beat the saturating one?
    gains, shifts, near_eq_sat = [], [], 0
    for r in done:
        l = r["ladder"]
        gains.append(l[r["n_sat"]]["error_tight_m"] / l[r["n_oracle"]]["error_tight_m"])
        shifts.append(r["n_sat"] - r["n_oracle"])
        near_eq_sat += int(r["n_near"] == r["n_sat"])
    # compare each fixed variant with the radial policy already propagated at beta=1
    traj = {d: json.loads((METRICS / f"r14_trajectory_{d}_beta_1.00.json"
                           ).read_text())["rows"] for d in ("A", "B")}
    rows_out = []
    for r in done:
        ref = next((x for x in traj[r["design"]]
                    if int(x["sobol_index"]) == r["index"]), None)
        if ref is None:
            continue
        e_at = ref["comparison"]["atallah_error_m"]
        l = r["ladder"]
        rows_out.append({
            "design": r["design"], "sobol_index": r["index"], "hp_km": r["hp_km"],
            "n_critical": r["n_critical"], "n_near": r["n_near"],
            "n_sat": r["n_sat"], "n_oracle": r["n_oracle"],
            "atallah_error_m": e_at,
            "error_near_m": l.get(r["n_near"], {}).get("error_tight_m"),
            "error_sat_m": l[r["n_sat"]]["error_tight_m"],
            "error_oracle_m": l[r["n_oracle"]]["error_tight_m"],
            "rho_vs_sat": l[r["n_sat"]]["error_tight_m"] / e_at if e_at else None,
            "rho_vs_oracle": l[r["n_oracle"]]["error_tight_m"] / e_at if e_at else None,
            "oracle_gain_over_sat": (l[r["n_sat"]]["error_tight_m"]
                                     / l[r["n_oracle"]]["error_tight_m"]),
            "ladder": {str(k): v["error_tight_m"] for k, v in sorted(l.items())}})
    summary = {
        "orbits": len(rows_out),
        "n_near_equals_n_sat": near_eq_sat,
        "oracle_gain_over_sat": stat(gains),
        "oracle_degree_shift_below_sat": stat([float(s) for s in shifts]),
        "orbits_where_oracle_is_sat": int(sum(1 for s in shifts if s == 0)),
        "rho_vs_sat": stat([r["rho_vs_sat"] for r in rows_out]),
        "rho_vs_oracle": stat([r["rho_vs_oracle"] for r in rows_out]),
        "fixed_beats_radial_vs_sat": int(sum(
            1 for r in rows_out if r["error_sat_m"] < r["atallah_error_m"])),
        "fixed_beats_radial_vs_oracle": int(sum(
            1 for r in rows_out if r["error_oracle_m"] < r["atallah_error_m"])),
    }
    payload = {"schema": "r15_fixed_oracle_v1", "created_utc": base.utc_now(),
               "beta": BETA, "offsets_below_saturating": OFFSETS,
               "note": ("F_oracle is a post-hoc lower envelope over the fixed "
                        "family under the budget, not a selectable policy"),
               "rows": rows_out, "raw": done, "failures": fails,
               "summary": summary, "source": prov}
    base.atomic_json(METRICS / "r15_fixed_oracle.json", payload)
    s = summary
    print(f"[fixed-oracle] orbits={s['orbits']}  N_near==N_sat on "
          f"{s['n_near_equals_n_sat']}  oracle==saturating on "
          f"{s['orbits_where_oracle_is_sat']}")
    print(f"  oracle gain over saturating: median {s['oracle_gain_over_sat']['median']:.2f}x "
          f"(max {s['oracle_gain_over_sat']['max']:.2f}x)")
    print(f"  radial vs saturating: rho median {s['rho_vs_sat']['median']:.3g}, "
          f"fixed better on {s['fixed_beats_radial_vs_sat']}/{s['orbits']}")
    print(f"  radial vs oracle:     rho median {s['rho_vs_oracle']['median']:.3g}, "
          f"fixed better on {s['fixed_beats_radial_vs_oracle']}/{s['orbits']}")
    print("[written] r15_fixed_oracle.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
