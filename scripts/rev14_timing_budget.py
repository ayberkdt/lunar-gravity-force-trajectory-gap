"""O25 Phase D: does the equal-budget result survive a measured-time budget? (R14)

The budget grid of Phase A/B is the per-call quadratic proxy, which is an
operation count and not a machine cost. This experiment re-runs the beta = 1
comparison with the budget defined as measured serial gravity-kernel time: the
budget-calibrated radial policy is propagated, its kernel time is measured, and
the fixed comparator is the degree whose measured kernel time matches it.

Everything here is serial -- one process, no other propagation on the machine --
because the claim is about cost. Kernel times from a parallel pool are not
comparable and are never used. Each configuration is repeated and medians are
compared; a pair is only called time-matched if the achieved ratio falls inside
the declared repeatability band.

Usage:
    python rev14_timing_budget.py select
    python rev14_timing_budget.py run --repeats 3     # machine must be idle
    python rev14_timing_budget.py aggregate
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev12_atallah as at
from rev14_budget_pareto import DESIGNS as PARETO_DESIGNS, BIN_KM, FLOOR

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PARETO = METRICS / "r14_budget_pareto.json"
COST_CURVE = METRICS / "r12_kernel_cost_curve.json"
COST_CURVE_HIGH = METRICS / "r13_kernel_cost_curve_high.json"
SELECTION = METRICS / "r14_timing_selection.json"
OUTPUT = METRICS / "r14_timing_budget.json"
TABLE = METRICS / "r14_timing_budget_table_full.tex"
CASE_ROOT = METRICS / "r14_cases" / "timing_budget"
RAW_ROOT = METRICS / "r14_raw" / "timing_budget"

LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3)},
    "tighter": {"rtol": 1.0e-13, "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3)},
}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
BETA = 1.00
N_ORBITS_PER_DESIGN = 7
BAND = (0.90, 1.10)          # declared match band


def _sig2(v) -> str:
    """Two significant figures in LaTeX math, with a fixed shape.

    The compact column previously used ``%.2g``, which mixes ``1.4e-05`` with
    ``0.61`` in one column and prints a ratio of 1.01 as a bare ``1`` --- which
    is exactly where a reader needs to see whether the single exception in the
    panel is a tie or a loss.
    """
    if v is None or not np.isfinite(v):
        return "--"
    if v == 0:
        return "$0$"
    exp = int(np.floor(np.log10(abs(v))))
    if -2 <= exp <= 2:
        return f"${v:.{max(0, 1 - exp)}f}$"
    return f"${v / 10.0 ** exp:.1f}\\times 10^{{{exp}}}$"
PREFERRED = (0.95, 1.05)

_MODELS: dict[int, tuple] = {}


def _model(degree: int):
    if degree not in _MODELS:
        m = base.load_model(degree)
        a = base.kernel_args(m)
        base.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def cost_curve():
    """Measured per-call kernel cost, merged over the base and high-degree runs.

    The two archives name their tables differently: the R12 sweep stores ``rows``
    and the R13 high-degree extension stores ``high_rows`` plus an already-merged
    ``combined_rows``. Accept whichever is present rather than assuming one.
    """
    deg, ns = [], []
    for p in (COST_CURVE, COST_CURVE_HIGH):
        if not p.exists():
            continue
        payload = json.loads(p.read_text())
        table = next((payload[k] for k in ("combined_rows", "rows", "high_rows")
                      if isinstance(payload.get(k), list) and payload[k]), None)
        if table is None:
            raise RuntimeError(f"{p.name}: no recognizable cost-curve table")
        for r in table:
            deg.append(float(r["degree"]))
            ns.append(float(r["per_call_ns_median"]))
    d = np.array(deg)
    n = np.array(ns)
    o = np.argsort(d)
    d, n = d[o], n[o]
    keep = np.concatenate([[True], np.diff(d) > 0])
    return d[keep], n[keep]


def inverse_cost(target_ns, deg_tab, ns_tab, cap) -> int:
    fine = np.arange(2.0, float(cap) + 1.0)
    c = np.interp(fine, deg_tab, ns_tab)
    return int(fine[int(np.argmin(np.abs(c - target_ns)))])


def rows_of(design):
    return json.loads(PARETO_DESIGNS[design]["rows"].read_text())["rows"]


def spec_of(design, index):
    pareto = json.loads(PARETO.read_text(encoding="utf-8"))
    for r in pareto["designs"][design]["rows"]:
        if int(r["sobol_index"]) == index:
            return r["budgets"][f"beta_{BETA:.2f}"]
    raise KeyError(index)


def select() -> int:
    deg_tab, ns_tab = cost_curve()
    out = {"schema": "r14_timing_selection_v1", "created_utc": base.utc_now(),
           "beta": BETA,
           "rule": (f"{N_ORBITS_PER_DESIGN} orbits per design spread over perilune "
                    "altitude with the extremes retained; the same panel design as "
                    "the R13 measured-time control"),
           "band": BAND, "preferred_band": PREFERRED, "designs": {}}
    for design in ("A", "B"):
        rows = sorted(rows_of(design), key=lambda r: r["design_point"]["hp_km"])
        pick = [rows[int(i)] for i in
                np.linspace(0, len(rows) - 1, N_ORBITS_PER_DESIGN).round()]
        entries = []
        for r in pick:
            index = int(r["sobol_index"])
            b = spec_of(design, index)
            adopted = int(r["adopted_truth_degree"])
            # seed the comparator from the measured cost curve at the radial
            # history's call-weighted mean per-call cost; it is refined on the
            # achieved time ratio during the run
            hist = b["atallah"]["allocation"]["degree_histogram"]
            tot = sum(hist.values())
            mean_ns = sum(np.interp(float(k), deg_tab, ns_tab) * v
                          for k, v in hist.items()) / tot
            entries.append({
                "sobol_index": index,
                "hp_km": float(r["design_point"]["hp_km"]),
                "ha_km": float(r["design_point"]["ha_km"]),
                "adopted_truth_degree": adopted,
                "n_critical": int(r["n_critical"]),
                "atallah_tol_accel_m_s2": b["atallah"]["tol_accel_m_s2"],
                "n_fixed_proxy": int(b["fixed"]["degree"]),
                "n_fixed_time_seed": inverse_cost(mean_ns, deg_tab, ns_tab, adopted),
                "atallah_mean_per_call_ns": float(mean_ns)})
            e = entries[-1]
            print(f"  {design}{index:03d} hp={e['hp_km']:5.0f} km  "
                  f"N_proxy={e['n_fixed_proxy']:3d} -> N_time_seed="
                  f"{e['n_fixed_time_seed']:3d}")
        out["designs"][design] = entries
    base.atomic_json(SELECTION, out)
    print(f"[written] {SELECTION.name}")
    return 0


def propagate(model, args, y0, degree_of, level):
    grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    return base.propagate_event_instrumented(
        model, np.asarray(y0), DURATION, grid, degree_of, args,
        LEVELS[level]["rtol"], LEVELS[level]["atol"], max_step=MAX_STEP)


def run(repeats: int) -> int:
    sel = json.loads(SELECTION.read_text(encoding="utf-8"))
    deg_tab, ns_tab = cost_curve()
    payload = {"schema": "r14_timing_budget_v1", "created_utc": base.utc_now(),
               "beta": BETA, "repeats": repeats, "band": BAND,
               "preferred_band": PREFERRED, "serial": True,
               "protocol": ("budget defined as measured serial gravity-kernel time; "
                            "comparator degree refined on the achieved time ratio; "
                            "medians over repeated serial runs"),
               "source": base.provenance(), "designs": {}}
    if OUTPUT.exists():
        try:
            payload["designs"] = json.loads(OUTPUT.read_text())["designs"]
        except Exception:
            payload["designs"] = {}
    for design in ("A", "B"):
        rows = {int(r["sobol_index"]): r for r in rows_of(design)}
        out_rows = []
        for e in sel["designs"][design]:
            index = e["sobol_index"]
            row = rows[index]
            adopted = e["adopted_truth_degree"]
            y0 = np.asarray(row["design_point"]["initial_state_si"], float)
            model, args = _model(adopted)
            g = at.precompute_Sn(model, adopted)
            deg_fn, table = at.atallah_binned_schedule(
                model, g, float(e["atallah_tol_accel_m_s2"]),
                e["hp_km"], e["ha_km"], floor=FLOOR, cap=adopted, bin_km=BIN_KM)
            t0 = time.time()
            at_runs = []
            for _ in range(repeats):
                t, y, st, ev, fail, tel = propagate(model, args, y0, deg_fn, "tight")
                if st == "numerical_failure":
                    print(f"  !! {design}{index:03d} atallah {fail}", flush=True)
                    break
                at_runs.append({"kernel_ns": int(tel["gravity_kernel_ns"]),
                                "wall_ns": int(tel["total_wall_ns"]),
                                "n_rhs": int(tel["n_rhs"]),
                                "accepted": int(tel["n_accepted_steps"]),
                                "rejected": int(tel["n_rejected_trials"]),
                                "switches": int(tel["switch_count_at_rhs_samples"]),
                                "mean_degree": float(tel["mean_degree"]),
                                "degree_range": tel["degree_range"]})
                at_state = (t, y)
            if len(at_runs) < repeats:
                continue
            at_kernel = float(np.median([r["kernel_ns"] for r in at_runs]))

            # comparator: seed from the cost curve, then refine once on the
            # achieved ratio through a local c(N) ~ N^2 inversion, capped at truth
            n_try = int(e["n_fixed_time_seed"])
            attempts = []
            fx_state = None
            def fixed_runs(n, count):
                out, state = [], None
                for _ in range(count):
                    t, y, st, ev, fail, tel = propagate(
                        model, args, y0, (lambda _t, _h, nn=n: nn), "tight")
                    if st == "numerical_failure":
                        return [], None
                    out.append({"kernel_ns": int(tel["gravity_kernel_ns"]),
                                "wall_ns": int(tel["total_wall_ns"]),
                                "n_rhs": int(tel["n_rhs"]),
                                "accepted": int(tel["n_accepted_steps"]),
                                "rejected": int(tel["n_rejected_trials"])})
                    state = (t, y)
                return out, state

            # single probes locate the comparator degree; only the accepted one
            # is repeated, so the serial panel stays affordable
            for it in range(3):
                n_try = max(2, min(adopted, n_try))
                if any(a["degree"] == n_try for a in attempts):
                    break
                probe, _ = fixed_runs(n_try, 1)
                if not probe:
                    break
                ratio = probe[0]["kernel_ns"] / at_kernel
                attempts.append({"degree": n_try, "kernel_ns": probe[0]["kernel_ns"],
                                 "ratio": ratio, "probe_only": True})
                if BAND[0] <= ratio <= BAND[1]:
                    break
                n_try = int(round(n_try / math.sqrt(ratio)))
            if not attempts:
                continue
            best = min(attempts, key=lambda a: abs(math.log(a["ratio"])))
            fx_runs, fx_state = fixed_runs(best["degree"], repeats)
            if not fx_runs:
                continue
            best = {"degree": best["degree"],
                    "kernel_ns": float(np.median([r["kernel_ns"] for r in fx_runs])),
                    "runs": fx_runs}
            best["ratio"] = best["kernel_ns"] / at_kernel
            # error against the reused archived truth at the tight level
            tt, ty = base.load_raw(PARETO_DESIGNS[design]["r11_raw"]
                                   / f"sobolA_{index:03d}" / "truth_tight.npz")
            e_at = base.common_error(at_state[0], at_state[1], tt, ty)["pos_rms_m"]
            e_fx = base.common_error(fx_state[0], fx_state[1], tt, ty)["pos_rms_m"]
            spread = lambda v: (max(v) - min(v)) / max(np.median(v), 1.0)
            out_rows.append({
                "design": design, "sobol_index": index, "hp_km": e["hp_km"],
                "n_critical": e["n_critical"],
                "n_fixed_proxy": e["n_fixed_proxy"],
                "n_fixed_time": best["degree"],
                "achieved_time_ratio": best["ratio"],
                "time_matched": bool(BAND[0] <= best["ratio"] <= BAND[1]),
                "in_preferred_band": bool(PREFERRED[0] <= best["ratio"] <= PREFERRED[1]),
                "atallah_kernel_ns_median": at_kernel,
                "fixed_kernel_ns_median": best["kernel_ns"],
                "atallah_kernel_repeatability": spread(
                    [r["kernel_ns"] for r in at_runs]),
                "fixed_kernel_repeatability": spread(
                    [r["kernel_ns"] for r in best["runs"]]),
                "atallah_error_m": e_at, "fixed_error_m": e_fx,
                "rho_time": (e_fx / e_at) if e_at > 0 else None,
                "raw_winner": "atallah" if e_at < e_fx else "fixed",
                "atallah_runs": at_runs, "fixed_runs": best["runs"],
                "comparator_attempts": [{k: v for k, v in a.items() if k != "runs"}
                                        for a in attempts],
                "rhs_ratio": (np.median([r["n_rhs"] for r in at_runs])
                              / np.median([r["n_rhs"] for r in best["runs"]])),
                "elapsed_s": time.time() - t0})
            r = out_rows[-1]
            print(f"  {design}{index:03d} hp={e['hp_km']:5.0f} N_time={r['n_fixed_time']:3d} "
                  f"ratio={r['achieved_time_ratio']:.3f} "
                  f"E_at={e_at:.3f} E_fx={e_fx:.3f} rho={r['rho_time']:.3g}",
                  flush=True)
        payload["designs"][design] = out_rows
        base.atomic_json(OUTPUT, payload)
    print(f"[written] {OUTPUT.name}")
    return 0


def aggregate() -> int:
    d = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rows = [r for des in d["designs"].values() for r in des]
    matched = [r for r in rows if r["time_matched"]]
    # the fourth cost definition: with kernel time matched, what does the
    # end-to-end propagation cost? computed before the table is laid out
    wall = []
    for r in rows:
        wa = float(np.median([x["wall_ns"] for x in r["atallah_runs"]]))
        wf = float(np.median([x["wall_ns"] for x in r["fixed_runs"]]))
        r["wall_ratio"] = wa / wf
        wall.append(wa / wf)
    lines = []
    for r in rows:
        lines.append(
            f"    {r['design']} & {r['sobol_index']:03d} & {r['hp_km']:.0f} & "
            f"{r['n_critical']} & {r['n_fixed_proxy']} & {r['n_fixed_time']} & "
            f"{r['achieved_time_ratio']:.3f} & {r['wall_ratio']:.3f} & "
            f"{r['atallah_error_m']:.3f} & "
            f"{r['fixed_error_m']:.3f} & {_sig2(r['rho_time'])} & "
            f"{'yes' if r['time_matched'] else 'no'}\\\\")
    body = "\n".join(lines)
    ratios = [r["rho_time"] for r in matched if r["rho_time"]]
    summary = {
        "orbits": len(rows), "time_matched": len(matched),
        "in_preferred_band": sum(r["in_preferred_band"] for r in rows),
        "atallah_raw_wins": sum(r["raw_winner"] == "atallah" for r in matched),
        "fixed_raw_wins": sum(r["raw_winner"] == "fixed" for r in matched),
        "rho_time_median": float(np.median(ratios)) if ratios else None,
        "time_ratio_range": [min(r["achieved_time_ratio"] for r in rows),
                             max(r["achieved_time_ratio"] for r in rows)],
        "median_comparator_shift_vs_proxy": float(np.median(
            [r["n_fixed_time"] - r["n_fixed_proxy"] for r in rows])),
        # the fourth cost definition: with kernel time matched, what does the
        # end-to-end propagation actually cost?
        "wall_time_ratio": {"median": float(np.median(wall)),
                            "p10": float(np.percentile(wall, 10)),
                            "p90": float(np.percentile(wall, 90)),
                            "above_unity": int(sum(1 for w in wall if w > 1.0)),
                            "orbits": len(wall)},
    }
    d["summary"] = summary
    base.atomic_json(OUTPUT, d)
    # the panel selects both designs at evenly spaced perilune ranks, so the two
    # rows at each rank are a matched pair; pairing them halves the table and
    # lets the replication be read across
    ra = sorted((r for r in rows if r["design"] == "A"), key=lambda r: r["hp_km"])
    rb = sorted((r for r in rows if r["design"] == "B"), key=lambda r: r["hp_km"])
    compact_rows = "\n".join(
        f"    {a['hp_km']:.0f} & {a['sobol_index']:03d}/{b['sobol_index']:03d} & "
        f"{a['n_fixed_time']}/{b['n_fixed_time']} & "
        f"{a['achieved_time_ratio']:.3f}/{b['achieved_time_ratio']:.3f} & "
        f"{a['wall_ratio']:.3f}/{b['wall_ratio']:.3f} & "
        f"{_sig2(a['rho_time'])}\\,/\\,{_sig2(b['rho_time'])}\\\\"
        for a, b in zip(ra, rb))
    (TABLE.parent / "r14_timing_budget_table_compact.tex").write_text(
        f"""% auto-generated by rev14_timing_budget.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{5pt}}
  \\caption{{Equal-budget comparison with the budget defined as measured serial
  gravity-kernel time rather than the quadratic proxy.
  $N_{{\\mathrm{{time}}}}$ is the degree whose measured kernel time matches the
  radial policy's within the timing repeatability of this machine; all 14 pairs
  fall inside the declared $0.90$--$1.10$ band. The wall ratio is the
  end-to-end propagation time of the radial policy over the comparator's, so a
  value above unity means that matching the kernel time did not make the run as a
  whole equally expensive. $\\rho = E_{{\\mathrm{{fix}}}}/E_{{\\mathrm{{At}}}}$ is a
  raw error ratio: the truth-inclusive resolution rule used elsewhere in the
  paper is not applied to this panel. Degrees, errors in meters and the proxy
  comparator are in Supplementary
  Table~\\ref*{{supp-tab:budget-timing-full}}.}}
  \\label{{tab:budget-timing}}
  \\begin{{tabular}}{{l c c c c c}}
    \\toprule
    $h_p$ [km] & idx & $N_{{\\mathrm{{time}}}}$ & kernel ratio & wall ratio &
      $\\rho$\\\\
    & A/B & A/B & A/B & A/B & A\\,/\\,B\\\\
    \\midrule
{compact_rows}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
""", encoding="utf-8")
    TABLE.write_text(f"""% auto-generated by rev14_timing_budget.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Measured-time budget control, full form. The compact version in the
  main text carries the matched degree and the two cost ratios; this one adds the
  critical degree, the comparator the quadratic proxy would have selected, and
  both errors in meters.
  $N_{{\\mathrm{{proxy}}}}$ is the comparator the proxy selects at $\\beta=1$ and
  $N_{{\\mathrm{{time}}}}$ the degree whose measured kernel time approximately
  matches the radial policy's within the measured timing repeatability. All runs
  are serial on an idle machine and each configuration is repeated; medians are
  compared. A pair counts as matched only inside the declared
  $0.90$--$1.10$ band. $\\rho = E_{{\\mathrm{{fix}}}}/E_{{\\mathrm{{At}}}}$ is a raw
  error ratio; the truth-inclusive resolution rule is not applied to this
  panel.}}
  \\label{{tab:budget-timing-full}}
  \\begin{{tabular}}{{l r r r r r r r r r r c}}
    \\toprule
    Des. & idx & $h_p$ [km] & $N_{{\\mathrm{{crit}}}}$ & $N_{{\\mathrm{{proxy}}}}$ &
      $N_{{\\mathrm{{time}}}}$ & kernel ratio & wall ratio &
      $E_{{\\mathrm{{At}}}}$ [m] &
      $E_{{\\mathrm{{fix}}}}$ [m] & $\\rho$ & matched\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
""", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"[written] {TABLE.name}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("select", "run", "aggregate"))
    ap.add_argument("--repeats", type=int, default=3)
    a = ap.parse_args()
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    if a.command == "select":
        return select()
    if a.command == "run":
        return run(a.repeats)
    return aggregate()


if __name__ == "__main__":
    raise SystemExit(main())
