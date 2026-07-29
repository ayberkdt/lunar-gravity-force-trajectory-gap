"""Bin-resolution control for the Atallah benchmark (R12).

The published rule is a function of the instantaneous radius. The benchmark
applies it on the same 10-km altitude bins as the paper's other schedules (the
original paper explicitly allows a stored lookup, and binning keeps the degree
constant inside an integration step). This control tests whether that
discretization changes the result: for a set of representative orbits it repeats
the benchmark with the degree evaluated at the exact current radius on every
right-hand-side call, and compares degree history, switch count, work, cost, and
seven-day trajectory error against the binned run.

Everything else is frozen to the campaign: the same per-orbit acceleration
tolerance (read from the archived binned sidecar, not recomputed), the same
degree floor and truth-degree cap, the same vector tolerances, arc, and output
grid, and the same reused truth trajectories.

Usage:
    python rev12_atallah_bincontrol.py pilot --workers 2
    python rev12_atallah_bincontrol.py run --workers 5
    python rev12_atallah_bincontrol.py table
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
import rev12_atallah as at
import rev12_atallah_campaign as camp

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
CASE_ROOT = METRICS / "r12_cases" / "atallah_bincontrol"
RAW_ROOT = METRICS / "r12_raw" / "atallah_bincontrol"
OUTPUT = METRICS / "r12_atallah_bincontrol.json"

# guard: a continuous-radius rule can collapse the adaptive step size, so cap the
# right-hand-side budget instead of letting a single orbit run unbounded.
RHS_BUDGET = 4_000_000

N_ORBITS = 10


def paths(index: int, level: str):
    return (CASE_ROOT / f"sobolA_{index:03d}" / f"atallah_continuous_{level}.json",
            RAW_ROOT / f"sobolA_{index:03d}" / f"atallah_continuous_{level}.npz")


def binned_sidecar(index: int, level: str) -> Path:
    return camp.CASE_ROOT / f"sobolA_{index:03d}" / f"atallah_{level}.json"


def select_orbits(rows: list[dict], n: int = N_ORBITS) -> list[dict]:
    """Representative orbits: evenly spaced in perilune altitude over the
    truth-surviving population, so the control spans the regimes where the bin
    width is a larger and a smaller fraction of the perilune--apolune span."""
    ok = sorted(rows, key=lambda r: r["design_point"]["hp_km"])
    if len(ok) <= n:
        return ok
    idx = np.linspace(0, len(ok) - 1, n).round().astype(int)
    return [ok[i] for i in sorted(set(int(i) for i in idx))]


def _budgeted(fn, budget: int = RHS_BUDGET):
    state = {"n": 0}

    def wrapped(t, h_m):
        state["n"] += 1
        if state["n"] > budget:
            raise RuntimeError(f"rhs budget {budget} exceeded (step collapse)")
        return fn(t, h_m)

    return wrapped


def worker(task: dict) -> dict:
    row = task["row"]
    index = int(row["sobol_index"])
    try:
        side = json.loads(binned_sidecar(index, "tight").read_text())
        cfg0 = side["config"]
        adopted = int(cfg0["adopted_truth_degree"])
        tol = float(cfg0["atallah_tol_accel_m_s2"])
        y0 = np.asarray(cfg0["initial_state_si"], dtype=float)
        model, args = camp._model(adopted)
        g = camp._g(adopted)
        cont = at.atallah_degree_fn(model, g, tol, floor=2, cap=adopted)

        out = {"index": index, "status": "complete", "tol": tol,
               "adopted_truth_degree": adopted,
               "n_critical": int(cfg0["n_critical"]),
               "design_point": {k: row["design_point"][k]
                                for k in ("hp_km", "ha_km", "incl_deg",
                                          "eccentricity")},
               "telemetry": {}}
        for level in ("tight", "tighter"):
            t, y, st, ev, fail, tel = camp._propagate(
                model, args, y0, _budgeted(cont), level)
            if st == "numerical_failure":
                return {"index": index, "status": "numerical_failure",
                        "where": f"continuous/{level}", "message": fail,
                        "telemetry": {level: tel}}
            cfg = {**{k: v for k, v in cfg0.items()
                      if k not in ("policy", "level", "policy_spec",
                                   "atallah_degree_table", "rtol",
                                   "atol_position_m", "atol_velocity_m_s")},
                   "policy": "atallah_continuous", "level": level,
                   "rtol": camp.LEVELS[level]["rtol"],
                   "atol_position_m": camp.LEVELS[level]["atol_position_m"],
                   "atol_velocity_m_s": camp.LEVELS[level]["atol_velocity_m_s"],
                   "bin_control": ("degree evaluated at the exact instantaneous "
                                   "radius on every RHS call; identical tolerance, "
                                   "floor and cap to the archived 10-km binned run"),
                   "policy_spec": {"kind": "atallah_radial_adaptive_continuous",
                                   "tol": tol}}
            sidecar, raw = paths(index, level)
            base.atomic_npz(raw, t_s=t, state_si=y)
            base.atomic_json(sidecar, {
                "schema": "r12_atallah_bincontrol_trajectory_v1",
                "created_utc": base.utc_now(), "config": cfg,
                "config_sha256": base.object_hash(cfg), "status": st,
                "event": ev, "telemetry": tel,
                "raw_path": str(raw.relative_to(ROOT)),
                "raw_sha256": base.file_hash(raw),
                "n_output_epochs": int(len(t)),
                "last_output_epoch_s": float(t[-1])})
            out["telemetry"][level] = tel
        return out
    except Exception as exc:
        return {"index": index, "status": "worker_error",
                "message": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc()}


# ------------------------------------------------------------------ aggregation
def _load(path: Path):
    return base.load_raw(path)


def orbit_summary(row: dict) -> dict:
    index = int(row["sobol_index"])
    truth = {lv: base.load_raw(camp.reuse_paths(index, "truth", lv)[1])
             for lv in camp.LEVELS}
    binned = {lv: base.load_raw(camp.paths(index, "atallah", lv)[1])
              for lv in camp.LEVELS}
    cont = {lv: base.load_raw(paths(index, lv)[1]) for lv in camp.LEVELS}

    def err(pol, lv):
        return base.common_error(pol[lv][0], pol[lv][1],
                                 truth[lv][0], truth[lv][1])["pos_rms_m"]

    def self_diff(pol):
        return base.common_error(pol["tight"][0], pol["tight"][1],
                                 pol["tighter"][0], pol["tighter"][1])["pos_rms_m"]

    truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                   truth["tighter"][0], truth["tighter"][1])["pos_rms_m"]
    tel_b = json.loads(binned_sidecar(index, "tight").read_text())["telemetry"]
    tel_c = json.loads(paths(index, "tight")[0].read_text())["telemetry"]
    e_b, e_c = err(binned, "tight"), err(cont, "tight")
    env_b = self_diff(binned) + truth_self
    env_c = self_diff(cont) + truth_self
    mutual = base.common_error(binned["tight"][0], binned["tight"][1],
                               cont["tight"][0], cont["tight"][1])["pos_rms_m"]

    def tel_view(t):
        return {"mean_degree": t["mean_degree"], "mean_degree_sq": t["mean_degree_sq"],
                "degree_range": t["degree_range"], "n_rhs": t["n_rhs"],
                "n_distinct_degrees": len(t["degree_counts"]),
                "switch_count_at_rhs_samples": t["switch_count_at_rhs_samples"],
                "gravity_kernel_ns": t["gravity_kernel_ns"],
                "n_accepted_steps": t["n_accepted_steps"],
                "n_rejected_trials": t["n_rejected_trials"]}

    return {
        "sobol_index": index, "name": row["name"],
        "design_point": {k: row["design_point"][k]
                         for k in ("hp_km", "ha_km", "incl_deg", "eccentricity")},
        "n_critical": int(row["n_critical"]),
        "adopted_truth_degree": int(row["adopted_truth_degree"]),
        "binned": {"error_m": e_b, "envelope_m": env_b, "telemetry": tel_view(tel_b)},
        "continuous": {"error_m": e_c, "envelope_m": env_c, "telemetry": tel_view(tel_c)},
        "binned_vs_continuous": {
            "mutual_difference_rms_m": mutual,
            "absolute_error_difference_m": abs(e_b - e_c),
            "resolution_threshold_m": env_b + env_c,
            "resolved": bool(abs(e_b - e_c) > env_b + env_c),
            "winner_if_resolved": (("continuous" if e_c < e_b else "binned")
                                   if abs(e_b - e_c) > env_b + env_c else None),
            "work_ratio_continuous_over_binned": (
                tel_c["mean_degree_sq"] / tel_b["mean_degree_sq"]),
            "mean_degree_difference": tel_c["mean_degree"] - tel_b["mean_degree"],
            "rhs_ratio_continuous_over_binned": tel_c["n_rhs"] / tel_b["n_rhs"],
            "kernel_time_ratio_continuous_over_binned": (
                tel_c["gravity_kernel_ns"] / tel_b["gravity_kernel_ns"])},
    }


def stat(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)])
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max())}


def aggregate(rows: list[dict], extra: dict) -> dict:
    cmpv = [r["binned_vs_continuous"] for r in rows]
    payload = {
        "schema": "r12_atallah_bincontrol_v1",
        "protocol": ("continuous-radius Atallah versus the archived 10-km binned "
                     "Atallah on representative design-A orbits; identical "
                     "tolerance, cap, integrator settings, arc and truth"),
        "orbits": len(rows), "rows": rows,
        "summary": {
            "resolved_error_differences": sum(c["resolved"] for c in cmpv),
            "unresolved": sum(not c["resolved"] for c in cmpv),
            "continuous_wins": sum(c["winner_if_resolved"] == "continuous" for c in cmpv),
            "binned_wins": sum(c["winner_if_resolved"] == "binned" for c in cmpv),
            "mutual_difference_rms_m": stat([c["mutual_difference_rms_m"] for c in cmpv]),
            "work_ratio_continuous_over_binned": stat(
                [c["work_ratio_continuous_over_binned"] for c in cmpv]),
            "rhs_ratio_continuous_over_binned": stat(
                [c["rhs_ratio_continuous_over_binned"] for c in cmpv]),
            "kernel_time_ratio_continuous_over_binned": stat(
                [c["kernel_time_ratio_continuous_over_binned"] for c in cmpv]),
            "mean_degree_difference": stat([c["mean_degree_difference"] for c in cmpv]),
            "binned_error_m": stat([r["binned"]["error_m"] for r in rows]),
            "continuous_error_m": stat([r["continuous"]["error_m"] for r in rows])},
        **extra}
    return payload


def build_table(payload: dict) -> str:
    s = payload["summary"]
    lines = []
    for r in payload["rows"]:
        c = r["binned_vs_continuous"]
        tb, tc = r["binned"]["telemetry"], r["continuous"]["telemetry"]
        lines.append(
            f"    {r['sobol_index']:d} & {r['design_point']['hp_km']:.0f} & "
            f"{tb['mean_degree']:.1f} & {tc['mean_degree']:.1f} & "
            f"{tb['n_distinct_degrees']:d} & {tc['n_distinct_degrees']:d} & "
            f"{tb['switch_count_at_rhs_samples']:d} & "
            f"{tc['switch_count_at_rhs_samples']:d} & "
            f"{c['work_ratio_continuous_over_binned']:.3f} & "
            f"{c['rhs_ratio_continuous_over_binned']:.2f} & "
            f"{r['binned']['error_m']:.3f} & {r['continuous']['error_m']:.3f} & "
            f"{'yes' if c['resolved'] else 'no'}\\\\")
    body = "\n".join(lines)
    md = s["mutual_difference_rms_m"]
    wr = s["work_ratio_continuous_over_binned"]
    return f"""% auto-generated by rev12_atallah_bincontrol.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\scriptsize
  \\setlength{{\\tabcolsep}}{{3pt}}
  \\caption{{Bin-resolution control for the Atallah benchmark. For
  {payload['orbits']} design-A orbits spanning the perilune range, the published
  rule evaluated at the exact instantaneous radius on every right-hand-side call
  (cont) is compared with the archived 10-km binned application (bin) at
  identical tolerance, degree cap, integrator settings, arc and truth.
  $\\bar N$ is the call-weighted mean degree, $n_{{\\mathrm{{lev}}}}$ the number of
  distinct degrees used, and $n_{{\\mathrm{{sw}}}}$ the degree changes seen at
  right-hand-side samples. $w$ is the ratio of $\\langle N^2\\rangle$ and
  $n_{{\\mathrm{{rhs}}}}$ the ratio of right-hand-side calls (cont/bin). $E$ is the
  seven-day position RMS against the same-tolerance truth; the last column
  reports whether the two error values are separated by the truth-inclusive
  envelope rule. Median mutual trajectory difference
  {md['median']:.3f}~m; median work ratio {wr['median']:.3f}.}}
  \\label{{tab:atallah-bincontrol}}
  \\begin{{tabular}}{{r r r r r r r r r r r r c}}
    \\toprule
    & $h_p$ & \\multicolumn{{2}}{{c}}{{$\\bar N$}} & \\multicolumn{{2}}{{c}}{{$n_{{\\mathrm{{lev}}}}$}} & \\multicolumn{{2}}{{c}}{{$n_{{\\mathrm{{sw}}}}$}} & & & \\multicolumn{{2}}{{c}}{{$E$ [m]}} & \\\\
    \\cmidrule(lr){{3-4}}\\cmidrule(lr){{5-6}}\\cmidrule(lr){{7-8}}\\cmidrule(lr){{11-12}}
    idx & [km] & bin & cont & bin & cont & bin & cont & $w$ & $n_{{\\mathrm{{rhs}}}}$ & bin & cont & sep.\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""


def command_table() -> int:
    payload = json.loads(OUTPUT.read_text())
    tex = build_table(payload)
    (METRICS / "r12_atallah_bincontrol_table.tex").write_text(tex, encoding="utf-8")
    s = payload["summary"]
    print("[written] r12_atallah_bincontrol_table.tex")
    print(json.dumps(s, indent=2))
    return 0


def run(rows, workers) -> int:
    for r in rows:
        for sub in (CASE_ROOT, RAW_ROOT):
            (sub / f"sobolA_{int(r['sobol_index']):03d}").mkdir(parents=True, exist_ok=True)
    started = base.utc_now()
    t0 = time.time()
    print(f"[bincontrol] orbits={len(rows)} workers={workers}", flush=True)
    results, failures = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(worker, {"row": r}): r for r in rows}
        for n, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            if rec["status"] != "complete":
                failures.append(rec)
                print(f"  !! {rec['index']:03d} {rec['status']} "
                      f"{rec.get('where', '')}: {rec.get('message', '')}", flush=True)
            else:
                results.append(rec)
            print(f"  [{n:2d}/{len(rows)}] idx={rec['index']:03d} {rec['status']} "
                  f"elapsed={(time.time() - t0) / 60:.1f}min", flush=True)
    summaries = []
    for r in sorted(rows, key=lambda x: int(x["sobol_index"])):
        idx = int(r["sobol_index"])
        if all(paths(idx, lv)[0].exists() for lv in camp.LEVELS):
            try:
                summaries.append(orbit_summary(r))
            except Exception as exc:
                failures.append({"index": idx, "status": "summary_error",
                                 "message": f"{type(exc).__name__}: {exc}"})
    payload = aggregate(summaries, {
        "started_utc": started, "ended_utc": base.utc_now(),
        "population": "A", "rhs_budget": RHS_BUDGET,
        "failures": failures, "session_wall_s": time.time() - t0,
        "complete": len(summaries) == len(rows) and not failures})
    base.atomic_json(OUTPUT, payload)
    print(f"[bincontrol] done orbits={len(summaries)}/{len(rows)} "
          f"failures={len(failures)}", flush=True)
    if summaries:
        print(json.dumps(payload["summary"], indent=2))
    return 0 if payload["complete"] else 3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("pilot", "run", "table"))
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--orbits", type=int, default=N_ORBITS)
    a = ap.parse_args()
    if a.command == "table":
        return command_table()
    rows = [r for r in camp.load_rows()
            if binned_sidecar(int(r["sobol_index"]), "tight").exists()]
    picks = select_orbits(rows, a.orbits)
    if a.command == "pilot":
        picks = [picks[0], picks[len(picks) // 2]]
    return run(picks, a.workers)


if __name__ == "__main__":
    raise SystemExit(main())
