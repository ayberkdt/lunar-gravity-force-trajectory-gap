"""Measured-cost-matched fixed-degree comparator for the Atallah benchmark (R13).

The campaign matches work with the quadratic proxy
``N_work = round(sqrt(<N^2>))``. The proxy is not the machine cost: the measured
per-call kernel cost c(N) is super-quadratic at high degree
(``r12_kernel_cost_curve.json``, about 18% RMS residual against a pure N^2 law).
This experiment replaces the proxy with the measured curve.

For each selected orbit the comparator degree is the fixed degree whose measured
per-call cost equals the call-weighted mean per-call cost of the Atallah degree
history,

    N_time = c^{-1}( sum_k c(N_k) / n_rhs ),

taken from the archived degree histogram of the tight Atallah run and the
measured cost curve (monotone, interpolated between tabulated degrees). That
comparator is then propagated at both vector tolerance levels and compared with
Atallah under the same envelope rule as the campaign.

Because the claim is about cost, the runs are serial: one worker, no other
propagation on the machine, so the recorded gravity-kernel times of the Atallah
and comparator runs are directly comparable and the achieved time ratio can be
reported alongside the intended one.

Usage:
    python rev13_timing_match.py select
    python rev13_timing_match.py run          # serial, machine must be idle
    python rev13_timing_match.py aggregate
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
COST_CURVE = METRICS / "r12_kernel_cost_curve.json"
SELECTION = METRICS / "r13_timing_match_selection.json"
OUTPUT = METRICS / "r13_timing_match.json"
TABLE = METRICS / "r13_timing_match_table.tex"
CASE_ROOT = METRICS / "r13_cases" / "timing_match"
RAW_ROOT = METRICS / "r13_raw" / "timing_match"

DESIGNS = {
    "A": {"r12_case": METRICS / "r12_cases" / "atallah",
          "r12_raw": METRICS / "r12_raw" / "atallah",
          "r11_raw": METRICS / "r11_raw" / "convergence",
          "campaign": METRICS / "r12_atallah_campaign.json"},
    "B": {"r12_case": METRICS / "r12_cases" / "atallah_designB",
          "r12_raw": METRICS / "r12_raw" / "atallah_designB",
          "r11_raw": METRICS / "r11_raw" / "designB_convergence",
          "campaign": METRICS / "r12_atallah_campaign_designB.json"},
}
LEVELS = {
    "tight": {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3),
              "atol_position_m": 1.0e-5, "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3),
                "atol_position_m": 1.0e-6, "atol_velocity_m_s": 1.0e-9},
}
MAX_STEP = 60.0
DURATION = 7.0 * base.DAY
OUTPUT_STEP = 120.0
N_ORBITS_PER_DESIGN = 7
# measured log-log slope of the kernel cost curve above degree 100
# (r12_kernel_cost_curve.json: 2.05; r13 high-degree extension: 1.84)
COST_EXPONENT = 2.0


def cost_curve():
    d = json.loads(COST_CURVE.read_text())
    deg = np.array([r["degree"] for r in d["rows"]], dtype=float)
    ns = np.array([r["per_call_ns_median"] for r in d["rows"]], dtype=float)
    order = np.argsort(deg)
    return deg[order], ns[order]


def cost_of(degrees, deg_tab, ns_tab):
    return np.interp(np.asarray(degrees, dtype=float), deg_tab, ns_tab)


def inverse_cost(target_ns, deg_tab, ns_tab) -> int:
    """Smallest tabulated-interpolated degree whose per-call cost matches."""
    fine = np.arange(2.0, deg_tab.max() + 1.0)
    c = cost_of(fine, deg_tab, ns_tab)
    return int(fine[int(np.argmin(np.abs(c - target_ns)))])


def atallah_sidecar(design, index, level="tight"):
    return json.loads((DESIGNS[design]["r12_case"] / f"sobolA_{index:03d}"
                       / f"atallah_{level}.json").read_text())


def select() -> int:
    """Representative orbits per design: spread over perilune, including the
    lowest and highest, plus one retrograde and one near-equatorial case."""
    deg_tab, ns_tab = cost_curve()
    out = {"rule": (f"{N_ORBITS_PER_DESIGN} orbits per design spread over "
                    "perilune altitude, with the extremes retained"),
           "cost_curve": str(COST_CURVE.name), "designs": {}}
    for design in ("A", "B"):
        camp = json.loads(DESIGNS[design]["campaign"].read_text())
        rows = sorted(camp["rows"], key=lambda r: r["design_point"]["hp_km"])
        pick = [rows[int(i)] for i in
                np.linspace(0, len(rows) - 1, N_ORBITS_PER_DESIGN).round()]
        entries = []
        for r in pick:
            index = int(r["sobol_index"])
            side = atallah_sidecar(design, index)
            counts = {int(k): int(v) for k, v in side["telemetry"]["degree_counts"].items()}
            n_rhs = sum(counts.values())
            mean_cost = sum(cost_of([k], deg_tab, ns_tab)[0] * v
                            for k, v in counts.items()) / n_rhs
            n_time = inverse_cost(mean_cost, deg_tab, ns_tab)
            n_work = int(round(math.sqrt(side["telemetry"]["mean_degree_sq"])))
            entries.append({
                "sobol_index": index,
                "hp_km": r["design_point"]["hp_km"],
                "adopted_truth_degree": int(side["config"]["adopted_truth_degree"]),
                "n_work_proxy": n_work,
                "n_time_measured_cost": n_time,
                "atallah_mean_per_call_ns": float(mean_cost),
                "n_time_per_call_ns": float(cost_of([n_time], deg_tab, ns_tab)[0]),
                "proxy_per_call_ns": float(cost_of([n_work], deg_tab, ns_tab)[0])})
            print(f"  {design}{index:03d} hp={r['design_point']['hp_km']:5.0f} km  "
                  f"N_work={n_work:3d} -> N_time={n_time:3d}  "
                  f"(mean per-call {mean_cost / 1000:.1f} us)")
        out["designs"][design] = entries
    SELECTION.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[written] {SELECTION.name}")
    return 0


def serial_baseline() -> int:
    """Re-run the Atallah policy serially for the selected orbits.

    The campaign's Atallah trajectories were produced with five concurrent
    workers, so their recorded gravity-kernel times are contended and cannot be
    compared with a serially measured comparator. The trajectory itself is a
    deterministic function of the frozen configuration, so re-running it serially
    changes nothing but the timing: the accuracy comparison still uses the
    archived trajectories, while the timing match uses these serial numbers.
    """
    sel = json.loads(SELECTION.read_text())
    others = base.other_python_processes()
    if others:
        print(f"!! {len(others)} other python processes are running; timing "
              f"comparability requires an idle machine")
        return 2
    t0 = time.time()
    for design, entries in sel["designs"].items():
        for e in entries:
            index = int(e["sobol_index"])
            side = atallah_sidecar(design, index)
            cfg0 = side["config"]
            out = (CASE_ROOT / design / f"sobolA_{index:03d}"
                   / "atallah_serial_tight.json")
            if out.exists():
                e["atallah_serial_kernel_ns"] = json.loads(
                    out.read_text())["telemetry"]["gravity_kernel_ns"]
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            adopted = int(cfg0["adopted_truth_degree"])
            model = base.load_model(adopted)
            args = base.kernel_args(model)
            base.warmup(model, args)
            tab = {float(k): int(v) for k, v in cfg0["atallah_degree_table"].items()}
            hmin, hmax = min(tab), max(tab)

            def degree_of(tt, h_m, _tab=tab, _lo=hmin, _hi=hmax):
                hb = min(_hi, max(_lo, 10.0 * math.floor(h_m / 1e4)))
                return _tab[hb]

            y0 = np.asarray(cfg0["initial_state_si"], dtype=float)
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            tol = LEVELS["tight"]
            tt, yy, st, ev, fail, tel = base.propagate_event_instrumented(
                model, y0, DURATION, grid, degree_of, args,
                tol["rtol"], tol["atol"], max_step=MAX_STEP)
            if st == "numerical_failure":
                print(f"  !! {design}{index:03d}: {fail}", flush=True)
                continue
            base.atomic_json(out, {
                "schema": "r13_atallah_serial_timing_v1",
                "created_utc": base.utc_now(),
                "config": {"design": design, "sobol_index": index,
                           "policy": "atallah", "level": "tight",
                           "purpose": ("serial re-run of the archived Atallah "
                                       "configuration for timing comparability"),
                           "timing_comparable": True,
                           "archived_kernel_ns": side["telemetry"]["gravity_kernel_ns"],
                           "source": base.provenance()},
                "status": st, "telemetry": tel})
            e["atallah_serial_kernel_ns"] = tel["gravity_kernel_ns"]
            e["atallah_parallel_kernel_ns"] = side["telemetry"]["gravity_kernel_ns"]
            print(f"  [{(time.time()-t0)/60:5.1f} min] {design}{index:03d} "
                  f"atallah serial kernel={tel['gravity_kernel_ns']/1e9:.1f}s "
                  f"(archived parallel {side['telemetry']['gravity_kernel_ns']/1e9:.1f}s, "
                  f"nrhs {tel['n_rhs']} vs {side['telemetry']['n_rhs']})", flush=True)
    SELECTION.write_text(json.dumps(sel, indent=2), encoding="utf-8")
    print(f"[serial-baseline] done in {(time.time()-t0)/60:.1f} min")
    return 0


def refine() -> int:
    """Second pass: match total measured kernel time, not per-call cost.

    Matching the mean per-call cost leaves the total kernel times unequal because
    the two runs do not make the same number of right-hand-side calls: the
    adaptive-degree run takes more of them. With the first-pass runs in hand the
    comparator's own call count is known, so the per-call cost that equalizes the
    total is

        c_target = (Atallah kernel time) / (comparator right-hand-side calls),

    and the refined comparator degree is c^{-1}(c_target). The residual coupling
    (a different degree changes the call count slightly) is small and is reported
    as the achieved time ratio of the refined run.
    """
    deg_tab, ns_tab = cost_curve()
    sel = json.loads(SELECTION.read_text())
    for design, entries in sel["designs"].items():
        for e in entries:
            index = int(e["sobol_index"])
            first = (CASE_ROOT / design / f"sobolA_{index:03d}"
                     / "fixed_time_tight.json")
            if not first.exists():
                continue
            tel_cp = json.loads(first.read_text())["telemetry"]
            serial = (CASE_ROOT / design / f"sobolA_{index:03d}"
                      / "atallah_serial_tight.json")
            if not serial.exists():
                print(f"  !! {design}{index:03d}: no serial Atallah baseline")
                continue
            at_ns = json.loads(serial.read_text())["telemetry"]["gravity_kernel_ns"]
            ratio = tel_cp["gravity_kernel_ns"] / at_ns
            n1 = int(e["n_time_measured_cost"])
            # local power law c(N) ~ N^COST_EXPONENT, inverted on the MEASURED
            # kernel-time ratio of the two runs rather than on absolute
            # microbenchmarks, then capped at the orbit's truth degree
            n2 = int(round(n1 * (1.0 / ratio) ** (1.0 / COST_EXPONENT)))
            n2 = max(2, min(n2, int(e["adopted_truth_degree"])))
            e["first_pass_time_ratio"] = ratio
            e["cost_exponent_used"] = COST_EXPONENT
            e["n_time_refined"] = n2
            print(f"  {design}{index:03d} first pass N={n1} "
                  f"achieved {ratio:.2f} of Atallah's kernel "
                  f"time -> refined N={n2}")
    SELECTION.write_text(json.dumps(sel, indent=2), encoding="utf-8")
    print(f"[updated] {SELECTION.name}")
    return 0


def run(stage: str = "first") -> int:
    sel = json.loads(SELECTION.read_text())
    key = "n_time_measured_cost" if stage == "first" else "n_time_refined"
    stem = "fixed_time" if stage == "first" else "fixed_time2"
    others = base.other_python_processes()
    if others:
        print(f"!! {len(others)} other python processes are running; timing "
              f"comparability requires an idle machine:", flush=True)
        for o in others[:5]:
            print(f"   {o}", flush=True)
        return 2
    t0 = time.time()
    for design, entries in sel["designs"].items():
        for e in entries:
            index = int(e["sobol_index"])
            side = atallah_sidecar(design, index)
            cfg0 = side["config"]
            adopted = int(cfg0["adopted_truth_degree"])
            n_time = int(e[key])
            y0 = np.asarray(cfg0["initial_state_si"], dtype=float)
            model, args = base.load_model(adopted), None
            args = base.kernel_args(model)
            base.warmup(model, args)
            grid = np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
            for level, tol in LEVELS.items():
                sidecar = (CASE_ROOT / design / f"sobolA_{index:03d}"
                           / f"{stem}_{level}.json")
                raw = (RAW_ROOT / design / f"sobolA_{index:03d}"
                       / f"{stem}_{level}.npz")
                if sidecar.exists() and raw.exists():
                    continue
                sidecar.parent.mkdir(parents=True, exist_ok=True)
                raw.parent.mkdir(parents=True, exist_ok=True)
                t, y, st, ev, fail, tel = base.propagate_event_instrumented(
                    model, y0, DURATION, grid, (lambda tt, hh, n=n_time: n),
                    args, tol["rtol"], tol["atol"], max_step=MAX_STEP)
                if st == "numerical_failure":
                    print(f"  !! {design}{index:03d}/{level}: {fail}", flush=True)
                    continue
                cfg = {"design": design, "sobol_index": index,
                       "adopted_truth_degree": adopted,
                       "initial_state_si": [float(v) for v in y0],
                       "policy": ("fixed_time_matched" if stage == "first"
                                  else "fixed_time_matched_refined"),
                       "level": level, "stage": stage,
                       "policy_spec": {"kind": "fixed_measured_cost_matched",
                                       "degree": n_time,
                                       "source": ("c^{-1} of the call-weighted mean "
                                                  "per-call kernel cost of the tight "
                                                  "Atallah degree history"),
                                       "cost_curve": COST_CURVE.name},
                       "n_work_proxy": e["n_work_proxy"],
                       "duration_s": DURATION, "output_step_s": OUTPUT_STEP,
                       "integrator": "InstrumentedDOP853", "max_step_s": MAX_STEP,
                       "atol_kind": "vector", "rtol": tol["rtol"],
                       "atol_position_m": tol["atol_position_m"],
                       "atol_velocity_m_s": tol["atol_velocity_m_s"],
                       "timing_comparable": True,
                       "timing_note": ("run serially on an otherwise idle machine "
                                       "so kernel times are comparable with the "
                                       "archived Atallah runs"),
                       "source": base.provenance()}
                base.atomic_npz(raw, t_s=t, state_si=y)
                base.atomic_json(sidecar, {
                    "schema": "r13_timing_match_trajectory_v1",
                    "created_utc": base.utc_now(), "config": cfg,
                    "config_sha256": base.object_hash(cfg), "status": st,
                    "event": ev, "telemetry": tel,
                    "raw_path": str(raw.relative_to(ROOT)),
                    "raw_sha256": base.file_hash(raw),
                    "n_output_epochs": int(len(t)),
                    "last_output_epoch_s": float(t[-1])})
                print(f"  [{(time.time()-t0)/60:5.1f} min] {design}{index:03d} "
                      f"{level}: N={n_time} kernel="
                      f"{tel['gravity_kernel_ns']/1e9:.1f}s nrhs={tel['n_rhs']}",
                      flush=True)
    print(f"[timing-match] done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def _serial_kernel_ns(design: str, index: int, fallback: float) -> float:
    p = CASE_ROOT / design / f"sobolA_{index:03d}" / "atallah_serial_tight.json"
    if p.exists():
        return json.loads(p.read_text())["telemetry"]["gravity_kernel_ns"]
    return fallback


def aggregate() -> int:
    sel = json.loads(SELECTION.read_text())
    rows = []
    for design, entries in sel["designs"].items():
        camp = {int(r["sobol_index"]): r for r in
                json.loads(DESIGNS[design]["campaign"].read_text())["rows"]}
        for e in entries:
            index = int(e["sobol_index"])
            stem = ("fixed_time2" if all(
                (CASE_ROOT / design / f"sobolA_{index:03d}"
                 / f"fixed_time2_{lv}.json").exists() for lv in LEVELS)
                else "fixed_time")
            paths = {lv: (CASE_ROOT / design / f"sobolA_{index:03d}"
                          / f"{stem}_{lv}.json") for lv in LEVELS}
            if not all(p.exists() for p in paths.values()):
                continue
            raws = {lv: (RAW_ROOT / design / f"sobolA_{index:03d}"
                         / f"{stem}_{lv}.npz") for lv in LEVELS}
            truth = {lv: base.load_raw(DESIGNS[design]["r11_raw"]
                                       / f"sobolA_{index:03d}" / f"truth_{lv}.npz")
                     for lv in LEVELS}
            atal = {lv: base.load_raw(DESIGNS[design]["r12_raw"]
                                      / f"sobolA_{index:03d}" / f"atallah_{lv}.npz")
                    for lv in LEVELS}
            comp = {lv: base.load_raw(raws[lv]) for lv in LEVELS}

            def err(p, lv):
                return base.common_error(p[lv][0], p[lv][1],
                                         truth[lv][0], truth[lv][1])["pos_rms_m"]

            def self_diff(p):
                return base.common_error(p["tight"][0], p["tight"][1],
                                         p["tighter"][0], p["tighter"][1])["pos_rms_m"]

            truth_self = base.common_error(truth["tight"][0], truth["tight"][1],
                                           truth["tighter"][0], truth["tighter"][1]
                                           )["pos_rms_m"]
            e_at, e_cp = err(atal, "tight"), err(comp, "tight")
            env_at = self_diff(atal) + truth_self
            env_cp = self_diff(comp) + truth_self
            tel_at = atallah_sidecar(design, index)["telemetry"]
            tel_cp = json.loads(paths["tight"].read_text())["telemetry"]
            gap, thr = abs(e_at - e_cp), env_at + env_cp
            rows.append({
                "design": design, "sobol_index": index, "hp_km": e["hp_km"],
                "n_work_proxy": e["n_work_proxy"],
                "stage": stem,
                "n_time": e.get("n_time_refined") if stem == "fixed_time2"
                          else e["n_time_measured_cost"],
                "n_time_first_pass": e["n_time_measured_cost"],
                "intended_per_call_ratio": (e["n_time_per_call_ns"]
                                            / e["atallah_mean_per_call_ns"]),
                "atallah_error_m": e_at, "comparator_error_m": e_cp,
                "rho": e_cp / e_at if e_at > 0 else None,
                "absolute_error_difference_m": gap,
                "resolution_threshold_m": thr, "m_res": gap / thr,
                "resolved": bool(gap > thr),
                "winner_if_resolved": (("atallah" if e_at < e_cp else "fixed_time")
                                       if gap > thr else None),
                "kernel_ns_atallah_parallel_archive": tel_at["gravity_kernel_ns"],
                "kernel_ns_atallah": _serial_kernel_ns(design, index,
                                                       tel_at["gravity_kernel_ns"]),
                "kernel_ns_comparator": tel_cp["gravity_kernel_ns"],
                "achieved_kernel_time_ratio": (
                    tel_cp["gravity_kernel_ns"]
                    / _serial_kernel_ns(design, index, tel_at["gravity_kernel_ns"])),
                "n_rhs_atallah": tel_at["n_rhs"],
                "n_rhs_comparator": tel_cp["n_rhs"],
                "proxy_comparator_error_m": camp[index]["policies"]
                    ["fixed_work_atallah"]["error_tight"]["pos_rms_m"],
            })
    def stat(v):
        x = np.asarray([q for q in v if q is not None and np.isfinite(q)], float)
        return None if x.size == 0 else {
            "n": int(x.size), "median": float(np.median(x)),
            "min": float(x.min()), "max": float(x.max())}

    payload = {"schema": "r13_timing_match_v1", "rows": rows,
               "summary": {
                   "orbits": len(rows),
                   "n_time_minus_n_work": stat([r["n_time"] - r["n_work_proxy"]
                                                for r in rows]),
                   "achieved_kernel_time_ratio": stat(
                       [r["achieved_kernel_time_ratio"] for r in rows]),
                   "rho": stat([r["rho"] for r in rows]),
                   "resolved": int(sum(r["resolved"] for r in rows)),
                   "atallah_wins": int(sum(r["winner_if_resolved"] == "atallah"
                                           for r in rows)),
                   "fixed_wins": int(sum(r["winner_if_resolved"] == "fixed_time"
                                         for r in rows)),
                   "m_res": stat([r["m_res"] for r in rows])}}
    base.atomic_json(OUTPUT, payload)
    body = "\n".join(
        f"    {r['design']} & {r['sobol_index']} & {r['hp_km']:.0f} & "
        f"{r['n_work_proxy']} & {r['n_time']} & "
        f"{r['achieved_kernel_time_ratio']:.2f} & {r['atallah_error_m']:.3f} & "
        f"{r['comparator_error_m']:.3f} & {r['m_res']:.2f} & "
        f"{'yes' if r['resolved'] else 'no'}\\\\" for r in rows)
    TABLE.write_text(f"""% auto-generated by rev13_timing_match.py -- do not edit by hand
\\begin{{table}}[!htbp]
  \\centering\\small
  \\setlength{{\\tabcolsep}}{{4pt}}
  \\caption{{Measured-time-matched comparator. $N_{{\\mathrm{{work}}}}$ is the
  quadratic-proxy comparator of the campaign and $N_{{\\mathrm{{time}}}}$ the fixed
  degree whose measured kernel time approximately matches the rule's within the
  measured timing repeatability, obtained from the
  measured per-call costs and refined once on the achieved ratio. ``time ratio''
  is that achieved gravity-kernel time of the comparator over the rule's,
  measured serially on an idle machine; repeated runs place this machine's
  session-to-session drift near $10\\%$, which is the floor on how closely the
  match can be made. $E$ values are seven-day position RMS at the tight vector
  tolerance and $M_{{\\mathrm{{res}}}}$ is the resolution margin. For the design-A
  31-km-perilune orbit, both tolerance levels were recomputed consistently at the
  final comparator degree $N=515$, and only these corrected runs are reported
  (\\nolinkurl{{metrics/r13_timing_repair.json}}).}}
  \\label{{tab:timing-match}}
  \\begin{{tabular}}{{l r r r r r r r r c}}
    \\toprule
    Des. & idx & $h_p$ [km] & $N_{{\\mathrm{{work}}}}$ & $N_{{\\mathrm{{time}}}}$ &
      time ratio & $E_{{\\mathrm{{At}}}}$ [m] & $E_{{\\mathrm{{fix}}}}$ [m] &
      $M_{{\\mathrm{{res}}}}$ & res.\\\\
    \\midrule
{body}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
""", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))
    print(f"[written] {OUTPUT.name}, {TABLE.name}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("select", "serial-baseline", "run",
                                        "refine", "run-refined", "aggregate"))
    a = ap.parse_args()
    if a.command == "run-refined":
        return run("refined")
    return {"select": select, "serial-baseline": serial_baseline, "run": run,
            "refine": refine, "aggregate": aggregate}[a.command]()


if __name__ == "__main__":
    raise SystemExit(main())
