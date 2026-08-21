"""Tables for the R18 span sweep.

Emits, into metrics/:
  r18_span_table.tex             main text, one row per k
  r18_span_detail_table_A/B.tex  supplement, one row per orbit
  r18_manuscript_descriptives.json

Usage:  python rev18_tables.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import table_design_block as tdb

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
ANCHOR_BETA = 1.00


def src_of(design: str):
    """The anchor-budget record this table reports.

    The span sweep gained a budget argument after this generator was written,
    and its records were renamed with a beta suffix at that point. The
    unsuffixed name this function used to return has not existed since, so the
    generator could no longer reproduce its own table; the file on disk was the
    one written before the rename. The budget is named explicitly here.
    """
    return METRICS / f"r18_span_sweep_{design}_beta_{ANCHOR_BETA:.2f}.json"
K_ALL = ("0.00", "0.25", "0.50", "0.75", "1.00")
ENDPOINT_LABEL = {"0.00": "constant", "1.00": "radial rule"}


def _stat(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return None
    a = np.asarray(v, dtype=float)
    return {"n": len(v), "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)),
            "p90": float(np.percentile(a, 90))}


def _resolved_vs_constant(row, k):
    """Does k beat the constant endpoint by more than the summed envelope?"""
    e0 = row["entries"]["0.00"]
    ek = row["entries"].get(k)
    if ek is None or ek.get("error_m") is None:
        return None
    thr = (e0.get("envelope_m") or 0.0) + (ek.get("envelope_m") or 0.0)
    diff = e0["error_m"] - ek["error_m"]
    if abs(diff) <= thr:
        return "unresolved"
    return "k" if diff > 0 else "constant"


def summary_table(designs: dict) -> str:
    """One block per design, so the replication is visible in the table
    rather than asserted in the caption."""
    lines = ["\\begin{tabular}{@{}c l r r r r r r@{}}", "\\toprule",
             "& $k$ & span & error [m] & work & best on & beats const. & "
             "loses to const. \\\\",
             " & & median & median & median & orbits & resolved & resolved "
             "\\\\",
             "\\midrule"]
    groups = []
    for design, d in sorted(designs.items()):
        rows = d["rows"]
        block: list[str] = []
        for k in K_ALL:
            errs = _stat([r["entries"].get(k, {}).get("error_m") for r in rows])
            spans = _stat([r["entries"].get(k, {}).get("span") for r in rows])
            works = _stat([r["entries"].get(k, {}).get(
                "total_work_ratio_vs_constant") for r in rows])
            best = sum(1 for r in rows if r["best_k"] == k)
            verdicts = [_resolved_vs_constant(r, k) for r in rows]
            wins = sum(1 for v in verdicts if v == "k")
            losses = sum(1 for v in verdicts if v == "constant")
            label = k + (f"~({ENDPOINT_LABEL[k]})" if k in ENDPOINT_LABEL
                         else "")
            span_s = f"{spans['median']:.2f}" if spans else "---"
            err_s = f"{errs['median']:.2f}" if errs else "---"
            work_s = (f"$\\times{works['median']:.3f}$" if works else "---")
            w = "---" if k == "0.00" else str(wins)
            l = "---" if k == "0.00" else str(losses)
            block.append(f"{label} & {span_s} & {err_s} & {work_s} & {best} & "
                         f"{w} & {l} \\\\")
        # the label is rotated into a five-row block, so it has to be short:
        # "Design A, 64 orbits" overfills the box. The orbit count is in the
        # caption instead.
        groups.append((f"Design {design}", block))
    lines += tdb.blocks(groups)
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def detail_table(d) -> str:
    rows = sorted(d["rows"], key=lambda r: r["hp_km"])
    lines = ["\\begin{tabular}{@{}l r r " + "r " * len(K_ALL) + "l@{}}",
             "\\toprule",
             "Orbit & $h_p$ & $N_{\\mathrm{crit}}$ & "
             + " & ".join(f"$k{{=}}{k}$" for k in K_ALL) + " & best \\\\",
             " & [km] & & " + " & ".join(["[m]"] * len(K_ALL)) + " & \\\\",
             "\\midrule"]
    for r in rows:
        cells = []
        for k in K_ALL:
            e = r["entries"].get(k, {}).get("error_m")
            cells.append(f"{e:.2f}" if e is not None else "---")
        lines.append(
            f"{r['sobol_index']:d} & {r['hp_km']:.0f} & {r['n_critical']} & "
            + " & ".join(cells) + f" & {r['best_k']} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def descriptives(d) -> dict:
    rows = d["rows"]
    out = {"orbits": len(rows), "beta": d["beta"], "by_k": {}}
    for k in K_ALL:
        verdicts = [_resolved_vs_constant(r, k) for r in rows]
        out["by_k"][k] = {
            "error_m": _stat([r["entries"].get(k, {}).get("error_m")
                              for r in rows]),
            "span": _stat([r["entries"].get(k, {}).get("span") for r in rows]),
            "best_on_orbits": sum(1 for r in rows if r["best_k"] == k),
            "resolved_beats_constant": sum(1 for v in verdicts if v == "k"),
            "resolved_loses_to_constant": sum(1 for v in verdicts
                                              if v == "constant"),
            "unresolved_vs_constant": sum(1 for v in verdicts
                                          if v == "unresolved"),
        }
    interior = [r for r in rows if r["interior_optimum"]]
    out["orbits_with_interior_best"] = len(interior)
    out["interior_best_resolved_against_constant"] = sum(
        1 for r in interior if r["best_beats_constant_resolved"])
    # work-match quality actually achieved by the bisection
    mism = []
    for r in rows:
        for k in ("0.25", "0.50", "0.75"):
            e = r["entries"].get(k, {})
            if e.get("work_mismatch") is not None:
                mism.append(abs(e["work_mismatch"]))
    out["abs_work_mismatch"] = _stat(mism)
    out["work_mismatch_max"] = max(mism) if mism else None
    return out


def main() -> int:
    designs = {}
    for design in ("A", "B"):
        p = src_of(design)
        if p.exists():
            designs[design] = json.loads(p.read_text(encoding="utf-8"))
    if not designs:
        raise SystemExit("no r18_span_sweep_*.json; run rev18_span_sweep.py")

    (METRICS / "r18_span_table.tex").write_text(summary_table(designs),
                                                encoding="utf-8")
    for design, d in designs.items():
        (METRICS / f"r18_span_detail_table_{design}.tex").write_text(
            detail_table(d), encoding="utf-8")
    desc = {design: descriptives(d) for design, d in designs.items()}
    (METRICS / "r18_manuscript_descriptives.json").write_text(
        json.dumps(desc, indent=2), encoding="utf-8")
    print("[written] r18_span_table.tex, per-design detail tables, "
          "r18_manuscript_descriptives.json")
    for design, dd in sorted(desc.items()):
        print(f"  design {design}: {dd['orbits']} orbits, interior best on "
              f"{dd['orbits_with_interior_best']}, resolved "
              f"{dd['interior_best_resolved_against_constant']}")
        for k in K_ALL:
            b = dd["by_k"][k]
            e, s = b["error_m"], b["span"]
            print(f"    k={k}  span={s['median']:.2f}  "
                  f"err={e['median']:.3f}  best={b['best_on_orbits']:2d}  "
                  f"vs const {b['resolved_beats_constant']}/"
                  f"{b['resolved_loses_to_constant']}/"
                  f"{b['unresolved_vs_constant']}")
        print(f"    |dW| median={dd['abs_work_mismatch']['median']:.5f} "
              f"max={dd['work_mismatch_max']:.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
