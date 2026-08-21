"""Resolution margins of the undecided comparisons in (O56) and (O57).

Diagnostic only. It decides where a tighter integration would be worth its
compute, not whether a comparison is resolved: the manuscript's rule stays
M_res > 1 (Eq. 9). Every quantity here is recomputed from the sealed case
trees rather than read from a summary, so the components add up to the
threshold the campaigns used.

    E_num,P = E_self,P + E_self,ref
    M_res   = |E_P - E_Q| / (E_num,P + E_num,Q)

Usage:  python diag_unresolved_margins.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14
import rev20_span_longarc as r20
import rev48_interior_timing as r48

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
K = "0.50"


def load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def rms(a, b):
    return base.common_error(a[0], a[1], b[0], b[1])["pos_rms_m"]


def o56_rows():
    rec = json.loads((METRICS / "r56_longarc_interior.json"
                      ).read_text(encoding="utf-8"))
    out = []
    for r in rec["rows"]:
        if r["resolved"]:
            continue
        e_int, e_fix = r["member_error_m"], r["comparator_error_m"]
        s_int = r["member_self_difference_m"]
        s_fix = r["comparator_self_difference_m"]
        s_ref = r["truth_self_difference_rms_m"]
        den = s_int + s_fix + 2 * s_ref
        out.append({"campaign": "O56", "orbit": f"A{r['sobol_index']:03d}",
                    "hp_km": r["hp_km"], "E_fix": e_fix, "E_int": e_int,
                    "gap": abs(e_fix - e_int), "E_self_fix": s_fix,
                    "E_self_int": s_int, "E_self_ref": s_ref,
                    "envelope": den, "M_res": abs(e_fix - e_int) / den,
                    "time_ratio": None})
    return out


def o57_rows():
    rec = json.loads((METRICS / "r64_interior_timing_tighter.json"
                      ).read_text(encoding="utf-8"))
    out = []
    for r in rec["rows"]:
        if r["resolved"]:
            continue
        design, idx = r["design"], int(r["sobol_index"])
        span = {int(x["sobol_index"]): x for x in r48.span_rows(design)}
        e_k = span[idx]["entries"][K]

        truth = {}
        for lv in ("tight", "tighter"):
            _, raw = r14.reuse_paths(design, idx, "truth", lv)
            truth[lv] = load(raw)
        s_ref = rms(truth["tight"], truth["tighter"])

        rd = (METRICS / "r64_raw" / "interior_timing_tighter" / design
              / f"sobol{design}_{idx:03d}")
        s_fix = rms(load(rd / "fixed_time2_tight.npz"),
                    load(rd / "fixed_time2_tighter.npz"))
        # The member's own envelope is archived as E_self,int + E_self,ref.
        s_int = max(0.0, (e_k.get("envelope_m") or 0.0) - s_ref)

        e_int, e_fix = r["member_error_m"], r["comparator_error_m"]
        den = s_int + s_fix + 2 * s_ref
        out.append({"campaign": "O57", "orbit": f"{design}{idx:03d}",
                    "hp_km": r["hp_km"], "E_fix": e_fix, "E_int": e_int,
                    "gap": abs(e_fix - e_int), "E_self_fix": s_fix,
                    "E_self_int": s_int, "E_self_ref": s_ref,
                    "envelope": den, "M_res": abs(e_fix - e_int) / den,
                    "time_ratio": r.get("achieved_time_ratio")})
    return out


def band(m):
    return "HIGH" if m > 0.7 else ("MID" if m > 0.3 else "LOW")


def main() -> int:
    rows = o56_rows() + o57_rows()
    rows.sort(key=lambda r: (-r["M_res"],))
    hdr = (f"{'camp':<6}{'orbit':<8}{'hp':>7}{'E_fix':>11}{'E_int':>11}"
           f"{'gap':>11}{'Eslf_fix':>10}{'Eslf_int':>10}{'Eslf_ref':>10}"
           f"{'envelope':>11}{'M_res':>8}{'T_f/T_i':>9}  band")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        tr = "--" if r["time_ratio"] is None else f"{r['time_ratio']:.2f}"
        print(f"{r['campaign']:<6}{r['orbit']:<8}{r['hp_km']:>7.1f}"
              f"{r['E_fix']:>11.4g}{r['E_int']:>11.4g}{r['gap']:>11.4g}"
              f"{r['E_self_fix']:>10.4g}{r['E_self_int']:>10.4g}"
              f"{r['E_self_ref']:>10.4g}{r['envelope']:>11.4g}"
              f"{r['M_res']:>8.3f}{tr:>9}  {band(r['M_res'])}")
    out = ROOT / "output" / "diag_unresolved_margins.json"
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    for b in ("HIGH", "MID", "LOW"):
        sel = [r for r in rows if band(r["M_res"]) == b]
        print(f"  {b:<5} {len(sel):2d}: "
              + ", ".join(f"{r['campaign']}/{r['orbit']}" for r in sel))
    print(f"[written] {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
