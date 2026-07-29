"""Figure + tables for the R8 alpha-margin robustness control (final).

Sources:
  metrics/r8_alpha_margin.json           primary (ladder) + no-ladder (exact)
  metrics/r8_alpha_margin_down.json      worst-case Kaula-down baseline
  metrics/r8_alpha_margin_workmatch.json per-alpha re-matched work
      comparators, policy self-envelopes, and cap audit
  metrics/r8_alpha_margin_guard*.json truth self-envelopes

Products:
  figures/fig_alpha_margin.pdf           two-panel main-text figure
  metrics/r8_alpha_margin_table.tex      compact per-alpha aggregate (main text)
  metrics/r8_alpha_margin_resolution_table.tex clear rank-resolution counts
  metrics/r8_alpha_margin_supplement.tex per-orbit ladder table
  metrics/r8_alpha_margin_resolution_conservative.json derived resolution audit
plus a stdout report of the outcome classification.

Quadratic-work-proxy ratios use the re-matched comparator of the workmatch stage,
rho_work^alpha = E_workmatched(alpha)/E_sched; rankings use the
threshold-free, truth-inclusive pairwise criterion
|E_A - E_B| > E_num,A + E_num,B, where
E_num,P = E_self,P + E_self,truth.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import paper_style as ps

ROOT = Path(__file__).resolve().parents[1]
ps.apply()

prim = json.loads((ROOT / "metrics/r8_alpha_margin.json").read_text())
down = json.loads((ROOT / "metrics/r8_alpha_margin_down.json").read_text())
wm_v2_path = ROOT / "metrics/r8_alpha_margin_workmatch_v2.json"
wm_v1_path = ROOT / "metrics/r8_alpha_margin_workmatch.json"
if wm_v2_path.exists():
    candidate = json.loads(wm_v2_path.read_text())
    wm = candidate if candidate.get("complete") else json.loads(wm_v1_path.read_text())
else:
    wm = json.loads(wm_v1_path.read_text())
for d in (prim, down, wm):
    if not d.get("complete"):
        raise SystemExit("incomplete source JSON")
stage2 = json.loads((ROOT / "metrics/r7_doe_matrix_stage2.json").read_text())
s2rows = {r["name"]: r for r in stage2["rows"]}

truth_env = {}
for filename in ("r8_alpha_margin_guard.json",
                 "r8_alpha_margin_down_guard.json"):
    guard = json.loads((ROOT / "metrics" / filename).read_text())
    if not guard.get("complete"):
        raise SystemExit(f"incomplete guard artifact: {filename}")
    for row in guard["rows"]:
        value = float(row["truth_envelope_rms_m"])
        previous = truth_env.get(row["name"])
        if previous is not None and not np.isclose(previous, value,
                                                    rtol=1e-12, atol=1e-9):
            raise SystemExit(f"inconsistent truth envelope for {row['name']}")
        truth_env[row["name"]] = value

runs = {}   # (family, alpha, name) -> run record from the alpha runs
for src in (prim, down):
    for r in src["rows"]:
        for k, st in r["runs"].items():
            runs[(st["family"], round(st["alpha"], 2), r["name"])] = st
wmrec = {}  # (family, alpha, name) -> workmatch record
caps = {}
resolution_audit = []
for r in wm["rows"]:
    for k, rec in r["families"].items():
        fam, a = k.rsplit("_a", 1)
        truth = float(r.get("truth_envelope_rms_m", truth_env[r["name"]]))
        sched_self = rec.get("E_num_sched_rms_m")
        source = rec.get("sched_envelope_source", "measured")
        if sched_self is None:
            if rec["sched_wins_work_alpha"] or rec["sched_wins_crit"]:
                raise SystemExit("an unmeasured schedule envelope affects a "
                                 f"schedule-favourable pair: {r['name']} {k}")
            sched_self = 100.0
            source = "conservative_100m_bound"
        work_threshold = (float(sched_self)
                          + float(rec["E_num_work_rms_m"]) + 2.0 * truth)
        crit_threshold = (float(sched_self)
                          + float(rec["E_num_crit_rms_m"]) + 2.0 * truth)
        rec["work_rank_resolved"] = bool(
            abs(rec["E_sched_rms_m"] - rec["E_work_alpha_rms_m"])
            > work_threshold)
        rec["crit_rank_resolved"] = bool(
            abs(rec["E_sched_rms_m"] - rec["E_crit_rms_m"])
            > crit_threshold)
        resolution_audit.append({
            "orbit": r["name"], "family_alpha": k,
            "E_self_sched_rms_m": float(sched_self),
            "schedule_envelope_source": source,
            "E_self_truth_rms_m": truth,
            "work_resolution_threshold_m": work_threshold,
            "crit_resolution_threshold_m": crit_threshold,
            "work_rank_resolved": rec["work_rank_resolved"],
            "crit_rank_resolved": rec["crit_rank_resolved"],
        })
        wmrec[(fam, round(float(a), 2), r["name"])] = rec
    for k, rec in r["cap_audit"].items():
        fam, a = k.rsplit("_a", 1)
        caps[(fam, round(float(a), 2), r["name"])] = rec

NAMES = [r["name"] for r in prim["rows"]]
FAMS = {"ladder": (1.0, 1.1, 1.2, 1.3, 1.5),
        "exact": (1.0, 1.2),
        "kdown": (1.0, 1.1, 1.2, 1.3, 1.5)}

(ROOT / "metrics/r8_alpha_margin_resolution_conservative.json").write_text(
    json.dumps({
        "schema": "r8_alpha_margin_resolution_conservative_v1",
        "source_workmatch": wm_v1_path.name,
        "truth_sources": ["r8_alpha_margin_guard.json",
                          "r8_alpha_margin_down_guard.json"],
        "resolution_rule": "E_num,P = E_self,P + E_self,truth; resolved iff "
                           "|E_A-E_B| > E_num,A + E_num,B",
        "unmeasured_schedule_rule": "100 m conservative self-difference bound; "
                                    "only raw schedule losses are eligible",
        "records": resolution_audit,
        "complete": len(resolution_audit) == 288,
    }, indent=2) + "\n", encoding="utf-8")
print("[written] metrics/r8_alpha_margin_resolution_conservative.json")


def col(fam, a, get):
    return np.array([get(runs[(fam, a, n)], wmrec[(fam, a, n)])
                     for n in NAMES])


def med_band(v):
    return (float(np.median(v)), float(np.percentile(v, 10)),
            float(np.percentile(v, 90)))


def wins(fam, a):
    """(raw_w, res_w_wins, res_w_n, raw_c, res_c_wins, res_c_n)"""
    rw = rww = rwn = rc = rcw = rcn = 0
    for n in NAMES:
        x = wmrec[(fam, a, n)]
        rw += x["sched_wins_work_alpha"]
        rc += x["sched_wins_crit"]
        if x["work_rank_resolved"]:
            rwn += 1
            rww += x["sched_wins_work_alpha"]
        if x["crit_rank_resolved"]:
            rcn += 1
            rcw += x["sched_wins_crit"]
    return rw, rww, rwn, rc, rcw, rcn


# ------------------------------------------------------------------ figure
up_rc = np.array([s2rows[n]["ratios"]["sched_up"]["rho_crit"]
                  for n in NAMES])
up_gr = np.array([s2rows[n]["policies"]["sched_up"]["grav_s"]
                  / s2rows[n]["policies"]["fixed_crit"]["grav_s"]
                  for n in NAMES])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.75),
                               constrained_layout=True)

ALP = FAMS["ladder"]
for get, color, label, reskey in (
    (lambda st, x: st["rho_crit"], ps.C1, r"$\rho_{\mathrm{crit}}$",
     "crit_rank_resolved"),
    (lambda st, x: x["rho_work_alpha"], ps.C2,
     r"$\rho_{\mathrm{work}}^{\alpha}$ (re-matched)",
     "work_rank_resolved"),
):
    med = [med_band(col("ladder", a, get)) for a in ALP]
    m, lo, hi = map(np.array, zip(*med))
    unres_major = np.array([
        sum(1 for n in NAMES
            if not wmrec[("ladder", a, n)][reskey]) > 12 for a in ALP])
    ax1.plot(ALP, m, "-", color=color, label=label)
    ax1.plot(np.array(ALP)[~unres_major], m[~unres_major], "o", color=color)
    ax1.plot(np.array(ALP)[unres_major], m[unres_major], "o", mfc="white",
             color=color)
    ax1.fill_between(ALP, lo, hi, color=color, alpha=0.18, lw=0)
    ex = [float(np.median(col("exact", a, get))) for a in FAMS["exact"]]
    ax1.plot(FAMS["exact"], ex, "s", mfc="none", color=color, ms=5,
             label=label + r" (no ladder)")
ax1.axhline(1.0, color="0.35", lw=0.7, ls=":")
ax1.plot([1.0], [float(np.median(up_rc))], "D", color=ps.C5, ms=5,
         label=r"archived up-quant.\ ($\rho_{\mathrm{crit}}$)")
ax1.set_yscale("log")
ax1.set_xlabel(r"margin factor $\alpha$")
ax1.set_ylabel(r"median error ratio ($>1$: sched.\ wins)")
ax1.legend(fontsize=6.0, loc="lower right")

med = [med_band(col("ladder", a,
                    lambda st, x: st["grav_time_ratio_vs_crit_rerun"]))
       for a in ALP]
m, lo, hi = map(np.array, zip(*med))
ax2.plot(ALP, m, "o-", color=ps.C3,
         label=r"$t_{\mathrm{grav}}(\alpha)/t_{\mathrm{grav}}(N_{\min})$")
ax2.fill_between(ALP, lo, hi, color=ps.C3, alpha=0.18, lw=0)
tsw = [float(np.median(col("ladder", a,
                           lambda st, x: x["grav_s_sched"]
                           / x["grav_s_work_alpha"]))) for a in ALP]
ax2.plot(ALP, tsw, "^--", color=ps.C4, ms=4,
         label=r"$t_{\mathrm{grav}}(\alpha)/t_{\mathrm{grav}}"
               r"(N_{\mathrm{work}}^{\alpha})$")
ax2.plot([1.0], [float(np.median(up_gr))], "D", color=ps.C5, ms=5,
         label=r"archived up-quant.")
ax2.axhline(1.0, color="0.35", lw=0.7, ls=":")
ax2.set_xlabel(r"margin factor $\alpha$")
ax2.set_ylabel(r"median gravity-time ratio")
ax2.legend(fontsize=6.0, loc="lower right")

if "--tables-only" not in sys.argv:
    fig.savefig(ROOT / "figures/fig_alpha_margin.pdf")
else:
    print("[skipped] figures/fig_alpha_margin.pdf (--tables-only)")
plt.close(fig)
if "--tables-only" not in sys.argv:
    print("[written] figures/fig_alpha_margin.pdf")

# ------------------------------------------------------------- main table
FAM_LABEL = {"ladder": "ladder (q10 up)", "exact": "no ladder",
             "kdown": "Kaula-down base"}
lines = [
    r"\setlength{\tabcolsep}{4pt}",
    r"\begin{tabular}{@{}l l r r r r@{}}",
    r"\toprule",
    r"Family & $\alpha$ & "
    r"$\tilde\rho_{\mathrm{work}}^{\alpha}$ & $\tilde\rho_{\mathrm{crit}}$ & "
    r"$t_{\mathrm{grav}}$ sav. & $t_{\mathrm{s}}/t_{\mathrm{w}}$ \\",
    r"\midrule",
]
for fam, alphas in FAMS.items():
    fam_cell = FAM_LABEL[fam]
    for a in alphas:
        rwa = col(fam, a, lambda st, x: x["rho_work_alpha"])
        rc = col(fam, a, lambda st, x: st["rho_crit"])
        sv = 1.0 - col(fam, a,
                       lambda st, x: st["grav_time_ratio_vs_crit_rerun"])
        tsw_ = col(fam, a, lambda st, x: x["grav_s_sched"]
                   / x["grav_s_work_alpha"])
        lines.append(
            f"{fam_cell} & {a:.2f} & "
            f"{np.median(rwa):.2f} & {np.median(rc):.2f} & "
            f"${100 * np.median(sv):.0f}$" + r"\,\% & "
            f"{np.median(tsw_):.2f}" + r" \\")
        fam_cell = ""
    lines.append(r"\midrule")
lines[-1] = r"\bottomrule"
lines.append(r"\end{tabular}")
(ROOT / "metrics/r8_alpha_margin_table.tex").write_text(
    "\n".join(lines), encoding="utf-8")
print("[written] metrics/r8_alpha_margin_table.tex")

# --------------------------------------------- supplement resolution table
resolution = [
    r"\begin{tabular}{@{}l r l r r r r@{}}",
    r"\toprule",
    r"Family & $\alpha$ & comparator & raw schedule wins & "
    r"resolved schedule wins & resolved fixed wins & unresolved \\",
    r"\midrule",
]
for fam, alphas in FAMS.items():
    fam_cell = FAM_LABEL[fam]
    for a in alphas:
        rw, rww, rwn, rc, rcw, rcn = wins(fam, a)
        resolution.append(
            f"{fam_cell} & {a:.2f} & work matched & {rw} & {rww} & "
            f"{rwn - rww} & {len(NAMES) - rwn}" + r" \\")
        resolution.append(
            f" &  & critical altitude & {rc} & {rcw} & "
            f"{rcn - rcw} & {len(NAMES) - rcn}" + r" \\")
        fam_cell = ""
    resolution.append(r"\midrule")
resolution[-1] = r"\bottomrule"
resolution.append(r"\end{tabular}")
(ROOT / "metrics/r8_alpha_margin_resolution_table.tex").write_text(
    "\n".join(resolution), encoding="utf-8")
print("[written] metrics/r8_alpha_margin_resolution_table.tex")

# ------------------------------------------------------- supplement table
sup = [
    r"\begin{longtable}{@{}l r r r r r r r r@{}}",
    r"\caption{Per-orbit uniform-margin ladder for the 24-orbit design. For "
    r"each orbit and margin factor $\alpha$, the inflated empirical schedule's "
    r"error $E_{\mathrm{sched}}$, the work-matched degree "
    r"$N_{\mathrm{work}}^{\alpha}$, the work- and critical-altitude ratios, the "
    r"gravity-time ratio against the critical-altitude comparator, the two "
    r"pairwise resolutions (work, critical; y resolved, n unresolved), and the "
    r"realized degree range. A ratio above unity favors the schedule.}"
    r"\label{tab:alpha-margin-per-orbit}\\",
    r"\toprule",
    r"Orbit & $\alpha$ & $E_{\mathrm{sched}}$ [m] & "
    r"$N_{\mathrm{work}}^{\alpha}$ & $\rho_{\mathrm{work}}^{\alpha}$ & "
    r"$\rho_{\mathrm{crit}}$ & $t_{\mathrm{grav}}/t_{\mathrm{crit}}$ & "
    r"res.\ (w,c) & $N$ range \\",
    r"\midrule",
    r"\endfirsthead",
    r"\caption[]{Per-orbit uniform-margin ladder (continued).}\\",
    r"\toprule",
    r"Orbit & $\alpha$ & $E_{\mathrm{sched}}$ [m] & "
    r"$N_{\mathrm{work}}^{\alpha}$ & $\rho_{\mathrm{work}}^{\alpha}$ & "
    r"$\rho_{\mathrm{crit}}$ & $t_{\mathrm{grav}}/t_{\mathrm{crit}}$ & "
    r"res.\ (w,c) & $N$ range \\",
    r"\midrule",
    r"\endhead",
]
for n in NAMES:
    first = True
    for a in FAMS["ladder"]:
        st = runs[("ladder", a, n)]
        x = wmrec[("ladder", a, n)]
        nm = n.replace("_", r"\_") if first else ""
        first = False
        res = (("y" if x["work_rank_resolved"] else "n") + ","
               + ("y" if x["crit_rank_resolved"] else "n"))
        sup.append(
            f"{nm} & {a:.2f} & {st['pos_rms_m']:.1f} & "
            f"{x['n_work_alpha']} & {x['rho_work_alpha']:.3f} & "
            f"{st['rho_crit']:.3f} & "
            f"{st['grav_time_ratio_vs_crit_rerun']:.2f} & {res} & "
            f"{st['degree_table_min']}--{st['degree_table_max']} \\\\")
    sup.append(r"\addlinespace")
sup += [r"\bottomrule", r"\end{longtable}"]
(ROOT / "metrics/r8_alpha_margin_supplement.tex").write_text(
    "\n".join(sup), encoding="utf-8")
print("[written] metrics/r8_alpha_margin_supplement.tex")

# ------------------------------------------------------------ outcome dump
print("\n--- outcome (re-matched comparators, pairwise resolution) ---")
for fam, alphas in FAMS.items():
    for a in alphas:
        rwa = col(fam, a, lambda st, x: x["rho_work_alpha"])
        rw, rww, rwn, rcr, rcw, rcn = wins(fam, a)
        cb = max(caps[(fam, a, n)]["frac_time_cap_binding"] for n in NAMES)
        print(f"{fam:6s} a={a}: med_rho_work_alpha={np.median(rwa):.3f} "
              f"res wins w {rww}/{rwn} c {rcw}/{rcn} cap_max={cb:.3f}")
