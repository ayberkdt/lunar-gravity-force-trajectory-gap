"""Figures promoted from the supplement to the main text (R19 editorial pass).

  fig_variational_parity.pdf   predicted vs measured equal-budget error ratio
  budget_pareto.pdf            regenerated with a third realized-work panel

Both draw only on frozen R14 records; nothing is recomputed here.

Usage:  python make_figures_r19.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

import paper_style as ps                  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
FIGS = ROOT / "figures"


def variational_parity() -> None:
    """Predicted against measured equal-budget ratio on the variational panel.

    The panel is the highest completed level of the R37 extension when that
    record exists, and the archived eight-orbit panel otherwise. Orbits whose
    propagated comparison the resolution rule leaves undecided are drawn open:
    they place a point on the plane but carry no verdict, so they are excluded
    from the sign tally (metrics/r37_scoring_amendment.json).
    """
    ext = METRICS / "r37_variational_extension.json"
    verdict = METRICS / "r37_panel_verdict.json"
    if ext.exists() and verdict.exists():
        d = json.loads(ext.read_text(encoding="utf-8"))
        v = json.loads(verdict.read_text(encoding="utf-8"))
        keep = {(c["design"], c["sobol_index"])
                for c in v["unresolved"]["orbits"]}
        n_per = v["levels"]["highest_complete_per_design"]
        pareto = json.loads((METRICS / "r14_budget_pareto.json"
                             ).read_text(encoding="utf-8"))
        members = set()
        for des in ("A", "B"):
            pr = [r for r in pareto["designs"][des]["rows"]
                  if not r["budgets"]["beta_1.00"]["censored"]]
            pr.sort(key=lambda r: r["hp_km"])
            idx = [int(i) for i in np.linspace(0, len(pr) - 1, n_per).round()]
            members |= {(des, int(pr[i]["sobol_index"])) for i in idx}
        rows = [r for r in d["rows"]
                if (r["design"], r["sobol_index"]) in members]
        panel_n = v["panel"]["orbits"]
    else:
        d = json.loads((METRICS / "r14_variational_budget.json"
                        ).read_text(encoding="utf-8"))
        rows = [r for r in d["rows"] if r.get("status") == "complete"]
        keep = set()
        panel_n = len(rows)

    pred = np.array([r["predicted_ratio_fixed_over_atallah"] for r in rows])
    meas = np.array([r["measured"]["fixed_budget"] / r["measured"]["atallah_budget"]
                     for r in rows])
    des = [r["design"] for r in rows]
    resolved = np.array([bool(r["measured"].get("resolved")) for r in rows])
    calib = np.array([r["calibration_ratio_fixed"] for r in rows
                      if r.get("calibration_ratio_fixed") is not None])

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lim = [min(pred.min(), meas.min()) * 0.45,
           max(pred.max(), meas.max()) * 2.2]
    ax.plot(lim, lim, color="0.55", lw=1.0, ls="--", zorder=1,
            label="exact agreement")
    ax.fill_between(lim, [x / 2 for x in lim], [x * 2 for x in lim],
                    color="0.88", zorder=0, label="within a factor of two")
    for design, marker in (("A", "o"), ("B", "s")):
        m = np.array([i for i, x in enumerate(des) if x == design], dtype=int)
        if m.size == 0:
            continue
        r = m[resolved[m]]
        u = m[~resolved[m]]
        if r.size:
            ax.scatter(pred[r], meas[r], marker=marker, s=26, zorder=3,
                       edgecolor="black", linewidth=0.4,
                       label=f"design {design}, resolved")
        if u.size:
            ax.scatter(pred[u], meas[u], marker=marker, s=26, zorder=2,
                       facecolor="none", edgecolor="0.55", linewidth=0.7,
                       label=f"design {design}, undecided")
    ax.axvline(1.0, color="0.75", lw=0.7, zorder=1)
    ax.axhline(1.0, color="0.75", lw=0.7, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"predicted $\rho$  (forced variational)")
    ax.set_ylabel(r"measured $\rho$  (propagated)")
    ax.legend(loc="upper left", fontsize=6, framealpha=0.9)

    ins = ax.inset_axes([0.62, 0.13, 0.34, 0.20])
    ins.axhline(1.0, color="0.55", lw=0.8, ls="--")
    ins.scatter(np.arange(len(calib)), calib, s=6, color="0.25")
    ins.set_yscale("log")
    ins.set_xticks([])
    ins.tick_params(labelsize=6)
    ins.set_title("comparator calibration", fontsize=6, pad=2)

    fig.tight_layout()
    out = FIGS / "fig_variational_parity.pdf"
    fig.savefig(out)
    plt.close(fig)
    agree = int(np.sum((pred[resolved] < 1.0) == (meas[resolved] < 1.0)))
    print(f"[written] {out.name}  (panel {panel_n} orbits, "
          f"{int(resolved.sum())} resolved, {agree}/{int(resolved.sum())} "
          f"sign agreement)")


def main() -> int:
    try:
        ps.apply()
    except Exception:                                          # noqa: BLE001
        pass
    FIGS.mkdir(exist_ok=True)
    variational_parity()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
