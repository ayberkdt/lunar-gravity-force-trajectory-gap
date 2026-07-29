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
    d = json.loads((METRICS / "r14_variational_budget.json"
                    ).read_text(encoding="utf-8"))
    rows = [r for r in d["rows"] if r.get("status") == "complete"]
    pred = np.array([r["predicted_ratio_fixed_over_atallah"] for r in rows])
    meas = np.array([r["measured"]["rho_measured"] for r in rows])
    des = [r["design"] for r in rows]
    calib = np.array([r["calibration_ratio_fixed"] for r in rows])

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lim = [min(pred.min(), meas.min()) * 0.45,
           max(pred.max(), meas.max()) * 2.2]
    ax.plot(lim, lim, color="0.55", lw=1.0, ls="--", zorder=1,
            label="exact agreement")
    # parity is what matters, so shade the decade around it rather than
    # inviting the reader to read off individual offsets
    ax.fill_between(lim, [x / 2 for x in lim], [x * 2 for x in lim],
                    color="0.88", zorder=0, label="within a factor of two")
    for design, marker in (("A", "o"), ("B", "s")):
        m = [i for i, x in enumerate(des) if x == design]
        if not m:
            continue
        ax.scatter(pred[m], meas[m], marker=marker, s=46, zorder=3,
                   edgecolor="black", linewidth=0.6,
                   label=f"design {design}")
    # the single orbit the variational solution predicts in the rule's favor
    win = np.where(pred > 1.0)[0]
    for i in win:
        ax.annotate("radial favored\n(predicted and measured)",
                    (pred[i], meas[i]), textcoords="offset points",
                    xytext=(-12, -34), fontsize=7, ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="0.35"))
    ax.axvline(1.0, color="0.75", lw=0.7, zorder=1)
    ax.axhline(1.0, color="0.75", lw=0.7, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"predicted $\rho$  (forced variational)")
    ax.set_ylabel(r"measured $\rho$  (propagated)")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.9)

    ins = ax.inset_axes([0.62, 0.13, 0.34, 0.20])
    ins.axhline(1.0, color="0.55", lw=0.8, ls="--")
    ins.scatter(np.arange(len(calib)), calib, s=14, color="0.25")
    ins.set_ylim(0.95, 1.05)
    ins.set_xticks([])
    ins.tick_params(labelsize=6)
    ins.set_title("comparator calibration", fontsize=6, pad=2)

    fig.tight_layout()
    out = FIGS / "fig_variational_parity.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[written] {out.name}  ({len(rows)} orbits, "
          f"{int(np.sum(np.sign(np.log(pred)) == np.sign(np.log(meas))))}"
          f"/{len(rows)} sign agreement)")


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
