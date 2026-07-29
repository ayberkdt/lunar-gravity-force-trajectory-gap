"""The interpolation curve: the paper's constructive result as a figure.

Upper panel: median seven-day position error against the member's median degree
span, one curve per coverage design, with the 25--75 percentile band. Lower
panel: realized total quadratic work relative to the constant endpoint, on the
same axis, which is what stops the upper panel from being read as a free lunch.

Both endpoints are marked, because the point of the figure is that neither is
the minimum.

Usage:  python make_figures_r20.py
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
K_ALL = ("0.00", "0.25", "0.50", "0.75", "1.00")


def series(design: str):
    p = METRICS / f"r18_span_sweep_{design}_beta_1.00.json"
    if not p.exists():
        return None
    rows = json.loads(p.read_text(encoding="utf-8"))["rows"]
    span, err_med, err_lo, err_hi, work = [], [], [], [], []
    for k in K_ALL:
        e = [r["entries"][k]["error_m"] for r in rows
             if k in r["entries"] and r["entries"][k].get("error_m")]
        s = [r["entries"][k]["span"] for r in rows
             if k in r["entries"] and r["entries"][k].get("span")]
        w = [r["entries"][k].get("total_work_ratio_vs_constant") for r in rows
             if k in r["entries"]
             and r["entries"][k].get("total_work_ratio_vs_constant")]
        if not e:
            return None
        span.append(float(np.median(s)) if s else 1.0)
        err_med.append(float(np.median(e)))
        err_lo.append(float(np.percentile(e, 25)))
        err_hi.append(float(np.percentile(e, 75)))
        work.append(float(np.median(w)) if w else 1.0)
    return (np.array(span), np.array(err_med), np.array(err_lo),
            np.array(err_hi), np.array(work))


def main() -> int:
    try:
        ps.apply()
    except Exception:                                          # noqa: BLE001
        pass
    data = {d: series(d) for d in ("A", "B")}
    data = {d: v for d, v in data.items() if v}
    if not data:
        raise SystemExit("no r18 span records found")

    fig, (ax, axw) = plt.subplots(
        2, 1, figsize=(5.0, 5.2), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.12})

    styles = {"A": dict(color="#1f4e79", marker="o", ls="-"),
              "B": dict(color="#a03a2f", marker="s", ls="--")}
    for d, (span, med, lo, hi, work) in data.items():
        st = styles[d]
        ax.fill_between(span, lo, hi, color=st["color"], alpha=0.13, lw=0)
        ax.plot(span, med, label=f"design {d}", ms=5, lw=1.5, **st)
        axw.plot(span, work, ms=4, lw=1.3, **st)
        # endpoints: the two policies the literature actually offers
        ax.scatter(span[[0, -1]], med[[0, -1]], s=95, facecolor="white",
                   edgecolor=st["color"], zorder=5, lw=1.4)

    ax.set_yscale("log")
    ax.set_ylabel("median 7-day position error [m]")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.annotate("constant\ndegree", xy=(1.0, data[next(iter(data))][1][0]),
                xytext=(1.03, 6.5), fontsize=7.5, ha="left",
                arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
    sp_last = data[next(iter(data))][0][-1]
    ax.annotate("budget-calibrated\nradial endpoint",
                xy=(sp_last, data[next(iter(data))][1][-1]),
                xytext=(3.0, 40.0), fontsize=7.5, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))

    axw.axhline(1.0, color="0.6", lw=0.8, ls=":")
    axw.set_ylabel("realized total quadratic\nwork / constant")
    axw.set_xlabel("median degree span (max/min) of the propagated schedule")
    axw.set_xscale("log")
    axw.set_xticks([1, 2, 3, 5])
    axw.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    # Setting major ticks does not stop matplotlib labelling the minor decade
    # subdivisions of a log axis; over this narrow span those labels run into
    # the major ones. The shared x-axis carries the upper panel too.
    axw.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

    fig.subplots_adjust(left=0.155, right=0.975, top=0.985, bottom=0.095)
    out = FIGS / "fig_span_curve.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[written] {out.name}  designs={list(data)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
