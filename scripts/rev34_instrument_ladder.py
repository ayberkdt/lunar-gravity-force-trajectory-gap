"""R34: the instrument ladder.

Recomputes, from archived records only, how each of the three force-level and
trajectory-level instruments orders the budget-calibrated radial policy against
its equal-budget constant degree at beta = 1:

  rung 1  defect RMS                    ||Delta a||        (unsigned, unweighted)
  rung 2  in-track displacement proxy   |int (T-tau) Delta a_I dtau|
                                        (signed, (T-tau)-weighted, no gradient
                                         coupling)
  rung 3  forced variational response   full Phi coupling  (8-orbit panel)
  rung 4  propagated seven-day arcs     full nonlinear dynamics along each
                                        policy's own path, plus the integrator
                                        (scored by the resolution rule)

No propagation is run here; every number is read from the frozen campaign
records of R14.  Emits metrics/r34_instrument_ladder.json and the LaTeX table
metrics/r34_instrument_ladder_table.tex.
"""

import json
import pathlib
import statistics

ROOT = pathlib.Path(__file__).resolve().parent.parent
METRICS = ROOT / "metrics"

BETA_KEY = "beta_1.00"
DESIGNS = ("A", "B")


def _median(xs):
    return statistics.median(xs) if xs else float("nan")


def force_rungs(design):
    """Rungs 1 and 2, over all 64 orbits of a design."""
    pareto = json.loads((METRICS / "r14_budget_pareto.json").read_text())
    rows = pareto["designs"][design]["rows"]

    defect_radial_better = 0
    proxy_radial_better = 0
    n = 0
    defect_ratios = []
    proxy_ratios = []

    for row in rows:
        entry = row.get("budgets", {}).get(BETA_KEY)
        if entry is None or entry.get("censored"):
            continue
        rad = entry.get("atallah", {}).get("defect")
        fix = entry.get("fixed", {}).get("defect")
        if not rad or not fix:
            continue
        n += 1

        d_rad = rad["defect_rms_m_s2"]
        d_fix = fix["defect_rms_m_s2"]
        defect_ratios.append(d_fix / d_rad)
        if d_rad < d_fix:
            defect_radial_better += 1

        # The record stores the signed integral; the proxy is its magnitude.
        p_rad = abs(rad["in_track_displacement_proxy_m"])
        p_fix = abs(fix["in_track_displacement_proxy_m"])
        proxy_ratios.append(p_fix / p_rad)
        if p_rad < p_fix:
            proxy_radial_better += 1

    return {
        "n_orbits": n,
        "defect_radial_better": defect_radial_better,
        "defect_median_ratio": _median(defect_ratios),
        "proxy_radial_better": proxy_radial_better,
        "proxy_median_ratio": _median(proxy_ratios),
    }


def propagated_rung(design):
    """Rung 4: the propagated comparison, read from the campaign's own summary.

    The verdict counts are taken from the frozen record rather than recomputed,
    so that this table cannot drift from the campaign that produced them.
    """
    traj = json.loads(
        (METRICS / f"r14_trajectory_{design}_beta_1.00.json").read_text()
    )
    s = traj["summary"]
    return {
        "n_orbits": s["orbits"],
        "resolved_fixed": s["resolved_fixed_wins"],
        "resolved_radial": s["resolved_atallah_wins"],
        "undecided": s["unresolved"],
        "raw_radial_better": s["raw_atallah_wins"],
        "median_ratio": s["rho_budget"]["median"],
    }


def variational_rung():
    """Rung 3: the forced-variational panel.

    The panel is the highest completed level of the R37 extension when that
    record exists, and the archived eight-orbit panel otherwise. Sign agreement
    is scored on resolved comparisons only, per r37_scoring_amendment.json; the
    count of orbits the instrument itself favors is over every panel orbit, so
    that this column means the same thing as the other rungs' columns.
    """
    verdict = METRICS / "r37_panel_verdict.json"
    if verdict.exists():
        v = json.loads(verdict.read_text())
        p = v["panel"]
        return {
            "source": "r37_panel_verdict.json",
            "n_orbits": p["orbits"],
            "level_per_design": v["levels"]["highest_complete_per_design"],
            "orbits_in_record": v["record"]["orbits"],
            "gradient_degree": 120,
            "resolved": p["resolved"],
            "unresolved": p["unresolved"],
            "sign_agreement": p["sign_agreement"],
            "sign_disagreement": p["sign_disagreement"],
            "radial_better_predicted": p["predicted_favors_radial_all"],
            "record_resolved": v["record"]["resolved"],
            "record_sign_agreement": v["record"]["sign_agreement"],
        }
    rec = json.loads((METRICS / "r14_variational_budget.json").read_text())
    n = sign_agree = radial_better_pred = 0
    for row in rec["rows"]:
        pred = row["predicted_ratio_fixed_over_atallah"]
        meas = row["measured"]
        meas_ratio = meas["fixed_budget"] / meas["atallah_budget"]
        n += 1
        if (pred < 1.0) == (meas_ratio < 1.0):
            sign_agree += 1
        if pred > 1.0:
            radial_better_pred += 1
    return {
        "source": "r14_variational_budget.json",
        "n_orbits": n,
        "gradient_degree": rec.get("gradient_degree"),
        "sign_agreement": sign_agree,
        "radial_better_predicted": radial_better_pred,
    }


def main():
    out = {
        "schema": "r34_instrument_ladder_v1",
        "beta": 1.0,
        "note": (
            "Derived from frozen R14 records; no propagation performed. "
            "Rungs 1, 2 and 4 are over all 64 orbits of each design; rung 3 "
            "is the forced-variational panel, whose size is read from "
            "r37_panel_verdict.json when the extension record exists."
        ),
        "designs": {},
        "variational_panel": variational_rung(),
    }
    for design in DESIGNS:
        out["designs"][design] = {
            **force_rungs(design),
            "propagated": propagated_rung(design),
        }

    (METRICS / "r34_instrument_ladder.json").write_text(
        json.dumps(out, indent=2) + "\n"
    )

    a = out["designs"]["A"]
    b = out["designs"]["B"]
    pa = a["propagated"]
    pb = b["propagated"]
    v = out["variational_panel"]
    v_rad = v["radial_better_predicted"]
    v_n = v["n_orbits"]

    lines = [
        r"\begin{tabular}{@{}l l cc cc@{}}",
        r"\toprule",
        r" & & \multicolumn{2}{c}{Orbits with lower metric under radial}"
        r" & \multicolumn{2}{c}{Median ratio} \\",
        r"\cmidrule(lr){3-4}\cmidrule(l){5-6}",
        r"Instrument & Carries & A & B & A & B \\",
        r"\midrule",
        r"Defect RMS & magnitude only & "
        f"{a['defect_radial_better']}/{a['n_orbits']} & "
        f"{b['defect_radial_better']}/{b['n_orbits']} & "
        f"{a['defect_median_ratio']:.2f} & {b['defect_median_ratio']:.2f}"
        r" \\",
        r"Proxy $D_I$ & $+$ sign, $(T-\tau)$ & "
        f"{a['proxy_radial_better']}/{a['n_orbits']} & "
        f"{b['proxy_radial_better']}/{b['n_orbits']} & "
        f"{a['proxy_median_ratio']:.2f} & {b['proxy_median_ratio']:.2f}"
        r" \\",
        r"Forced variational & $+$ gradient $\mathbf G$ & "
        rf"\multicolumn{{2}}{{c}}{{{v_rad}/{v_n}}} & "
        r"\multicolumn{2}{c}{---} \\",
        r"Propagated arcs & $+$ nonlinear dynamics, integrator & "
        f"{pa['raw_radial_better']}/{pa['n_orbits']} & "
        f"{pb['raw_radial_better']}/{pb['n_orbits']} & "
        f"{pa['median_ratio']:.2f} & {pb['median_ratio']:.2f}"
        r" \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    (METRICS / "r34_instrument_ladder_table.tex").write_text(
        "\n".join(lines) + "\n"
    )

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
