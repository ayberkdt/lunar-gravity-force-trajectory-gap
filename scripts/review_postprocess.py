"""Low-cost post-processing requested by the submission review.

This script deliberately uses archived machine-readable results.  It does not
rerun an orbit propagation or require the external Lunaris checkout.

Outputs
-------
metrics/review_postprocess.json
    Fit/safe exponent results and the corrected 28-day checkpoint audit.
metrics/review_static_selection_table.tex
    Static-selection table used by the manuscript.
metrics/review_stage3_checkpoints.tex
    Weekly checkpoint table used by the supplement.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parents[1]
METRICS = BASE / "metrics"
R_REF_M = 1_738_000.0


def _load(name: str) -> dict:
    return json.loads((METRICS / name).read_text(encoding="utf-8"))


def _write(name: str, text: str) -> None:
    (METRICS / name).write_text(text, encoding="utf-8")


def _spectrum_power() -> list[float]:
    archived = _load("r1_spectrum_pfit.json")["spectrum_arrays"]
    power = [0.0] * (max(archived["n"]) + 1)
    for degree, sigma_coeff in zip(
        archived["n"], archived["sigma_coeff_rms"], strict=True
    ):
        power[degree] = float(sigma_coeff) ** 2 * (2.0 * degree + 1.0)
    return power


def _empirical_nmin(power: list[float], altitude_km: float, epsilon: float) -> int:
    radius = R_REF_M + altitude_km * 1_000.0
    log_ratio = math.log(R_REF_M / radius)
    terms = []
    for degree in range(2, len(power)):
        factor = math.sqrt((degree + 1.0) * (2.0 * degree + 1.0))
        amplitude = factor * math.exp(degree * log_ratio) * math.sqrt(power[degree])
        terms.append(amplitude * amplitude)
    total = math.fsum(terms)
    budget = epsilon * epsilon * total
    tail = total
    for index, term in enumerate(terms):
        if tail <= budget:
            return index + 1
        tail -= term
    return len(power) - 1


def _proxy_nmin_many(
    altitude_km: float, exponents: np.ndarray, epsilon: float
) -> np.ndarray:
    radius = R_REF_M + altitude_km * 1_000.0
    log_ratio = math.log(R_REF_M / radius)
    # Determine one safe finite degree range from the slowest-decaying
    # exponent, then evaluate the full p-grid in one vectorized pass.
    degree = np.arange(2.0, 10_001.0)
    log_factor = 0.5 * np.log(degree + 1.0) + np.log(2.0 * degree + 1.0)
    slow_log_terms = 2.0 * (
        log_factor + degree * log_ratio - float(exponents.min()) * np.log(degree)
    )
    slow_terms = np.exp(slow_log_terms - slow_log_terms.max())
    below = np.flatnonzero(slow_terms < 1.0e-32)
    stop = int(below[0] + 1) if below.size else len(degree)
    degree = degree[:stop]
    log_factor = log_factor[:stop]
    log_degree = np.log(degree)
    log_terms = 2.0 * (
        log_factor[None, :]
        + degree[None, :] * log_ratio
        - exponents[:, None] * log_degree[None, :]
    )
    row_max = log_terms.max(axis=1, keepdims=True)
    terms = np.exp(log_terms - row_max)
    tail = np.cumsum(terms[:, ::-1], axis=1)[:, ::-1]
    threshold = (epsilon * epsilon * tail[:, :1])
    first = np.argmax(tail <= threshold, axis=1)
    if np.any(first == 0):
        raise RuntimeError("proxy degree range was too short")
    return first + 1


def _calibrate(
    power: list[float], altitude_start: int, altitude_stop: int, epsilon: float
) -> dict:
    altitudes = list(range(altitude_start, altitude_stop + 1, 5))
    empirical = [_empirical_nmin(power, h, epsilon) for h in altitudes]
    exponents_np = np.round(np.arange(1.5, 2.001, 0.001), 3)
    exponents = exponents_np.tolist()
    predictions_np = np.column_stack(
        [_proxy_nmin_many(h, exponents_np, epsilon) for h in altitudes]
    )
    predictions = predictions_np.tolist()
    scores = [
        sum((pred - obs) ** 2 for pred, obs in zip(row, empirical, strict=True))
        for row in predictions
    ]
    best_index = min(range(len(scores)), key=scores.__getitem__)
    safe_indices = [
        i
        for i, row in enumerate(predictions)
        if all(pred >= obs for pred, obs in zip(row, empirical, strict=True))
    ]
    if not safe_indices:
        raise RuntimeError("p-grid contains no conservative exponent")
    safe_index = max(safe_indices)
    return {
        "altitude_domain_km": [altitude_start, altitude_stop],
        "altitude_step_km": 5,
        "epsilon": epsilon,
        "p_fit": exponents[best_index],
        "p_fit_rms_degree_mismatch": math.sqrt(scores[best_index] / len(altitudes)),
        "p_safe": exponents[safe_index],
        "p_safe_min_degree_margin": min(
            pred - obs
            for pred, obs in zip(predictions[safe_index], empirical, strict=True)
        ),
        "altitudes_km": altitudes,
        "empirical_nmin": empirical,
        "fit_nmin": predictions[best_index],
        "safe_nmin": predictions[safe_index],
    }


def _stage3_audit() -> dict:
    archived = _load("r7_doe_stage3_longarc.json")
    rows = []
    for orbit in archived["rows"]:
        if orbit.get("impacted"):
            rows.append(
                {
                    "name": orbit["name"],
                    "impacted": True,
                    "impact_day": orbit["impact_day"],
                }
            )
            continue
        schedules = ("sched_emp", "sched_down", "sched_up")
        best_name = min(
            schedules,
            key=lambda name: orbit["policies"][name]["pos_rms_m"],
        )
        best = orbit["policies"][best_name]
        empirical = orbit["policies"]["sched_emp"]

        def checkpoints(policy: dict) -> list[dict]:
            return [
                {
                    "day": day,
                    "cumulative_pos_rms_m": policy["checkpoints"][f"d{day}"][
                        "pos_rms_m"
                    ],
                    "instantaneous_in_track_m": policy["checkpoints"][f"d{day}"][
                        "in_track_at_checkpoint_m"
                    ],
                }
                for day in (7, 14, 21, 28)
            ]

        rows.append(
            {
                "name": orbit["name"],
                "impacted": False,
                "best_schedule": best_name,
                "best_schedule_28d_pos_rms_m": best["pos_rms_m"],
                "best_schedule_checkpoints": checkpoints(best),
                "empirical_schedule_28d_pos_rms_m": empirical["pos_rms_m"],
                "empirical_schedule_checkpoints": checkpoints(empirical),
            }
        )
    return {
        "source": "metrics/r7_doe_stage3_longarc.json",
        "finding": (
            "The previous prose mixed empirical-schedule in-track checkpoints "
            "with best-schedule full-arc RMS values."
        ),
        "rows": rows,
    }


def _static_table(calibrations: list[dict]) -> str:
    full = next(c for c in calibrations if c["altitude_domain_km"] == [20, 300])
    fit_domain = next(
        c for c in calibrations if c["altitude_domain_km"] == [50, 300]
        and c["epsilon"] == 1.0e-2
    )
    operational = next(
        c for c in calibrations if c["altitude_domain_km"] == [80, 300]
    )
    return "\n".join(
        [
            r"\begin{tabular}{@{}lll@{}}",
            r"\toprule",
            r"Use & Selection & Scope \\",
            r"\midrule",
            r"Offline high-accuracy propagation & exact empirical lookup & adopted coefficient product \\",
            rf"Compact analytic description & $p_{{\mathrm{{fit}}}}={fit_domain['p_fit']:.3f}$ & 50--300~km fit \\",
            rf"Conservative one-parameter rule & $p_{{\mathrm{{safe}}}}={full['p_safe']:.3f}$ & 20--300~km \\",
            rf"Conservative operational rule & $p_{{\mathrm{{safe}}}}={operational['p_safe']:.3f}$ & 80--300~km \\",
            r"Localized high-fidelity analysis & directional $p_{95}/p_{99}$ lookup & mission geometry \\",
            r"\bottomrule",
            r"\end{tabular}",
            "",
        ]
    )


def _stage3_table(audit: dict) -> str:
    labels = {
        "c2_50x300_polar": r"$50\times300$ polar",
        "c3_50x300_i60": r"$50\times300$, $i=60^\circ$",
        "c6_lro_30x216": r"LRO-like $30\times216$",
    }
    lines = [
        r"\begin{tabular}{@{}llrrrrrrrr@{}}",
        r"\toprule",
        r"& & \multicolumn{4}{c}{cumulative position RMS [m]} & \multicolumn{4}{c}{instantaneous in-track [m]} \\",
        r"\cmidrule(lr){3-6}\cmidrule(l){7-10}",
        r"Geometry & schedule & d7 & d14 & d21 & d28 & d7 & d14 & d21 & d28 \\",
        r"\midrule",
    ]
    for row in audit["rows"]:
        if row.get("impacted"):
            continue
        cps = row["best_schedule_checkpoints"]
        rms = [point["cumulative_pos_rms_m"] for point in cps]
        intrack = [point["instantaneous_in_track_m"] for point in cps]
        schedule = row["best_schedule"].removeprefix("sched_")
        lines.append(
            f"{labels[row['name']]} & {schedule} & "
            + " & ".join(f"{value:.1f}" for value in rms)
            + " & "
            + " & ".join(f"{value:.1f}" for value in intrack)
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    return "\n".join(lines)


def main() -> None:
    power = _spectrum_power()
    calibrations = [
        _calibrate(power, 20, 300, 1.0e-2),
        _calibrate(power, 50, 300, 1.0e-2),
        _calibrate(power, 80, 300, 1.0e-2),
        _calibrate(power, 50, 300, 1.0e-3),
    ]
    audit = _stage3_audit()
    payload = {
        "method": (
            "Post-processing of archived coefficient-spectrum and stage-3 JSON; "
            "no orbit propagation"
        ),
        "calibrations": calibrations,
        "stage3_checkpoint_audit": audit,
    }
    _write("review_postprocess.json", json.dumps(payload, indent=2) + "\n")
    _write("review_static_selection_table.tex", _static_table(calibrations))
    _write("review_stage3_checkpoints.tex", _stage3_table(audit))
    print(json.dumps({"calibrations": calibrations, "stage3": audit}, indent=2))


if __name__ == "__main__":
    main()
