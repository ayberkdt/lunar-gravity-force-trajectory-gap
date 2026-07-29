"""Dense-altitude and bootstrap sensitivity analysis for p_tail*."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rev3_common import SEED, commit_sha, degree_power, load_model, working_tree_clean
from rev_field import _nmin_from_sigma, _nmin_proxy, sigma_a


OUT = Path(__file__).resolve().parents[1] / "metrics" / "supplemental_pstar_uncertainty.json"


def _best_p(model, power, altitudes, p_grid, weights=None):
    """Return the integer-degree SSE optimum for a deterministic altitude design."""
    empirical = np.array([
        _nmin_from_sigma(sigma_a(model, model.r_ref + h * 1e3, "vector", power), 1e-2)
        for h in altitudes
    ], dtype=int)
    predicted = np.array([
        [_nmin_proxy(model.r_ref, model.r_ref + h * 1e3, float(p), 1e-2, "vector")
         for h in altitudes]
        for p in p_grid
    ], dtype=int)
    squared = (predicted - empirical[None, :]) ** 2
    if weights is None:
        weights = np.ones(len(altitudes), dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights /= weights.sum()
    scores = squared @ weights
    best_index = int(np.argmin(scores))
    exact = p_grid[np.isclose(scores, scores[best_index], rtol=0.0, atol=1e-14)]
    return {
        "best_p": float(p_grid[best_index]),
        "exact_minimum_p_interval": [float(exact[0]), float(exact[-1])],
        "weighted_rms_degree_mismatch": float(np.sqrt(scores[best_index])),
    }


def _orbit_dwell_weights(altitudes, hp_km, ha_km, moon_radius_m):
    """Deterministic nearest-grid occupancy for uniform mean-anomaly samples."""
    a = moon_radius_m + 0.5 * (hp_km + ha_km) * 1e3
    e = (ha_km - hp_km) * 1e3 / (2.0 * a)
    mean_anomaly = (np.arange(200_000, dtype=float) + 0.5) * (2.0 * np.pi / 200_000)
    eccentric_anomaly = mean_anomaly.copy()
    for _ in range(12):
        eccentric_anomaly -= (
            eccentric_anomaly - e * np.sin(eccentric_anomaly) - mean_anomaly
        ) / (1.0 - e * np.cos(eccentric_anomaly))
    sampled_h = (a * (1.0 - e * np.cos(eccentric_anomaly)) - moon_radius_m) / 1e3
    idx = np.abs(sampled_h[:, None] - altitudes[None, :]).argmin(axis=1)
    return np.bincount(idx, minlength=len(altitudes)).astype(float)


def main() -> None:
    model = load_model(1800)
    power = degree_power(model)
    altitudes = np.arange(50.0, 300.0 + 0.1, 5.0)
    p_grid = np.round(np.arange(1.50, 2.001, 0.001), 3)
    empirical = np.array([
        _nmin_from_sigma(sigma_a(model, model.r_ref + h * 1e3, "vector", power), 1e-2)
        for h in altitudes
    ], dtype=int)
    predicted = np.array([
        [_nmin_proxy(model.r_ref, model.r_ref + h * 1e3, float(p), 1e-2, "vector")
         for h in altitudes]
        for p in p_grid
    ], dtype=int)
    squared = (predicted - empirical[None, :]) ** 2
    scores = squared.sum(axis=1)
    best_index = int(np.argmin(scores))
    min_score = int(scores[best_index])
    exact = p_grid[scores == min_score]

    rng = np.random.default_rng(SEED)
    n_boot = 2000
    counts = np.zeros((len(altitudes), n_boot), dtype=np.int16)
    for column in range(n_boot):
        counts[:, column] = np.bincount(
            rng.integers(0, len(altitudes), size=len(altitudes)),
            minlength=len(altitudes),
        )
    boot_scores = squared @ counts
    boot_p = p_grid[np.argmin(boot_scores, axis=0)]
    q = np.quantile(boot_p, [0.025, 0.5, 0.975])

    deterministic = []
    for lo, hi, step in [
        (50.0, 300.0, 2.5),
        (50.0, 300.0, 5.0),
        (50.0, 300.0, 10.0),
        (50.0, 250.0, 5.0),
        (80.0, 300.0, 5.0),
    ]:
        hs = np.arange(lo, hi + 0.1, step)
        result = _best_p(model, power, hs, p_grid)
        deterministic.append({
            "label": f"uniform_{lo:g}_{hi:g}_step_{step:g}_km",
            "altitude_range_km": [lo, hi],
            "step_km": step,
            "weighting": "uniform altitude-grid weights",
            **result,
        })

    for hp, ha in [(50.0, 300.0), (30.0, 216.0)]:
        hs = np.arange(hp, ha + 0.1, 5.0)
        weights = _orbit_dwell_weights(hs, hp, ha, model.r_ref)
        result = _best_p(model, power, hs, p_grid, weights)
        deterministic.append({
            "label": f"orbit_dwell_{hp:g}x{ha:g}_km",
            "altitude_range_km": [hp, ha],
            "step_km": 5.0,
            "weighting": "nearest-grid occupancy from 200000 uniform mean-anomaly samples",
            **result,
        })

    # Leave-one-contiguous-block-out sensitivity on the baseline 51-point grid.
    lobo = []
    for block_index, held in enumerate(np.array_split(np.arange(len(altitudes)), 5)):
        keep = np.ones(len(altitudes), dtype=bool)
        keep[held] = False
        result = _best_p(model, power, altitudes[keep], p_grid)
        lobo.append({
            "block_index": block_index,
            "held_out_altitude_km": [float(altitudes[held[0]]), float(altitudes[held[-1]])],
            **result,
        })

    sensitivity_values = [x["best_p"] for x in deterministic + lobo]
    delta_sse = len(altitudes)  # one additional mean-square degree over the grid
    accepted = p_grid[scores <= min_score + delta_sse]

    payload = {
        "method": "Dense 5-km altitude grid; integer-degree SSE; nonparametric altitude-grid bootstrap",
        "altitudes_km": altitudes.tolist(),
        "epsilon": 1e-2,
        "p_grid": p_grid.tolist(),
        "sse": scores.astype(int).tolist(),
        "empirical_nmin": empirical.tolist(),
        "predicted_nmin_at_best_p": predicted[best_index].tolist(),
        "best_p": float(p_grid[best_index]),
        "minimum_sse": min_score,
        "exact_minimum_p_interval": [float(exact[0]), float(exact[-1])],
        "rms_degree_mismatch": float(np.sqrt(min_score / len(altitudes))),
        "bootstrap": {
            "seed": SEED,
            "resamples": n_boot,
            "interpretation": "sensitivity to altitude-grid sampling, not a gravity-model covariance interval",
            "p_quantiles_2p5_50_97p5": q.tolist(),
            "p_mean": float(np.mean(boot_p)),
            "p_std": float(np.std(boot_p, ddof=1)),
        },
        "deterministic_altitude_weighting_sensitivity": {
            "interpretation": "design sensitivity, not a sampling-confidence interval",
            "configurations": deterministic,
            "leave_one_contiguous_block_out": lobo,
            "best_p_envelope": [float(min(sensitivity_values)),
                                float(max(sensitivity_values))],
            "near_optimal_objective_region": {
                "definition": "J(p) <= J_min + number_of_altitudes",
                "delta_J": int(delta_sse),
                "interpretation": "at most one additional mean-square degree of mismatch",
                "p_interval": [float(accepted[0]), float(accepted[-1])],
            },
        },
        "repo_commit_sha": commit_sha(),
        "repo_working_tree_clean": working_tree_clean(),
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: payload[k] for k in (
        "best_p", "minimum_sse", "exact_minimum_p_interval", "rms_degree_mismatch")}, indent=2))
    print("bootstrap", payload["bootstrap"])


if __name__ == "__main__":
    main()
