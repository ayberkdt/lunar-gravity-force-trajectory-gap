"""P0-4: dense fixed-degree sweep with phase and inclination controls.

Seven-day arcs against the N=300 reference at the tight tolerance
(rtol 1e-12, atol 1e-5) for a dense degree ladder, run for three
configurations of the 50 x 300 km orbit:

  phaseA : polar, perilune start (archived r1 geometry),
  phaseB : polar, apolune start (phase-shifted control),
  inc60  : 60-degree inclination, perilune start.

For each configuration the omitted-band acceleration (band (N, 300]) is
also projected on the instantaneous along-track direction of the reference
trajectory, giving the coherence metric that the in-track drift mechanism
predicts (nonzero mean along-track omitted force).
"""

from __future__ import annotations

import math

import numpy as np

from rev3_common import (DAY, OMEGA_MOON, dump, err_stats, kernel_args,
                         load_model, eccentric_state, propagate, warmup)

import sys
from pathlib import Path

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import sh_potential_accel_fixed  # noqa: E402

DEGREES = [80, 90, 100, 106, 110, 120, 127, 130, 138, 140, 150, 180, 220]
RTOL, ATOL = 1e-12, 1e-5


def band_track_projection(model, truth, degrees, stride=5):
    """Mean/RMS of the omitted-band (N,300] acceleration projected on the
    along-track direction of the reference trajectory."""
    t = truth.t[::stride]
    pos = truth.y[:3, ::stride].T
    vel = truth.y[3:, ::stride].T
    th = OMEGA_MOON * t
    c, s = np.cos(th), np.sin(th)
    xb = c * pos[:, 0] + s * pos[:, 1]
    yb = -s * pos[:, 0] + c * pos[:, 1]
    body = np.column_stack([xb, yb, pos[:, 2]])
    vhat = vel / np.linalg.norm(vel, axis=1, keepdims=True)
    out = {}
    for N in degrees:
        _, ab = sh_potential_accel_fixed(
            body, model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
            degree_max=300, degree_min=N)
        ai = np.column_stack([c * ab[:, 0] - s * ab[:, 1],
                              s * ab[:, 0] + c * ab[:, 1], ab[:, 2]])
        proj = np.sum(ai * vhat, axis=1)
        out[str(N)] = {
            "mean_m_s2": float(np.mean(proj)),
            "rms_m_s2": float(np.sqrt(np.mean(proj ** 2))),
            "coherence_abs_mean_over_rms": float(
                abs(np.mean(proj)) / np.sqrt(np.mean(proj ** 2))),
            "band_total_rms_m_s2": float(
                np.sqrt(np.mean(np.sum(ai ** 2, axis=1)))),
        }
    return out


def daily_rms_slope(sol_y, truth_y, t_grid):
    dp = np.linalg.norm(sol_y[:3] - truth_y[:3], axis=0)
    days = np.arange(1, 8)
    rms = []
    for k in days:
        m = (t_grid > (k - 1) * DAY) & (t_grid <= k * DAY)
        rms.append(float(np.sqrt(np.mean(dp[m] ** 2))))
    lr = np.polyfit(np.log(days), np.log(np.maximum(rms, 1e-12)), 1)
    return float(lr[0]), rms


def main() -> int:
    model = load_model(300)
    args = kernel_args(model)
    warmup(model, args)
    dur = 7.0 * DAY
    t_grid = np.arange(0.0, dur + 1.0, 120.0)

    configs = {
        "phaseA": eccentric_state(model, 50.0, 300.0),
        "phaseB": eccentric_state(model, 50.0, 300.0, at_apolune=True),
        "inc60": eccentric_state(model, 50.0, 300.0, incl_deg=60.0),
    }

    payload_rows = {}
    projections = {}
    for cname, y0 in configs.items():
        print(f"== config {cname} ==")
        truth, rhs_t, wall_t = propagate(model, y0, dur, t_grid,
                                         lambda t, h: 300, args, RTOL, ATOL)
        rows = []
        for N in DEGREES:
            sol, rhs, wall = propagate(model, y0, dur, t_grid,
                                       lambda t, h, N=N: N, args, RTOL, ATOL)
            st = err_stats(sol.y, truth.y)
            slope, daily = daily_rms_slope(sol.y, truth.y, t_grid)
            st.update({"degree": N, "n_rhs": rhs.n_calls, "wall_s": wall,
                       "grav_s": rhs.grav_ns / 1e9,
                       "daily_rms_m": daily, "growth_slope_loglog": slope})
            rows.append(st)
            print(f"  N={N:3d}: RMS {st['pos_rms_m']:9.2f} m  "
                  f"final {st['pos_final_m']:9.2f} m  slope {slope:.2f}")
        payload_rows[cname] = {
            "truth": {"n_rhs": rhs_t.n_calls, "wall_s": wall_t,
                      "grav_s": rhs_t.grav_ns / 1e9},
            "rows": rows,
        }
        projections[cname] = band_track_projection(model, truth, DEGREES)

    dump("r3_degree_sweep.json", {
        "scenario": {"perilune_km": 50.0, "apolune_km": 300.0,
                     "duration_s": dur, "truth_degree": 300,
                     "integrator": "DOP853", "rtol": RTOL, "atol": ATOL,
                     "output_step_s": 120.0,
                     "rotation": "uniform sidereal about polar axis",
                     "configs": {"phaseA": "polar, perilune start",
                                 "phaseB": "polar, apolune start",
                                 "inc60": "inclination 60 deg, perilune start"}},
        "degrees": DEGREES,
        "results": payload_rows,
        "omitted_band_track_projection": projections,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
