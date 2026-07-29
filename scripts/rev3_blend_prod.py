"""P0-11: potential-blend strengthening — curl convergence, analytic curl,
and production-kernel timing.

1. Finite-difference curl of each transition policy at a mid-band point for
   a ladder of FD steps (0.125 m ... 8 m), demonstrating that the corrected
   potential blend's measured curl is FD-limited while the acceleration
   blend's converges to its nonzero analytic value
   curl a = grad(w) x (a_hi - a_lo).
2. Production per-call timing: the fixed/switch/acceleration-blend policies
   as production sh_accel_fixed_numba calls with a preallocated workspace,
   the corrected blend as two single-point potential+acceleration calls,
   plus the potential-path fixed control that isolates the path overhead.

Config matches r2_potential_blend.json (N_lo=30, N_hi=120, smoothstep band
between 50 and 200 km altitude).
"""

from __future__ import annotations

import math
import time

import numpy as np

from rev3_common import REPO, dump, load_model

import sys

sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    sh_accel_fixed_numba, sh_potential_accel_fixed)

N_LO, N_HI = 30, 120
ALT_NEAR, ALT_FAR = 50e3, 200e3
FD_STEPS = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
REPS = 3000


def weight_and_deriv(r, R):
    alt = r - R
    if alt <= ALT_NEAR:
        return 1.0, 0.0
    if alt >= ALT_FAR:
        return 0.0, 0.0
    t = (ALT_FAR - alt) / (ALT_FAR - ALT_NEAR)
    s = t * t * (3.0 - 2.0 * t)
    dsdt = 6.0 * t * (1.0 - t)
    return s, dsdt * (-1.0 / (ALT_FAR - ALT_NEAR))


def main() -> int:
    model = load_model(120)
    R, mu = model.r_ref, model.mu
    C, S = model.c_coeffs, model.s_coeffs

    def U_and_a(pos, degree):
        V, a = sh_potential_accel_fixed(
            np.asarray(pos, float).reshape(1, 3), C, S, mu, R, degree, -1)
        return float(V[0]), a[0]

    def fixed(pos):
        return U_and_a(pos, N_HI)[1]

    def switch(pos):
        alt = float(np.linalg.norm(pos)) - R
        return U_and_a(pos, N_HI if alt <= 0.5 * (ALT_NEAR + ALT_FAR)
                       else N_LO)[1]

    def blend_accel(pos):
        r = float(np.linalg.norm(pos))
        w, _ = weight_and_deriv(r, R)
        return (1.0 - w) * U_and_a(pos, N_LO)[1] + w * U_and_a(pos, N_HI)[1]

    def blend_pot(pos):
        r = float(np.linalg.norm(pos))
        w, dwdr = weight_and_deriv(r, R)
        U_lo, a_lo = U_and_a(pos, N_LO)
        U_hi, a_hi = U_and_a(pos, N_HI)
        a = (1.0 - w) * a_lo + w * a_hi
        if dwdr != 0.0:
            a = a + (U_hi - U_lo) * dwdr * (np.asarray(pos) / r)
        return a

    fields = {"fixed": fixed, "switch": switch, "blend_accel": blend_accel,
              "blend_pot": blend_pot}

    lat, lon = math.radians(25.0), math.radians(40.0)
    u = np.array([math.cos(lat) * math.cos(lon),
                  math.cos(lat) * math.sin(lon), math.sin(lat)])
    pos_mid = (R + 118e3) * u

    def curl_at(field, pos, h):
        J = np.zeros((3, 3))
        for j in range(3):
            e = np.zeros(3)
            e[j] = h
            J[:, j] = (field(pos + e) - field(pos - e)) / (2.0 * h)
        c = np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0],
                      J[1, 0] - J[0, 1]])
        return float(np.linalg.norm(c))

    print("== curl FD-step convergence at 118 km mid-band ==")
    curl_rows = []
    for h in FD_STEPS:
        row = {"h_fd_m": h}
        for name, f in fields.items():
            row[name] = curl_at(f, pos_mid, h)
        curl_rows.append(row)
        print("  h=%.3f m: " % h +
              "  ".join(f"{k} {row[k]:.3e}" for k in fields))

    # analytic curl of the acceleration blend: grad(w) x (a_hi - a_lo)
    r = float(np.linalg.norm(pos_mid))
    _, dwdr = weight_and_deriv(r, R)
    _, a_lo = U_and_a(pos_mid, N_LO)
    _, a_hi = U_and_a(pos_mid, N_HI)
    rhat = pos_mid / r
    curl_analytic_accel = float(np.linalg.norm(
        dwdr * np.cross(rhat, a_hi - a_lo)))
    print(f"analytic |curl| of blend_accel: {curl_analytic_accel:.3e} 1/s^2 "
          "(exact 0 for fixed/blend_pot as gradients of scalars)")

    # ---- production timing (preallocated workspace, single-point kernel)
    model300 = load_model(300)
    ws = model300.make_workspace()
    args = (model300.r_ref, model300.mu, model300.c_coeffs, model300.s_coeffs,
            model300.diag_coeffs, model300.subdiag_coeffs,
            model300.a_coeffs, model300.b_coeffs, model300.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)
    sh_accel_fixed_numba(*pos_mid, 120, *args)  # warm-up

    def timed(fn, reps=REPS):
        samples = []
        for _ in range(7):
            t0 = time.perf_counter_ns()
            for _ in range(reps):
                fn()
            samples.append((time.perf_counter_ns() - t0) / reps / 1000.0)
        samples.sort()
        return {"median_us": samples[3], "q1_us": samples[1],
                "q3_us": samples[5], "reps": reps, "blocks": 7}

    x, y, z = pos_mid

    timing = {
        "prod_accel_N120_fixed": timed(
            lambda: sh_accel_fixed_numba(x, y, z, 120, *args)),
        "prod_accel_N30_switch_far_side": timed(
            lambda: sh_accel_fixed_numba(x, y, z, 30, *args)),
        "prod_accel_blend_two_calls": timed(
            lambda: (sh_accel_fixed_numba(x, y, z, 30, *args),
                     sh_accel_fixed_numba(x, y, z, 120, *args))),
        "potential_path_fixed_N120_control": timed(
            lambda: U_and_a(pos_mid, 120), reps=300),
        "potential_blend_two_calls": timed(
            lambda: (U_and_a(pos_mid, 30), U_and_a(pos_mid, 120)), reps=300),
    }
    for k, v in timing.items():
        print(f"  {k}: {v['median_us']:.1f} us median")

    dump("r3_blend_prod.json", {
        "config": {"N_lo": N_LO, "N_hi": N_HI, "alt_near_m": ALT_NEAR,
                   "alt_far_m": ALT_FAR,
                   "probe": "118 km altitude, lat 25 deg, lon 40 deg"},
        "curl_fd_convergence": curl_rows,
        "curl_analytic_blend_accel_1_s2": curl_analytic_accel,
        "curl_analytic_note": "fixed, switch (piecewise), and the corrected "
            "potential blend are gradients of scalar potentials, so their "
            "curl is identically zero; FD values measure stencil error only",
        "production_timing_us": timing,
        "timing_note": "prod_* rows use the production single-point numba "
            "kernel with preallocated workspace; potential_* rows use the "
            "batch potential+acceleration path (one call returns both), "
            "whose per-call overhead is measured by the fixed control row",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
