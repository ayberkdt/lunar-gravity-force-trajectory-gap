"""P0-7: deterministic quadrature verification of the degree-RMS identity.

The per-degree acceleration RMS identity (Eq. 3) is checked with
Gauss-Legendre x uniform-longitude product quadrature, which integrates the
squared isolated-degree-n acceleration exactly (the integrand is a
trigonometric polynomial of bounded degree), and with a Sobol sample-size
convergence sequence 1500 -> 10000 -> 50000 that brackets the archived
1500-direction Monte Carlo check.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss

from rev3_common import OUT, SEED, dump, load_model, degree_power  # noqa: F401

import sys
from pathlib import Path

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import sh_potential_accel_fixed  # noqa: E402

DEGREES = (10, 60, 120, 200, 300)
ALTS_KM = (50.0, 100.0)
SOBOL_SIZES = (1500, 10000, 50000)


def sigma_analytic(model, r, power, n, kind):
    ratio_n = math.exp(n * math.log(model.r_ref / r))
    base = (model.mu / r**2) * ratio_n * math.sqrt(power[n])
    if kind == "radial":
        return (n + 1.0) * base
    return math.sqrt((n + 1.0) * (2.0 * n + 1.0)) * base


def isolated_rms(model, xyz, w, n):
    """Weighted sphere RMS of the isolated degree-n acceleration field.

    w are quadrature weights summing to 1 (uniform sphere measure)."""
    _, a = sh_potential_accel_fixed(
        xyz, model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
        degree_max=n, degree_min=n - 1)
    u = xyz / np.linalg.norm(xyz, axis=1, keepdims=True)
    vec2 = np.sum(a**2, axis=1)
    rad2 = np.sum(a * u, axis=1) ** 2
    return (math.sqrt(float(np.sum(w * vec2))),
            math.sqrt(float(np.sum(w * rad2))),
            vec2, rad2)


def gl_grid(r, n_gl, n_phi):
    """Product Gauss-Legendre (colatitude) x uniform (longitude) grid with
    normalized weights for the uniform sphere measure."""
    u_nodes, u_w = leggauss(n_gl)  # integrates du over [-1, 1]
    phi = 2.0 * math.pi * np.arange(n_phi) / n_phi
    su = np.sqrt(1.0 - u_nodes**2)
    x = r * np.outer(su, np.cos(phi)).ravel()
    y = r * np.outer(su, np.sin(phi)).ravel()
    z = r * np.repeat(u_nodes, n_phi)
    w = np.repeat(u_w, n_phi) / (2.0 * n_phi)  # sums to 1
    return np.column_stack([x, y, z]), w


def sobol_dirs(n, seed):
    from scipy.stats import qmc
    from scipy.special import ndtri
    s = qmc.Sobol(d=3, scramble=True, seed=seed)
    u = s.random(n)
    v = ndtri(np.clip(u, 1e-12, 1 - 1e-12))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def main() -> int:
    model = load_model(requested_degree=300)
    power = degree_power(model)

    gl_rows = []
    for h_km in ALTS_KM:
        r = model.r_ref + h_km * 1e3
        for n in DEGREES:
            row = {"altitude_km": h_km, "degree": n}
            for tag, n_gl, n_phi in (("exact", n + 10, 2 * n + 12),
                                     ("double", 2 * n + 10, 4 * n + 12)):
                xyz, w = gl_grid(r, n_gl, n_phi)
                rv, rr, _, _ = isolated_rms(model, xyz, w, n)
                av = sigma_analytic(model, r, power, n, "vector")
                ar = sigma_analytic(model, r, power, n, "radial")
                row[f"{tag}_n_points"] = int(len(w))
                row[f"{tag}_vector_ratio"] = rv / av
                row[f"{tag}_radial_ratio"] = rr / ar
            row["quadrature_self_consistency_vector"] = abs(
                row["exact_vector_ratio"] - row["double_vector_ratio"])
            row["quadrature_self_consistency_radial"] = abs(
                row["exact_radial_ratio"] - row["double_radial_ratio"])
            gl_rows.append(row)
            print(f"GL  h={h_km:5.0f} n={n:3d}: vec ratio "
                  f"{row['exact_vector_ratio']:.10f} rad ratio "
                  f"{row['exact_radial_ratio']:.10f} "
                  f"(self-cons {row['quadrature_self_consistency_vector']:.2e})")

    sob_rows = []
    for h_km in ALTS_KM:
        r = model.r_ref + h_km * 1e3
        for n in DEGREES:
            av = sigma_analytic(model, r, power, n, "vector")
            ar = sigma_analytic(model, r, power, n, "radial")
            for m in SOBOL_SIZES:
                xyz = r * sobol_dirs(m, SEED)
                w = np.full(m, 1.0 / m)
                rv, rr, vec2, rad2 = isolated_rms(model, xyz, w, n)
                # standard error of the mean-square estimate, propagated to RMS
                se_vec = float(np.std(vec2, ddof=1) / math.sqrt(m) / (2.0 * rv))
                se_rad = float(np.std(rad2, ddof=1) / math.sqrt(m) / (2.0 * rr))
                sob_rows.append({
                    "altitude_km": h_km, "degree": n, "n_dirs": m,
                    "vector_ratio": rv / av, "radial_ratio": rr / ar,
                    "vector_ratio_se": se_vec / av,
                    "radial_ratio_se": se_rad / ar,
                })
                print(f"Sob h={h_km:5.0f} n={n:3d} m={m:6d}: vec "
                      f"{rv/av:.5f}+-{se_vec/av:.5f} rad {rr/ar:.5f}+-{se_rad/ar:.5f}")

    dump("r3_quadrature_verification.json", {
        "seed": SEED,
        "note": "isolated degree-n via batch kernel band (n-1, n]; GL exact "
                "quadrature (nodes n+10 x 2n+12, doubled control) vs analytic "
                "radial (n+1) and vector sqrt((n+1)(2n+1)) forms; Sobol "
                "sample-size convergence with standard errors",
        "gl_rows": gl_rows,
        "sobol_rows": sob_rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
