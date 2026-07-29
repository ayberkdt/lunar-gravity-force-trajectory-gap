"""Faithful implementation of the Atallah et al. (2022) radial-adaptive degree
selection, plus the paper's own verification, then the perilune-fidelity-matched
JGGRX benchmark helpers (R12).

Reference
---------
A. M. Atallah, A. Bani Younes, R. M. Woollands, J. L. Junkins,
"Analytical Radial Adaptive Method for Spherical Harmonic Gravity Models,"
J. Astronaut. Sci. 69(3):745-766, 2022.

The acceleration method (their Section 3.2). The per-degree dominating term of
the truncation-error series is their Eq. (28)/(36):

  a_hat_n(r) = (mu/R^2) (R/r)^{n+2} Lbar_{n+1} * S_n,

  S_n = ( 2 N_{n,0}/N_{n+1,1} + (n+1) N_{n,0}/N_{n+1,0} ) Chat_{n,0}
        + sum_{m=1}^n [ N_{nm}/N_{n+1,m+1} ( 1 + (n-m+2)(n-m+1) N_{nm}/N_{n+1,m-1} )
                        + (n-m+1) N_{nm}/N_{n+1,m} ] Chat_{n,m},

with
  Lbar_{n+1} = (2/sqrt(pi)) sqrt(6(n+1))                         (their Eq. 17)
  N_{nm}     = sqrt( (n-m)! (2n+1) (2-delta_{0m}) / (n+m)! )     (their Eq. 4)
  Chat_{nm}  = |Cbar_{nm}| + |Sbar_{nm}|   (normalized coefficients, Eq. 18)

The required degree (their Eq. 20/30) is

  N_req(r, tol) = min N such that  sum_{n=N+1}^{Nmax} a_hat_n(r) < tol,

where tol is an absolute acceleration tolerance [m/s^2] and Nmax is the model's
maximum degree. The method uses the ACTUAL normalized coefficients of the field,
not an analytic Kaula rule; porting to another body is only a change of
coefficients, mu, and R (their Conclusion).

S_n is radius-independent, so it is precomputed once per model; N_req(r,tol) is
then a cheap tail sum of a geometric-like sequence.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln

from rev3_common import load_model, kernel_args, warmup
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba


def _logN(n: int, m: int) -> float:
    """log of the Eq.(4) normalization factor N_{nm}."""
    if m < 0 or m > n:
        return -np.inf
    delta = 1.0 if m == 0 else 0.0
    return 0.5 * (gammaln(n - m + 1) + math.log(2 * n + 1) + math.log(2 - delta)
                  - gammaln(n + m + 1))


def precompute_Sn(model, n_max: int | None = None) -> np.ndarray:
    """Radius-independent per-degree factor  Lbar_{n+1} * S_n  for n=0..n_max.

    Returns an array ``g[n]`` such that a_hat_n(r) = (mu/R^2)(R/r)^{n+2} g[n].

    Vectorized over order m per degree n, using a precomputed log-factorial table
    so it is fast (and correct) up to the field's full degree (needed because at
    low perilune the bound decays slowly and high-degree terms still matter).
    """
    N = int(model.max_degree if n_max is None else min(n_max, model.max_degree))
    C, S = model.c_coeffs, model.s_coeffs
    lg = gammaln(np.arange(2 * N + 4) + 1.0)  # lg[k] = log(k!)

    def logN_vec(nn, mm):
        # Eq.(4): 0.5*(log((n-m)!) + log(2n+1) + log(2-delta_0m) - log((n+m)!))
        mm = np.asarray(mm)
        delta = np.where(mm == 0, 1.0, 0.0)
        return 0.5 * (lg[nn - mm] + math.log(2 * nn + 1) + np.log(2.0 - delta)
                      - lg[nn + mm])

    g = np.zeros(N + 1)
    for n in range(2, N + 1):
        Lbar = (2.0 / math.sqrt(math.pi)) * math.sqrt(6.0 * (n + 1))
        m = np.arange(0, n + 1)
        Chat = np.abs(C[n, m]) + np.abs(S[n, m])
        logNnm = logN_vec(n, m)
        # ratios to degree n+1 normalization factors
        r_mp1 = np.exp(logNnm - logN_vec(n + 1, m + 1))          # N_nm / N_{n+1,m+1}
        r_m = np.exp(logNnm - logN_vec(n + 1, m))                # N_nm / N_{n+1,m}
        # N_nm / N_{n+1,m-1}: only defined for m>=1 (m-1>=0)
        r_mm1 = np.zeros(n + 1)
        if n >= 1:
            mm = m[1:]
            r_mm1[1:] = np.exp(logNnm[1:] - logN_vec(n + 1, mm - 1))
        fac = (n - m + 2) * (n - m + 1)
        weight = r_mp1 * (1.0 + fac * r_mm1) + (n - m + 1) * r_m
        # m = 0 uses the special first-term form (Eq. 28 first line)
        weight[0] = (2.0 * np.exp(logNnm[0] - logN_vec(n + 1, 1))
                     + (n + 1) * np.exp(logNnm[0] - logN_vec(n + 1, 0)))
        g[n] = Lbar * float(np.sum(weight * Chat))
    return g


def a_hat_series(r_m: float, model, g: np.ndarray) -> np.ndarray:
    """a_hat_n(r) for n=0..len(g)-1 [m/s^2]."""
    n = np.arange(len(g), dtype=float)
    ratio = model.r_ref / r_m
    return (model.mu / model.r_ref**2) * ratio**(n + 2.0) * g


def n_req(r_m: float, tol: float, model, g: np.ndarray,
          floor: int = 2, cap: int | None = None) -> int:
    """min N such that sum_{n>N} a_hat_n(r) < tol  (their Eq. 20/30)."""
    a = a_hat_series(r_m, model, g)
    # tail[N] = sum_{n>N} a_hat_n
    tail = np.concatenate([np.cumsum(a[::-1])[::-1][1:], [0.0]])
    idx = np.argmax(tail < tol)  # first N with tail < tol
    if not (tail < tol).any():
        N = len(g) - 1
    else:
        N = int(idx)
    hi = (len(g) - 1) if cap is None else min(cap, len(g) - 1)
    return int(min(hi, max(floor, N)))


def tol_for_degree(r_m: float, target_N: int, model, g: np.ndarray) -> float:
    """tol that makes n_req(r, tol) == target_N (perilune-fidelity matching).

    N_req = target_N requires sum_{n>target_N} a_hat_n < tol <= sum_{n>=target_N}.
    Return a tol just above the lower bracket so N_req resolves to target_N.
    """
    a = a_hat_series(r_m, model, g)
    tail_above = float(np.sum(a[target_N + 1:]))
    return tail_above * (1.0 + 1e-9)


def atallah_degree_fn(model, g: np.ndarray, tol: float, floor: int = 2,
                      cap: int | None = None):
    """Return degree_of(t, h_m) evaluating N_req at the current radius (exact)."""
    r_ref = model.r_ref

    def degree_of(t, h_m, _tol=tol):
        return n_req(r_ref + h_m, _tol, model, g, floor, cap)

    return degree_of


def atallah_binned_schedule(model, g: np.ndarray, tol: float, hp_km: float,
                            ha_km: float, floor: int = 2, cap: int | None = None,
                            bin_km: float = 10.0):
    """Atallah N_req(r,tol) discretized onto ``bin_km`` altitude bins, matching
    the binning convention of the paper's other altitude schedules.

    Binning (rather than per-RHS-call evaluation at the exact radius) keeps the
    truncation degree constant within an integration step, which avoids the
    within-step degree inconsistency that otherwise collapses the adaptive step
    size. Returns (degree_of(t, h_m), table) where ``table`` maps bin-floor
    altitude [km] to the selected degree.
    """
    lo = math.floor((hp_km - 1.0) / bin_km) * bin_km
    hi = math.ceil((ha_km + 1.0) / bin_km) * bin_km
    bins = np.arange(max(0.0, lo), hi + bin_km, bin_km)
    table = {}
    for hk in bins:
        table[float(hk)] = n_req(model.r_ref + hk * 1e3, tol, model, g, floor, cap)
    hmax = max(table)
    hmin = min(table)

    def degree_of(t, h_m):
        hb = min(hmax, max(hmin, bin_km * math.floor(h_m / 1e3 / bin_km)))
        return table[hb]

    return degree_of, table


# ------------------------------------------------------------------ verification
def actual_degree_contribution_max(model, args, r_m: float, n: int,
                                   n_lat: int = 19, n_lon: int = 24) -> float:
    """max over a lat/lon grid of the degree-n acceleration contribution
    magnitude |a_n| = |a_{<=n} - a_{<=n-1}| on the real field."""
    lats = np.linspace(-math.pi / 2 + 1e-3, math.pi / 2 - 1e-3, n_lat)
    lons = np.linspace(0.0, 2 * math.pi, n_lon, endpoint=False)
    worst = 0.0
    for phi in lats:
        for lam in lons:
            x = r_m * math.cos(phi) * math.cos(lam)
            y = r_m * math.cos(phi) * math.sin(lam)
            z = r_m * math.sin(phi)
            an = np.array(sh_accel_fixed_numba(x, y, z, n, *args))
            anm1 = np.array(sh_accel_fixed_numba(x, y, z, n - 1, *args))
            worst = max(worst, float(np.linalg.norm(an - anm1)))
    return worst


def actual_truncation_error_max(model, args, r_m: float, N: int, N_max: int,
                                n_lat: int = 19, n_lon: int = 24) -> float:
    """max over lat/lon of |a_{N_max} - a_N| on the real field [m/s^2]."""
    lats = np.linspace(-math.pi / 2 + 1e-3, math.pi / 2 - 1e-3, n_lat)
    lons = np.linspace(0.0, 2 * math.pi, n_lon, endpoint=False)
    worst = 0.0
    for phi in lats:
        for lam in lons:
            x = r_m * math.cos(phi) * math.cos(lam)
            y = r_m * math.cos(phi) * math.sin(lam)
            z = r_m * math.sin(phi)
            aN = np.array(sh_accel_fixed_numba(x, y, z, N, *args))
            aM = np.array(sh_accel_fixed_numba(x, y, z, N_max, *args))
            worst = max(worst, float(np.linalg.norm(aM - aN)))
    return worst
