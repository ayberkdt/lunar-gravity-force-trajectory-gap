"""R6-A: calibration transfer to a spectrally distinct, pre-GRAIL field.

Lunar Prospector LP150Q (Konopliv et al., 150x150, PDS jgl150q1.sha,
SHA-256 61dd4f9b...3514) predates GRAIL, carries a different GM, and is
Kaula-constrained at high degree, so its spectrum genuinely differs from
the GRAIL family. The procedure of Eqs. (3)-(7) is repeated verbatim:

  1. per-degree RMS spectrum and OLS spectral-slope fits over candidate
     observationally dominated windows;
  2. empirical minimum degrees N_emp(h) from the tail criterion
     (eps = 1e-2), with the model's 150-degree cap handled explicitly;
  3. dense-grid effective tail-budget exponent p* on the altitude range
     where N_emp is safely below the cap;
  4. cross-application: the JGGRX-calibrated exponent (1.759) applied to
     LP150Q versus LP150Q's own fit, classical p = 2, and the
     conventional attenuation threshold;
  5. a 24-hour trajectory-level fixed-degree mapping on LP150Q
     (100 km circular polar, truth N = 150) as an orbit-level check.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from rev3_common import (DAY, dump, err_stats, kernel_args, propagate,
                         warmup, REPO)

import sys

sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import GravityModel  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "data" / "jgl150q1.sha"
P_JGGRX = 1.759
EPS = 1e-2
CAP = 150


def load_lp150q():
    header = np.array(open(DATA).readline().split(","), dtype=float)
    R = header[0] * 1e3
    mu = header[1] * 1e9
    nmax = int(header[3])
    rows = np.loadtxt(DATA, delimiter=",", skiprows=1,
                      usecols=(0, 1, 2, 3))
    P = np.zeros(nmax + 1)
    for n, m, C, S in rows:
        P[int(n)] += C * C + S * S
    return R, mu, nmax, P


def sigma_n(P):
    n = np.arange(len(P), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(P / (2.0 * n + 1.0))


def ols_slope(sig, n1, n2):
    n = np.arange(len(sig), dtype=float)
    m = (n >= n1) & (n <= n2) & (sig > 0)
    x = np.log10(n[m])
    y = np.log10(sig[m])
    A = np.vstack([np.ones_like(x), -x]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    return float(coef[1]), float(np.sqrt(np.mean(resid**2)))


def n_emp(R, mu, P, h_m, eps=EPS, nmax=CAP):
    r = R + h_m
    n = np.arange(len(P), dtype=float)
    ratio = np.exp(n * math.log(R / r))
    sig2 = (n + 1.0) * (2.0 * n + 1.0) * ratio**2 * P * (mu / r**2) ** 2
    total = float(np.sum(sig2[2:nmax + 1]))
    csum = np.cumsum(sig2[::-1])[::-1]
    for N in range(2, nmax):
        if float(csum[N + 1]) <= eps * eps * total:
            return N
    return nmax


def n_spec(R, h_m, p, eps=EPS, nmax=2000):
    """One-parameter Kaula-type proxy recommendation (no coefficient data)."""
    r = R + h_m
    n = np.arange(2, nmax + 1, dtype=float)
    ratio = np.exp(n * math.log(R / r))
    a2 = (n + 1.0) * (2.0 * n + 1.0)**2 * ratio**2 * n**(-2.0 * p)
    total = float(np.sum(a2))
    csum = np.cumsum(a2[::-1])[::-1]
    for i, N in enumerate(range(2, nmax)):
        if float(csum[i + 1]) <= eps * eps * total:
            return N
    return nmax


def n_att(R, h_m, f=1e-3):
    r = R + h_m
    return int(math.ceil(math.log(f) / math.log(R / r)))


def main() -> int:
    R, mu, nmax, P = load_lp150q()
    print(f"LP150Q: R={R/1e3:.1f} km, GM={mu:.6e} m3/s2, nmax={nmax}")
    sig = sigma_n(P)

    fits = {}
    for (a, b) in ((10, 60), (10, 90), (10, 150), (30, 90)):
        p, rms = ols_slope(sig, a, b)
        fits[f"{a}_{b}"] = {"p_spec": p, "resid_rms_dex": rms}
        print(f"  OLS [{a},{b}]: p_spec={p:.3f} (resid {rms:.3f} dex)")

    alts_km = np.arange(50.0, 300.1, 5.0)
    emp = {f"{h:.0f}": n_emp(R, mu, P, h * 1e3) for h in alts_km}
    # cap safety: use only altitudes where N_emp <= 140
    fit_alts = [h for h in alts_km if emp[f"{h:.0f}"] <= 140]
    print(f"  N_emp(50km)={emp['50']}, N_emp(100km)={emp['100']}, "
          f"altitudes below cap: {len(fit_alts)}/{len(alts_km)} "
          f"(from {fit_alts[0]:.0f} km)")

    grid = np.arange(1.40, 2.2001, 0.001)
    sse = []
    for p in grid:
        s = sum((n_spec(R, h * 1e3, p) - emp[f"{h:.0f}"])**2 for h in fit_alts)
        sse.append(s)
    p_star = float(grid[int(np.argmin(sse))])
    rms_mismatch = math.sqrt(min(sse) / len(fit_alts))
    print(f"  LP150Q p* = {p_star:.3f} (RMS mismatch {rms_mismatch:.2f} deg "
          f"over {len(fit_alts)} altitudes)")

    comp_rows = []
    for h in (50.0, 80.0, 100.0, 150.0, 200.0, 300.0):
        ne = emp[f"{h:.0f}"]
        comp_rows.append({
            "altitude_km": h, "N_emp": ne,
            "proxy_p_lp150q": n_spec(R, h * 1e3, p_star),
            "proxy_p_jggrx_1.759": n_spec(R, h * 1e3, P_JGGRX),
            "proxy_p_2.0": n_spec(R, h * 1e3, 2.0),
            "atten_f1e-3": n_att(R, h * 1e3),
        })
        r = comp_rows[-1]
        print(f"  h={h:5.0f}: emp {ne:3d}  own-p* {r['proxy_p_lp150q']:3d}  "
              f"JGGRX-p {r['proxy_p_jggrx_1.759']:3d}  p2 {r['proxy_p_2.0']:3d}  "
              f"att {r['atten_f1e-3']:3d}")

    # ---- 24 h trajectory-level mapping on LP150Q
    traj = None
    try:
        model = GravityModel.from_file(str(DATA), requested_degree=150)
        args = kernel_args(model)
        warmup(model, args)
        r0 = model.r_ref + 100e3
        v0 = math.sqrt(model.mu / r0)
        y0 = np.array([r0, 0.0, 0.0, 0.0, 0.0, v0])
        t_grid = np.arange(0.0, DAY + 1.0, 120.0)
        truth, _, _ = propagate(model, y0, DAY, t_grid, lambda t, h: 150,
                                args, 1e-12, 1e-5)
        traj = []
        for N in (40, 50, 60, 70, 80, 90, 100, 120, 140):
            sol, rhs, wall = propagate(model, y0, DAY, t_grid,
                                       lambda t, h, N=N: N, args, 1e-12, 1e-5)
            st = err_stats(sol.y, truth.y)
            traj.append({"degree": N, "pos_rms_m": st["pos_rms_m"],
                         "pos_final_m": st["pos_final_m"]})
            print(f"  traj N={N:3d}: RMS {st['pos_rms_m']:8.3f} m")
    except Exception as exc:  # noqa: BLE001
        traj = {"error": f"kernel load/propagation failed: {exc}"}
        print("  trajectory mapping skipped:", exc)

    dump("r6_lp150q_transfer.json", {
        "model": {"name": "LP150Q (Lunar Prospector, Konopliv et al.)",
                  "source": "NASA PDS lp-l-rss-5-gravity-v1 jgl150q1.sha",
                  "sha256": "61dd4f9bcbfe5552171395d56d91f6e3c8bf63cf"
                            "c2f57d387f203357e2003514",
                  "R_m": R, "GM_m3_s2": mu, "nmax": nmax,
                  "note": "Kaula-constrained at high degree; pre-GRAIL"},
        "spectral_fits": fits,
        "N_emp_eps1e-2": emp,
        "p_star": {"value": p_star, "grid": [1.40, 2.20, 0.001],
                   "eps": EPS, "altitudes_used_km": [float(a) for a in fit_alts],
                   "rms_mismatch_deg": rms_mismatch,
                   "cap_policy": "altitudes with N_emp > 140 excluded "
                                 "(regularization/cap contamination)"},
        "comparison_rows": comp_rows,
        "trajectory_mapping_24h_100km": traj,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
