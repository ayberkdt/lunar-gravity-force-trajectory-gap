"""Round-2 item 3: calibration-transfer test on a second lunar field.

Applies the identical calibration procedure (Eqs. degree-RMS, tail budget,
spectral OLS, effective-exponent argmin) to the GSFC GRAIL solution
GGGRX_1200L, an independent lunar gravity model (different analysis center
and GM estimate than the JPL JGGRX_1800F used in the paper), and compares:

  * spectral slope p_spec over the observationally dominated band,
  * effective tail-budget exponent p*,
  * empirical minimum degree N_emp(h) vs. the p*-proxy.

The claim under test is procedure transfer, not number transfer: the same
recipe should yield a coherent effective exponent and a proxy that tracks
N_emp, even if the exponent value differs from JGGRX's.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.physics.spherical_harmonics import GravityModel  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "metrics"


def load_gggrx_packed(path: Path, n_max: int):
    """Parse the packed SHADR layout of GGGRX_1200L.

    Header line 1: GM [m^3/s^2], R [m], source-url.
    Coefficient lines: n, m, Cnm, Snm, sigmaC, sigmaS (fully normalized).
    Returns (mu, R, C, S) with C, S shaped (n_max+1, n_max+1).
    """
    with path.open("r") as fh:
        header = fh.readline().split(",")
        mu = float(header[0])
        R = float(header[1])
        C = np.zeros((n_max + 1, n_max + 1))
        S = np.zeros((n_max + 1, n_max + 1))
        for line in fh:
            parts = line.split(",")
            if len(parts) < 4:
                continue
            n = int(parts[0]); m = int(parts[1])
            if n > n_max:
                break
            C[n, m] = float(parts[2])
            S[n, m] = float(parts[3])
    return mu, R, C, S
GGGRX = Path(r"C:\Users\ayber\.tudat\resource\gravity_models\Moon\gggrx_1200l_sha.tab")
FIT_LO, FIT_HI = 10, 600
CAL_ALTS = np.arange(50.0, 300.0 + 0.1, 5.0)
DIAGNOSTIC_ALTS = [20.0, 30.0]


def commit_sha() -> str:
    r = subprocess.run(["git", "-c", "safe.directory=D:/Masaustu/LUNAR_SIMULATION",
                        "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return r.stdout.strip()


def degree_power(model):
    N = model.max_degree
    C, S = model.c_coeffs, model.s_coeffs
    return np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                     for n in range(N + 1)])


def sigma_a_vec(model, r, power):
    n = np.arange(len(power), dtype=np.float64)
    ratio_n = np.exp(n * math.log(model.r_ref / r))
    return (model.mu / r ** 2) * np.sqrt((n + 1.0) * (2.0 * n + 1.0)) * ratio_n * np.sqrt(power)


def nmin_emp(model, r, power, eps=1e-2):
    sig = sigma_a_vec(model, r, power)
    sq = sig ** 2
    total = float(np.sum(sq[2:]))
    if total <= 0:
        return 0
    budget = eps * eps * total
    tail = total
    for n in range(2, len(sq)):
        if tail <= budget:
            return n - 1
        tail -= sq[n]
    return len(sq) - 1


def nmin_proxy(R, r, p, eps=1e-2):
    ratio = R / r
    log_r = math.log(ratio)
    terms = []
    peak = 0.0
    n = 2
    while n <= 100_000:
        fac = math.sqrt(n + 1.0) * (2.0 * n + 1.0)
        a = fac * math.exp(n * log_r) / float(n) ** p
        t = a * a
        terms.append(t)
        peak = max(peak, t)
        if t < peak * 1e-32:
            break
        n += 1
    total = math.fsum(terms)
    budget = eps * eps * total
    tail = total
    for i, t in enumerate(terms):
        if tail <= budget:
            return max(i + 1, 1)
        tail -= t
    return n


def ols_slope(n_arr, sigma_c, lo, hi):
    m = (n_arr >= lo) & (n_arr <= hi) & (sigma_c > 0)
    x = np.log10(n_arr[m].astype(float))
    y = np.log10(sigma_c[m])
    A = np.vstack([np.ones_like(x), -x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    logK, p = coef
    resid = y - A @ coef
    dof = len(x) - 2
    s2 = float(np.sum(resid ** 2) / dof)
    cov = s2 * np.linalg.inv(A.T @ A)
    return {"p": float(p), "logK": float(logK), "p_se": float(np.sqrt(cov[1, 1])),
            "resid_rms_dex": float(np.sqrt(np.mean(resid ** 2)))}


def main() -> int:
    commit = commit_sha()
    print("commit:", commit)
    if not GGGRX.exists():
        raise SystemExit(f"second model not found: {GGGRX}")
    OUT.mkdir(exist_ok=True)

    mu, R, C, S = load_gggrx_packed(GGGRX, 700)
    model = GravityModel.from_arrays(degree_max=700, r_ref=R, mu=mu,
                                     c_coeffs_full=C, s_coeffs_full=S)
    print(f"GGGRX_1200L loaded: max_degree={model.max_degree}, "
          f"mu={model.mu:.1f}, R={model.r_ref:.0f}")
    power = degree_power(model)
    n_arr = np.arange(len(power))
    sigma_c = np.sqrt(power / (2.0 * n_arr.clip(1) + 1.0))

    fit = ols_slope(n_arr, sigma_c, FIT_LO, FIT_HI)
    fit_300 = ols_slope(n_arr, sigma_c, 10, 300)
    print(f"  spectral slope p_spec[10,600] = {fit['p']:.3f} +- {fit['p_se']:.3f} "
          f"(resid {fit['resid_rms_dex']:.3f} dex)")

    all_alts = DIAGNOSTIC_ALTS + CAL_ALTS.tolist()
    emp = {h: nmin_emp(model, model.r_ref + h * 1e3, power) for h in all_alts}
    grid = np.round(np.arange(1.40, 2.201, 0.001), 3)
    band = CAL_ALTS.tolist()

    def sse(p, hs):
        return sum((nmin_proxy(model.r_ref, model.r_ref + h * 1e3, float(p)) - emp[h]) ** 2
                   for h in hs)

    pstar = float(min(grid, key=lambda p: sse(p, band)))
    # Contiguous-band holdout on the same dense altitude grid.
    cal = CAL_ALTS[CAL_ALTS <= 100.0].tolist()
    val = CAL_ALTS[CAL_ALTS >= 150.0].tolist()
    pstar_cal = float(min(grid, key=lambda p: sse(p, cal)))
    print(f"  effective p* (50-300 km) = {pstar:.2f}; "
          f"holdout p*(50-100) = {pstar_cal:.2f}, val SSE = {sse(pstar_cal, val)}")

    rows = []
    for h in all_alts:
        r = model.r_ref + h * 1e3
        rows.append({
            "altitude_km": h,
            "emp": emp[h],
            "proxy_pstar": nmin_proxy(model.r_ref, r, pstar),
            "proxy_p1_7": nmin_proxy(model.r_ref, r, 1.7),
            "proxy_p2_0": nmin_proxy(model.r_ref, r, 2.0),
        })
        print(f"  h={h:5.0f}  emp={emp[h]:4d}  p*={rows[-1]['proxy_pstar']:4d}  "
              f"p1.7={rows[-1]['proxy_p1_7']:4d}  p2.0={rows[-1]['proxy_p2_0']:4d}")

    payload = {
        "repo_commit_sha": commit,
        "model": {"name": "GGGRX_1200L", "source": "GSFC GRAIL solution",
                  "file": str(GGGRX), "max_degree_loaded": int(model.max_degree),
                  "mu_m3_s2": float(model.mu), "reference_radius_m": float(model.r_ref)},
        "reference_model": {"name": "JGGRX_1800F", "note": "JPL, used in the paper"},
        "spectral_slope_10_600": fit,
        "spectral_slope_10_300": fit_300,
        "effective_pstar_50_300": pstar,
        "calibration_design": {"altitude_range_km": [50.0, 300.0],
                               "step_km": 5.0, "altitude_count": len(band),
                               "p_grid_step": 0.001,
                               "objective": "integer-degree SSE"},
        "effective_pstar_holdout": {"cal_band_km": cal, "val_band_km": val,
                                    "pstar_cal": pstar_cal,
                                    "sse_val_pstar": int(sse(pstar_cal, val)),
                                    "sse_val_p2_0": int(sse(2.0, val))},
        "criteria_rows": rows,
        "spectrum_arrays": {"n": n_arr[2:].tolist(),
                            "sigma_coeff_rms": sigma_c[2:].tolist()},
    }
    (OUT / "r2_gggrx_transfer.json").write_text(json.dumps(payload, indent=2),
                                                encoding="utf-8")
    print("[written] r2_gggrx_transfer.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
