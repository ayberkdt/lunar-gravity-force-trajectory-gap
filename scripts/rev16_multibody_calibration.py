"""Round-16: cross-model and cross-body transfer of the calibration procedure.

Two questions are separated here.

  (1) *Cross-solution reproducibility.* The same recipe is applied to three
      independent GRAIL solutions of the same body -- JPL JGGRX_1800F (used
      throughout the paper), GSFC GRGM1200A, and GSFC GGGRX_1200L -- which
      share the underlying tracking data. Agreement is expected; disagreement
      would indicate that the procedure is sensitive to solution details it
      should not see.

  (2) *Cross-body transfer.* The same recipe is applied to five fields of four
      other bodies (Earth GOCO05c and EGM96, Mars JGMRO120D, Mercury
      JGMESS_160A, Venus SHGJ180U). Here the claim under test is that the
      *procedure* transfers while the *numbers* do not: each adopted model has
      its own effective tail-budget exponent, and the classical p = 2 is not a
      universal value.

To make the comparison geometric rather than accidental, every body is
evaluated over the same dimensionless altitude band, expressed as the
attenuation ratio R/r. The band is the one used for the lunar calibration in
the main text (50-300 km over R = 1738 km), so the geometric attenuation
ladder is identical across bodies and only the spectrum differs.

No orbit propagation is involved: this is a field-level calculation on the
coefficient sets, using the same degree-RMS identity, tail-budget criterion,
and dense-grid effective-exponent argmin as the main text.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "metrics"

# --- field inventory -------------------------------------------------------
# fmt tags describe the header layout only; the coefficient table is
# "n, m, C, S, [sigmaC, sigmaS]" in every case, fully (4*pi) normalized with no
# Condon-Shortley phase, which is the convention the paper's kernel assumes.
#   pds_km   : "R[km], GM[km^3/s^2], sigmaGM, nmax, mmax, norm, lon0, lat0"
#   pds_si   : "GM[m^3/s^2], R[m], J2, nmax, mmax, norm, lon0, lat0"
#   tudat_si : "GM[m^3/s^2], R[m], source-url"
#   plain    : "GM[m^3/s^2] R[m]"
FIELDS = [
    dict(key="JGGRX_1800F", body="Moon", center="JPL", nmax_file=1800,
         fmt="pds_km", role="reference",
         path=r"C:\Users\ayber\Desktop\Makale\lunaris_repo\data\gravity_models\jggrx_1800f_sha.tab"),
    dict(key="GRGM1200A", body="Moon", center="GSFC", nmax_file=1200,
         fmt="pds_km", role="cross_solution",
         path=r"C:\Users\ayber\Desktop\lunaris external validation\gravity_models\gggrx_1200a_sha.tab"),
    dict(key="GGGRX_1200L", body="Moon", center="GSFC", nmax_file=1200,
         fmt="tudat_si", role="cross_solution",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Moon\gggrx_1200l_sha.tab"),
    dict(key="GOCO05c", body="Earth", center="GOCO consortium", nmax_file=720,
         fmt="plain", role="cross_body",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Earth\GOCO05c.txt"),
    dict(key="EGM96", body="Earth", center="NIMA/GSFC", nmax_file=360,
         fmt="plain", role="cross_body",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Earth\egm96.txt"),
    dict(key="JGMRO120D", body="Mars", center="JPL", nmax_file=120,
         fmt="plain", role="cross_body",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Mars\jgmro120d.txt"),
    dict(key="JGMESS_160A", body="Mercury", center="JPL", nmax_file=160,
         fmt="pds_si", role="cross_body",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Mercury\jgmess_160a_sha.tab"),
    dict(key="SHGJ180U", body="Venus", center="JPL", nmax_file=180,
         fmt="pds_si", role="cross_body",
         path=r"C:\Users\ayber\.tudat\resource\gravity_models\Venus\shgj180u.a01"),
]

# Dimensionless band: the main-text lunar calibration grid, 50-300 km at 5 km
# steps over R = 1738 km, carried to every body as the attenuation ratio R/r.
LUNAR_R_KM = 1738.0
CAL_ALT_KM_MOON = np.arange(50.0, 300.0 + 0.1, 5.0)
CAL_RATIOS = LUNAR_R_KM / (LUNAR_R_KM + CAL_ALT_KM_MOON)
# Contiguous holdout: calibrate low, validate high (same split as the main text).
CAL_MASK_LOW = CAL_ALT_KM_MOON <= 100.0
CAL_MASK_HIGH = CAL_ALT_KM_MOON >= 150.0

P_GRID = np.round(np.arange(1.200, 3.2001, 0.001), 3)
EPS = 1e-2
SPEC_FIT_LO = 10
# Degrees above this are excluded from every spectral fit so that the fit band
# is set by the procedure, not by whichever model happens to extend furthest.
SPEC_FIT_HI_CAP = 600
# A model's own tail is only trustworthy well below its truncation, so the fit
# band and the reported empirical degrees are checked against this fraction.
CAP_GUARD = 0.70


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path, fmt: str, nmax: int):
    """Read a fully normalized coefficient table into (mu, R, C, S)."""
    C = np.zeros((nmax + 1, nmax + 1))
    S = np.zeros((nmax + 1, nmax + 1))
    with path.open("r") as fh:
        head = fh.readline()
        tok = [x for x in head.replace(",", " ").split() if x]
        if fmt == "pds_km":
            R, mu = float(tok[0]) * 1e3, float(tok[1]) * 1e9
        elif fmt in ("pds_si", "tudat_si", "plain"):
            mu, R = float(tok[0]), float(tok[1])
        else:
            raise ValueError(f"unknown header format {fmt!r}")
        n_seen = 0
        for line in fh:
            parts = [x for x in line.replace(",", " ").split() if x]
            if len(parts) < 4:
                continue
            n, m = int(float(parts[0])), int(float(parts[1]))
            if n > nmax:
                break
            C[n, m] = float(parts[2])
            S[n, m] = float(parts[3])
            n_seen = max(n_seen, n)
    return mu, R, C, S, n_seen


def degree_power(C, S):
    """Sigma_n^2 = sum_m (Cnm^2 + Snm^2); requires 4*pi normalization."""
    return np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                     for n in range(C.shape[0])])


def sigma_a_vec(mu, R, r, power):
    """Per-degree acceleration RMS, total-vector form (Eq. 5 of the paper)."""
    n = np.arange(len(power), dtype=np.float64)
    return (mu / r ** 2) * np.sqrt((n + 1.0) * (2.0 * n + 1.0)) \
        * np.exp(n * math.log(R / r)) * np.sqrt(power)


def nmin_emp(mu, R, r, power, eps=EPS):
    sq = sigma_a_vec(mu, R, r, power) ** 2
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


def nmin_proxy(ratio, p, eps=EPS):
    """Compact power-law tail rule at attenuation ratio R/r and exponent p."""
    log_r = math.log(ratio)
    terms, peak, n = [], 0.0, 2
    while n <= 100_000:
        a = math.sqrt(n + 1.0) * (2.0 * n + 1.0) * math.exp(n * log_r) / float(n) ** p
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
    resid = y - A @ coef
    dof = max(len(x) - 2, 1)
    cov = float(np.sum(resid ** 2) / dof) * np.linalg.inv(A.T @ A)
    return {"p": float(coef[1]), "logK": float(coef[0]),
            "p_se": float(np.sqrt(cov[1, 1])),
            "resid_rms_dex": float(np.sqrt(np.mean(resid ** 2))),
            "band": [int(lo), int(hi)], "n_points": int(m.sum())}


# Proxy degrees depend only on the ratio and the exponent, never on the body,
# so they are computed once and shared across all fields.
_PROXY_CACHE: dict[tuple[float, float], int] = {}


def proxy_cached(ratio, p):
    key = (round(ratio, 12), float(p))
    if key not in _PROXY_CACHE:
        _PROXY_CACHE[key] = nmin_proxy(ratio, p)
    return _PROXY_CACHE[key]


def calibrate(emp_by_ratio, ratios):
    """Dense-grid argmin of integer-degree SSE, as in the main text."""
    best_p, best_sse = None, math.inf
    for p in P_GRID:
        sse = sum((proxy_cached(q, float(p)) - emp_by_ratio[q]) ** 2 for q in ratios)
        if sse < best_sse:
            best_p, best_sse = float(p), sse
    return best_p, best_sse


def low_degree_dominance(mu, R, power, ratio=None):
    """Why a relative tail budget can be uninformative.

    The criterion spends a budget defined as a fraction of the *total*
    perturbing power. When one low degree carries almost all of that total, the
    budget is exhausted before the tail begins and the compact rule has nothing
    left to describe. This reports the share carried by degree 2 alone and by
    degrees 2-4, and -- as a diagnostic, not as an adopted criterion -- the
    empirical degree that results when degrees 2-4 are removed from the total
    the budget is taken relative to.
    """
    if ratio is None:
        ratio = float(CAL_RATIOS[0])
    r = R / ratio
    sq = sigma_a_vec(mu, R, r, power) ** 2
    total = float(np.sum(sq[2:]))
    if total <= 0:
        return None
    share_j2 = float(sq[2] / total)
    share_2_4 = float(np.sum(sq[2:5]) / total)

    tail_total = float(np.sum(sq[5:]))
    budget = EPS * EPS * tail_total
    tail = tail_total
    n_excl = None
    for n in range(5, len(sq)):
        if tail <= budget:
            n_excl = n - 1
            break
        tail -= sq[n]
    return {
        "ratio_R_over_r": ratio,
        "degree2_power_share": share_j2,
        "degree2_amplitude_share": float(math.sqrt(share_j2)),
        "residual_amplitude_share_above_degree2": float(math.sqrt(1.0 - share_j2)),
        "degrees_2_to_4_power_share": share_2_4,
        "n_emp_excluding_degrees_2_to_4": n_excl,
        "note": ("diagnostic only; the paper does not adopt an "
                 "oblateness-removed budget, which would be a different "
                 "criterion needing its own calibration and validation"),
    }


def analyze(field):
    path = Path(field["path"])
    if not path.exists():
        raise SystemExit(f"missing field file: {path}")
    nmax = field["nmax_file"]
    mu, R, C, S, n_seen = load(path, field["fmt"], nmax)
    power = degree_power(C, S)
    n_arr = np.arange(len(power))
    # Per-coefficient RMS: sigma_n / sqrt(2n+1).
    sigma_c = np.sqrt(power / (2.0 * n_arr.clip(1) + 1.0))

    fit_hi = min(n_seen, SPEC_FIT_HI_CAP)
    spec = ols_slope(n_arr, sigma_c, SPEC_FIT_LO, fit_hi)

    emp = {float(q): nmin_emp(mu, R, R / q, power) for q in CAL_RATIOS}
    ratios_all = [float(q) for q in CAL_RATIOS]
    ratios_low = [float(q) for q in CAL_RATIOS[CAL_MASK_LOW]]
    ratios_high = [float(q) for q in CAL_RATIOS[CAL_MASK_HIGH]]

    p_fit, sse_fit = calibrate(emp, ratios_all)
    p_low, _ = calibrate(emp, ratios_low)
    sse_hold_fit = sum((proxy_cached(q, p_low) - emp[q]) ** 2 for q in ratios_high)
    sse_hold_p2 = sum((proxy_cached(q, 2.0) - emp[q]) ** 2 for q in ratios_high)
    sse_p2 = sum((proxy_cached(q, 2.0) - emp[q]) ** 2 for q in ratios_all)

    # p_safe: smallest exponent on the grid that never underselects. A smaller
    # exponent gives a larger degree, so the search runs downward from p_fit.
    p_safe = None
    for p in sorted(P_GRID[P_GRID <= p_fit], reverse=True):
        if all(proxy_cached(q, float(p)) >= emp[q] for q in ratios_all):
            p_safe = float(p)
            break

    emp_max = max(emp.values())
    cap_ok = emp_max <= CAP_GUARD * n_seen

    rows = []
    for alt, q in zip(CAL_ALT_KM_MOON, CAL_RATIOS):
        q = float(q)
        rows.append({
            "moon_equivalent_altitude_km": float(alt),
            "ratio_R_over_r": q,
            "altitude_km": float((R / q - R) / 1e3),
            "emp": emp[q],
            "proxy_p_fit": proxy_cached(q, p_fit),
            "proxy_p2": proxy_cached(q, 2.0),
        })

    return {
        **{k: field[k] for k in ("key", "body", "center", "role", "fmt")},
        "file": str(path),
        "file_sha256": sha256(path),
        "max_degree_in_file": int(n_seen),
        "mu_m3_s2": float(mu),
        "reference_radius_m": float(R),
        "spectral_slope": spec,
        "p_fit": p_fit,
        "sse_p_fit": int(sse_fit),
        "p_safe": p_safe,
        "sse_p2": int(sse_p2),
        "holdout": {"p_fit_low_band": p_low,
                    "sse_high_band_p_fit_low": int(sse_hold_fit),
                    "sse_high_band_p2": int(sse_hold_p2)},
        "rms_mismatch_p_fit": float(math.sqrt(sse_fit / len(ratios_all))),
        "rms_mismatch_p2": float(math.sqrt(sse_p2 / len(ratios_all))),
        "emp_range": [int(min(emp.values())), int(emp_max)],
        "cap_guard_ok": bool(cap_ok),
        "cap_guard_note": (
            "" if cap_ok else
            f"max N_emp {emp_max} exceeds {CAP_GUARD:.2f} x model cap {n_seen}; "
            "the model's own tail is truncated within the criterion's band"),
        "low_degree_dominance": low_degree_dominance(mu, R, power),
        "degenerate": bool(emp_max <= 10),
        "degenerate_note": (
            "empirical degrees collapse to single digits: the degree variance is "
            "dominated by one low-degree term, so the relative tail budget is met "
            "almost immediately and the compact power-law rule is not informative"
            if emp_max <= 10 else ""),
        "criteria_rows": rows,
        "spectrum_arrays": {"n": n_arr[2:].tolist(),
                            "sigma_coeff_rms": sigma_c[2:].tolist()},
    }


def main() -> int:
    OUT.mkdir(exist_ok=True)
    results = []
    for field in FIELDS:
        print(f"--- {field['key']} ({field['body']}, {field['center']})")
        res = analyze(field)
        results.append(res)
        s = res["spectral_slope"]
        print(f"    p_spec[{s['band'][0]},{s['band'][1]}] = {s['p']:.3f} +- {s['p_se']:.3f} "
              f"(resid {s['resid_rms_dex']:.3f} dex)")
        print(f"    p_fit = {res['p_fit']:.3f}  p_safe = {res['p_safe']}  "
              f"N_emp range {res['emp_range']}  cap_ok={res['cap_guard_ok']}")
        print(f"    RMS degree mismatch: p_fit {res['rms_mismatch_p_fit']:.2f}  "
              f"p=2 {res['rms_mismatch_p2']:.2f}")
        if res["degenerate"]:
            print(f"    NOTE: {res['degenerate_note']}")
        if not res["cap_guard_ok"]:
            print(f"    NOTE: {res['cap_guard_note']}")

    payload = {
        "description": "Cross-solution and cross-body transfer of the "
                       "spectrum-calibrated truncation procedure.",
        "criterion": {"omission_budget_eps": EPS,
                      "degree_rms_form": "total-vector, sqrt((n+1)(2n+1))"},
        "band": {"definition": "dimensionless attenuation ratio R/r, matched to "
                               "the main-text lunar 50-300 km grid over R=1738 km",
                 "moon_equivalent_altitudes_km": CAL_ALT_KM_MOON.tolist(),
                 "ratios_R_over_r": [float(q) for q in CAL_RATIOS]},
        "exponent_grid": {"lo": float(P_GRID[0]), "hi": float(P_GRID[-1]),
                          "step": 0.001},
        "spectral_fit_band": {"lo": SPEC_FIT_LO, "hi_cap": SPEC_FIT_HI_CAP},
        "environment": {"python": sys.version.split()[0],
                        "numpy": np.__version__,
                        "platform": platform.platform()},
        "fields": results,
    }
    (OUT / "r16_multibody_calibration.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    print("\n[written] metrics/r16_multibody_calibration.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
