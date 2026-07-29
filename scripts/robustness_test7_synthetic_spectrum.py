"""Robustness Test 7: synthetic-spectrum recovery of the calibration pipeline.

The central spectral claim is that the descriptive spectral slope p_spec and
the effective tail-budget exponent p*_tail are *distinct* for the real GRAIL
field.  A referee could argue this is a numerical coincidence of JGGRX.  This
test drives the identical pipeline (per-degree RMS sigma_n, OLS slope, tail
criterion N_emp(h), and the one-parameter proxy fit for p*) with synthetic
coefficient spectra whose exponents are known by construction:

  7a  pure power law  sigma_n = n^-p0
        -> both p_spec and p*_tail must return p0 (they coincide).
  7b  broken power law (steep core p1, shallow tail p2, continuous at n_b)
        -> p_spec over a low-degree window returns ~p1, while p*_tail is set
           by the high-degree tail and returns ~p2: the two separate
           *by construction*, not as a JGGRX artifact.
  7c  regularized tail  sigma_n = n^-p0 * exp(-(n/n_c)^2)  (GRAIL-like damping)
        -> p*_tail rises with lowering altitude because the damped tail
           contributes differently to the altitude-resolved budget.

Pure spectral arithmetic; no propagation.  Uses the Moon reference radius so
the (R/r)^n altitude weighting matches the production calibration.

Run: ``.venv\\Scripts\\python.exe robustness_test7_synthetic_spectrum.py``
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from rev6_lp150q import ols_slope, n_emp, n_spec

OUT = Path(__file__).resolve().parents[1] / "metrics"
R_MOON = 1738000.0
MU_MOON = 4902800306330.2
NMAX = 900
EPS = 1.0e-2
FIT_ALT_KM = np.arange(50.0, 300.1, 5.0)
P_GRID = np.arange(1.00, 3.0001, 0.001)


def sigma_to_power(sigma_n_vals):
    """Invert sigma_n = sqrt(P_n/(2n+1)) to the degree-power array P_n."""
    n = np.arange(len(sigma_n_vals), dtype=float)
    return (2.0 * n + 1.0) * sigma_n_vals ** 2


def fit_p_star(P, eps=EPS, cap=NMAX, alt_cap_degree=None):
    """Grid-search the tail-budget exponent p* that best reproduces the
    empirical N_emp(h) altitude profile with the one-parameter proxy."""
    emp = {h: n_emp(R_MOON, MU_MOON, P, h * 1e3, eps=eps, nmax=cap)
           for h in FIT_ALT_KM}
    limit = alt_cap_degree if alt_cap_degree is not None else int(0.93 * cap)
    fit_alts = [h for h in FIT_ALT_KM if emp[h] <= limit]
    if len(fit_alts) < 5:
        fit_alts = list(FIT_ALT_KM)
    sse = [sum((n_spec(R_MOON, h * 1e3, p, eps=eps, nmax=2 * cap) - emp[h]) ** 2
               for h in fit_alts) for p in P_GRID]
    j = int(np.argmin(sse))
    p_star = float(P_GRID[j])
    rms = math.sqrt(sse[j] / len(fit_alts))
    return p_star, rms, {f"{h:.0f}": int(emp[h]) for h in fit_alts}


def make_sigma(kind, **kw):
    n = np.arange(NMAX + 1, dtype=float)
    sig = np.zeros(NMAX + 1)
    nz = n >= 2
    if kind == "pure":
        p0 = kw["p0"]
        sig[nz] = n[nz] ** (-p0)
    elif kind == "broken":
        p1, p2, nb = kw["p1"], kw["p2"], kw["nb"]
        lo = nz & (n <= nb)
        hi = n > nb
        sig[lo] = n[lo] ** (-p1)
        # continuity at nb: K2 * nb^-p2 = nb^-p1  ->  K2 = nb^(p2-p1)
        sig[hi] = (nb ** (p2 - p1)) * n[hi] ** (-p2)
    elif kind == "regularized":
        p0, nc = kw["p0"], kw["nc"]
        sig[nz] = n[nz] ** (-p0) * np.exp(-(n[nz] / nc) ** 2)
    else:
        raise ValueError(kind)
    return sig


def analyze(sigma, windows):
    P = sigma_to_power(sigma)
    spec = {f"{a}_{b}": {"p_spec": (lambda r: r[0])(ols_slope(sigma, a, b)),
                         "resid_dex": ols_slope(sigma, a, b)[1]}
            for (a, b) in windows}
    p_star, rms, emp = fit_p_star(P)
    return {"p_spec_windows": spec, "p_star_tail": p_star,
            "p_star_rms_degrees": rms, "n_emp_fit_sample": emp}


def main() -> int:
    windows = [(10, 60), (10, 120), (30, 200), (200, 600)]
    payload = {
        "schema": "robustness_test7_synthetic_spectrum_v1",
        "reference_radius_m": R_MOON, "mu_m3_s2": MU_MOON,
        "nmax": NMAX, "eps_tail_fraction": EPS,
        "ols_windows": [list(w) for w in windows],
        "cases": {},
    }

    # 7a pure power laws: p_spec and p*_tail must both recover p0.
    pure = {}
    for p0 in (1.5, 1.759, 2.0, 2.5):
        res = analyze(make_sigma("pure", p0=p0), windows)
        res["injected_p0"] = p0
        pure[f"p0_{p0}"] = res
        print(f"[7a pure p0={p0}] p_spec[10,120]="
              f"{res['p_spec_windows']['10_120']['p_spec']:.3f}  "
              f"p*_tail={res['p_star_tail']:.3f}", flush=True)
    payload["cases"]["pure_power_law"] = pure

    # 7b broken power laws: p_spec (core window) != p*_tail (tail-driven).
    broken = {}
    for (p1, p2, nb) in ((2.5, 1.5, 120), (2.0, 1.6, 150), (1.6, 2.4, 150)):
        res = analyze(make_sigma("broken", p1=p1, p2=p2, nb=nb), windows)
        res["injected"] = {"p1_core": p1, "p2_tail": p2, "n_break": nb}
        broken[f"p1_{p1}_p2_{p2}_nb{nb}"] = res
        print(f"[7b broken p1={p1} p2={p2} nb={nb}] "
              f"p_spec[10,60]={res['p_spec_windows']['10_60']['p_spec']:.3f}  "
              f"p_spec[200,600]={res['p_spec_windows']['200_600']['p_spec']:.3f}  "
              f"p*_tail={res['p_star_tail']:.3f}", flush=True)
    payload["cases"]["broken_power_law"] = broken

    # 7c regularized tail: altitude dependence of p*_tail.
    reg = {}
    for nc in (300, 450, 600):
        sigma = make_sigma("regularized", p0=1.8, nc=nc)
        P = sigma_to_power(sigma)
        # altitude-resolved effective exponent: fit p* on low vs high bands
        low = [h for h in FIT_ALT_KM if h <= 120.0]
        high = [h for h in FIT_ALT_KM if h >= 180.0]

        def band_pstar(alts):
            emp = {h: n_emp(R_MOON, MU_MOON, P, h * 1e3, eps=EPS, nmax=NMAX)
                   for h in alts}
            sse = [sum((n_spec(R_MOON, h * 1e3, p, eps=EPS, nmax=2 * NMAX)
                        - emp[h]) ** 2 for h in alts) for p in P_GRID]
            return float(P_GRID[int(np.argmin(sse))])

        full = analyze(sigma, windows)
        entry = {"injected": {"p0": 1.8, "n_c": nc},
                 "p_spec_windows": full["p_spec_windows"],
                 "p_star_tail_all": full["p_star_tail"],
                 "p_star_low_alt_50_120km": band_pstar(low),
                 "p_star_high_alt_180_300km": band_pstar(high)}
        reg[f"nc_{nc}"] = entry
        print(f"[7c reg nc={nc}] p*_all={entry['p_star_tail_all']:.3f}  "
              f"p*_low={entry['p_star_low_alt_50_120km']:.3f}  "
              f"p*_high={entry['p_star_high_alt_180_300km']:.3f}", flush=True)
    payload["cases"]["regularized_tail"] = reg

    out = OUT / "robustness_test7_synthetic_spectrum.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {out.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
