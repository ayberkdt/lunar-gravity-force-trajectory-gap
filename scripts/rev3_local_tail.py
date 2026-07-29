"""P0-13: local/percentile complement to the sphere-RMS truncation rule.

At each altitude the discarded-tail acceleration a_tail(u) = a(600) - a(N)
is evaluated over 4096 scrambled-Sobol directions. Its distribution
(median, RMS, p95, p99, max), normalized by the sphere-RMS total
perturbation, quantifies how far mascon-region local errors exceed the
sphere average, and defines percentile-based degree recommendations
N_min(p95) alongside the sphere-RMS N_min for eps = 1e-2 and 1e-3.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rev3_common import REPO, SEED, dump, load_model, degree_power

import sys

sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.spherical_harmonics import sh_potential_accel_fixed  # noqa: E402

ALTS_KM = (30.0, 50.0, 80.0, 100.0, 150.0, 200.0)
N_LADDER = list(range(20, 290, 10))
N_MAX = 600
N_DIRS = 4096
EPSES = (1e-2, 1e-3)


def sobol_dirs(n, seed):
    from scipy.stats import qmc
    from scipy.special import ndtri
    s = qmc.Sobol(d=3, scramble=True, seed=seed)
    u = s.random(n)
    v = ndtri(np.clip(u, 1e-12, 1 - 1e-12))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def pert_total_rms(model, power, r, n_max):
    n = np.arange(len(power), dtype=np.float64)
    ratio_n = np.exp(n * math.log(model.r_ref / r))
    sig2 = (n + 1.0) * (2.0 * n + 1.0) * ratio_n**2 * power * \
        (model.mu / r**2) ** 2
    return math.sqrt(float(np.sum(sig2[2:n_max + 1])))


def main() -> int:
    model = load_model(N_MAX)
    power = degree_power(model)
    dirs = sobol_dirs(N_DIRS, SEED)

    rows = []
    recs = []
    for h_km in ALTS_KM:
        r = model.r_ref + h_km * 1e3
        xyz = r * dirs
        scale = pert_total_rms(model, power, r, N_MAX)
        _, a_full = sh_potential_accel_fixed(
            xyz, model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
            degree_max=N_MAX, degree_min=-1)
        frac = {}
        for N in N_LADDER:
            _, a_N = sh_potential_accel_fixed(
                xyz, model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
                degree_max=N, degree_min=-1)
            mag = np.linalg.norm(a_full - a_N, axis=1) / scale
            stats = {
                "altitude_km": h_km, "degree": N,
                "median": float(np.median(mag)),
                "rms": float(np.sqrt(np.mean(mag**2))),
                "p95": float(np.percentile(mag, 95)),
                "p99": float(np.percentile(mag, 99)),
                "max": float(np.max(mag)),
            }
            stats["p95_over_rms"] = stats["p95"] / stats["rms"]
            stats["max_over_rms"] = stats["max"] / stats["rms"]
            rows.append(stats)
            frac[N] = stats

        for eps in EPSES:
            def first_below(key):
                for N in N_LADDER:
                    if frac[N][key] <= eps:
                        return N
                return None
            # sphere-RMS recommendation from the sampled tail (consistent
            # estimator) and from the analytic spectrum lookup
            n_arr = np.arange(len(power), dtype=np.float64)
            ratio_n = np.exp(n_arr * math.log(model.r_ref / r))
            sig2 = (n_arr + 1.0) * (2.0 * n_arr + 1.0) * ratio_n**2 * power * \
                (model.mu / r**2) ** 2
            total = float(np.sum(sig2[2:N_MAX + 1]))
            csum = np.cumsum(sig2[::-1])[::-1]
            n_emp = None
            for N in range(2, N_MAX):
                if float(csum[N + 1]) <= eps * eps * total:
                    n_emp = N
                    break
            recs.append({
                "altitude_km": h_km, "eps": eps,
                "N_min_sampled_rms": first_below("rms"),
                "N_min_p95": first_below("p95"),
                "N_min_p99": first_below("p99"),
                "N_min_empirical_spectrum": n_emp,
            })
            print(f"h={h_km:5.0f} eps={eps:g}: N_rms(sampled)="
                  f"{recs[-1]['N_min_sampled_rms']} N_p95="
                  f"{recs[-1]['N_min_p95']} N_p99={recs[-1]['N_min_p99']} "
                  f"N_emp(spec)={n_emp}")

    dump("r3_local_tail.json", {
        "seed": SEED, "n_dirs": N_DIRS, "n_max_reference": N_MAX,
        "normalization": "per-direction |a(600)-a(N)| divided by the "
                         "analytic sphere-RMS total perturbation (n>=2..600)",
        "rows": rows,
        "recommendations": recs,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
