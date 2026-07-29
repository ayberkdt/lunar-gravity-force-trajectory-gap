"""Revision-1 field-level experiments (commit-pinned, seeded).

R1-B  Degree-RMS formulation check: analytic radial-only and total-vector
      per-degree RMS versus direct seeded sampling of the isolated degree-n
      acceleration field (production batch kernel, degree_min feature).
R1-C  Formal p calibration: OLS fit of log10 sigma_n = log10 K - p log10 n
      over restricted and full degree ranges, block-bootstrap CIs, holdout
      validation, usage-based argmin, N_min sensitivity to p and epsilon.
R1-E2 Truncation criteria recomputed with radial-only AND total-vector
      weightings on the current commit.
R1-E1 Band shares with bootstrap 95% CIs, a second seed, and a Sobol
      low-discrepancy sampling control at 80 km.
R1-E5 Switching jump with bootstrap 95% CIs.
R1-E3 Kernel timing: median and IQR over repeated blocks (not best-only).
R1-E4 Blend analysis rerun (curl, omitted term) on the current commit.

Outputs are written as r1_*.json into the paper's metrics/ directory.
"""

from __future__ import annotations

import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.common.lunar_data import resolve_lunar_gravity_path  # noqa: E402
from lunaris.common.math_utils import recommended_sh_degree  # noqa: E402
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel,
    sh_accel_adaptive_blend_numba,
    sh_accel_fixed_numba,
    sh_potential_accel_fixed,
    _apply_smoothstep,
    _compute_sh_acceleration_dual_numba,
)

OUT = Path(__file__).resolve().parents[1] / "metrics"
SEED = 20260719
SEED_ALT = 990417


def commit_sha() -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return r.stdout.strip()


def _dump(name: str, payload: dict) -> None:
    payload["repo_commit_sha"] = COMMIT
    (OUT / name).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[written] {name}")


def _dirs(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _sobol_dirs(n: int, seed: int) -> np.ndarray:
    from scipy.stats import qmc
    from scipy.special import ndtri
    s = qmc.Sobol(d=3, scramble=True, seed=seed)
    u = s.random(n)
    v = ndtri(np.clip(u, 1e-12, 1 - 1e-12))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _kernel_args(model: GravityModel):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


# ------------------------------------------------------------------ R1-B
def degree_power(model: GravityModel) -> np.ndarray:
    """Sum_m (Cnm^2 + Snm^2) per degree from the coefficient file."""
    N = model.max_degree
    C, S = model.c_coeffs, model.s_coeffs
    return np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                     for n in range(N + 1)])


def sigma_a(model: GravityModel, r: float, kind: str, power: np.ndarray) -> np.ndarray:
    """Analytic per-degree acceleration RMS over the sphere of radius r.

    kind='radial': (mu/r^2)(n+1)(R/r)^n sqrt(P_n)
    kind='vector': (mu/r^2) sqrt((n+1)(2n+1)) (R/r)^n sqrt(P_n)
    """
    n = np.arange(len(power), dtype=np.float64)
    ratio_n = np.exp(n * math.log(model.r_ref / r))
    base = (model.mu / r ** 2) * ratio_n * np.sqrt(power)
    if kind == "radial":
        return (n + 1.0) * base
    if kind == "vector":
        return np.sqrt((n + 1.0) * (2.0 * n + 1.0)) * base
    raise ValueError(kind)


def r1_verification(model: GravityModel, power: np.ndarray) -> None:
    print("== R1-B: isolated-degree sampling verification ==")
    n_dirs = 1500
    rows = []
    for h_km in (50.0, 100.0):
        r = model.r_ref + h_km * 1e3
        xyz = r * _dirs(n_dirs, SEED)
        u = xyz / r
        for n in (10, 60, 120, 200, 300):
            _, a = sh_potential_accel_fixed(
                xyz, model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
                degree_max=n, degree_min=n - 1)
            rms_vec = float(np.sqrt(np.mean(np.sum(a ** 2, axis=1))))
            a_rad = np.sum(a * u, axis=1)
            rms_rad = float(np.sqrt(np.mean(a_rad ** 2)))
            an_vec = float(sigma_a(model, r, "vector", power)[n])
            an_rad = float(sigma_a(model, r, "radial", power)[n])
            rows.append({
                "altitude_km": h_km, "degree": n, "n_dirs": n_dirs,
                "sampled_vector_rms": rms_vec, "analytic_vector_rms": an_vec,
                "vector_ratio": rms_vec / an_vec,
                "sampled_radial_rms": rms_rad, "analytic_radial_rms": an_rad,
                "radial_ratio": rms_rad / an_rad,
            })
            print(f"  h={h_km:5.0f} n={n:3d}  vec ratio {rms_vec/an_vec:.4f}  "
                  f"rad ratio {rms_rad/an_rad:.4f}")
    _dump("r1_degree_rms_verification.json", {
        "seed": SEED, "n_dirs": n_dirs,
        "note": "isolated degree-n field via batch kernel degree_min=n-1; "
                "sampled sphere RMS vs analytic radial/vector formulas",
        "rows": rows,
    })


# ------------------------------------------------------------------ R1-C
def _nmin_from_sigma(sig: np.ndarray, eps: float) -> int:
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


def _nmin_proxy(R: float, r: float, p: float, eps: float, kind: str,
                n_max: int = 100_000) -> int:
    """Kaula-proxy N_min with radial or vector degree factor."""
    ratio = R / r
    terms = []
    peak = 0.0
    log_r = math.log(ratio)
    n = 2
    while n <= n_max:
        if kind == "radial":
            fac = (n + 1.0) * math.sqrt(2.0 * n + 1.0)
        else:
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


def r1_pfit_and_criteria(model: GravityModel, power: np.ndarray) -> None:
    print("== R1-C: formal p calibration ==")
    n_arr = np.arange(len(power))
    S_n = np.sqrt(power)
    sigma_c = S_n / np.sqrt(2.0 * n_arr.clip(1) + 1.0)  # per-coefficient RMS

    def ols(lo: int, hi: int):
        m = (n_arr >= lo) & (n_arr <= hi) & (sigma_c > 0)
        x = np.log10(n_arr[m].astype(float))
        y = np.log10(sigma_c[m])
        A = np.vstack([np.ones_like(x), -x]).T
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        logK, p = coef
        yhat = A @ coef
        resid = y - yhat
        dof = len(x) - 2
        s2 = float(np.sum(resid ** 2) / dof)
        cov = s2 * np.linalg.inv(A.T @ A)
        return {"range": [lo, hi], "p": float(p), "logK": float(logK),
                "p_se_ols": float(np.sqrt(cov[1, 1])),
                "resid_rms_dex": float(np.sqrt(np.mean(resid ** 2))),
                "x": x, "y": y, "resid": resid}

    def block_bootstrap_p(lo: int, hi: int, n_boot: int = 2000,
                          block: int = 25, seed: int = SEED):
        m = (n_arr >= lo) & (n_arr <= hi) & (sigma_c > 0)
        x = np.log10(n_arr[m].astype(float))
        y = np.log10(sigma_c[m])
        nblk = len(x) // block
        rng = np.random.default_rng(seed)
        ps = []
        for _ in range(n_boot):
            idx = []
            for b in rng.integers(0, nblk, size=nblk):
                idx.extend(range(b * block, min((b + 1) * block, len(x))))
            idx = np.array(idx)
            A = np.vstack([np.ones(len(idx)), -x[idx]]).T
            coef, *_ = np.linalg.lstsq(A, y[idx], rcond=None)
            ps.append(float(coef[1]))
        return [float(np.percentile(ps, 2.5)), float(np.percentile(ps, 97.5))]

    fit_600 = ols(10, 600)
    fit_full = ols(2, 1800)
    fit_hold = ols(10, 300)
    ci_600 = block_bootstrap_p(10, 600)
    # holdout residuals on 301..600 using the 10..300 fit
    m = (n_arr >= 301) & (n_arr <= 600)
    x_h = np.log10(n_arr[m].astype(float))
    y_h = np.log10(sigma_c[m])
    pred = fit_hold["logK"] - fit_hold["p"] * x_h
    holdout_rms = float(np.sqrt(np.mean((y_h - pred) ** 2)))
    print(f"  p(10-600) = {fit_600['p']:.3f} +- {fit_600['p_se_ols']:.3f} "
          f"(bootstrap CI {ci_600[0]:.3f}..{ci_600[1]:.3f})")
    print(f"  p(2-1800) = {fit_full['p']:.3f};  p(10-300) = {fit_hold['p']:.3f}, "
          f"holdout(301-600) resid {holdout_rms:.3f} dex vs "
          f"in-fit {fit_hold['resid_rms_dex']:.3f} dex")

    # empirical N_min, both weightings, and proxies
    alts = [20.0, 30.0, 50.0, 80.0, 100.0, 150.0, 200.0, 300.0]
    crit_rows = []
    for h in alts:
        r = model.r_ref + h * 1e3
        sig_vec = sigma_a(model, r, "vector", power)
        sig_rad = sigma_a(model, r, "radial", power)
        row = {
            "altitude_km": h,
            "emp_vector_eps1e2": _nmin_from_sigma(sig_vec, 1e-2),
            "emp_radial_eps1e2": _nmin_from_sigma(sig_rad, 1e-2),
            "emp_vector_eps1e3": _nmin_from_sigma(sig_vec, 1e-3),
            "proxy_vec_p1_7": _nmin_proxy(model.r_ref, r, 1.7, 1e-2, "vector"),
            "proxy_vec_p2_0": _nmin_proxy(model.r_ref, r, 2.0, 1e-2, "vector"),
            "proxy_vec_pfit": _nmin_proxy(model.r_ref, r, fit_600["p"], 1e-2, "vector"),
            "proxy_rad_p1_7": _nmin_proxy(model.r_ref, r, 1.7, 1e-2, "radial"),
            "atten_1e3": recommended_sh_degree(h, model.r_ref, attenuation_floor=1e-3),
        }
        crit_rows.append(row)
        print(f"  h={h:5.0f}  emp_vec={row['emp_vector_eps1e2']:4d} "
              f"emp_rad={row['emp_radial_eps1e2']:4d}  "
              f"p1.7={row['proxy_vec_p1_7']:4d}  p2.0={row['proxy_vec_p2_0']:4d} "
              f"pfit={row['proxy_vec_pfit']:4d}  atten={row['atten_1e3']:4d}")

    # usage-based p*: argmin over grid of sum (N_spec - N_emp)^2, h in 50..300
    grid = np.round(np.arange(1.40, 2.101, 0.02), 3)
    best = None
    scores = []
    for p in grid:
        s = 0
        for row in crit_rows:
            if row["altitude_km"] < 50.0:
                continue
            npred = _nmin_proxy(model.r_ref, model.r_ref + row["altitude_km"] * 1e3,
                                float(p), 1e-2, "vector")
            s += (npred - row["emp_vector_eps1e2"]) ** 2
        scores.append([float(p), int(s)])
        if best is None or s < best[1]:
            best = (float(p), int(s))
    print(f"  usage-based p* (vector, 50-300 km) = {best[0]:.2f} (SSE {best[1]})")

    # sensitivity: N_min vs p and eps at selected altitudes
    sens = []
    for h in (50.0, 100.0, 200.0):
        r = model.r_ref + h * 1e3
        for p in (1.5, 1.6, 1.7, 1.8, 1.9, 2.0):
            for eps in (3e-2, 1e-2, 3e-3, 1e-3):
                sens.append({"altitude_km": h, "p": p, "eps": eps,
                             "nmin": _nmin_proxy(model.r_ref, r, p, eps, "vector")})

    _dump("r1_spectrum_pfit.json", {
        "seed": SEED,
        "sigma_model": "per-coefficient RMS sigma_n = sqrt(P_n/(2n+1)); "
                       "P_n = sum_m Cnm^2+Snm^2 from the coefficient file",
        "fit_10_600": {k: v for k, v in fit_600.items() if k not in ("x", "y", "resid")},
        "fit_10_600_bootstrap_ci95": ci_600,
        "fit_2_1800": {k: v for k, v in fit_full.items() if k not in ("x", "y", "resid")},
        "fit_10_300": {k: v for k, v in fit_hold.items() if k not in ("x", "y", "resid")},
        "holdout_301_600_resid_rms_dex": holdout_rms,
        "usage_based_pstar": {"p": best[0], "sse": best[1],
                              "grid_scores": scores,
                              "band_km": [50.0, 300.0], "eps": 1e-2},
        "criteria_rows": crit_rows,
        "sensitivity": sens,
        "spectrum_arrays": {
            "n": n_arr[2:].tolist(),
            "sigma_coeff_rms": sigma_c[2:].tolist(),
        },
    })


# ------------------------------------------------------------------ R1-E1
def r1_band_shares(model300: GravityModel) -> None:
    print("== R1-E1: band shares with bootstrap CIs ==")
    m = model300
    ws = m.make_workspace()

    def shares_at(h_km: float, dirs: np.ndarray):
        r = m.r_ref + h_km * 1e3
        per_dir = {"b2_60": [], "b61_100": [], "tail": [], "pert": []}
        for u in dirs:
            pos = r * u
            a60 = m.accel_fixed(pos, degree=60, workspace=ws)
            a100 = m.accel_fixed(pos, degree=100, workspace=ws)
            a300 = m.accel_fixed(pos, degree=300, workspace=ws)
            a_pm = -(m.mu / (r * r)) * u
            per_dir["pert"].append(float(np.sum((a300 - a_pm) ** 2)))
            per_dir["b2_60"].append(float(np.sum((a60 - a_pm) ** 2)))
            per_dir["b61_100"].append(float(np.sum((a100 - a60) ** 2)))
            per_dir["tail"].append(float(np.sum((a300 - a100) ** 2)))
        return {k: np.array(v) for k, v in per_dir.items()}

    def summarize(sq: dict, n_boot: int = 2000, seed: int = SEED):
        n = len(sq["pert"])
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, n, size=(n_boot, n))
        out = {}
        for key in ("b2_60", "b61_100", "tail"):
            share = math.sqrt(sq[key].mean()) / math.sqrt(sq["pert"].mean())
            boots = np.sqrt(sq[key][idx].mean(axis=1)) / np.sqrt(sq["pert"][idx].mean(axis=1))
            out[key] = {"share": share,
                        "ci95": [float(np.percentile(boots, 2.5)),
                                 float(np.percentile(boots, 97.5))]}
        out["pert_rms"] = math.sqrt(sq["pert"].mean())
        return out

    rows = []
    for h in (30.0, 50.0, 80.0, 100.0, 150.0, 200.0):
        sq = shares_at(h, _dirs(1000, SEED))
        s = summarize(sq)
        s["altitude_km"] = h
        rows.append(s)
        print(f"  h={h:5.0f}  61-100 {100*s['b61_100']['share']:.2f}% "
              f"[{100*s['b61_100']['ci95'][0]:.2f},{100*s['b61_100']['ci95'][1]:.2f}]")

    # controls at 80 km: alternative seed and Sobol sampling
    ctrl = {}
    for name, dirs in (("seed_alt", _dirs(1000, SEED_ALT)),
                       ("sobol", _sobol_dirs(1000, SEED))):
        s = summarize(shares_at(80.0, dirs), seed=SEED_ALT)
        ctrl[name] = {k: s[k] for k in ("b2_60", "b61_100", "tail")}
        print(f"  80 km control {name}: 61-100 {100*s['b61_100']['share']:.2f}%")

    _dump("r1_band_shares.json", {
        "seed": SEED, "seed_alt": SEED_ALT, "n_dirs": 1000,
        "quantity": "normalized band-difference RMS (non-additive)",
        "rows": rows, "controls_80km": ctrl,
    })


# ------------------------------------------------------------------ R1-E5
def r1_switch_jump(model300: GravityModel) -> None:
    print("== R1-E5: switching jump with CIs ==")
    m = model300
    ws = m.make_workspace()
    dirs = _dirs(400, SEED)
    rng = np.random.default_rng(SEED)
    rows = []
    for h_km in (50.0, 100.0):
        r = m.r_ref + h_km * 1e3
        pert_sq = []
        jump_sq = {q: [] for q in (5, 10, 25, 50)}
        for u in dirs:
            pos = r * u
            a120 = m.accel_fixed(pos, degree=120, workspace=ws)
            a_pm = -(m.mu / (r * r)) * u
            pert_sq.append(float(np.sum((a120 - a_pm) ** 2)))
            for q in jump_sq:
                a_lo = m.accel_fixed(pos, degree=120 - q, workspace=ws)
                jump_sq[q].append(float(np.sum((a120 - a_lo) ** 2)))
        pert_sq = np.array(pert_sq)
        idx = rng.integers(0, len(dirs), size=(2000, len(dirs)))
        for q, sq in jump_sq.items():
            sq = np.array(sq)
            ratio = math.sqrt(sq.mean()) / math.sqrt(pert_sq.mean())
            boots = np.sqrt(sq[idx].mean(axis=1)) / np.sqrt(pert_sq[idx].mean(axis=1))
            rows.append({"altitude_km": h_km, "step": q,
                         "rms_jump_m_s2": math.sqrt(sq.mean()),
                         "jump_over_pert": ratio,
                         "ci95": [float(np.percentile(boots, 2.5)),
                                  float(np.percentile(boots, 97.5))]})
            print(f"  h={h_km:5.0f} q={q:2d}: {100*ratio:.3f}% "
                  f"[{100*np.percentile(boots,2.5):.3f},{100*np.percentile(boots,97.5):.3f}]")
    _dump("r1_switch_jump.json", {"seed": SEED, "n_dirs": 400, "rows": rows})


# ------------------------------------------------------------------ R1-E3
def r1_timing(model1800: GravityModel) -> None:
    print("== R1-E3: kernel timing (median/IQR) ==")
    args = _kernel_args(model1800)
    r = model1800.r_ref + 80e3
    lat, lon = math.radians(25.0), math.radians(40.0)
    x = r * math.cos(lat) * math.cos(lon)
    y = r * math.cos(lat) * math.sin(lon)
    z = r * math.sin(lat)
    sh_accel_fixed_numba(x, y, z, 10, *args)  # warm-up

    rows = []
    for n in (10, 20, 40, 60, 80, 120, 160, 200, 300, 400, 600, 800, 1200, 1800):
        reps = max(5, min(400, int(4_000_000 / (n * n + 1))))
        blocks = []
        for _ in range(9):
            t0 = time.perf_counter_ns()
            for _ in range(reps):
                sh_accel_fixed_numba(x, y, z, n, *args)
            blocks.append((time.perf_counter_ns() - t0) / reps / 1000.0)
        blocks.sort()
        rows.append({"degree": n, "reps": reps,
                     "median_us": statistics.median(blocks),
                     "q1_us": blocks[2], "q3_us": blocks[6],
                     "best_us": blocks[0]})
        print(f"  N={n:5d}  median {rows[-1]['median_us']:10.2f} us "
              f"IQR [{blocks[2]:.2f},{blocks[6]:.2f}]")

    _compute_sh_acceleration_dual_numba(x, y, z, 60, 120, *args)
    dual_blocks, two_blocks = [], []
    for _ in range(9):
        t0 = time.perf_counter_ns()
        for _ in range(1000):
            _compute_sh_acceleration_dual_numba(x, y, z, 60, 120, *args)
        dual_blocks.append((time.perf_counter_ns() - t0) / 1000 / 1000.0)
        t0 = time.perf_counter_ns()
        for _ in range(1000):
            sh_accel_fixed_numba(x, y, z, 60, *args)
            sh_accel_fixed_numba(x, y, z, 120, *args)
        two_blocks.append((time.perf_counter_ns() - t0) / 1000 / 1000.0)
    _dump("r1_kernel_timing.json", {
        "timer": "perf_counter_ns; 9 blocks per degree; median and IQR reported",
        "degree_sweep": rows,
        "dual_pass": {"pair": [60, 120],
                      "dual_median_us": statistics.median(dual_blocks),
                      "two_calls_median_us": statistics.median(two_blocks)},
    })


# ------------------------------------------------------------------ R1-E4
def r1_blend(model300: GravityModel) -> None:
    print("== R1-E4: blend analysis rerun ==")
    m = model300
    args = _kernel_args(m)
    n_far, n_near = 30, 120
    alt_far, alt_near = 200e3, 50e3
    step = 10
    lat, lon = math.radians(25.0), math.radians(40.0)
    u = np.array([math.cos(lat) * math.cos(lon),
                  math.cos(lat) * math.sin(lon), math.sin(lat)])

    def blend_at(pos):
        ax, ay, az = sh_accel_adaptive_blend_numba(
            pos[0], pos[1], pos[2], n_far, n_near, alt_far, alt_near, step, *args)
        return np.array([ax, ay, az])

    def fixed_at(pos, n):
        ax, ay, az = sh_accel_fixed_numba(pos[0], pos[1], pos[2], n, *args)
        return np.array([ax, ay, az])

    def curl_of(field, pos, h_fd=0.5):
        J = np.zeros((3, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = h_fd
            J[:, j] = (field(pos + e) - field(pos - e)) / (2.0 * h_fd)
        c = np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]])
        return float(np.linalg.norm(c))

    pos_mid = (m.r_ref + 118e3) * u
    curl_blend = curl_of(blend_at, pos_mid)
    curl_fixed = curl_of(lambda p: fixed_at(p, 120), pos_mid)

    def wchain(rr):
        tt = (alt_far - (rr - m.r_ref)) / (alt_far - alt_near)
        ss = _apply_smoothstep(tt)
        dd = n_far + ss * (n_near - n_far)
        kk = int((dd - n_far) // step)
        lo = n_far + kk * step
        hi = min(lo + step, n_near)
        if hi == lo:
            return lo, hi, 0.0
        return lo, hi, min(max((dd - lo) / (hi - lo), 0.0), 1.0)

    omitted_max, at_km = 0.0, float("nan")
    for h in np.linspace(alt_near + 500.0, alt_far - 500.0, 300):
        r = m.r_ref + h
        lo, hi, w = wchain(r)
        lo2, hi2, w2 = wchain(r + 25.0)
        if (lo2, hi2) != (lo, hi) or hi == lo:
            continue
        dw = abs(w2 - w) / 25.0
        pos = r * u
        U_lo = m.potential_fixed(pos, degree=lo)
        U_hi = m.potential_fixed(pos, degree=hi)
        om = abs(U_hi - U_lo) * dw
        if om > omitted_max:
            omitted_max, at_km = om, h / 1e3
    r_at = m.r_ref + at_km * 1e3
    pos_at = r_at * u
    a_pert_at = float(np.linalg.norm(
        fixed_at(pos_at, 120) + (m.mu / r_at ** 2) * u))
    print(f"  curl blend {curl_blend:.3e} vs fixed {curl_fixed:.3e}; "
          f"omitted max {omitted_max:.3e} at {at_km:.1f} km "
          f"(local pert {a_pert_at:.3e})")
    _dump("r1_blend.json", {
        "config": {"degree_far": n_far, "degree_near": n_near,
                   "alt_far_m": alt_far, "alt_near_m": alt_near, "step": step},
        "curl_1_s2": {"blend": curl_blend, "fixed": curl_fixed, "fd_step_m": 0.5,
                      "at_altitude_km": 118.0},
        "omitted_term": {"max_m_s2": omitted_max, "at_altitude_km": at_km,
                         "local_pert_m_s2_at_same_altitude": a_pert_at},
    })


def main() -> int:
    global COMMIT
    COMMIT = commit_sha()
    print("commit:", COMMIT)
    OUT.mkdir(exist_ok=True)
    gf = resolve_lunar_gravity_path(None)
    model1800 = GravityModel.from_file(str(gf), requested_degree=1800)
    model300 = GravityModel.from_file(str(gf), requested_degree=300)
    power = degree_power(model1800)

    r1_verification(model1800, power)
    r1_pfit_and_criteria(model1800, power)
    r1_band_shares(model300)
    r1_switch_jump(model300)
    r1_timing(model1800)
    r1_blend(model300)
    print("field-level revision experiments complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
