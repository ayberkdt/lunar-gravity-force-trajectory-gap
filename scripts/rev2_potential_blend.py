"""Round-2 item 2: corrected potential-level blend vs. alternatives.

Four gravity policies are compared on the same transition band
(N_lo=30 above 200 km, N_hi=120 below 50 km, smoothstep weight w(r)):

  fixed        a = a(N_hi)                         (reference field)
  switch       a = a(N_lo) or a(N_hi), discrete    (discontinuous)
  blend_accel  a = (1-w) a_lo + w a_hi             (C0, non-conservative)
  blend_pot    a = (1-w) a_lo + w a_hi
                   + (U_hi - U_lo) (dw/dr) r_hat    (a = grad U_w, conservative)

Measurements
------------
1. Curl at a mid-band point (central differences) for each policy.
2. Specific-energy drift over a NON-rotating body-fixed eccentric orbit
   whose perilune sits below the band and apolune above it, so it crosses
   the transition twice per revolution. In a static conservative field the
   specific energy E = |v|^2/2 - U is an exact integral; a non-conservative
   field breaks it. This isolates the conservation property from any
   rotating-frame effect.
3. Position error of each policy against the fixed N_hi=120 field over the
   same arc (a truncation/scheduling error, not a conservation metric).
4. Per-call kernel wall time for each policy.

The corrected blend uses the production scalar-potential kernel
(sh_potential_accel_fixed, which returns U and a=grad U together).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.common.lunar_data import resolve_lunar_gravity_path  # noqa: E402
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel,
    sh_potential_accel_fixed,
)

OUT = Path(__file__).resolve().parents[1] / "metrics"
N_LO, N_HI = 30, 120
ALT_NEAR, ALT_FAR = 50e3, 200e3  # w=1 below near, w=0 above far


def commit_sha() -> str:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                       text=True, check=True)
    return r.stdout.strip()


def weight_and_deriv(r: float, R: float):
    """Smoothstep weight w(alt) and dw/dr; w=1 (high degree) at low altitude."""
    alt = r - R
    if alt <= ALT_NEAR:
        return 1.0, 0.0
    if alt >= ALT_FAR:
        return 0.0, 0.0
    t = (ALT_FAR - alt) / (ALT_FAR - ALT_NEAR)          # 0 at far, 1 at near
    s = t * t * (3.0 - 2.0 * t)                          # smoothstep
    dsdt = 6.0 * t * (1.0 - t)
    dtdr = -1.0 / (ALT_FAR - ALT_NEAR)
    return s, dsdt * dtdr


def make_evaluators(model: GravityModel):
    R = model.r_ref
    C, S, mu = model.c_coeffs, model.s_coeffs, model.mu

    def U_and_a(pos, degree):
        V, a = sh_potential_accel_fixed(
            np.asarray(pos, float).reshape(1, 3), C, S, mu, R, degree, -1)
        return float(V[0]), a[0]

    def fixed(pos):
        _, a = U_and_a(pos, N_HI)
        return a

    def switch(pos):
        alt = float(np.linalg.norm(pos)) - R
        deg = N_HI if alt <= 0.5 * (ALT_NEAR + ALT_FAR) else N_LO
        _, a = U_and_a(pos, deg)
        return a

    def blend_accel(pos):
        r = float(np.linalg.norm(pos))
        w, _ = weight_and_deriv(r, R)
        _, a_lo = U_and_a(pos, N_LO)
        _, a_hi = U_and_a(pos, N_HI)
        return (1.0 - w) * a_lo + w * a_hi

    def blend_pot(pos):
        r = float(np.linalg.norm(pos))
        w, dwdr = weight_and_deriv(r, R)
        U_lo, a_lo = U_and_a(pos, N_LO)
        U_hi, a_hi = U_and_a(pos, N_HI)
        a = (1.0 - w) * a_lo + w * a_hi
        if dwdr != 0.0:
            a = a + (U_hi - U_lo) * dwdr * (np.asarray(pos) / r)
        return a

    def potential_of_policy(pos, policy):
        """U_w for the blended-potential policy (for the energy integral)."""
        r = float(np.linalg.norm(pos))
        w, _ = weight_and_deriv(r, R)
        U_lo, _ = U_and_a(pos, N_LO)
        U_hi, _ = U_and_a(pos, N_HI)
        if policy == "blend_pot":
            return (1.0 - w) * U_lo + w * U_hi
        if policy == "fixed":
            return U_hi
        if policy == "switch":
            alt = r - R
            return U_hi if alt <= 0.5 * (ALT_NEAR + ALT_FAR) else U_lo
        return None  # blend_accel has no potential

    return {"fixed": fixed, "switch": switch, "blend_accel": blend_accel,
            "blend_pot": blend_pot}, potential_of_policy


def curl_at(field, pos, h_fd=0.5):
    J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h_fd
        J[:, j] = (field(pos + e) - field(pos - e)) / (2.0 * h_fd)
    c = np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]])
    return float(np.linalg.norm(c))


def propagate_bodyfixed(field, y0, dur, t_eval):
    """Non-rotating body-fixed propagation: a = field(pos) directly."""
    def rhs(t, y):
        a = field(y[:3])
        return (y[3], y[4], y[5], a[0], a[1], a[2])
    sol = solve_ivp(rhs, (0.0, dur), y0, method="DOP853", rtol=1e-11, atol=1e-5,
                    t_eval=t_eval, dense_output=False)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def main() -> int:
    commit = commit_sha()
    print("commit:", commit)
    OUT.mkdir(exist_ok=True)
    model = GravityModel.from_file(str(resolve_lunar_gravity_path(None)),
                                   requested_degree=120)
    R, mu = model.r_ref, model.mu
    fields, pot_of = make_evaluators(model)

    # geometry in a fixed direction for the curl probe
    lat, lon = math.radians(25.0), math.radians(40.0)
    u = np.array([math.cos(lat) * math.cos(lon),
                  math.cos(lat) * math.sin(lon), math.sin(lat)])
    pos_mid = (R + 118e3) * u

    print("== curl at 118 km mid-band ==")
    curl = {name: curl_at(f, pos_mid) for name, f in fields.items()}
    for name, c in curl.items():
        print(f"  {name:12s} |curl| = {c:.3e} 1/s^2")

    # eccentric body-fixed orbit crossing the band: perilune 30 km, apolune 260 km
    rp, ra = R + 30e3, R + 260e3
    a_sma = 0.5 * (rp + ra)
    vp = math.sqrt(mu * (2.0 / rp - 1.0 / a_sma))
    y0 = np.array([rp, 0.0, 0.0, 0.0, vp * math.cos(math.radians(20.0)),
                   vp * math.sin(math.radians(20.0))])  # inclined in body frame
    period = 2.0 * math.pi * math.sqrt(a_sma ** 3 / mu)
    dur = 6.0 * period
    t_eval = np.linspace(0.0, dur, 4001)
    print(f"  orbit period {period:.1f} s, {dur/period:.0f} revolutions")

    truth = propagate_bodyfixed(fields["fixed"], y0, dur, t_eval)

    results = {}
    for name, field in fields.items():
        t0 = time.perf_counter()
        sol = propagate_bodyfixed(field, y0, dur, t_eval)
        wall = time.perf_counter() - t0
        # per-call kernel time
        reps = 300
        tt = time.perf_counter_ns()
        for _ in range(reps):
            field(pos_mid)
        per_call_us = (time.perf_counter_ns() - tt) / reps / 1000.0
        # specific energy using the policy's own potential (where defined)
        E_series = None
        if name != "blend_accel":
            E = []
            for k in range(sol.y.shape[1]):
                pos = sol.y[:3, k]
                vel = sol.y[3:, k]
                U = pot_of(pos, name)
                E.append(0.5 * float(vel @ vel) - U)
            E = np.array(E)
            E_series = E
        # For blend_accel there is no scalar potential of its own; measure the
        # energy against the SAME blended potential U_w the corrected policy
        # integrates. The drift of |v|^2/2 - U_w then equals exactly the work
        # done by the omitted (U_hi - U_lo) grad(w) term, i.e. the physical
        # non-conservation of the acceleration blend.
        else:
            E = []
            for k in range(sol.y.shape[1]):
                pos = sol.y[:3, k]
                vel = sol.y[3:, k]
                U = pot_of(pos, "blend_pot")
                E.append(0.5 * float(vel @ vel) - U)
            E = np.array(E)
            E_series = E
        E0 = E_series[0]
        e_drift_rel = float(np.max(np.abs(E_series - E0)) / abs(E0))
        e_final_rel = float(abs(E_series[-1] - E0) / abs(E0))
        # secular fit slope of energy vs time (rel per revolution)
        rev = t_eval / period
        slope = float(np.polyfit(rev, (E_series - E0) / abs(E0), 1)[0])
        # position error vs fixed truth
        dp = np.linalg.norm(sol.y[:3] - truth.y[:3], axis=0)
        results[name] = {
            "curl_1_s2": curl[name],
            "per_call_us": per_call_us,
            "wall_s_12rev": wall,
            "energy_peak_rel_drift": e_drift_rel,
            "energy_final_rel_drift": e_final_rel,
            "energy_secular_rel_per_rev": slope,
            "pos_rms_vs_fixed_m": float(np.sqrt(np.mean(dp ** 2))),
            "pos_max_vs_fixed_m": float(np.max(dp)),
            "energy_rel_series": ((E_series - E0) / abs(E0))[::50].tolist(),
        }
        print(f"  {name:12s} E peak drift {e_drift_rel:.2e} rel, "
              f"secular {slope:+.2e}/rev, pos-vs-fixed "
              f"{results[name]['pos_rms_vs_fixed_m']:.2f} m RMS, "
              f"{per_call_us:.1f} us/call")

    payload = {
        "repo_commit_sha": commit,
        "config": {"N_lo": N_LO, "N_hi": N_HI, "alt_near_m": ALT_NEAR,
                   "alt_far_m": ALT_FAR, "weight": "smoothstep single band",
                   "test_frame": "non-rotating body-fixed (conservative test)",
                   "orbit": {"perilune_km": 30.0, "apolune_km": 260.0,
                             "revolutions": round(dur / period),
                             "period_s": period},
                   "energy_note": "E = |v|^2/2 - U; U is the policy's own "
                                  "blended potential where defined; for "
                                  "blend_accel (no potential) the fixed-field "
                                  "U is used as a diagnostic"},
        "t_rev_series": (t_eval / period)[::50].tolist(),
        "results": results,
    }
    (OUT / "r2_potential_blend.json").write_text(json.dumps(payload, indent=2),
                                                 encoding="utf-8")
    print("[written] r2_potential_blend.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
