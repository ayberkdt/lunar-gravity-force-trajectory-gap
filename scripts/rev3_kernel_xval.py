"""P0-8: high-degree external validation of the production SH kernel.

Two independent references:
  1. SHTOOLS (pyshtools MakeGravGridPoint, Fortran backend) for the
     gravity vector at N = 120, 300, 600, 1200, 1800;
  2. the repository-independent classical-recurrence implementation
     (validation/independent/independent_sh.py, different ALF source,
     acceleration via high-order numerical gradient) for the potential and
     acceleration at moderate degrees.

Points: seeded random off-pole and near-pole directions at 30, 100 and
500 km altitude. Reported: relative acceleration error (max/RMS), radial
and tangential component errors, relative potential error, and the
(r, theta, phi) -> Cartesian convention resolution.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rev3_common import REPO, SEED, dump, load_model

import sys

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO))
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    sh_accel_fixed_numba, sh_potential_accel_fixed)
from validation.independent import independent_sh  # noqa: E402

import pyshtools  # noqa: E402
from pyshtools.backends.shtools import MakeGravGridPoint  # noqa: E402

DEGREES = (120, 300, 600, 1200, 1800)
INDEP_DEGREES = (120, 300, 600)
ALTS_KM = (30.0, 100.0, 500.0)


def kernel_args(model):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


def make_points(rng, n_offpole=128, n_nearpole=16):
    pts = []
    for _ in range(n_offpole):
        while True:
            v = rng.normal(size=3)
            v /= np.linalg.norm(v)
            if abs(math.degrees(math.asin(v[2]))) < 80.0:
                pts.append(v)
                break
    for k in range(n_nearpole):
        lat = math.radians(89.6 + 0.09 * k) * (1 if k % 2 == 0 else -1)
        lon = rng.uniform(0.0, 2.0 * math.pi)
        pts.append(np.array([math.cos(lat) * math.cos(lon),
                             math.cos(lat) * math.sin(lon),
                             math.sin(lat)]))
    return np.array(pts)


def sph_to_cart(gr, gth, gphi, lat_rad, lon_rad, theta_south: bool):
    cl, sl = math.cos(lon_rad), math.sin(lon_rad)
    cf, sf = math.cos(lat_rad), math.sin(lat_rad)
    rhat = np.array([cf * cl, cf * sl, sf])
    that = np.array([sf * cl, sf * sl, -cf])  # increasing colatitude (south)
    if not theta_south:
        that = -that
    phat = np.array([-sl, cl, 0.0])
    return gr * rhat + gth * that + gphi * phat


def main() -> int:
    model = load_model(1800)
    args = kernel_args(model)
    sh_accel_fixed_numba(model.r_ref + 100e3, 0.0, 0.0, 60, *args)

    L = model.max_degree
    cilm = np.zeros((2, L + 1, L + 1))
    cilm[0] = model.c_coeffs
    cilm[1] = model.s_coeffs
    cilm[0, 0, 0] = 1.0  # structural monopole (GRAIL SHADR omits degree 0)

    rng = np.random.default_rng(SEED)
    units_by_alt = {h: make_points(rng) for h in ALTS_KM}
    units = units_by_alt[ALTS_KM[0]]

    # resolve the (r,theta,phi) convention once, at full degree, first point
    conv_votes = []
    for h_km in ALTS_KM[:1]:
        r = model.r_ref + h_km * 1e3
        for u in units[:4]:
            xyz = r * u
            lat = math.asin(u[2])
            lon = math.atan2(u[1], u[0])
            a_lun = np.array(sh_accel_fixed_numba(*xyz, 1800, *args))
            g = MakeGravGridPoint(cilm, model.mu, model.r_ref, r,
                                  math.degrees(lat), math.degrees(lon),
                                  lmax=1800)
            for ts in (True, False):
                a_ref = sph_to_cart(g[0], g[1], g[2], lat, lon, ts)
                conv_votes.append((ts, float(np.linalg.norm(a_ref - a_lun) /
                                             np.linalg.norm(a_lun))))
    best = {}
    for ts, e in conv_votes:
        best.setdefault(ts, []).append(e)
    theta_south = (np.median(best[True]) < np.median(best[False]))
    print(f"theta convention: {'south (colatitude)' if theta_south else 'north (latitude)'}")

    rows = []
    for h_km in ALTS_KM:
        r = model.r_ref + h_km * 1e3
        units = units_by_alt[h_km]
        for N in DEGREES:
            rel, rel_rad, rel_tan = [], [], []
            for i, u in enumerate(units):
                xyz = r * u
                lat = math.asin(u[2])
                lon = math.atan2(u[1], u[0])
                a_lun = np.array(sh_accel_fixed_numba(*xyz, N, *args))
                g = MakeGravGridPoint(cilm, model.mu, model.r_ref, r,
                                      math.degrees(lat), math.degrees(lon),
                                      lmax=N)
                a_ref = sph_to_cart(g[0], g[1], g[2], lat, lon, theta_south)
                scale = np.linalg.norm(a_ref)
                d = a_lun - a_ref
                rel.append(float(np.linalg.norm(d) / scale))
                dr = float(np.dot(d, u))
                rel_rad.append(abs(dr) / scale)
                rel_tan.append(float(np.linalg.norm(d - dr * u)) / scale)
            rows.append({
                "altitude_km": h_km, "degree": N, "n_points": len(units),
                "n_near_pole": 16,
                "accel_rel_max": float(np.max(rel)),
                "accel_rel_rms": float(np.sqrt(np.mean(np.square(rel)))),
                "accel_rel_radial_max": float(np.max(rel_rad)),
                "accel_rel_tangential_max": float(np.max(rel_tan)),
            })
            print(f"SHTOOLS h={h_km:5.0f} N={N:4d}: rel max "
                  f"{rows[-1]['accel_rel_max']:.3e} rms "
                  f"{rows[-1]['accel_rel_rms']:.3e}")

    indep_rows = []
    units = units_by_alt[ALTS_KM[0]]
    for h_km in (30.0, 100.0):
        r = model.r_ref + h_km * 1e3
        for N in INDEP_DEGREES:
            relU, relA = [], []
            for u in units[:6]:
                xyz = (r * u).reshape(1, 3)
                V_lun, a_lun = sh_potential_accel_fixed(
                    xyz, model.c_coeffs, model.s_coeffs, model.mu,
                    model.r_ref, degree_max=N)
                U_ind = independent_sh.geopotential(
                    xyz[0], mu=model.mu, r_ref=model.r_ref,
                    c_coeffs=model.c_coeffs, s_coeffs=model.s_coeffs,
                    degree=N)
                a_ind = independent_sh.acceleration(
                    xyz[0], mu=model.mu, r_ref=model.r_ref,
                    c_coeffs=model.c_coeffs, s_coeffs=model.s_coeffs,
                    degree=N)
                relU.append(abs(float(V_lun[0]) - float(U_ind)) /
                            abs(float(U_ind)))
                relA.append(float(np.linalg.norm(a_lun[0] - a_ind) /
                                  np.linalg.norm(a_ind)))
            indep_rows.append({
                "altitude_km": h_km, "degree": N,
                "potential_rel_max": float(np.max(relU)),
                "accel_rel_max": float(np.max(relA)),
                "note": "acceleration reference is a numerical gradient; "
                        "its stencil error dominates at tight levels",
            })
            print(f"indep h={h_km:5.0f} N={N:4d}: U rel {np.max(relU):.3e} "
                  f"a rel {np.max(relA):.3e}")

    dump("r3_kernel_xval.json", {
        "seed": SEED,
        "references": {
            "shtools": f"pyshtools {pyshtools.__version__} MakeGravGridPoint",
            "independent": "validation/independent/independent_sh.py "
                           "(classical recurrence + numerical gradient)",
        },
        "theta_convention_resolved": "south/colatitude" if theta_south
                                     else "north/latitude",
        "sampling": {"off_pole_per_altitude": 128,
                     "near_pole_per_altitude": 16,
                     "total_unique_points": 432,
                     "altitudes_km": list(ALTS_KM)},
        "shtools_rows": rows,
        "independent_rows": indep_rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
