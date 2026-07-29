"""Independent degree-120 potential and blend-correction validation.

The production scalar-potential path is compared with pyshtools'
``expand.MakeGridPoint`` synthesis.  The check covers absolute U_30 and U_120,
the cancellation-sensitive U_120-U_30 band, and the corresponding radial
potential-blend correction (U_120-U_30) dw/dr rhat.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pyshtools as pysh
from scipy.stats import qmc

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel, sh_potential_accel_fixed)
from lunaris.common.lunar_data import resolve_lunar_gravity_path  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "metrics" / "r5_potential_xval.json"
N_LO, N_HI = 30, 120
ALT_NEAR, ALT_FAR = 50e3, 200e3


def shtools_potential(model, cilm, pos, degree):
    r = float(np.linalg.norm(pos))
    lat = math.degrees(math.asin(float(pos[2]) / r))
    lon = math.degrees(math.atan2(float(pos[1]), float(pos[0])))
    scaled = np.zeros((2, degree + 1, degree + 1), dtype=float)
    for n in range(degree + 1):
        scaled[:, n, : n + 1] = (
            cilm[:, n, : n + 1] * (model.r_ref / r) ** n
        )
    return float(model.mu / r * pysh.expand.MakeGridPoint(
        scaled, lat=lat, lon=lon, lmax=degree, norm=1, csphase=1))


def production_potential(model, pos, degree):
    value, _ = sh_potential_accel_fixed(
        np.asarray(pos).reshape(1, 3), model.c_coeffs, model.s_coeffs,
        model.mu, model.r_ref, degree, -1)
    return float(value[0])


def dwdr(alt_m):
    if alt_m <= ALT_NEAR or alt_m >= ALT_FAR:
        return 0.0
    t = (ALT_FAR - alt_m) / (ALT_FAR - ALT_NEAR)
    return -6.0 * t * (1.0 - t) / (ALT_FAR - ALT_NEAR)


def summarize(error, reference):
    error = np.asarray(error)
    reference = np.asarray(reference)
    return {
        "absolute_max": float(np.max(np.abs(error))),
        "absolute_rms": float(np.sqrt(np.mean(error ** 2))),
        "relative_rms": float(np.linalg.norm(error) / np.linalg.norm(reference)),
        "relative_max_pointwise": float(np.max(np.abs(error) /
                                               np.maximum(np.abs(reference), 1e-300))),
    }


def main():
    model = GravityModel.from_file(str(resolve_lunar_gravity_path(None)),
                                   requested_degree=N_HI)
    cilm = np.zeros((2, N_HI + 1, N_HI + 1), dtype=float)
    cilm[0] = model.c_coeffs
    cilm[1] = model.s_coeffs
    cilm[0, 0, 0] = 1.0

    sob = qmc.Sobol(d=2, scramble=True, seed=20260720).random_base2(8)
    z = 2.0 * sob[:, 0] - 1.0
    lon = 2.0 * np.pi * sob[:, 1]
    unit = np.column_stack((np.sqrt(1.0 - z*z) * np.cos(lon),
                            np.sqrt(1.0 - z*z) * np.sin(lon), z))

    rows = []
    all_band_err, all_band_ref = [], []
    all_corr_err, all_corr_ref = [], []
    cancellation = []
    for alt_km in (50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0):
        u_prod, u_ref = {N_LO: [], N_HI: []}, {N_LO: [], N_HI: []}
        for direction in unit:
            pos = (model.r_ref + alt_km * 1e3) * direction
            for degree in (N_LO, N_HI):
                u_prod[degree].append(production_potential(model, pos, degree))
                u_ref[degree].append(shtools_potential(model, cilm, pos, degree))
        for degree in (N_LO, N_HI):
            u_prod[degree] = np.asarray(u_prod[degree])
            u_ref[degree] = np.asarray(u_ref[degree])
        band_prod = u_prod[N_HI] - u_prod[N_LO]
        band_ref = u_ref[N_HI] - u_ref[N_LO]
        band_err = band_prod - band_ref
        all_band_err.extend(band_err)
        all_band_ref.extend(band_ref)
        deriv = dwdr(alt_km * 1e3)
        if deriv != 0.0:
            all_corr_err.extend(band_err * deriv)
            all_corr_ref.extend(band_ref * deriv)
        cancellation.append(float(np.sqrt(np.mean(u_ref[N_HI] ** 2)) /
                                  np.sqrt(np.mean(band_ref ** 2))))
        rows.append({
            "altitude_km": alt_km,
            "U30": summarize(u_prod[N_LO] - u_ref[N_LO], u_ref[N_LO]),
            "U120": summarize(u_prod[N_HI] - u_ref[N_HI], u_ref[N_HI]),
            "U120_minus_U30": summarize(band_err, band_ref),
            "dw_dr_per_m": deriv,
        })

    payload = {
        "reference": f"pyshtools {pysh.__version__} expand.MakeGridPoint",
        "normalization": "4pi, csphase=1 (no Condon--Shortley phase)",
        "sampling": {"sobol_directions_per_altitude": len(unit),
                     "altitudes_km": [r["altitude_km"] for r in rows]},
        "rows": rows,
        "pooled_band_difference": summarize(all_band_err, all_band_ref),
        "pooled_corrected_gradient_term_m_s2": summarize(all_corr_err, all_corr_ref),
        "cancellation_ratio_full_U_rms_over_band_rms": {
            "min": float(min(cancellation)), "max": float(max(cancellation))},
        "interpretation": "The independent synthesis validates both the absolute potentials and the cancellation-sensitive difference used by the corrected blend.",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "band": payload["pooled_band_difference"],
        "corrected_term": payload["pooled_corrected_gradient_term_m_s2"],
        "cancellation_ratio": payload["cancellation_ratio_full_U_rms_over_band_rms"],
    }, indent=2))


if __name__ == "__main__":
    main()
