"""P0-9: off-grid parity of the Lunaris quaternion-table rotation vs SPICE.

The Lunaris SH path samples J2000->MOON_PA as a SPICE-derived quaternion
table (cadence = ephemeris_step_s; normalized, sign-continuous) with SLERP
between nodes at arbitrary integrator-stage times. The archived rotation
probes were taken at exact table nodes; this experiment quantifies parity
on off-grid time sets:

  nodes      : exact table nodes,
  offgrid    : uniform-random inter-node times,
  stagelike  : DOP853 stage times c_i * h from random step sizes 40-300 s.

Windows: 5-day and 30-day; cadences: 5 s (tight contract) and 30 s
(baseline contract). Metrics: Frobenius matrix difference, principal
rotation angle of R_lun R_spice^T, and the induced SH acceleration
difference (degree 120) at 50-km-altitude states.
"""

from __future__ import annotations

import math

import numpy as np

from rev3_common import REPO, SEED, dump, load_model

import sys

sys.path.insert(0, str(REPO / "src"))
from lunaris.physics.ephemeris import EphemerisManager, build_tables  # noqa: E402
from lunaris.physics.spherical_harmonics import sh_accel_fixed_numba  # noqa: E402

import spiceypy as sp  # noqa: E402

KERNEL_DIR = r"C:\Users\ayber\Desktop\lunaris\data\ephemeris_models"
KERNELS = [f"{KERNEL_DIR}\\naif0012.tls", f"{KERNEL_DIR}\\pck00011.tpc",
           f"{KERNEL_DIR}\\gm_de440.tpc", f"{KERNEL_DIR}\\de440s.bsp",
           f"{KERNEL_DIR}\\moon_de440_250416.tf",
           f"{KERNEL_DIR}\\moon_pa_de440_200625.bpc"]
START = "2025-01-01 00:00:00 TDB"
DOP853_C = np.array([0.0, 0.0526001519587677, 0.0789002279381515,
                     0.118350341907227, 0.281649658092773, 0.333333333333333,
                     0.25, 0.307692307692308, 0.651282051282051, 0.6,
                     0.857142857142857, 1.0])


def quat_to_matrix(q):
    w, x, y, z = (float(v) for v in q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


def principal_angle(R_a, R_b):
    """Return the relative rotation angle without the small-angle loss of
    precision of acos((trace(R)-1)/2).

    The atan2 form uses the antisymmetric part for sin(theta) and remains
    well conditioned when the two matrices agree to machine precision.
    """
    D = R_a @ R_b.T
    sin_theta = 0.5 * np.linalg.norm(
        [D[2, 1] - D[1, 2], D[0, 2] - D[2, 0], D[1, 0] - D[0, 1]])
    cos_theta = 0.5 * (np.trace(D) - 1.0)
    return float(math.atan2(sin_theta, cos_theta))


def kernel_args(model):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


def probe_sets(rng, dur_s, dt_s, n=400):
    nodes = (rng.integers(0, int(dur_s / dt_s) + 1, size=n) * dt_s).astype(float)
    offgrid = np.clip((rng.integers(0, int(dur_s / dt_s), size=n)
                       + rng.uniform(0.05, 0.95, size=n)) * dt_s, 0, dur_s)
    base = rng.uniform(0.0, dur_s - 400.0, size=n)
    h = rng.uniform(40.0, 300.0, size=n)
    ci = rng.choice(DOP853_C, size=n)
    stagelike = np.clip(base + ci * h, 0.0, dur_s)
    return {"nodes": nodes, "offgrid": offgrid, "stagelike": stagelike}


def main() -> int:
    model = load_model(120)
    args = kernel_args(model)
    sh_accel_fixed_numba(model.r_ref + 50e3, 0.0, 0.0, 120, *args)
    rng = np.random.default_rng(SEED)

    for k in KERNELS:
        sp.furnsh(k)
    et0 = sp.str2et(START)

    test_dirs = rng.normal(size=(20, 3))
    test_dirs /= np.linalg.norm(test_dirs, axis=1, keepdims=True)
    r_test = model.r_ref + 50e3

    results = []
    for dur_days, dt_s in ((5.0, 5.0), (5.0, 30.0), (30.0, 5.0), (30.0, 30.0)):
        dur = dur_days * 86400.0
        tables = build_tables(
            start_utc=START, duration_s=dur, output_dt_s=dt_s,
            kernels=tuple(KERNELS), inertial_frame="J2000",
            fixed_frame="MOON_PA", observer="MOON",
            include_third_body=False, clear_kernels_after=False,
            clean_kernels_before=False, auto_fix_kernel_paths=False,
            need_moon_fixed_rotation=True)
        mgr = EphemerisManager.from_tables(tables)

        for set_name, times in probe_sets(rng, dur, dt_s).items():
            fro, ang, dacc = [], [], []
            for t in times:
                q = mgr.get_inertial_to_fixed_rotation(float(t))
                R_lun = quat_to_matrix(q)
                R_sp = sp.pxform("J2000", "MOON_PA", et0 + float(t))
                fro.append(float(np.linalg.norm(R_lun - R_sp, ord="fro")))
                ang.append(principal_angle(R_lun, R_sp))
            # induced acceleration difference at the worst probe
            iworst = int(np.argmax(fro))
            t = float(times[iworst])
            R_lun = quat_to_matrix(mgr.get_inertial_to_fixed_rotation(t))
            R_sp = sp.pxform("J2000", "MOON_PA", et0 + t)
            for u in test_dirs:
                xi = r_test * u
                a1 = np.array(sh_accel_fixed_numba(*(R_lun @ xi), 120, *args))
                a2 = np.array(sh_accel_fixed_numba(*(R_sp @ xi), 120, *args))
                ai1 = R_lun.T @ a1
                ai2 = R_sp.T @ a2
                dacc.append(float(np.linalg.norm(ai1 - ai2) /
                                  np.linalg.norm(ai2)))
            fro = np.array(fro)
            ang = np.array(ang)
            results.append({
                "duration_days": dur_days, "table_dt_s": dt_s,
                "probe_set": set_name, "n_probes": len(times),
                "frobenius_max": float(fro.max()),
                "frobenius_rms": float(np.sqrt(np.mean(fro**2))),
                "principal_angle_max_rad": float(ang.max()),
                "principal_angle_rms_rad": float(np.sqrt(np.mean(ang**2))),
                "induced_pos_max_m_at_r1788km": float(ang.max() * r_test),
                "induced_accel_rel_max_deg120_at_worst_probe":
                    float(np.max(dacc)),
            })
            print(f"{dur_days:4.0f} d dt={dt_s:4.0f} s {set_name:9s}: "
                  f"fro max {fro.max():.3e}, angle max {ang.max():.3e} rad, "
                  f"induced accel rel {np.max(dacc):.3e}")

    dump("r3_rotation_parity.json", {
        "seed": SEED,
        "policy": {
            "table": "SPICE-sampled J2000->MOON_PA quaternion table",
            "interpolation": "SLERP between adjacent nodes "
                             "(interp_quat_slerp)",
            "sign_continuity": "enforced at table build: q[k] flipped if "
                               "dot(q[k-1], q[k]) < 0",
            "normalization": "per-node unit normalization at table build",
            "cadences_s": [5.0, 30.0],
            "epoch": START,
        },
        "rows": results,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
