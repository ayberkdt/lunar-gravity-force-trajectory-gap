"""Self-contained check that the measurement instrument works.

Loads the Lunar Prospector LP150Q product shipped in ``data/`` (1.4 MB, no
external download needed), evaluates the spherical-harmonic acceleration at a
few fixed points and degrees, and compares against values recorded from the
source tree that produced the published results.

    python verify_snapshot.py

Exit status is 0 when every value matches to the last bit.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel,
    sh_accel_fixed_numba,
)

GRAVITY_FILE = ROOT / "data" / "jgl150q1.sha"
POINTS = [(1.9e6, 3.0e5, 7.0e5), (1.2e6, -1.1e6, 1.0e6), (1.0e5, 2.0e5, 1.85e6)]
DEGREES = (20, 50, 80)

# Recorded from the source tree that produced the published results.
EXPECTED = [
    [-1.086496168467129, -0.17146674522897332, -0.4003022727720085],
    [-0.8433169411912707, 0.7735600592964192, -0.703384692735618],
    [-0.07537696176016051, -0.1514498237439886, -1.4012006846967908],
    [-1.0864925903361602, -0.17147076228739955, -0.4002991534046025],
    [-0.8433143912778865, 0.7735674111419282, -0.703392058838559],
    [-0.07539755682304185, -0.1514368452456295, -1.4011948208892162],
    [-1.0864925908595566, -0.17147074219382977, -0.40029916869967513],
    [-0.8433144576227641, 0.7735674393306156, -0.703391738223319],
    [-0.07540168688200513, -0.15143644150163138, -1.4011989460683283],
]
EXPECTED_MU = 4902801076000.0
EXPECTED_R_REF = 1738000.0


def main() -> int:
    if not GRAVITY_FILE.exists():
        print(f"[fail] missing {GRAVITY_FILE}")
        return 1

    model = GravityModel.from_file(str(GRAVITY_FILE), requested_degree=80)
    ws = model.make_workspace()
    args = (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs, model.a_coeffs,
            model.b_coeffs, model.scale_m_table, ws.P, ws.dP, ws.cos_m,
            ws.sin_m)

    ok = True
    if float(model.mu) != EXPECTED_MU or float(model.r_ref) != EXPECTED_R_REF:
        print(f"[fail] header mismatch: mu={model.mu} r_ref={model.r_ref}")
        ok = False

    got = []
    for degree in DEGREES:
        for x, y, z in POINTS:
            got.append([float(v) for v in
                        sh_accel_fixed_numba(x, y, z, degree, *args)])

    for i, (g, e) in enumerate(zip(got, EXPECTED)):
        if g != e:
            print(f"[fail] point {i}: {g} != {e}")
            ok = False

    if ok:
        print(f"[ok] {len(got)} accelerations bitwise-identical to the "
              f"published snapshot")
        print(f"[ok] kernel: {sh_accel_fixed_numba.__module__}")
        return 0
    print(json.dumps(got, indent=2))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
