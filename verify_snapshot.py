"""Self-contained check that the measurement instrument works.

Loads the Lunar Prospector LP150Q product shipped in ``data/`` (1.4 MB, no
external download needed), evaluates the spherical-harmonic acceleration at a
few fixed points and degrees, and compares against values recorded from the
archived snapshot used for the submitted manuscript.

    python verify_snapshot.py                   # cross-platform tolerance check
    python verify_snapshot.py --strict-bitwise  # exact equality

The default compares against a declared relative tolerance. Exact reproduction
to the last bit depends on the LLVM version Numba compiles through, the CPU's
fused-multiply-add behaviour and the platform's libm, so it is expected on the
archived Windows 11 x64 environment and not guaranteed elsewhere. A difference
at the tolerance level is a portability artifact, not a broken archive; a
difference above it is a real problem.

Exit status is 0 when every value passes the selected check.
"""

from __future__ import annotations

import argparse
import math
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

# Tolerance for the default check, relative to the component magnitude. Set
# well above float64 rounding noise and far below any difference that would
# indicate a genuine change in the kernel.
RELATIVE_TOLERANCE = 1e-12

# Reference values recorded from source snapshot
# 27e9ab86ed61d623f78c453ea2054348f1044c23, release tag paper-truncation-v1.0,
# which is the snapshot the campaign manifests pin across R10-R23.
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


def evaluate() -> tuple[GravityModel, list[list[float]]]:
    model = GravityModel.from_file(str(GRAVITY_FILE), requested_degree=80)
    workspace = model.make_workspace()
    args = (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs, model.a_coeffs,
            model.b_coeffs, model.scale_m_table, workspace.P, workspace.dP,
            workspace.cos_m, workspace.sin_m)
    values = []
    for degree in DEGREES:
        for x, y, z in POINTS:
            values.append([float(v) for v in
                           sh_accel_fixed_numba(x, y, z, degree, *args)])
    return model, values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-bitwise", action="store_true",
                        help="require exact equality instead of the default "
                             "relative tolerance")
    args = parser.parse_args()

    if not GRAVITY_FILE.exists():
        print(f"[fail] missing {GRAVITY_FILE}")
        return 1

    model, got = evaluate()

    ok = True
    if float(model.mu) != EXPECTED_MU or float(model.r_ref) != EXPECTED_R_REF:
        print(f"[fail] header mismatch: mu={model.mu} r_ref={model.r_ref}")
        ok = False

    worst = 0.0
    exact = 0
    for index, (actual, expected) in enumerate(zip(got, EXPECTED)):
        for component, (a, e) in enumerate(zip(actual, expected)):
            if a == e:
                exact += 1
                continue
            relative = abs(a - e) / max(abs(e), 1e-300)
            worst = max(worst, relative)
            if args.strict_bitwise:
                print(f"[fail] point {index} component {component}: "
                      f"{a!r} != {e!r}")
                ok = False
            elif not (relative <= RELATIVE_TOLERANCE
                      or math.isclose(a, e, rel_tol=RELATIVE_TOLERANCE)):
                print(f"[fail] point {index} component {component}: "
                      f"{a!r} vs {e!r}, relative {relative:.3e} exceeds "
                      f"{RELATIVE_TOLERANCE:.0e}")
                ok = False

    total = len(EXPECTED) * 3
    if not ok:
        return 1

    if args.strict_bitwise:
        print(f"[ok] {total} components bitwise-identical to the archived "
              f"manuscript snapshot")
    elif exact == total:
        print(f"[ok] {total} components bitwise-identical to the archived "
              f"manuscript snapshot")
    else:
        print(f"[ok] {total} components within {RELATIVE_TOLERANCE:.0e} "
              f"relative tolerance of the archived manuscript snapshot "
              f"({exact} exact, worst relative difference {worst:.3e})")
    print(f"[ok] kernel: {sh_accel_fixed_numba.__module__}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
