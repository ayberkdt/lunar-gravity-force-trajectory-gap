"""Shared harness for the JGCD pre-submission campaigns (J1, J2, J3).

Three questions are asked with one piece of machinery, because the point of all
three is that nothing except the named variable changes.

  J1  Does the force--trajectory reversal survive a *different* GRAIL gravity
      solution? Only the coefficient set changes (GSFC GRGM1200A in place of
      JPL JGGRX_1800F), and the recipe is recalibrated on the new spectrum
      rather than inheriting the primary field's numbers.

  J2  Does it survive *realistic* lunar dynamics? Only the force model changes
      (DE440 ephemerides, MOON_PA orientation, Earth and Sun third-body
      gravity, cannonball SRP with lunar eclipse). The extra accelerations are
      common to reference and policy, so what is still being measured is the
      truncation policy difference.

  J3  Does the ranking survive the numerical-resolution criterion? Only the
      integration tolerance changes, over three levels rather than two.

Everything downstream of those three variables -- the orbit box, the arc, the
policy definitions, the error metric, the resolution rule -- is the code the
manuscript already uses, imported rather than restated.

The field is selected through the environment, not through a function argument,
because ProcessPoolExecutor on Windows re-imports this module in every child
and a parent-only patch would leave the children propagating the wrong field.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

import rev3_common as rc                                          # noqa: E402
from rev3_common import DAY, OMEGA_MOON, REPO, commit_sha, working_tree_clean
from lunaris.physics.spherical_harmonics import (                 # noqa: E402
    GravityModel, sh_accel_fixed_numba)

# --------------------------------------------------------------- field select
FIELD_ENV = "JCAMP_GRAVITY_FILE"
FIELD_KEY_ENV = "JCAMP_GRAVITY_KEY"

FIELDS = {
    "JGGRX_1800F": {
        "center": "JPL", "max_degree_in_file": 1800,
        "path": r"D:\Masaustu\LUNAR_SIMULATION\data\gravity_models"
                r"\jggrx_1800f_sha.tab.txt",
        "sha256": "d2a552067a78bf1d2755807ae14ee1d6843a8f6a4228e01ce59a66551"
                  "6738fec"},
    "GRGM1200A": {
        "center": "GSFC", "max_degree_in_file": 1200,
        "path": r"C:\Users\ayber\Desktop\lunaris external validation"
                r"\gravity_models\gggrx_1200a_sha.tab",
        "sha256": "fa04c3dce9376948ad243f3df74144e2602f12d183ea4d179604ed0a79"
                  "da7ded"},
}


def select_field(key: str) -> None:
    """Declare the field for this process *and every child it spawns*."""
    if key not in FIELDS:
        raise SystemExit(f"unknown field {key!r}: {sorted(FIELDS)}")
    os.environ[FIELD_KEY_ENV] = key
    os.environ[FIELD_ENV] = FIELDS[key]["path"]


def field_key() -> str:
    return os.environ.get(FIELD_KEY_ENV, "JGGRX_1800F")


def field_path() -> Path:
    return Path(os.environ.get(FIELD_ENV, FIELDS["JGGRX_1800F"]["path"]))


def install_field() -> None:
    """Point every model factory in the imported stack at the declared field.

    ``rev10_sobol_confirmatory`` binds ``load_model`` by name at import, and
    the R11/R14 workers reach the model through *that* binding, so patching
    ``rev3_common`` alone is not enough. Both are rebound here, together with
    the path resolver the provenance block reads, so that a record can never
    claim one field while the kernel used another.
    """
    path = field_path()
    if not path.is_file():
        raise SystemExit(f"declared gravity file is missing: {path}")

    def load_model(requested_degree: int = 300):
        return GravityModel.from_file(str(path),
                                      requested_degree=requested_degree)

    rc.load_model = load_model
    import rev10_sobol_confirmatory as base
    base.load_model = load_model
    base.resolve_lunar_gravity_path = lambda _p=None: path
    return None


# ------------------------------------------------------------------ contract
LEVELS = {
    "loose": {"rtol": 1.0e-11, "atol": np.array([1.0e-4] * 3 + [1.0e-7] * 3),
              "atol_position_m": 1.0e-4, "atol_velocity_m_s": 1.0e-7},
    "tight": {"rtol": 1.0e-12, "atol": np.array([1.0e-5] * 3 + [1.0e-8] * 3),
              "atol_position_m": 1.0e-5, "atol_velocity_m_s": 1.0e-8},
    "tighter": {"rtol": 1.0e-13, "atol": np.array([1.0e-6] * 3 + [1.0e-9] * 3),
                "atol_position_m": 1.0e-6, "atol_velocity_m_s": 1.0e-9},
    # SciPy clamps DOP853 to rtol >= 2.220446049250313e-14 and warns. Declaring
    # anything below that would put a number in the record that the integrator
    # never used, so the tightest level is set just above the clamp: it is the
    # tightest relative tolerance this integrator admits, which is itself the
    # honest statement about how far the numerical envelope can be pushed.
    "tightest": {"rtol": 2.5e-14,
                 "atol": np.array([1.0e-7] * 3 + [1.0e-10] * 3),
                 "atol_position_m": 1.0e-7, "atol_velocity_m_s": 1.0e-10,
                 "note": "rtol is at the DOP853 relative-tolerance floor"},
}
MAX_STEP = 60.0
DURATION = 7.0 * DAY
OUTPUT_STEP = 120.0
BIN_KM = 10.0
FLOOR = 2
EPS_TAIL = 1.0e-3          # empirical tail criterion, as in the main text
CAP_CRIT = 250             # the manuscript's cap on the critical degree


def out_grid() -> np.ndarray:
    return np.arange(0.0, DURATION + 0.5 * OUTPUT_STEP, OUTPUT_STEP)


# --------------------------------------------------------------------- worker
_MODELS: dict[int, tuple] = {}
_GCACHE: dict[int, np.ndarray] = {}
_POWER: dict[int, np.ndarray] = {}


def model_for(degree: int):
    if degree not in _MODELS:
        install_field()
        m = rc.load_model(degree)
        a = rc.kernel_args(m)
        rc.warmup(m, a)
        _MODELS[degree] = (m, a)
    return _MODELS[degree]


def power_for(degree: int) -> np.ndarray:
    if degree not in _POWER:
        m, _ = model_for(degree)
        _POWER[degree] = rc.degree_power(m)
    return _POWER[degree]


def atallah_g(degree: int) -> np.ndarray:
    if degree not in _GCACHE:
        import rev12_atallah as at
        m, _ = model_for(degree)
        _GCACHE[degree] = at.precompute_Sn(m, degree)
    return _GCACHE[degree]


# ------------------------------------------------------------------- policies
def emp_nmin_exact(power: np.ndarray, r_ref: float, h_m: float,
                   eps: float = EPS_TAIL) -> int:
    """Unquantized empirical tail-criterion degree at one altitude."""
    n = np.arange(len(power), dtype=np.float64)
    r = r_ref + h_m
    sig = (np.sqrt((n + 1.0) * (2.0 * n + 1.0))
           * np.exp(n * math.log(r_ref / r)) * np.sqrt(power))
    sq = sig ** 2
    total = float(np.sum(sq[2:]))
    budget = eps * eps * total
    tail = total
    for k in range(2, len(sq)):
        if tail <= budget:
            return k - 1
        tail -= sq[k]
    return len(sq) - 1


def critical_degree(power: np.ndarray, r_ref: float, hp_km: float) -> int:
    return int(min(CAP_CRIT, emp_nmin_exact(power, r_ref, hp_km * 1e3)))


def empirical_table(power: np.ndarray, r_ref: float) -> dict:
    """The 10-km-quantized lookup schedule the manuscript's prepass uses."""
    table = {}
    for hk in np.arange(20.0, 561.0, 10.0):
        nmin = emp_nmin_exact(power, r_ref, hk * 1e3)
        table[float(hk)] = max(60, min(CAP_CRIT, (nmin // 10) * 10))
    return table


def table_sched(table: dict):
    hmax, hmin = max(table), min(table)

    def f(t, h_m):
        hb = min(hmax, max(hmin, 10.0 * math.floor(h_m / 1e3 / 10.0)))
        return table[hb]

    return f


def fixed_degree_for(beta: float, n_crit: int, cap: int) -> tuple[int, bool]:
    """N_F(beta) = argmin_N |N^2 - beta N_crit^2|, no interpolation."""
    target = beta * n_crit ** 2
    n_real = math.sqrt(target)
    cands = [max(1, int(math.floor(n_real))), max(1, int(math.ceil(n_real)))]
    n = min(cands, key=lambda k: abs(k ** 2 - target))
    return (min(n, cap), n > cap)


def degrees_from_table(table: dict, h_km: np.ndarray) -> np.ndarray:
    keys = np.array(sorted(table))
    vals = np.array([table[k] for k in keys], dtype=int)
    hb = np.clip(BIN_KM * np.floor(h_km / BIN_KM), min(table), max(table))
    idx = np.clip(np.searchsorted(keys, hb - 1e-9), 0, len(keys) - 1)
    return vals[idx]


def calibrate_radial(degree: int, hp_km: float, ha_km: float, cap: int,
                     h_km: np.ndarray, target_work: float,
                     work_tolerance: float = 0.01) -> dict:
    """Bisect the Atallah accuracy parameter so the binned degree history
    spends exactly the declared budget. This is the R14 Phase-A calibration,
    re-derived for whichever field and altitude history it is handed."""
    import rev12_atallah as at
    model, _ = model_for(degree)
    g = atallah_g(degree)

    def work_of(tol):
        _, table = at.atallah_binned_schedule(model, g, tol, hp_km, ha_km,
                                              floor=FLOOR, cap=cap,
                                              bin_km=BIN_KM)
        table = {float(k): int(v) for k, v in table.items()}
        deg = degrees_from_table(table, h_km)
        return float(np.mean(deg.astype(float) ** 2)), table, deg

    lo_log, hi_log = -18.0, 2.0
    w_hi, t_hi, d_hi = work_of(10.0 ** hi_log)
    w_lo, t_lo, d_lo = work_of(10.0 ** lo_log)
    if target_work > w_lo:
        return {"tol": 10.0 ** lo_log, "table": t_lo, "degrees": d_lo,
                "work": w_lo, "attainable": False, "limit": "cap",
                "mismatch": abs(w_lo / target_work - 1.0)}
    if target_work < w_hi:
        return {"tol": 10.0 ** hi_log, "table": t_hi, "degrees": d_hi,
                "work": w_hi, "attainable": False, "limit": "floor",
                "mismatch": abs(w_hi / target_work - 1.0)}
    best = None
    for _ in range(200):
        mid = 0.5 * (lo_log + hi_log)
        w, table, deg = work_of(10.0 ** mid)
        err = abs(w / target_work - 1.0)
        if best is None or err < best["mismatch"]:
            best = {"tol": 10.0 ** mid, "table": table, "degrees": deg,
                    "work": w, "mismatch": err}
        if w > target_work:
            lo_log = mid
        else:
            hi_log = mid
        if hi_log - lo_log < 1e-12:
            break
    best["attainable"] = bool(best["mismatch"] < work_tolerance)
    best["limit"] = None if best["attainable"] else "integer_degree_discretization"
    return best


# ------------------------------------------------------------- force metrics
def accel_inertial(r_vec, t, n, args):
    """Inertial acceleration of the uniformly rotating body-fixed N-field,
    identical to the propagator's right-hand-side convention."""
    x, y, z = r_vec
    th = OMEGA_MOON * t
    c, s = math.cos(th), math.sin(th)
    axb, ayb, azb = sh_accel_fixed_numba(c * x + s * y, -s * x + c * y, z,
                                         int(n), *args)
    return np.array([c * axb - s * ayb, s * axb + c * ayb, azb])


def _defect_summary(mags: np.ndarray) -> dict:
    """``J_F`` is the time average of the defect magnitude, which is the
    quantity the campaign registrations name; the RMS is carried alongside
    because the manuscript's Pareto records use it and the two must not be
    silently interchanged."""
    return {"J_force_mean_m_s2": float(np.mean(mags)),
            "J_force_rms_m_s2": float(np.sqrt(np.mean(mags ** 2))),
            "J_force_max_m_s2": float(np.max(mags)),
            "n_epochs": int(len(mags))}


def force_defects(t: np.ndarray, R: np.ndarray, degree_sets: dict,
                  truth_degree: int, args,
                  to_fixed=None) -> dict:
    """Deterministic truncation defect along the reference, for several
    policies at once.

    The reference acceleration is the same for every policy at a given epoch,
    so it is evaluated once per epoch rather than once per policy; on the
    degree-900 orbits that is the difference between one and two full kernel
    evaluations per epoch per policy.

    ``to_fixed`` converts an inertial position to body-fixed coordinates. With
    it omitted the uniformly rotating convention of the gravity-only system is
    used; the full-dynamics campaign passes the MOON_PA transform instead. The
    defect is a difference of two accelerations in the same frame, so the
    common third-body and radiation terms cancel exactly and never enter.
    """
    keys = list(degree_sets)
    mags = {k: np.empty(len(t)) for k in keys}
    for i in range(len(t)):
        ti = float(t[i])
        if to_fixed is None:
            a_ref = accel_inertial(R[:, i], ti, truth_degree, args)
            for k in keys:
                d = accel_inertial(R[:, i], ti, int(degree_sets[k][i]),
                                   args) - a_ref
                mags[k][i] = float(np.linalg.norm(d))
        else:
            rb = to_fixed(ti, R[:, i])
            a_ref = np.array(sh_accel_fixed_numba(float(rb[0]), float(rb[1]),
                                                  float(rb[2]), truth_degree,
                                                  *args))
            for k in keys:
                a_p = np.array(sh_accel_fixed_numba(
                    float(rb[0]), float(rb[1]), float(rb[2]),
                    int(degree_sets[k][i]), *args))
                mags[k][i] = float(np.linalg.norm(a_p - a_ref))
    return {k: _defect_summary(mags[k]) for k in keys}


def force_defect(t: np.ndarray, R: np.ndarray, degrees: np.ndarray,
                 truth_degree: int, args) -> dict:
    return force_defects(t, R, {"p": degrees}, truth_degree, args)["p"]


def trajectory_error(Y: np.ndarray, T: np.ndarray) -> dict:
    """Position error against the reference on the common output grid."""
    stats = rc.err_stats(Y, T)
    return {"J_traj_rms_m": stats["pos_rms_m"],
            "pos_max_m": stats["pos_max_m"],
            "pos_final_m": stats["pos_final_m"],
            "ric_rms_m": stats["ric_rms_m"]}


# ------------------------------------------------------------------------ IO
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def object_hash(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":"),
                                     ensure_ascii=False,
                                     allow_nan=False).encode("utf-8")
                          ).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              allow_nan=False), encoding="utf-8")
    os.replace(tmp, path)


def atomic_npz(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with tmp.open("wb") as fh:
        np.savez_compressed(fh, **arrays)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def load_states(path) -> np.ndarray:
    """States from a trajectory archive, whichever naming it uses.

    The manuscript's own raw files store ``t_s`` / ``state_si``; the J-campaign
    files store ``t`` / ``y``. Reading both here keeps the two conventions from
    leaking into the campaign code, and keeps the archived trajectories usable
    unchanged rather than rewritten into a new convention.
    """
    with np.load(path) as z:
        for key in ("state_si", "y"):
            if key in z:
                return z[key]
    raise KeyError(f"no state array in {path}")


def load_times(path) -> np.ndarray:
    with np.load(path) as z:
        for key in ("t_s", "t"):
            if key in z:
                return z[key]
    raise KeyError(f"no time array in {path}")


def provenance() -> dict:
    kernel = (REPO / "src" / "lunaris" / "physics" /
              "spherical_harmonics.py").resolve()
    gpath = field_path()
    return {
        "lunaris_repo": str(REPO.resolve()),
        "lunaris_commit": commit_sha(),
        "lunaris_working_tree_clean": working_tree_clean(),
        "kernel_path": str(kernel),
        "kernel_sha256": sha256_file(kernel),
        "gravity_key": field_key(),
        "gravity_path": str(gpath),
        "gravity_sha256": sha256_file(gpath),
        "harness": str(Path(__file__).resolve()),
        "harness_sha256": sha256_file(Path(__file__).resolve()),
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
    }


# ------------------------------------------------------------ resolution rule
def self_difference(Y_a: np.ndarray, Y_b: np.ndarray) -> float:
    """Position RMS between the two tolerance levels of the same policy."""
    d = np.linalg.norm(Y_a[:3] - Y_b[:3], axis=0)
    return float(np.sqrt(np.mean(d ** 2)))


def resolved(err_a: float, err_b: float, env_a: float, env_b: float) -> bool:
    """The R10 rule, unchanged: a comparison counts only when the gap exceeds
    the two policies' numerical envelopes, each of which already includes the
    reference's own self-difference."""
    return abs(err_a - err_b) > (env_a + env_b)


def log_line(path: Path, msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")
