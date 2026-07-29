"""Robustness controls requested after the publication-readiness review.

The experiment set contains four deliberately scoped controls:

1. A PEFRL symplectic integration of the autonomous, body-fixed field.  A
   fixed-degree Hamiltonian is compared with a discontinuous 30/120-degree
   switch over a band-crossing orbit, including step refinement.
2. A seven-day solver-independence comparison: DOP853 and Radau each compare
   fixed/scheduled policies with their own same-solver N=300 reference.
3. A seven-day expanded-force comparison with DE440/MOON_PA orientation,
   Earth and Sun differential third-body gravity, and cannonball SRP with
   lunar eclipse.  These deterministic forces are common to policy and
   reference; they are not treated as a stochastic noise floor.
4. A seven-day LRO-design-geometry control initialized at 30 x 216 km,
   i=90 deg, argument of periapsis=270 deg.  This is called LRO-like or
   quasi-frozen geometry, not a newly differential-corrected frozen orbit.

Outputs:
  metrics/r4_robustness_controls.json
  figures/fig_symplectic_switch.pdf
  figures/fig_robustness_controls.pdf
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from paper_style import C1, C2, C3, C4, C5, C6, apply as apply_style
from rev3_common import (
    DAY,
    Rhs,
    alt_sched,
    commit_sha,
    dump,
    eccentric_state,
    err_stats,
    kernel_args,
    load_model,
    make_p_table,
    warmup,
    working_tree_clean,
)


BASE = Path(__file__).resolve().parents[1]
FIGS = BASE / "figures"
METRICS = BASE / "metrics"
REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))

from lunaris.common.math_utils import coe_to_rv  # noqa: E402
from lunaris.common.type_defs import SpacecraftProps  # noqa: E402
from lunaris.core.propagation.integrators.symplectic import _pefrl_step  # noqa: E402
from lunaris.physics.ephemeris import EphemerisManager, build_tables  # noqa: E402
from lunaris.physics.solar_effects import SRPConfig, accel_srp  # noqa: E402
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    sh_accel_fixed_numba,
    sh_potential_accel_fixed,
)
from lunaris.physics.third_body_effects import calc_3rd_body_accel  # noqa: E402


KERNEL_DIR = Path(r"C:\Users\ayber\Desktop\lunaris\data\ephemeris_models")
SPICE_KERNELS = (
    str(KERNEL_DIR / "naif0012.tls"),
    str(KERNEL_DIR / "pck00011.tpc"),
    str(KERNEL_DIR / "gm_de440.tpc"),
    str(KERNEL_DIR / "de440s.bsp"),
    str(KERNEL_DIR / "moon_de440_250416.tf"),
    str(KERNEL_DIR / "moon_pa_de440_200625.bpc"),
)
START_UTC = "2025-01-01T00:00:00Z"
RTOL = 1.0e-11
ATOL = np.array([1.0e-4] * 3 + [1.0e-7] * 3)
MAX_STEP_S = 60.0
OUT_STEP_S = 120.0


def _run_ivp(model, args, y0, degree_of, *, method, duration_s,
             rtol=RTOL, atol=ATOL, max_step=MAX_STEP_S):
    rhs = Rhs(model, degree_of, args)
    grid = np.arange(0.0, duration_s + 0.1, OUT_STEP_S)
    t0 = time.perf_counter()
    sol = solve_ivp(
        rhs,
        (0.0, duration_s),
        y0,
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        t_eval=grid,
    )
    if not sol.success:
        raise RuntimeError(f"{method} failed: {sol.message}")
    return sol, {
        "method": method,
        "n_rhs": int(rhs.n_calls),
        "wall_s": float(time.perf_counter() - t0),
        "gravity_wall_s": float(rhs.grav_ns / 1.0e9),
        "mean_degree_squared": float(rhs.sum_deg_sq / rhs.n_calls),
    }


def solver_control(model, args, duration_s: float) -> dict:
    print("== solver independence ==", flush=True)
    y0 = eccentric_state(model, 50.0, 300.0)
    schedule = alt_sched(make_p_table(model, 1.0e-3, 60, policy="down"))
    policies = {
        "reference_N300": lambda t, h: 300,
        "fixed_N138": lambda t, h: 138,
        "schedule_down": schedule,
    }
    out: dict[str, dict] = {}
    states: dict[str, dict[str, np.ndarray]] = {}
    for method in ("DOP853", "Radau"):
        out[method] = {"runs": {}, "errors_vs_same_solver_N300": {}}
        states[method] = {}
        for name, degree_of in policies.items():
            sol, info = _run_ivp(
                model, args, y0, degree_of, method=method,
                duration_s=duration_s,
            )
            states[method][name] = sol.y
            out[method]["runs"][name] = info
            print(
                f"  {method:7s} {name:16s}: rhs={info['n_rhs']:7d} "
                f"wall={info['wall_s']:7.1f}s",
                flush=True,
            )
        ref = states[method]["reference_N300"]
        for name in ("fixed_N138", "schedule_down"):
            out[method]["errors_vs_same_solver_N300"][name] = err_stats(
                states[method][name], ref
            )
    out["cross_solver_state_differences"] = {
        name: err_stats(states["Radau"][name], states["DOP853"][name])
        for name in policies
    }
    out["scenario"] = {
        "duration_s": duration_s,
        "orbit": "50 x 300 km polar, perilune start, uniform rotation",
        "rtol": RTOL,
        "atol_position_m": float(ATOL[0]),
        "atol_velocity_m_s": float(ATOL[3]),
        "max_step_s": MAX_STEP_S,
        "output_step_s": OUT_STEP_S,
    }
    return out


class ExpandedForceRhs:
    """MOON_PA SH plus Earth/Sun differential gravity and eclipsed SRP."""

    def __init__(self, model, args, ephem, degree_of, *, use_third_body=True,
                 use_srp=True, track_components=False):
        self.model = model
        self.args = args
        self.ephem = ephem
        self.degree_of = degree_of
        self.use_third_body = bool(use_third_body)
        self.use_srp = bool(use_srp)
        self.track_components = bool(track_components)
        self.sc = SpacecraftProps(mass_kg=100.0, area_m2=1.0, cr=1.3)
        self.srp = SRPConfig(enable_moon_eclipse=True, enable_earth_eclipse=False)
        self.n_calls = 0
        self.sum_degree = 0.0
        self.sum_deg_sq = 0.0
        self.min_degree = 10 ** 9
        self.max_degree = -1
        self.degree_counts: dict[int, int] = {}
        self.grav_ns = 0
        self.component_sum_sq = {"earth_3b": 0.0, "sun_3b": 0.0, "srp": 0.0}
        self._rb = np.empty(3)
        self._ai = np.empty(3)
        self._earth = np.empty(3)
        self._sun = np.empty(3)

    def __call__(self, t, y):
        self.n_calls += 1
        r_i = np.asarray(y[:3], dtype=float)
        self.ephem.transform_inertial_to_fixed(t, r_i, out=self._rb)
        altitude_m = float(np.linalg.norm(r_i) - self.model.r_ref)
        degree = int(self.degree_of(t, altitude_m))
        self.sum_degree += float(degree)
        self.sum_deg_sq += float(degree * degree)
        self.min_degree = min(self.min_degree, degree)
        self.max_degree = max(self.max_degree, degree)
        self.degree_counts[degree] = self.degree_counts.get(degree, 0) + 1
        t0 = time.perf_counter_ns()
        axb, ayb, azb = sh_accel_fixed_numba(
            float(self._rb[0]), float(self._rb[1]), float(self._rb[2]),
            degree, *self.args,
        )
        self.grav_ns += time.perf_counter_ns() - t0
        self.ephem.transform_fixed_to_inertial(
            t, np.array([axb, ayb, azb]), out=self._ai
        )
        acc = self._ai.copy()
        self.ephem.get_earth_position(t, out=self._earth)
        self.ephem.get_sun_position(t, out=self._sun)
        a_earth = np.zeros(3)
        a_sun = np.zeros(3)
        a_srp = np.zeros(3)
        if self.use_third_body:
            a_earth = calc_3rd_body_accel(r_i, self._earth, self.ephem.tables.mu_earth_m3s2)
            a_sun = calc_3rd_body_accel(r_i, self._sun, self.ephem.tables.mu_sun_m3s2)
            acc += a_earth + a_sun
        if self.use_srp:
            sx, sy, sz = accel_srp(
                float(r_i[0]), float(r_i[1]), float(r_i[2]),
                float(self._sun[0]), float(self._sun[1]), float(self._sun[2]),
                float(self._earth[0]), float(self._earth[1]), float(self._earth[2]),
                float(self.srp.R_moon_m), float(self.srp.R_earth_m),
                float(self.srp.AU_m), float(self.srp.P0),
                float(self.sc.cr), float(self.sc.area_m2), float(self.sc.mass_kg),
                bool(self.srp.enable_moon_eclipse),
                bool(self.srp.enable_earth_eclipse),
            )
            a_srp[:] = (sx, sy, sz)
            acc += a_srp
        if self.track_components:
            self.component_sum_sq["earth_3b"] += float(a_earth @ a_earth)
            self.component_sum_sq["sun_3b"] += float(a_sun @ a_sun)
            self.component_sum_sq["srp"] += float(a_srp @ a_srp)
        return np.array([y[3], y[4], y[5], acc[0], acc[1], acc[2]])

    def info(self, wall_s):
        out = {
            "n_rhs": int(self.n_calls),
            "wall_s": float(wall_s),
            "gravity_wall_s": float(self.grav_ns / 1.0e9),
            "mean_active_degree": float(self.sum_degree / self.n_calls),
            "mean_degree_squared": float(self.sum_deg_sq / self.n_calls),
            "rms_active_degree": float(math.sqrt(self.sum_deg_sq / self.n_calls)),
            "min_active_degree": int(self.min_degree),
            "max_active_degree": int(self.max_degree),
            "degree_counts": {str(k): int(v) for k, v in sorted(self.degree_counts.items())},
        }
        if self.track_components:
            out["component_rms_acceleration_m_s2"] = {
                key: math.sqrt(value / self.n_calls)
                for key, value in self.component_sum_sq.items()
            }
        return out


def _run_expanded(model, args, ephem, y0, degree_of, duration_s, *,
                  use_third_body=True, use_srp=True, track_components=False):
    rhs = ExpandedForceRhs(
        model, args, ephem, degree_of,
        use_third_body=use_third_body,
        use_srp=use_srp,
        track_components=track_components,
    )
    grid = np.arange(0.0, duration_s + 0.1, OUT_STEP_S)
    t0 = time.perf_counter()
    sol = solve_ivp(
        rhs, (0.0, duration_s), y0, method="DOP853", rtol=RTOL,
        atol=ATOL, max_step=MAX_STEP_S, t_eval=grid,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, rhs.info(time.perf_counter() - t0)


def build_ephemeris(duration_s):
    print("== building DE440/MOON_PA ephemeris tables ==", flush=True)
    tables = build_tables(
        start_utc=START_UTC,
        duration_s=duration_s,
        output_dt_s=60.0,
        kernels=SPICE_KERNELS,
        inertial_frame="J2000",
        fixed_frame="MOON_PA",
        observer="MOON",
        include_third_body=True,
        clear_kernels_after=True,
        need_moon_fixed_rotation=True,
    )
    return EphemerisManager.from_tables(tables)


def expanded_force_control(model, args, ephem, duration_s) -> dict:
    print("== expanded force control ==", flush=True)
    y0 = eccentric_state(model, 50.0, 300.0)
    schedule_down = alt_sched(make_p_table(model, 1.0e-3, 60, policy="down"))
    schedule_up = alt_sched(make_p_table(model, 1.0e-3, 60, policy="up"))
    policies = {
        "reference_N300": lambda t, h: 300,
        "fixed_N138": lambda t, h: 138,
        "schedule_down": schedule_down,
        "schedule_up": schedule_up,
    }
    states = {}
    runs = {}
    for name, degree_of in policies.items():
        sol, info = _run_expanded(
            model, args, ephem, y0, degree_of, duration_s,
            track_components=name == "reference_N300",
        )
        states[name] = sol.y
        runs[name] = info
        print(
            f"  full {name:16s}: rhs={info['n_rhs']:7d} "
            f"wall={info['wall_s']:7.1f}s",
            flush=True,
        )
    no_srp, no_srp_info = _run_expanded(
        model, args, ephem, y0, policies["reference_N300"], duration_s,
        use_third_body=True, use_srp=False,
    )
    errors = {
        name: err_stats(states[name], states["reference_N300"])
        for name in policies if name != "reference_N300"
    }
    return {
        "scenario": {
            "duration_s": duration_s,
            "orbit": "50 x 300 km polar, perilune start",
            "orientation": "DE440/MOON_PA, 60 s quaternion table",
            "forces": ["JGGRX SH", "Earth point-mass differential gravity",
                       "Sun point-mass differential gravity",
                       "cannonball SRP with lunar eclipse"],
            "spacecraft": {"mass_kg": 100.0, "area_m2": 1.0, "cr": 1.3},
            "rtol": RTOL,
            "atol_position_m": float(ATOL[0]),
            "atol_velocity_m_s": float(ATOL[3]),
            "max_step_s": MAX_STEP_S,
        },
        "runs": runs,
        "errors_vs_full_force_N300": errors,
        "srp_on_minus_off_N300": err_stats(states["reference_N300"], no_srp.y),
        "third_body_only_N300_run": no_srp_info,
        "interpretation": "Earth/Sun/SRP are common deterministic forces, not an error floor; policy errors are measured against an N=300 run with the same forces.",
    }


def _rotation_matrix_i2f(ephem, t):
    eye = np.eye(3)
    return np.column_stack([
        ephem.transform_inertial_to_fixed(t, eye[:, j]) for j in range(3)
    ])


def lro_like_state(model, ephem):
    rp = model.r_ref + 30.0e3
    ra = model.r_ref + 216.0e3
    a = 0.5 * (rp + ra)
    e = (ra - rp) / (ra + rp)
    r_f, v_f = coe_to_rv(
        a, e, math.radians(90.0), math.radians(0.0),
        math.radians(270.0), 0.0, model.mu,
    )
    M = _rotation_matrix_i2f(ephem, 0.0)
    # The ephemeris table clamps negative query times to its first sample, so a
    # centred difference at the epoch would underestimate the frame rate.  Use
    # the second-order forward derivative instead.
    dt = 1.0
    Mdot = (
        -3.0 * M
        + 4.0 * _rotation_matrix_i2f(ephem, dt)
        - _rotation_matrix_i2f(ephem, 2.0 * dt)
    ) / (2.0 * dt)
    r_i = M.T @ r_f
    v_i = M.T @ (v_f - Mdot @ r_i)
    return np.concatenate((r_i, v_i)), {
        "perilune_altitude_km": 30.0,
        "apolune_altitude_km": 216.0,
        "semi_major_axis_km": a / 1.0e3,
        "eccentricity_derived_from_apses": e,
        "inclination_deg_moon_fixed": 90.0,
        "raan_deg_moon_fixed_declared": 0.0,
        "argument_of_periapsis_deg_moon_fixed": 270.0,
        "true_anomaly_deg": 0.0,
        "source_scope": "LRO commissioning-orbit design geometry; not a new differential correction",
    }


def lro_geometry_control(model, args, ephem, duration_s) -> dict:
    print("== LRO-like quasi-frozen geometry control ==", flush=True)
    y0, orbit_meta = lro_like_state(model, ephem)
    schedule_down = alt_sched(make_p_table(model, 1.0e-3, 60, policy="down"))
    schedule_up = alt_sched(make_p_table(model, 1.0e-3, 60, policy="up"))
    policies = {
        "reference_N300": lambda t, h: 300,
        "fixed_N138": lambda t, h: 138,
        "fixed_N194_empirical": lambda t, h: 194,
        "fixed_N219_calibrated": lambda t, h: 219,
        # Nearest integer degrees to sqrt(<N^2>) of the archived schedules.
        # These are quadratic-work matched, not wall-time matched.
        "fixed_N104_cost_down": lambda t, h: 104,
        "fixed_N112_cost_up": lambda t, h: 112,
        "schedule_down": schedule_down,
        "schedule_up": schedule_up,
    }
    states, runs = {}, {}
    for name, degree_of in policies.items():
        sol, info = _run_expanded(
            model, args, ephem, y0, degree_of, duration_s,
            use_third_body=False, use_srp=False,
        )
        states[name] = sol.y
        runs[name] = info
        print(
            f"  LRO  {name:16s}: rhs={info['n_rhs']:7d} "
            f"wall={info['wall_s']:7.1f}s",
            flush=True,
        )
    ref_alt = np.linalg.norm(states["reference_N300"][:3], axis=0) - model.r_ref
    return {
        "scenario": {
            "duration_s": duration_s,
            "orientation": "DE440/MOON_PA, 60 s quaternion table",
            "forces": ["JGGRX SH only"],
            "orbit": orbit_meta,
            "reference_altitude_range_km": [
                float(np.min(ref_alt) / 1.0e3), float(np.max(ref_alt) / 1.0e3)
            ],
            "qualification": "LRO-like quasi-frozen design geometry control, not a validated frozen-family discovery",
            "fixed_comparators": {
                "legacy": 138,
                "empirical_critical_altitude": 194,
                "calibrated_critical_altitude": 219,
                "quadratic_work_matched_to_schedule_down": 104,
                "quadratic_work_matched_to_schedule_up": 112,
            },
            "cost_match_definition": "nearest integer fixed degree to sqrt(mean N^2) of the archived matching schedule; wall time is reported separately",
        },
        "runs": runs,
        "errors_vs_same_force_N300": {
            name: err_stats(states[name], states["reference_N300"])
            for name in policies if name != "reference_N300"
        },
    }


def _potential(model, pos, degree):
    V, _ = sh_potential_accel_fixed(
        np.asarray(pos, float).reshape(1, 3),
        model.c_coeffs, model.s_coeffs, model.mu, model.r_ref,
        int(degree), -1,
    )
    return float(V[0])


def _symplectic_run(model, args, y0, period, *, step_s, switched):
    threshold_m = 145.0e3

    def degree(pos):
        if not switched:
            return 120
        return 120 if np.linalg.norm(pos) - model.r_ref <= threshold_m else 30

    def accel(t, y):
        n = degree(y[:3])
        return np.asarray(sh_accel_fixed_numba(
            float(y[0]), float(y[1]), float(y[2]), n, *args
        ))

    n_steps = int(math.ceil(period / step_s))
    y = np.asarray(y0, float).copy()
    t = 0.0
    sample_every = max(1, int(round(10.0 / step_s)))
    rows = []
    switches = []
    last_n = degree(y[:3])
    for k in range(n_steps + 1):
        if k % sample_every == 0 or k == n_steps:
            n = degree(y[:3])
            U = _potential(model, y[:3], n)
            E = 0.5 * float(y[3:] @ y[3:]) - U
            rows.append((t, E, n))
        if k == n_steps:
            break
        h = min(step_s, period - t)
        y = _pefrl_step(accel, t, y, h)
        t += h
        new_n = degree(y[:3])
        if new_n != last_n:
            U_before = _potential(model, y[:3], last_n)
            U_after = _potential(model, y[:3], new_n)
            switches.append({
                "t_s": float(t), "from_degree": int(last_n),
                "to_degree": int(new_n),
                "specific_energy_jump_J_kg": float(-(U_after - U_before)),
            })
            last_n = new_n
    arr = np.asarray(rows, float)
    rel = (arr[:, 1] - arr[0, 1]) / abs(arr[0, 1])
    return {
        "step_s": float(step_s),
        "switched": bool(switched),
        "max_abs_relative_energy_change": float(np.max(np.abs(rel))),
        "final_relative_energy_change": float(rel[-1]),
        "switches": switches,
        "series": {
            "t_over_period": (arr[:, 0] / period).tolist(),
            "relative_energy_change": rel.tolist(),
            "degree": arr[:, 2].astype(int).tolist(),
        },
    }


def symplectic_control() -> dict:
    print("== symplectic switching control ==", flush=True)
    model = load_model(120)
    args = kernel_args(model)
    warmup(model, args)
    rp = model.r_ref + 30.0e3
    ra = model.r_ref + 260.0e3
    a = 0.5 * (rp + ra)
    vp = math.sqrt(model.mu * (2.0 / rp - 1.0 / a))
    y0 = np.array([rp, 0.0, 0.0, 0.0, vp, 0.0])
    period = 2.0 * math.pi * math.sqrt(a ** 3 / model.mu)
    runs = {
        "fixed_step2": _symplectic_run(
            model, args, y0, period, step_s=2.0, switched=False
        ),
        "switch_step2": _symplectic_run(
            model, args, y0, period, step_s=2.0, switched=True
        ),
        "switch_step1": _symplectic_run(
            model, args, y0, period, step_s=1.0, switched=True
        ),
    }
    print(
        "  max |dE/E0| fixed2={:.3e}, switch2={:.3e}, switch1={:.3e}".format(
            runs["fixed_step2"]["max_abs_relative_energy_change"],
            runs["switch_step2"]["max_abs_relative_energy_change"],
            runs["switch_step1"]["max_abs_relative_energy_change"],
        ),
        flush=True,
    )
    return {
        "scenario": {
            "integrator": "PEFRL fourth-order symplectic split",
            "frame": "autonomous non-rotating body-fixed control",
            "orbit": "30 x 260 km equatorial, one revolution",
            "period_s": period,
            "fixed_degree": 120,
            "switch_degrees": [30, 120],
            "switch_altitude_km": 145.0,
            "energy_definition": "specific E = |v|^2/2 - U_N; piecewise U_N for the switch",
        },
        "runs": runs,
    }


def make_figures(payload):
    apply_style()
    FIGS.mkdir(exist_ok=True)

    sym = payload["symplectic_control"]
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(5.3, 4.4), sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.12},
    )
    for key, color, ls, label in (
        ("fixed_step2", C1, "-", "fixed $N=120$, PEFRL 2 s"),
        ("switch_step2", C5, "-", "switch $30/120$, PEFRL 2 s"),
        ("switch_step1", C2, "--", "switch $30/120$, PEFRL 1 s"),
    ):
        s = sym["runs"][key]["series"]
        x = np.asarray(s["t_over_period"])
        y = np.asarray(s["relative_energy_change"])
        ax.plot(x, y, color=color, ls=ls, label=label)
        ax2.semilogy(x, np.maximum(np.abs(y), 1.0e-16), color=color, ls=ls)
    for event in sym["runs"]["switch_step1"]["switches"]:
        x = event["t_s"] / sym["scenario"]["period_s"]
        ax.axvline(x, color="0.55", lw=0.7, ls=":")
        ax2.axvline(x, color="0.55", lw=0.7, ls=":")
    ax.set_ylabel(r"$(E-E_0)/|E_0|$")
    ax.legend(fontsize=7.3, loc="best")
    ax2.set_ylabel(r"absolute $|\Delta E/E_0|$")
    ax2.set_xlabel("Time [revolutions]")
    fig.savefig(FIGS / "fig_symplectic_switch.pdf")
    plt.close(fig)

    sol = payload["solver_independence"]
    full = payload["expanded_force_control"]
    lro = payload["lro_geometry_control"]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.75), sharey=True)
    categories = ["fixed 138", "schedule down", "schedule up"]
    x = np.arange(3)
    solver_values = {
        "DOP853": [
            sol["DOP853"]["errors_vs_same_solver_N300"]["fixed_N138"]["pos_rms_m"],
            sol["DOP853"]["errors_vs_same_solver_N300"]["schedule_down"]["pos_rms_m"],
            np.nan,
        ],
        "Radau": [
            sol["Radau"]["errors_vs_same_solver_N300"]["fixed_N138"]["pos_rms_m"],
            sol["Radau"]["errors_vs_same_solver_N300"]["schedule_down"]["pos_rms_m"],
            np.nan,
        ],
    }
    width = 0.34
    solver_x = np.arange(2)
    axs[0].bar(solver_x - width / 2, solver_values["DOP853"][:2], width, color=C1, label="DOP853")
    axs[0].bar(solver_x + width / 2, solver_values["Radau"][:2], width, color=C2, hatch="//", label="Radau")
    axs[0].set_title("Solver control")
    axs[0].legend(fontsize=7.2)
    axs[0].set_xticks(solver_x, categories[:2], rotation=28, ha="right")

    grav = json.loads((METRICS / "r3_longarc_matrix.json").read_text(encoding="utf-8"))
    g_rows = {r["run"]: r for r in grav["cases"]["M_moonpa"]["rows"]}
    gravity_values = [
        g_rows["fixed_138"]["pos_rms_m"], g_rows["sched_down"]["pos_rms_m"],
        g_rows["sched_up"]["pos_rms_m"],
    ]
    full_values = [
        full["errors_vs_full_force_N300"]["fixed_N138"]["pos_rms_m"],
        full["errors_vs_full_force_N300"]["schedule_down"]["pos_rms_m"],
        full["errors_vs_full_force_N300"]["schedule_up"]["pos_rms_m"],
    ]
    axs[1].bar(x - width / 2, gravity_values, width, color=C3, label="gravity only")
    axs[1].bar(x + width / 2, full_values, width, color=C4, hatch="//", label="+ Earth/Sun/SRP")
    axs[1].set_title("Expanded-force control")
    axs[1].legend(fontsize=7.0)

    lro_keys = [
        "fixed_N138", "fixed_N194_empirical", "fixed_N219_calibrated",
        "schedule_down", "schedule_up",
    ]
    lro_labels = ["fixed 138", "fixed 194", "fixed 219", "sched. down", "sched. up"]
    lro_values = [
        lro["errors_vs_same_force_N300"][key]["pos_rms_m"] for key in lro_keys
    ]
    lro_x = np.arange(len(lro_keys))
    bars = axs[2].bar(lro_x, lro_values, 0.68, color=[C1, C3, C4, C5, C2])
    axs[2].bar_label(bars, labels=[f"{v:.0f}" for v in lro_values], padding=2, fontsize=6.2)
    axs[2].set_title("LRO-like geometry")
    for axx in axs:
        axx.set_yscale("log")
        axx.grid(axis="x", visible=False)
    axs[1].set_xticks(x, categories, rotation=28, ha="right")
    axs[2].set_xticks(lro_x, lro_labels, rotation=32, ha="right")
    axs[0].set_ylabel("7-day RMS position error vs. same-case $N=300$ [m]")
    fig.savefig(FIGS / "fig_robustness_controls.pdf", pad_inches=0.08)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true", help="one-day solver-only timing pilot")
    parser.add_argument("--figures-only", action="store_true", help="regenerate figures from saved metrics")
    parser.add_argument("--lro-only", action="store_true", help="rerun corrected LRO-like control and regenerate figures")
    args_cli = parser.parse_args()
    if args_cli.figures_only:
        payload = json.loads((METRICS / "r4_robustness_controls.json").read_text(encoding="utf-8"))
        make_figures(payload)
        print("[written] robustness figures", flush=True)
        return 0
    model = load_model(300)
    kargs = kernel_args(model)
    warmup(model, kargs)
    if args_cli.pilot:
        pilot = solver_control(model, kargs, DAY)
        print(json.dumps(pilot, indent=2))
        return 0

    if args_cli.lro_only:
        payload = json.loads((METRICS / "r4_robustness_controls.json").read_text(encoding="utf-8"))
        duration = float(payload["lro_geometry_control"]["scenario"]["duration_s"])
        ephem = build_ephemeris(duration)
        payload["lro_geometry_control"] = lro_geometry_control(
            model, kargs, ephem, duration
        )
        (METRICS / "r4_robustness_controls.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
        make_figures(payload)
        print("[written] corrected LRO-like control and robustness figures", flush=True)
        return 0

    duration = 7.0 * DAY
    payload = {
        "schema": "r4_robustness_controls_v1",
        "solver_independence": solver_control(model, kargs, duration),
    }
    ephem = build_ephemeris(duration)
    payload["expanded_force_control"] = expanded_force_control(
        model, kargs, ephem, duration
    )
    payload["lro_geometry_control"] = lro_geometry_control(
        model, kargs, ephem, duration
    )
    payload["symplectic_control"] = symplectic_control()
    payload["source_notes"] = {
        "LRO_geometry": "Mesarch et al., AAS 23-238 / NASA NTRS 20230010945: 30 x 216 km, i approximately 90 deg, AOP approximately 270 deg.",
        "scope": "The four controls test numerical robustness and one published quasi-frozen design geometry; they do not establish a universal scheduling theorem.",
    }
    payload["repo_commit_sha"] = commit_sha()
    payload["repo_working_tree_clean"] = working_tree_clean()
    (METRICS / "r4_robustness_controls.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    make_figures(payload)
    print("[written] r4_robustness_controls.json and robustness figures", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
