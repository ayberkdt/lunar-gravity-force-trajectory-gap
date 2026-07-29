"""Independent, non-propagating correctness audit for the R10 campaign.

The audit deliberately does not call the R10 propagation or post-processing
helpers. It reconstructs designs and metrics from frozen JSON/NPZ artifacts,
checks all formal sidecar/hash/config contracts, and evaluates a small set of
instantaneous analytic and metamorphic gravity identities.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import qmc


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
PYTHON_CODES = ROOT / "python_codes"
LUNARIS = Path(r"D:\Masaustu\LUNAR_SIMULATION")
OUTPUT = METRICS / "r10_code_correctness_audit.json"

sys.path.insert(0, str(PYTHON_CODES))
sys.path.insert(0, str(LUNARIS / "src"))

import rev3_common as common  # noqa: E402
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    sh_accel_fixed_numba,
    sh_potential_accel_fixed,
)
from rev9_potential_blend_longarc import BlendRhs  # noqa: E402


def load_json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value: object) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def local_path(value: str) -> Path:
    return ROOT / Path(value.replace("\\", "/"))


checks: list[dict] = []


def check(identifier: str, passed: bool, details: dict, severity: str = "critical") -> None:
    checks.append(
        {
            "id": identifier,
            "severity_if_failed": severity,
            "passed": bool(passed),
            "details": details,
        }
    )


raw_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}


def load_raw(path: Path) -> tuple[np.ndarray, np.ndarray]:
    key = path.resolve().as_posix()
    if key not in raw_cache:
        with np.load(path) as data:
            raw_cache[key] = (
                np.asarray(data["t_s"], dtype=float),
                np.asarray(data["state_si"], dtype=float),
            )
    return raw_cache[key]


def independent_metrics(solution: np.ndarray, truth: np.ndarray) -> dict:
    if solution.shape != truth.shape or solution.shape[0] != 6:
        raise ValueError(f"incompatible state arrays: {solution.shape}, {truth.shape}")
    delta_position = (solution[:3] - truth[:3]).T
    delta_velocity = (solution[3:] - truth[3:]).T
    position_norm = np.sqrt(np.einsum("ij,ij->i", delta_position, delta_position))
    velocity_norm = np.sqrt(np.einsum("ij,ij->i", delta_velocity, delta_velocity))

    position = truth[:3].T
    velocity = truth[3:].T
    radial = position / np.sqrt(
        np.einsum("ij,ij->i", position, position)
    )[:, None]
    cross_track = np.cross(position, velocity)
    cross_track /= np.sqrt(
        np.einsum("ij,ij->i", cross_track, cross_track)
    )[:, None]
    in_track = np.cross(cross_track, radial)

    ric = {
        "radial": np.einsum("ij,ij->i", delta_position, radial),
        "in_track": np.einsum("ij,ij->i", delta_position, in_track),
        "cross_track": np.einsum("ij,ij->i", delta_position, cross_track),
    }
    return {
        "pos_rms_m": float(np.sqrt(np.mean(position_norm * position_norm))),
        "pos_max_m": float(np.max(position_norm)),
        "pos_final_m": float(position_norm[-1]),
        "vel_rms_m_s": float(np.sqrt(np.mean(velocity_norm * velocity_norm))),
        "vel_max_m_s": float(np.max(velocity_norm)),
        "vel_final_m_s": float(velocity_norm[-1]),
        "ric_rms_m": {
            name: float(np.sqrt(np.mean(values * values)))
            for name, values in ric.items()
        },
        "ric_max_m": {
            name: float(np.max(np.abs(values))) for name, values in ric.items()
        },
        "ric_final_m": {name: float(values[-1]) for name, values in ric.items()},
    }


def numeric_leaves(value: object, prefix: str = "") -> dict[str, float]:
    result: dict[str, float] = {}
    if isinstance(value, dict):
        for key, child in value.items():
            result.update(numeric_leaves(child, f"{prefix}.{key}" if prefix else key))
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        result[prefix] = float(value)
    return result


def metric_deviation(recomputed: dict, recorded: dict) -> tuple[float, float, list[str]]:
    a = numeric_leaves(recomputed)
    b = numeric_leaves(recorded)
    common_keys = sorted(set(a) & set(b))
    maximum_absolute = 0.0
    maximum_relative = 0.0
    failures: list[str] = []
    for key in common_keys:
        absolute = abs(a[key] - b[key])
        relative = absolute / max(abs(b[key]), 1.0e-30)
        maximum_absolute = max(maximum_absolute, absolute)
        maximum_relative = max(maximum_relative, relative)
        if not math.isclose(a[key], b[key], rel_tol=5.0e-12, abs_tol=1.0e-9):
            failures.append(key)
    return maximum_absolute, maximum_relative, failures


def rotation_z(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotation_x(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def independent_state(orbit: dict, mu: float) -> np.ndarray:
    semimajor = float(orbit["semimajor_axis_m"])
    eccentricity = float(orbit["eccentricity"])
    parameter = semimajor * (1.0 - eccentricity * eccentricity)
    radius_pf = np.array([parameter / (1.0 + eccentricity), 0.0, 0.0])
    velocity_pf = np.array(
        [0.0, math.sqrt(mu / parameter) * (1.0 + eccentricity), 0.0]
    )
    transform = (
        rotation_z(math.radians(orbit["raan_deg"]))
        @ rotation_x(math.radians(orbit["incl_deg"]))
        @ rotation_z(math.radians(orbit["argp_deg"]))
    )
    return np.concatenate((transform @ radius_pf, transform @ velocity_pf))


baseline = load_json("metrics/r10_sobolA_baseline.json")
corrected = load_json("metrics/r10_sobolA_baseline_truth_corrected.json")
truth_audit = load_json("metrics/r10_sobolA_truth_audit.json")
truth_audit_extended = load_json("metrics/r10_sobolA_truth_audit_extended.json")
convergence = load_json("metrics/r10_sobolA_convergence.json")
blend = load_json("metrics/r10_blend_lro_convergence.json")
aggregate = load_json("metrics/r10_aggregate_summary.json")
descriptives = load_json("metrics/r10_manuscript_descriptives.json")
corrected_extended = load_json(
    "metrics/r10_sobolA_baseline_truth_corrected_extended.json"
)
extended_aggregate = load_json("metrics/r10_truth_audit_extended_aggregate.json")
design_a = load_json("metrics/r10_sobolA_design.json")
design_b = load_json("metrics/r10_sobolB_design_frozen.json")


# 1. Every formal R10 trajectory is tied to a valid sidecar, config, and hash.
formal_directories = (
    METRICS / "r10_cases" / "baseline",
    METRICS / "r10_cases" / "truth_audit",
    METRICS / "r10_cases" / "truth_audit_extended",
    METRICS / "r10_cases" / "convergence",
    METRICS / "r10_cases" / "blend_lro_convergence",
)
sidecars = sorted(
    path
    for directory in formal_directories
    for path in directory.rglob("*.json")
    if "smoke" not in path.name.lower()
    and ".failure." not in path.name.lower()
    and ".invalid." not in path.name.lower()
)
sidecar_failures: list[str] = []
kernel_hashes: set[str] = set()
gravity_hashes: set[str] = set()
for sidecar in sidecars:
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    raw = local_path(meta["raw_path"])
    config = meta["config"]
    source = config.get("source", {})
    if source.get("kernel_sha256"):
        kernel_hashes.add(source["kernel_sha256"])
    if source.get("gravity_sha256"):
        gravity_hashes.add(source["gravity_sha256"])
    if meta.get("status") != "complete":
        sidecar_failures.append(f"{sidecar.relative_to(ROOT)}: status")
        continue
    if not raw.exists() or sha256(raw) != meta.get("raw_sha256"):
        sidecar_failures.append(f"{sidecar.relative_to(ROOT)}: raw hash")
        continue
    if canonical_hash(config) != meta.get("config_sha256"):
        sidecar_failures.append(f"{sidecar.relative_to(ROOT)}: config hash")
    times, states = load_raw(raw)
    if (
        states.shape != (6, len(times))
        or not np.all(np.isfinite(states))
        or not np.all(np.diff(times) > 0.0)
        or not math.isclose(times[0], 0.0, abs_tol=0.0)
        or not math.isclose(times[-1], config["duration_s"], abs_tol=1.0e-9)
        or not np.array_equal(states[:, 0], np.asarray(config["initial_state_si"]))
    ):
        sidecar_failures.append(f"{sidecar.relative_to(ROOT)}: raw contract")
check(
    "formal_trajectory_sidecar_hash_and_shape_contract",
    len(sidecars) == 536 and not sidecar_failures,
    {
        "formal_sidecars": len(sidecars),
        "expected": 536,
        "failures": sidecar_failures,
        "unique_kernel_hashes": sorted(kernel_hashes),
        "unique_gravity_hashes": sorted(gravity_hashes),
    },
)


# 2. Regenerate both Sobol designs and reconstruct every Cartesian state.
model_300 = common.load_model(300)
design_failures: list[str] = []
maximum_state_error = 0.0
maximum_longitude_error = 0.0
for design, seed in ((design_a, 20260723), (design_b, 20260724)):
    regenerated = qmc.Sobol(d=5, scramble=True, seed=seed).random_base2(m=6)
    stored = np.asarray([orbit["u"] for orbit in design["orbits"]], dtype=float)
    if not np.array_equal(regenerated, stored):
        design_failures.append(f"{design['family']}: Sobol coordinates")
    for orbit, coordinates in zip(design["orbits"], regenerated):
        hp = 30.0 + 120.0 * coordinates[0]
        ha_min = max(180.0, hp + 30.0)
        ha = ha_min + coordinates[1] * (600.0 - ha_min)
        inclination = 180.0 * coordinates[2]
        argument = 360.0 * coordinates[3]
        longitude = 360.0 * coordinates[4]
        offset = math.degrees(
            math.atan2(
                math.cos(math.radians(inclination))
                * math.sin(math.radians(argument)),
                math.cos(math.radians(argument)),
            )
        )
        raan = (longitude - offset) % 360.0
        expected = {
            "hp_km": hp,
            "ha_km": ha,
            "incl_deg": inclination,
            "argp_deg": argument,
            "requested_perilune_lon_deg_bodyfixed_t0": longitude,
            "raan_deg": raan,
        }
        for key, value in expected.items():
            if not math.isclose(float(orbit[key]), float(value), abs_tol=1.0e-12):
                design_failures.append(f"{orbit['name']}: {key}")
        state = independent_state(orbit, model_300.mu)
        maximum_state_error = max(
            maximum_state_error,
            float(np.max(np.abs(state - np.asarray(orbit["initial_state_si"])))),
        )
        actual_longitude = math.degrees(math.atan2(state[1], state[0])) % 360.0
        wrapped = (actual_longitude - longitude + 180.0) % 360.0 - 180.0
        maximum_longitude_error = max(maximum_longitude_error, abs(wrapped))
check(
    "independent_sobol_and_initial_state_reconstruction",
    not design_failures
    and maximum_state_error <= 1.0e-9
    and maximum_longitude_error <= 1.0e-10,
    {
        "orbits_reconstructed": 128,
        "maximum_state_component_error_si": maximum_state_error,
        "maximum_perilune_longitude_error_deg": maximum_longitude_error,
        "failures": design_failures,
    },
)


# 3. Recompute degree telemetry, work matching, and critical degrees.
models = {300: model_300, 600: common.load_model(600)}
power = {
    degree: np.asarray(
        [
            np.sum(
                model.c_coeffs[n, : n + 1] ** 2
                + model.s_coeffs[n, : n + 1] ** 2
            )
            for n in range(model.max_degree + 1)
        ],
        dtype=float,
    )
    for degree, model in models.items()
}


def critical_degree(model, degree_power: np.ndarray, altitude_m: float) -> int:
    indices = np.arange(len(degree_power), dtype=float)
    radius = model.r_ref + altitude_m
    attenuation = np.exp(indices * math.log(model.r_ref / radius))
    sigma = (
        np.sqrt((indices + 1.0) * (2.0 * indices + 1.0))
        * attenuation
        * np.sqrt(degree_power)
    )
    squared = sigma * sigma
    total = float(np.sum(squared[2:]))
    budget = 1.0e-6 * total
    tail = total
    selected = len(squared) - 1
    for n in range(2, len(squared)):
        if tail <= budget:
            selected = n - 1
            break
        tail -= squared[n]
    return min(250, max(60, selected))


degree_failures: list[str] = []
maximum_mean_n2_error = 0.0
for row in baseline["rows"]:
    telemetry = row["policies"]["schedule_empirical"]["telemetry"]
    counts = {int(key): int(value) for key, value in telemetry["degree_counts"].items()}
    count_sum = sum(counts.values())
    mean_n2 = sum(degree * degree * count for degree, count in counts.items()) / count_sum
    maximum_mean_n2_error = max(
        maximum_mean_n2_error, abs(mean_n2 - telemetry["mean_degree_sq"])
    )
    expected_work = int(round(math.sqrt(mean_n2)))
    expected_critical = critical_degree(
        models[row["truth_degree"]],
        power[row["truth_degree"]],
        row["design_point"]["hp_km"] * 1000.0,
    )
    if count_sum != telemetry["n_rhs"]:
        degree_failures.append(f"{row['name']}: degree-count sum")
    if expected_work != row["n_work"]:
        degree_failures.append(f"{row['name']}: N_work")
    if expected_critical != row["n_critical"]:
        degree_failures.append(f"{row['name']}: N_critical")
check(
    "independent_degree_telemetry_and_comparator_reconstruction",
    not degree_failures,
    {
        "orbits": 64,
        "maximum_mean_N2_absolute_error": maximum_mean_n2_error,
        "failures": degree_failures,
    },
)


# 4. Independently recompute baseline and truth-corrected policy metrics.
maximum_metric_absolute = 0.0
maximum_metric_relative = 0.0
metric_failures: list[str] = []
corrected_by_index = {row["sobol_index"]: row for row in corrected["rows"]}
audit_by_index = {row["sobol_index"]: row for row in truth_audit["rows"]}
recomputed_rho_work: list[float] = []
recomputed_rho_critical: list[float] = []
recomputed_savings: list[float] = []
recomputed_corrected_errors_by_index: dict[int, dict[str, dict]] = {}
recomputed_ratios_by_index: dict[int, tuple[float, float]] = {}

for row in baseline["rows"]:
    index = row["sobol_index"]
    truth_meta = json.loads(
        (
            METRICS
            / "r10_cases"
            / "baseline"
            / f"sobolA_{index:03d}"
            / "truth_baseline.json"
        ).read_text(encoding="utf-8")
    )
    _, baseline_truth = load_raw(local_path(truth_meta["raw_path"]))
    adopted = corrected_by_index[index]["adopted_truth_degree"]
    if adopted == 900:
        audit_meta = json.loads(
            (
                METRICS
                / "r10_cases"
                / "truth_audit"
                / f"sobolA_{index:03d}"
                / "truth_N900_baseline.json"
            ).read_text(encoding="utf-8")
        )
        _, adopted_truth = load_raw(local_path(audit_meta["raw_path"]))
    else:
        adopted_truth = baseline_truth

    recomputed_errors: dict[str, dict] = {}
    for policy in (
        "schedule_empirical",
        "fixed_work",
        "fixed_critical",
        "schedule_up",
        "schedule_down",
    ):
        sidecar = json.loads(
            (
                METRICS
                / "r10_cases"
                / "baseline"
                / f"sobolA_{index:03d}"
                / f"{policy}_baseline.json"
            ).read_text(encoding="utf-8")
        )
        _, states = load_raw(local_path(sidecar["raw_path"]))
        baseline_stats = independent_metrics(states, baseline_truth)
        absolute, relative, failures = metric_deviation(
            baseline_stats, sidecar["error_against_truth"]
        )
        maximum_metric_absolute = max(maximum_metric_absolute, absolute)
        maximum_metric_relative = max(maximum_metric_relative, relative)
        metric_failures.extend(
            f"{row['name']}/{policy}/baseline/{key}" for key in failures
        )
        recomputed_errors[policy] = independent_metrics(states, adopted_truth)
        absolute, relative, failures = metric_deviation(
            recomputed_errors[policy],
            corrected_by_index[index]["policy_errors"][policy],
        )
        maximum_metric_absolute = max(maximum_metric_absolute, absolute)
        maximum_metric_relative = max(maximum_metric_relative, relative)
        metric_failures.extend(
            f"{row['name']}/{policy}/corrected/{key}" for key in failures
        )

    empirical = recomputed_errors["schedule_empirical"]["pos_rms_m"]
    rho_work = recomputed_errors["fixed_work"]["pos_rms_m"] / empirical
    rho_critical = recomputed_errors["fixed_critical"]["pos_rms_m"] / empirical
    recomputed_corrected_errors_by_index[index] = recomputed_errors
    recomputed_ratios_by_index[index] = (rho_work, rho_critical)
    recomputed_rho_work.append(rho_work)
    recomputed_rho_critical.append(rho_critical)
    recomputed_savings.append(
        1.0
        - row["policies"]["schedule_empirical"]["telemetry"]["gravity_kernel_ns"]
        / row["policies"]["fixed_critical"]["telemetry"]["gravity_kernel_ns"]
    )
    if not math.isclose(rho_work, corrected_by_index[index]["rho_work"], rel_tol=5e-12):
        metric_failures.append(f"{row['name']}: rho_work")
    if not math.isclose(
        rho_critical, corrected_by_index[index]["rho_crit"], rel_tol=5e-12
    ):
        metric_failures.append(f"{row['name']}: rho_crit")

summary_failures: list[str] = []
summary_expectations = {
    "rho_work_median": (
        float(np.median(recomputed_rho_work)),
        aggregate["baseline"]["truth_corrected_rho_work"]["median"],
    ),
    "rho_work_p10": (
        float(np.percentile(recomputed_rho_work, 10)),
        aggregate["baseline"]["truth_corrected_rho_work"]["p10"],
    ),
    "rho_work_p90": (
        float(np.percentile(recomputed_rho_work, 90)),
        aggregate["baseline"]["truth_corrected_rho_work"]["p90"],
    ),
    "rho_critical_median": (
        float(np.median(recomputed_rho_critical)),
        aggregate["baseline"]["truth_corrected_rho_crit"]["median"],
    ),
    "saving_median": (
        float(np.median(recomputed_savings)),
        descriptives["median_gravity_time_saving_vs_critical"],
    ),
}
for name, (actual, archived) in summary_expectations.items():
    if not math.isclose(actual, archived, rel_tol=5.0e-12, abs_tol=1.0e-12):
        summary_failures.append(name)
if int(np.sum(np.asarray(recomputed_rho_work) > 1.0)) != 7:
    summary_failures.append("raw work win count")
if int(np.sum(np.asarray(recomputed_rho_critical) > 1.0)) != 2:
    summary_failures.append("raw critical win count")
check(
    "independent_raw_metric_and_aggregate_recomputation",
    not metric_failures and not summary_failures,
    {
        "policy_truth_metric_pairs": 640,
        "maximum_metric_absolute_error": maximum_metric_absolute,
        "maximum_metric_relative_error": maximum_metric_relative,
        "metric_failures": metric_failures,
        "summary_failures": summary_failures,
        "recomputed": {
            "raw_work_wins": int(np.sum(np.asarray(recomputed_rho_work) > 1.0)),
            "raw_critical_wins": int(
                np.sum(np.asarray(recomputed_rho_critical) > 1.0)
            ),
            "rho_work_median": float(np.median(recomputed_rho_work)),
            "rho_critical_median": float(np.median(recomputed_rho_critical)),
            "gravity_saving_median": float(np.median(recomputed_savings)),
        },
    },
)


# 5. Recompute the four-case extended truth audit and merged population.
extended_failures: list[str] = []
extended_by_index = {
    row["sobol_index"]: row for row in truth_audit_extended["rows"]
}
for index, row in extended_by_index.items():
    baseline_truth = load_raw(
        METRICS
        / "r10_raw"
        / "baseline"
        / f"sobolA_{index:03d}"
        / "truth_baseline.npz"
    )[1]
    n900_truth = load_raw(
        METRICS
        / "r10_raw"
        / "truth_audit_extended"
        / f"sobolA_{index:03d}"
        / "truth_N900_baseline.npz"
    )[1]
    truth_difference = independent_metrics(baseline_truth, n900_truth)
    _, _, failures = metric_deviation(truth_difference, row["N600_to_N900"])
    extended_failures.extend(f"{index}/truth/{key}" for key in failures)
    policy_errors: dict[str, dict] = {}
    for policy in (
        "schedule_empirical",
        "fixed_work",
        "fixed_critical",
        "schedule_up",
        "schedule_down",
    ):
        states = load_raw(
            METRICS
            / "r10_raw"
            / "baseline"
            / f"sobolA_{index:03d}"
            / f"{policy}_baseline.npz"
        )[1]
        policy_errors[policy] = independent_metrics(states, n900_truth)
        _, _, failures = metric_deviation(
            policy_errors[policy], row["policy_errors_against_N900"][policy]
        )
        extended_failures.extend(f"{index}/{policy}/{key}" for key in failures)
    threshold = min(
        5.0,
        0.05 * min(error["pos_rms_m"] for error in policy_errors.values()),
    )
    passed = truth_difference["pos_rms_m"] < threshold
    if (
        not math.isclose(
            threshold, row["acceptance_threshold_m"], rel_tol=5e-12, abs_tol=1e-9
        )
        or passed != row["passes"]
        or row["adopted_truth_degree"] != (600 if passed else 900)
    ):
        extended_failures.append(f"{index}/acceptance")
    empirical = policy_errors["schedule_empirical"]["pos_rms_m"]
    recomputed_corrected_errors_by_index[index] = policy_errors
    recomputed_ratios_by_index[index] = (
        policy_errors["fixed_work"]["pos_rms_m"] / empirical,
        policy_errors["fixed_critical"]["pos_rms_m"] / empirical,
    )

merged_rows = {
    row["sobol_index"]: row for row in corrected_extended["rows"]
}
for index, row in merged_rows.items():
    rho_work, rho_critical = recomputed_ratios_by_index[index]
    if not math.isclose(rho_work, row["rho_work"], rel_tol=5e-12, abs_tol=1e-12):
        extended_failures.append(f"{index}/merged/rho_work")
    if not math.isclose(
        rho_critical, row["rho_crit"], rel_tol=5e-12, abs_tol=1e-12
    ):
        extended_failures.append(f"{index}/merged/rho_crit")
    for policy, metrics in recomputed_corrected_errors_by_index[index].items():
        _, _, failures = metric_deviation(metrics, row["policy_errors"][policy])
        extended_failures.extend(
            f"{index}/merged/{policy}/{key}" for key in failures
        )

merged_work = np.asarray(
    [recomputed_ratios_by_index[index][0] for index in sorted(merged_rows)]
)
merged_critical = np.asarray(
    [recomputed_ratios_by_index[index][1] for index in sorted(merged_rows)]
)
archived_merged = extended_aggregate["aggregate_with_extension"]
merged_expectations = (
    (np.median(merged_work), archived_merged["rho_work"]["median"], "work median"),
    (
        np.percentile(merged_work, 10),
        archived_merged["rho_work"]["p10"],
        "work p10",
    ),
    (
        np.percentile(merged_work, 90),
        archived_merged["rho_work"]["p90"],
        "work p90",
    ),
    (
        np.median(merged_critical),
        archived_merged["rho_crit"]["median"],
        "critical median",
    ),
    (
        np.percentile(merged_critical, 10),
        archived_merged["rho_crit"]["p10"],
        "critical p10",
    ),
    (
        np.percentile(merged_critical, 90),
        archived_merged["rho_crit"]["p90"],
        "critical p90",
    ),
)
for actual, archived, label in merged_expectations:
    if not math.isclose(float(actual), float(archived), rel_tol=5e-12, abs_tol=1e-12):
        extended_failures.append(f"aggregate/{label}")
if int(np.sum(merged_work > 1.0)) != archived_merged["rho_work"]["raw_schedule_wins"]:
    extended_failures.append("aggregate/work wins")
if (
    int(np.sum(merged_critical > 1.0))
    != archived_merged["rho_crit"]["raw_schedule_wins"]
):
    extended_failures.append("aggregate/critical wins")
if (
    extended_aggregate["new_contested_indices"]
    or extended_aggregate["requires_additional_convergence"]
):
    extended_failures.append("unexpected new contested case")
check(
    "independent_extended_truth_audit_and_merged_population_recomputation",
    not extended_failures,
    {
        "extended_indices": sorted(extended_by_index),
        "extended_adopted_N900": [
            index
            for index, row in extended_by_index.items()
            if row["adopted_truth_degree"] == 900
        ],
        "merged_raw_work_wins": int(np.sum(merged_work > 1.0)),
        "merged_raw_critical_wins": int(np.sum(merged_critical > 1.0)),
        "merged_rho_work_median": float(np.median(merged_work)),
        "merged_rho_critical_median": float(np.median(merged_critical)),
        "new_contested_indices": extended_aggregate["new_contested_indices"],
        "failures": extended_failures,
    },
)


# 6. Recompute every selective-convergence envelope and decision from raw data.
convergence_failures: list[str] = []
maximum_convergence_absolute = 0.0
recomputed_counts = {
    "fixed_work": {"schedule": 0, "fixed": 0, "unresolved": 0},
    "fixed_critical": {"schedule": 0, "fixed": 0, "unresolved": 0},
}
for row in convergence["rows"]:
    index = row["sobol_index"]

    def convergence_states(policy: str, level: str) -> np.ndarray:
        path = (
            METRICS
            / "r10_raw"
            / "convergence"
            / f"sobolA_{index:03d}"
            / f"{policy}_{level}.npz"
        )
        return load_raw(path)[1]

    truth_tight = convergence_states("truth", "tight")
    truth_tighter = convergence_states("truth", "tighter")
    truth_self = independent_metrics(truth_tight, truth_tighter)["pos_rms_m"]
    maximum_convergence_absolute = max(
        maximum_convergence_absolute,
        abs(truth_self - row["truth_self_difference_rms_m"]),
    )
    policy_values: dict[str, dict] = {}
    for policy in ("schedule_empirical", "fixed_work", "fixed_critical"):
        tight = convergence_states(policy, "tight")
        tighter = convergence_states(policy, "tighter")
        tight_error = independent_metrics(tight, truth_tight)
        tighter_error = independent_metrics(tighter, truth_tighter)
        self_difference = independent_metrics(tight, tighter)["pos_rms_m"]
        envelope = self_difference + truth_self
        recorded = row["policies"][policy]
        for actual, archived, label in (
            (
                tight_error,
                recorded["errors_against_same_tolerance_truth"]["tight"],
                "tight",
            ),
            (
                tighter_error,
                recorded["errors_against_same_tolerance_truth"]["tighter"],
                "tighter",
            ),
        ):
            absolute, _, failures = metric_deviation(actual, archived)
            maximum_convergence_absolute = max(maximum_convergence_absolute, absolute)
            convergence_failures.extend(
                f"{index}/{policy}/{label}/{key}" for key in failures
            )
        if not math.isclose(
            self_difference,
            recorded["self_difference_rms_m"],
            rel_tol=5e-12,
            abs_tol=1e-9,
        ):
            convergence_failures.append(f"{index}/{policy}/self")
        if not math.isclose(
            envelope,
            recorded["truth_inclusive_envelope_m"],
            rel_tol=5e-12,
            abs_tol=1e-9,
        ):
            convergence_failures.append(f"{index}/{policy}/envelope")
        policy_values[policy] = {
            "tight_error": tight_error["pos_rms_m"],
            "envelope": envelope,
        }

    schedule = policy_values["schedule_empirical"]
    for comparator in ("fixed_work", "fixed_critical"):
        fixed = policy_values[comparator]
        gap = abs(schedule["tight_error"] - fixed["tight_error"])
        threshold = schedule["envelope"] + fixed["envelope"]
        resolved = gap > threshold
        if not resolved:
            outcome = "unresolved"
        elif schedule["tight_error"] < fixed["tight_error"]:
            outcome = "schedule"
        else:
            outcome = "fixed"
        recomputed_counts[comparator][outcome] += 1
        recorded = row["comparisons"][comparator]
        recorded_outcome = (
            "unresolved"
            if not recorded["resolved"]
            else (
                "schedule"
                if recorded["winner_if_resolved"] == "schedule_empirical"
                else "fixed"
            )
        )
        if (
            not math.isclose(
                gap,
                recorded["absolute_error_difference_m"],
                rel_tol=5e-12,
                abs_tol=1e-9,
            )
            or not math.isclose(
                threshold,
                recorded["resolution_threshold_m"],
                rel_tol=5e-12,
                abs_tol=1e-9,
            )
            or outcome != recorded_outcome
        ):
            convergence_failures.append(f"{index}/{comparator}/decision")
check(
    "independent_selective_convergence_recomputation",
    not convergence_failures
    and recomputed_counts["fixed_work"]
    == {"schedule": 7, "fixed": 9, "unresolved": 1}
    and recomputed_counts["fixed_critical"]
    == {"schedule": 0, "fixed": 17, "unresolved": 0},
    {
        "orbits": 17,
        "comparisons": 34,
        "maximum_absolute_metric_error": maximum_convergence_absolute,
        "counts": recomputed_counts,
        "failures": convergence_failures,
    },
)


# 7. Recompute the complete LRO-like blend convergence summary from raw data.
blend_failures: list[str] = []


def blend_states(policy: str, level: str) -> np.ndarray:
    return load_raw(
        METRICS / "r10_raw" / "blend_lro_convergence" / f"{policy}_{level}.npz"
    )[1]


truth_baseline = blend_states("truth_N600", "baseline")
truth_tighter = blend_states("truth_N600", "tighter")
blend_truth_self = independent_metrics(truth_baseline, truth_tighter)["pos_rms_m"]
blend_recomputed: dict[str, dict] = {}
for policy in ("fixed_N120", "corrected_blend"):
    baseline_state = blend_states(policy, "baseline")
    tighter_state = blend_states(policy, "tighter")
    baseline_error = independent_metrics(baseline_state, truth_baseline)
    tighter_error = independent_metrics(tighter_state, truth_tighter)
    self_difference = independent_metrics(baseline_state, tighter_state)["pos_rms_m"]
    envelope = self_difference + blend_truth_self
    blend_recomputed[policy] = {
        "baseline": baseline_error,
        "tighter": tighter_error,
        "self": self_difference,
        "envelope": envelope,
    }
    recorded = blend["summary"]["policies"][policy]
    for actual, archived, label in (
        (
            baseline_error,
            recorded["error_against_same_tolerance_truth"]["baseline"],
            "baseline",
        ),
        (
            tighter_error,
            recorded["error_against_same_tolerance_truth"]["tighter"],
            "tighter",
        ),
    ):
        _, _, failures = metric_deviation(actual, archived)
        blend_failures.extend(f"{policy}/{label}/{key}" for key in failures)
    if not math.isclose(
        self_difference,
        recorded["self_difference_rms_m"],
        rel_tol=5e-12,
        abs_tol=1e-9,
    ):
        blend_failures.append(f"{policy}/self")
    if not math.isclose(
        envelope,
        recorded["truth_inclusive_envelope_m"],
        rel_tol=5e-12,
        abs_tol=1e-9,
    ):
        blend_failures.append(f"{policy}/envelope")
blend_gap = abs(
    blend_recomputed["corrected_blend"]["baseline"]["pos_rms_m"]
    - blend_recomputed["fixed_N120"]["baseline"]["pos_rms_m"]
)
blend_threshold = (
    blend_recomputed["corrected_blend"]["envelope"]
    + blend_recomputed["fixed_N120"]["envelope"]
)
blend_resolved = blend_gap > blend_threshold
if (
    not math.isclose(
        blend_truth_self,
        blend["summary"]["truth_self_difference_rms_m"],
        rel_tol=5e-12,
        abs_tol=1e-9,
    )
    or not math.isclose(
        blend_gap,
        blend["summary"]["comparison"]["absolute_baseline_error_difference_m"],
        rel_tol=5e-12,
        abs_tol=1e-9,
    )
    or not math.isclose(
        blend_threshold,
        blend["summary"]["comparison"]["resolution_threshold_m"],
        rel_tol=5e-12,
        abs_tol=1e-9,
    )
    or blend_resolved != blend["summary"]["comparison"]["resolved"]
):
    blend_failures.append("comparison")
check(
    "independent_blend_convergence_recomputation",
    not blend_failures and not blend_resolved,
    {
        "truth_self_difference_m": blend_truth_self,
        "baseline_gap_m": blend_gap,
        "resolution_threshold_m": blend_threshold,
        "resolved": blend_resolved,
        "failures": blend_failures,
    },
)


# 8. Instantaneous analytic and metamorphic gravity checks.
model_600 = models[600]
args_600 = common.kernel_args(model_600)
directions = np.asarray(
    [[1.0, 0.2, -0.1], [0.3, -0.8, 0.5], [-0.6, 0.4, 0.7]], dtype=float
)
directions /= np.linalg.norm(directions, axis=1)[:, None]

central_relative_errors: list[float] = []
kernel_potential_relative_errors: list[float] = []
potential_gradient_relative_errors: list[float] = []
for direction, altitude in zip(directions, (30.0e3, 100.0e3, 300.0e3)):
    position = direction * (model_600.r_ref + altitude)
    acceleration_n0 = np.asarray(
        sh_accel_fixed_numba(*position, 0, *args_600), dtype=float
    )
    point_mass = -model_600.mu * position / np.linalg.norm(position) ** 3
    central_relative_errors.append(
        float(np.linalg.norm(acceleration_n0 - point_mass) / np.linalg.norm(point_mass))
    )
    for degree in (30, 120):
        _, acceleration_potential = sh_potential_accel_fixed(
            position.reshape(1, 3),
            model_600.c_coeffs,
            model_600.s_coeffs,
            model_600.mu,
            model_600.r_ref,
            degree,
            -1,
        )
        acceleration_kernel = np.asarray(
            sh_accel_fixed_numba(*position, degree, *args_600), dtype=float
        )
        kernel_potential_relative_errors.append(
            float(
                np.linalg.norm(acceleration_kernel - acceleration_potential[0])
                / np.linalg.norm(acceleration_kernel)
            )
        )
        if degree == 120:
            gradient = np.zeros(3)
            step = 3.0
            for axis in range(3):
                plus = position.copy()
                minus = position.copy()
                plus[axis] += step
                minus[axis] -= step
                potential_plus, _ = sh_potential_accel_fixed(
                    plus.reshape(1, 3),
                    model_600.c_coeffs,
                    model_600.s_coeffs,
                    model_600.mu,
                    model_600.r_ref,
                    degree,
                    -1,
                )
                potential_minus, _ = sh_potential_accel_fixed(
                    minus.reshape(1, 3),
                    model_600.c_coeffs,
                    model_600.s_coeffs,
                    model_600.mu,
                    model_600.r_ref,
                    degree,
                    -1,
                )
                gradient[axis] = (potential_plus[0] - potential_minus[0]) / (
                    2.0 * step
                )
            potential_gradient_relative_errors.append(
                float(
                    np.linalg.norm(gradient - acceleration_kernel)
                    / np.linalg.norm(acceleration_kernel)
                )
            )
check(
    "instantaneous_point_mass_kernel_and_potential_identities",
    max(central_relative_errors) < 5.0e-15
    and max(kernel_potential_relative_errors) < 5.0e-12
    and max(potential_gradient_relative_errors) < 1.0e-7,
    {
        "maximum_N0_point_mass_relative_error": max(central_relative_errors),
        "maximum_fixed_kernel_vs_potential_acceleration_relative_error": max(
            kernel_potential_relative_errors
        ),
        "maximum_finite_difference_potential_gradient_relative_error": max(
            potential_gradient_relative_errors
        ),
    },
)


def independent_blended_potential(position: np.ndarray) -> float:
    radius = float(np.linalg.norm(position))
    altitude = radius - model_600.r_ref
    coordinate = min(1.0, max(0.0, (200.0e3 - altitude) / 150.0e3))
    weight = coordinate * coordinate * (3.0 - 2.0 * coordinate)
    low, _ = sh_potential_accel_fixed(
        position.reshape(1, 3),
        model_600.c_coeffs,
        model_600.s_coeffs,
        model_600.mu,
        model_600.r_ref,
        30,
        -1,
    )
    high, _ = sh_potential_accel_fixed(
        position.reshape(1, 3),
        model_600.c_coeffs,
        model_600.s_coeffs,
        model_600.mu,
        model_600.r_ref,
        120,
        -1,
    )
    return float((1.0 - weight) * low[0] + weight * high[0])


blend_rhs = BlendRhs(model_600, args_600, "blend_potential_corrected")
corrected_blend_gradient_errors: list[float] = []
for direction, altitude in zip(directions, (60.0e3, 125.0e3, 190.0e3)):
    position = direction * (model_600.r_ref + altitude)
    state = np.concatenate((position, np.zeros(3)))
    acceleration = np.asarray(blend_rhs(0.0, state)[3:], dtype=float)
    gradient = np.zeros(3)
    step = 3.0
    for axis in range(3):
        plus = position.copy()
        minus = position.copy()
        plus[axis] += step
        minus[axis] -= step
        gradient[axis] = (
            independent_blended_potential(plus)
            - independent_blended_potential(minus)
        ) / (2.0 * step)
    corrected_blend_gradient_errors.append(
        float(np.linalg.norm(acceleration - gradient) / np.linalg.norm(acceleration))
    )
check(
    "corrected_blend_is_gradient_of_blended_potential",
    max(corrected_blend_gradient_errors) < 1.0e-7,
    {
        "transition_altitudes_tested_km": [60.0, 125.0, 190.0],
        "maximum_relative_error": max(corrected_blend_gradient_errors),
    },
)


# 9. Frame rotation and RIC orthonormality checks.
frame_errors: list[float] = []
rhs = common.Rhs(model_600, lambda t, h: 120, args_600)
body_position = directions[1] * (model_600.r_ref + 100.0e3)
body_acceleration = np.asarray(
    sh_accel_fixed_numba(*body_position, 120, *args_600), dtype=float
)
for epoch in (0.0, 0.25 * common.DAY, 3.75 * common.DAY):
    angle = common.OMEGA_MOON * epoch
    body_to_inertial = rotation_z(angle)
    inertial_position = body_to_inertial @ body_position
    state = np.concatenate((inertial_position, np.zeros(3)))
    rhs_acceleration = np.asarray(rhs(epoch, state)[3:], dtype=float)
    expected = body_to_inertial @ body_acceleration
    frame_errors.append(
        float(np.linalg.norm(rhs_acceleration - expected) / np.linalg.norm(expected))
    )

sample_truth = load_raw(
    METRICS / "r10_raw" / "baseline" / "sobolA_000" / "truth_baseline.npz"
)[1]
positions = sample_truth[:3].T
velocities = sample_truth[3:].T
radial = positions / np.linalg.norm(positions, axis=1)[:, None]
cross_track = np.cross(positions, velocities)
cross_track /= np.linalg.norm(cross_track, axis=1)[:, None]
in_track = np.cross(cross_track, radial)
orthogonality = np.stack(
    (
        np.einsum("ij,ij->i", radial, in_track),
        np.einsum("ij,ij->i", radial, cross_track),
        np.einsum("ij,ij->i", in_track, cross_track),
        np.linalg.norm(radial, axis=1) - 1.0,
        np.linalg.norm(in_track, axis=1) - 1.0,
        np.linalg.norm(cross_track, axis=1) - 1.0,
    )
)
maximum_orthonormality_error = float(np.max(np.abs(orthogonality)))
check(
    "body_inertial_rotation_and_RIC_orthonormality",
    max(frame_errors) < 5.0e-14 and maximum_orthonormality_error < 5.0e-15,
    {
        "maximum_frame_covariance_relative_error": max(frame_errors),
        "maximum_RIC_orthonormality_error": maximum_orthonormality_error,
    },
)


# 10. Carry-forward external evidence, restricted to the unchanged kernel.
kernel_xval = load_json("metrics/r3_kernel_xval.json")
tudat = load_json("metrics/external_tudat_evidence_matrix_2026_07_19.json")
shtools_maximum = max(
    row["accel_rel_max"] for row in kernel_xval["shtools_rows"]
)
tudat_position_maximum = max(row["sh_position_max_m"] for row in tudat["runs"])
tudat_acceleration_maximum = max(
    row["sh_acceleration_relative_max"] for row in tudat["runs"]
)
check(
    "archived_external_kernel_and_trajectory_cross_validation",
    tudat["status"] == "PASS"
    and all(row["all_checks_pass"] for row in tudat["runs"])
    and shtools_maximum < 5.0e-13
    and tudat_position_maximum < 0.001,
    {
        "scope_note": (
            "The audit separately confirmed by git diff that "
            "spherical_harmonics.py is unchanged between the R3 validation "
            "commit and the R10 release commit."
        ),
        "SHTOOLS_maximum_acceleration_relative_error": shtools_maximum,
        "Tudat_run_count": len(tudat["runs"]),
        "Tudat_maximum_position_difference_m": tudat_position_maximum,
        "Tudat_maximum_acceleration_relative_error": tudat_acceleration_maximum,
    },
)


critical_failures = [
    item["id"]
    for item in checks
    if not item["passed"] and item["severity_if_failed"] == "critical"
]
advisory_failures = [
    item["id"]
    for item in checks
    if not item["passed"] and item["severity_if_failed"] != "critical"
]
payload = {
    "schema": "r10_code_correctness_audit_v1",
    "audit_type": "non-propagating independent recomputation and metamorphic audit",
    "new_trajectory_propagations": 0,
    "formal_trajectory_artifacts_checked": len(sidecars),
    "checks_passed": sum(item["passed"] for item in checks),
    "checks_total": len(checks),
    "critical_failures": critical_failures,
    "advisory_failures": advisory_failures,
    "overall_status": "PASS" if not critical_failures else "FAIL",
    "checks": checks,
    "limitations": [
        "Instantaneous metamorphic checks exercise the production gravity kernel; "
        "they are not a second full propagator.",
        "The archived Tudat and SHTOOLS comparisons are independent implementations "
        "but use the same public gravity coefficients and physical contract.",
        "No finite validation suite can prove the absence of every software defect.",
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(
    f"[written] {OUTPUT.relative_to(ROOT)} "
    f"status={payload['overall_status']} "
    f"checks={payload['checks_passed']}/{payload['checks_total']}"
)
if critical_failures:
    print("[critical failures] " + ", ".join(critical_failures))
    raise SystemExit(1)
