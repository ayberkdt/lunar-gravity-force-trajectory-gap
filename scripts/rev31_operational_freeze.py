"""R31: freeze one population of operational elliptical lunar orbits.

Every orbit designs A, B and C sample is a near-circular low lunar orbit. The
frozen factor box runs from 30 x 180 km to 150 x 600 km, so its widest arc has
an eccentricity of about 0.135, and the R30 high_apolune stratum -- the most
eccentric sub-box the coverage designs contain -- only reaches 0.12. The paper
therefore has no evidence at all about the geometry a relay or a
communications-orbit designer actually flies.

That geometry exists and has been flown. Kaguya released two subsatellites into
elliptical lunar orbits: Ouna at 100 x 800 km and Okina at 100 x 2400 km, both
near-polar. This population brackets both of them -- perilune 80-120 km,
apolune 700-2500 km -- which spans eccentricities of roughly 0.15 to 0.40, three
times anything the coverage designs contain.

This is the one population in the campaign that leaves the frozen factor box,
and it is registered separately for exactly that reason. An R30 stratum is a
sub-box of a range that was already sampled and can only locate an existing
claim inside it. This extends the range, so whatever it finds is reported as an
extension of scope and is never folded into a population-level statement about
the coverage designs.

The orbit map is the pinned one in everything except the apolune range it is
allowed to draw from, which is the single thing that has to change for the box
to widen. The perilune-longitude solve, the RAAN construction, the truth-degree
rule and the initial-state construction are the archived functions, called
unchanged, so a point here differs from a coverage-design point in its factor
values and in nothing else.

Cost note, since it decides whether this runs at all: perilune stays above
80 km, so every adopted truth degree is 300 rather than 900, and the integrator
is step-capped at 60 s regardless of period. This is one of the cheap
populations despite being the widest.

Usage:  python rev31_operational_freeze.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SEED = 20260731              # the archive's sequence, continued past the strata
NAME = "operational_elliptical"
DESIGN_KEY = "OE"
FAMILY = "sobol_operational_elliptical"

HP_KM = (80.0, 120.0)        # brackets the 100 km perilune both relays flew
HA_KM = (700.0, 2500.0)      # brackets Ouna (800 km) and Okina (2400 km)

DESIGN_OUT = METRICS / f"r31_{NAME}_design_frozen.json"
PREREG_OUT = METRICS / "r31_preregistration.json"

FLOWN = [
    {"spacecraft": "Ouna (VRAD/Vstar), SELENE/Kaguya subsatellite",
     "perilune_km": 100, "apolune_km": 800, "note": "near-polar"},
    {"spacecraft": "Okina (Rstar), SELENE/Kaguya subsatellite",
     "perilune_km": 100, "apolune_km": 2400, "note": "near-polar"},
]

OUTCOMES = {
    "N_extends": (
        "The realized-work tally at the declared budget favours the interior "
        "member here as it does on the coverage designs. The allocation result "
        "is then reported as holding on a flown elliptical geometry outside the "
        "sampled box, stated as an extension of scope and kept separate from "
        "the population-level statements, which remain about the box."),
    "O_boundary": (
        "The tally favours the constant degree here. The coverage-design result "
        "is then reported as bounded by geometry: it holds on near-circular low "
        "lunar orbits and does not carry to the wide elliptical arcs a relay "
        "flies. This is a limitation of the paper's claim and is written as "
        "one, in the main text and not only in the supplement."),
    "P_undecided": (
        "Too few comparisons resolve to produce a tally. Reported as a geometry "
        "on which the measurement does not decide, with its resolved and "
        "unresolved counts."),
}


def orbit_from_u(index: int, row: np.ndarray, model) -> dict:
    """rev10_sobol_confirmatory.orbit_from_u with one range widened.

    Reproduced rather than imported because the archived map hardcodes an
    apolune ceiling of 600 km, which is the ceiling this population exists to
    cross. Everything else below is the archived body line for line, and the
    initial state comes from the archived constructor.
    """
    hp_km = HP_KM[0] + (HP_KM[1] - HP_KM[0]) * float(row[0])
    ha_min_km = max(HA_KM[0], hp_km + 30.0)
    ha_km = ha_min_km + float(row[1]) * (HA_KM[1] - ha_min_km)
    incl_deg = 180.0 * float(row[2])
    argp_deg = 360.0 * float(row[3])
    requested_lon = 360.0 * float(row[4])
    inclination = math.radians(incl_deg)
    argument = math.radians(argp_deg)
    longitude_offset = math.degrees(
        math.atan2(math.cos(inclination) * math.sin(argument),
                   math.cos(argument))
    )
    raan_deg = (requested_lon - longitude_offset) % 360.0
    orbit = {
        "name": f"{FAMILY}_{index:03d}",
        "family": FAMILY,
        "sobol_index": index,
        "u": [float(x) for x in row],
        "hp_km": hp_km,
        "ha_km": ha_km,
        "incl_deg": incl_deg,
        "argp_deg": argp_deg,
        "requested_perilune_lon_deg_bodyfixed_t0": requested_lon,
        "raan_deg": raan_deg,
        "nu0_deg": 0.0,
    }
    rp = model.r_ref + hp_km * 1000.0
    ra = model.r_ref + ha_km * 1000.0
    orbit["semimajor_axis_m"] = 0.5 * (rp + ra)
    orbit["eccentricity"] = (ra - rp) / (ra + rp)
    state = base.initial_state(model, orbit)
    actual_lon = math.degrees(math.atan2(state[1], state[0])) % 360.0
    actual_lat = math.degrees(math.asin(state[2] / np.linalg.norm(state[:3])))
    orbit["initial_state_si"] = [float(x) for x in state]
    orbit["reconstructed_perilune_lon_deg_bodyfixed_t0"] = actual_lon
    orbit["reconstructed_perilune_lat_deg_bodyfixed_t0"] = actual_lat
    orbit["longitude_roundtrip_error_deg"] = base.wrapped_delta_deg(
        actual_lon, requested_lon)
    orbit["truth_degree"] = 600 if hp_km < 50.0 else 300
    return orbit


def main() -> int:
    protocol_sha = base.object_hash(
        {"rule": "the archive's seed sequence, continued past the R30 strata",
         "seed": SEED, "hp_km": list(HP_KM), "ha_km": list(HA_KM)})
    model = base.load_model(300)
    points = qmc.Sobol(d=5, scramble=True, seed=SEED).random_base2(m=6)
    orbits = [orbit_from_u(i, row, model) for i, row in enumerate(points)]
    counts = {k: 0 for k in ("prograde", "high_inclination", "retrograde")}
    for o in orbits:
        counts[base.inclination_regime(o["incl_deg"])] += 1

    design = {
        "schema": "r31_operational_design_v1",
        "family": FAMILY,
        "population": NAME,
        "design_key": DESIGN_KEY,
        "role": ("an elliptical geometry outside the frozen factor box, "
                 "bracketing flown lunar relay orbits"),
        "seed": SEED,
        "seed_rule": ("the archive's seed sequence continued past the R30 "
                      "strata; arithmetic, not selected"),
        "sample_count": len(orbits),
        "dimension": 5,
        "generator": "scrambled Sobol random_base2(m=6)",
        "factor_box": {"hp_km": list(HP_KM), "ha_km": list(HA_KM),
                       "incl_deg": [0.0, 180.0], "argp_deg": [0.0, 360.0]},
        "outside_frozen_box": ("the coverage designs cap apolune at 600 km; "
                               "this population runs to 2500 km. Perilune and "
                               "inclination stay inside their sampled ranges."),
        "orbit_map": ("rev10_sobol_confirmatory.orbit_from_u reproduced with "
                      "the apolune range widened; the perilune-longitude "
                      "solve, RAAN construction, truth-degree rule and "
                      "initial-state constructor are the archived functions"),
        "flown_reference_orbits": FLOWN,
        "protocol_sha256": protocol_sha,
        "propagation_status": "frozen_pending_base_generation",
        "realized_perilune_km": [min(o["hp_km"] for o in orbits),
                                 max(o["hp_km"] for o in orbits)],
        "realized_apolune_km": [min(o["ha_km"] for o in orbits),
                                max(o["ha_km"] for o in orbits)],
        "realized_inclination_range_deg": [min(o["incl_deg"] for o in orbits),
                                           max(o["incl_deg"] for o in orbits)],
        "realized_eccentricity": [min(o["eccentricity"] for o in orbits),
                                  max(o["eccentricity"] for o in orbits)],
        "inclination_regime_counts": counts,
        "orbits": orbits,
    }
    design["design_sha256"] = base.object_hash(design)
    DESIGN_OUT.write_text(json.dumps(design, indent=2), encoding="utf-8")

    prereg = {
        "schema": "r31_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R31: one 64-orbit population of operational elliptical "
                     "lunar orbits, propagated to test whether the allocation "
                     "result survives outside the near-circular box the paper "
                     "samples"),
        "question": ("every orbit in designs A, B and C is near-circular "
                     "(eccentricity below 0.14). Relay and communications "
                     "orbits are not. Does the declared-budget result hold on "
                     "the geometry that is actually flown for those roles?"),
        "why_separate_from_R30": (
            "an R30 stratum is a sub-box of a range already sampled and can "
            "only locate an existing claim inside it; the R30 registration says "
            "so explicitly. This population widens the range, so it is a scope "
            "extension and is registered, reported and bounded as one. It is "
            "never pooled with the strata or with the coverage designs."),
        "not_blind": ("the coverage-design results are known. The protections "
                      "are the arithmetic seed, the single-draw commitment, and "
                      "outcomes fixed before propagation."),
        "single_draw_commitment": ("one draw. If the numbers disagree with the "
                                   "coverage designs, the disagreement is the "
                                   "result; the seed is not changed, the box is "
                                   "not redrawn, and no orbit is dropped after "
                                   "its numbers are known."),
        "protocol": ("identical to the R30 strata: same truth-degree rule, "
                     "tolerance levels, seven-day arc, output grid, resolution "
                     "rule, calibration and ladder, and the same two-policy "
                     "base scope, which is declared in r30_preregistration.json "
                     "and inherited here"),
        "budget": "the declared budget beta = 1 first; further budgets only after",
        "population": NAME,
        "design": {"seed": SEED, "design_key": DESIGN_KEY,
                   "design_sha256": design["design_sha256"],
                   "file": DESIGN_OUT.name,
                   "factor_box": design["factor_box"],
                   "flown_reference_orbits": FLOWN},
        "outcomes": OUTCOMES,
        "verdict_rule": ("the R19 realized-work tally, as everywhere else: "
                         "interior if resolved_interior_wins exceeds "
                         "resolved_fixed_wins, constant if the reverse, split "
                         "if equal, undecided if nothing resolves"),
        "reporting_commitment": (
            "the phrase 'outside the sampled factor box' travels with every "
            "number from this population, in the table caption and in the "
            "sentence that quotes it. A favourable result here does not widen "
            "any abstract claim; an unfavourable one does bound it, and that "
            "bound goes in the main text rather than only in the supplement."),
        "prohibited": [
            "re-seeding or redrawing the box",
            "pooling this population with the strata or the coverage designs",
            "quoting it without its out-of-box qualifier",
            "widening an abstract or conclusion claim on the strength of it",
        ],
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")

    print(f"[r31] {DESIGN_OUT.name}  seed={SEED}  "
          f"sha={design['design_sha256'][:16]}")
    print(f"      hp {design['realized_perilune_km'][0]:.0f}-"
          f"{design['realized_perilune_km'][1]:.0f} km, "
          f"ha {design['realized_apolune_km'][0]:.0f}-"
          f"{design['realized_apolune_km'][1]:.0f} km, "
          f"i {design['realized_inclination_range_deg'][0]:.1f}-"
          f"{design['realized_inclination_range_deg'][1]:.1f} deg, "
          f"e {design['realized_eccentricity'][0]:.3f}-"
          f"{design['realized_eccentricity'][1]:.3f}")
    print(f"      regimes {counts}, all adopted truth degree 300 "
          f"(perilune min {design['realized_perilune_km'][0]:.0f} km)")
    print(f"[r31] {PREREG_OUT.name} "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
