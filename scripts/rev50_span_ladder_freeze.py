"""R50: freeze a paired radial-span ladder, and the registration that governs it.

The paper names radial span as the clearest geometric discriminator it has, and
then says, honestly, that it never measured it as one factor. The evidence is a
single sub-box that restricts the span (SH) and one step out to an elliptical
population (OE) that moves perilune, eccentricity, dwell and degree demand
together. A referee reading that reads a direction, not a dependence.

This population measures the dependence. Sixteen orbit identities are drawn once
by the pinned map -- perilune, inclination, argument of perilune and the
body-fixed perilune longitude -- and each identity is then flown at four
apolunes:

    ha = 300, 600, 1200, 2400 km

Everything else about an identity is held bit-identical across its four members:
the same perilune, the same inclination, the same argument of perilune, the same
requested perilune longitude, the same true anomaly at epoch. The only quantity
that changes is the apolune, so the only geometric quantity that changes is the
radial span. Each identity is therefore its own control and the comparison is
paired rather than between populations.

The four levels are chosen to cross the boundary the paper's two results sit on
either side of. At hp = 80-120 km the coverage box allows apolunes up to 600 km,
so levels 300 and 600 are *inside* the sampled factor box and are the box's
narrow and widest arcs; 1200 and 2400 are *outside* it, and 2400 with a 100 km
perilune is the geometry Kaguya's Okina subsatellite flew and the upper end of
the R31 operational population. One ladder therefore runs from the geometry
where the paper reports the constant degree winning to the geometry where it
reports the radial endpoint winning, along a single controlled axis.

Perilune stays in 80-120 km, the same band R31 drew from, so every adopted truth
degree is 300 and this is one of the cheap populations: its parent measured a
72-minute base and 35-46 minute ladders on the same sixty-four orbits.

Block B is a second sixteen identities at the same four levels, drawn with the
next seed. It doubles the panel behind every level and is run only if the clock
allows. It does not change, add or move a level.

This script propagates nothing.

Usage:  python rev50_span_ladder_freeze.py
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

SEED_A = 20260732            # the archive's sequence, continued past R31's 20260731
SEED_B = 20260733

HP_KM = (80.0, 120.0)        # the R31 perilune band: every truth degree is 300
APOLUNE_LEVELS_KM = (300.0, 600.0, 1200.0, 2400.0)
IN_BOX_LEVELS_KM = (300.0, 600.0)      # the coverage box caps apolune at 600 km
IDENTITIES = 16

BLOCKS = {
    "span_ladder_a": {"seed": SEED_A, "design_key": "RS1", "role": "primary"},
    "span_ladder_b": {"seed": SEED_B, "design_key": "RS2",
                      "role": "second panel, conditional on the clock"},
}

PREREG_OUT = METRICS / "r50_preregistration.json"

FLOWN = [
    {"spacecraft": "Okina (Rstar), SELENE/Kaguya subsatellite",
     "perilune_km": 100, "apolune_km": 2400,
     "note": "the top level of this ladder is the geometry it flew"},
]

OUTCOMES = {
    "T_span_dependence": (
        "The per-level tally moves with radial span in one direction, and at "
        "every budget that decides. The paper then states a controlled "
        "dependence: holding perilune, inclination, argument of perilune and "
        "perilune longitude fixed, increasing the radial span moves the budget "
        "crossing. The discussion sentence that currently reads the span as a "
        "direction is replaced by the measured dependence, and the OE result "
        "stops being the only wide-span evidence."),
    "U_threshold": (
        "The verdict does not move level by level but flips once, between two "
        "adjacent levels, and stays flipped. Reported as a threshold in radial "
        "span rather than as a gradient, with the two levels it sits between "
        "named. The paper does not describe a threshold as a dependence."),
    "V_no_dependence": (
        "The verdict is the same at every level at a given budget, within the "
        "resolved counts. The radial-span reading is then withdrawn from the "
        "discussion and the conclusion: the difference between the coverage "
        "designs and the elliptical population is attributed to something the "
        "ladder holds fixed -- perilune, dwell or degree demand -- and that "
        "correction goes in the main text, not only in the supplement. This "
        "outcome contradicts a sentence the paper currently prints, and it is "
        "written before any orbit of this population propagates."),
    "W_undecided": (
        "Too few comparisons resolve per level to produce level tallies. "
        "Reported as undecided with its resolved and unresolved counts per "
        "level, and read as no evidence either way rather than as agreement "
        "with the existing sentence."),
}


def identity_rows(seed: int) -> np.ndarray:
    """Sixteen draws of the pinned five-column unit cube.

    Column u1 is drawn and then not used: it is the apolune coordinate, and the
    apolune is prescribed by the level rather than sampled. It is kept in the
    record so the draw is the archived generator called unchanged.
    """
    return qmc.Sobol(d=5, scramble=True, seed=seed).random_base2(m=4)


def orbit_from_identity(index: int, identity: int, row: np.ndarray,
                        ha_km: float, family: str, model) -> dict:
    """rev10_sobol_confirmatory.orbit_from_u with the apolune prescribed.

    Reproduced rather than imported for the same reason R31 reproduced it: the
    archived map draws the apolune from a coordinate and caps it at 600 km, and
    a prescribed apolune above that cap is the one thing this population needs.
    Every other line is the archived body, and the perilune-longitude solve, the
    RAAN construction, the truth-degree rule and the initial-state constructor
    are the archived functions called unchanged.
    """
    hp_km = HP_KM[0] + (HP_KM[1] - HP_KM[0]) * float(row[0])
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
        "name": f"{family}_{index:03d}",
        "family": family,
        "sobol_index": index,
        "u": [float(x) for x in row],
        "identity_index": identity,
        "apolune_level_km": ha_km,
        "apolune_level_inside_factor_box": ha_km in IN_BOX_LEVELS_KM,
        "hp_km": hp_km,
        "ha_km": ha_km,
        "radial_span_km": ha_km - hp_km,
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


def build(name: str, spec: dict, model, protocol_sha: str) -> dict:
    family = f"sobol_{name}"
    rows = identity_rows(spec["seed"])
    orbits = []
    index = 0
    for identity, row in enumerate(rows):
        for ha_km in APOLUNE_LEVELS_KM:
            orbits.append(orbit_from_identity(index, identity, row, ha_km,
                                              family, model))
            index += 1

    # the pairing is a claim about the numbers, so it is checked here rather
    # than asserted in prose: the four members of an identity must agree to the
    # bit on everything except the apolune.
    held = ("hp_km", "incl_deg", "argp_deg",
            "requested_perilune_lon_deg_bodyfixed_t0", "raan_deg", "nu0_deg")
    for identity in range(len(rows)):
        group = [o for o in orbits if o["identity_index"] == identity]
        if len(group) != len(APOLUNE_LEVELS_KM):
            raise SystemExit(f"identity {identity} has {len(group)} members")
        for field in held:
            if len({o[field] for o in group}) != 1:
                raise SystemExit(f"identity {identity} does not hold {field} "
                                 f"fixed across its apolune levels")

    counts = {k: 0 for k in ("prograde", "high_inclination", "retrograde")}
    for o in orbits:
        counts[base.inclination_regime(o["incl_deg"])] += 1

    payload = {
        "schema": "r50_span_ladder_design_v1",
        "family": family,
        "population": name,
        "design_key": spec["design_key"],
        "role": ("a paired apolune ladder: one draw of orbit identities, each "
                 "flown at four apolunes, so radial span varies alone"),
        "block_role": spec["role"],
        "seed": spec["seed"],
        "seed_rule": ("the archive's seed sequence continued past R31's "
                      f"20260731: A {SEED_A}, B {SEED_B}. Arithmetic, not "
                      "selected."),
        "sample_count": len(orbits),
        "identities": len(rows),
        "apolune_levels_km": list(APOLUNE_LEVELS_KM),
        "levels_inside_factor_box_km": list(IN_BOX_LEVELS_KM),
        "dimension": 5,
        "generator": "scrambled Sobol random_base2(m=4), one draw per identity",
        "unused_coordinate": ("u1, the apolune coordinate of the pinned map. "
                              "The apolune is prescribed by the level, so the "
                              "coordinate is drawn and recorded but not read."),
        "factor_box": {"hp_km": list(HP_KM),
                       "ha_km": [min(APOLUNE_LEVELS_KM),
                                 max(APOLUNE_LEVELS_KM)],
                       "incl_deg": [0.0, 180.0], "argp_deg": [0.0, 360.0]},
        "held_fixed_within_identity": list(held),
        "outside_frozen_box": ("levels 300 and 600 km are inside the coverage "
                               "box, which caps apolune at 600 km; levels 1200 "
                               "and 2400 km are outside it and carry that "
                               "qualifier wherever they are quoted"),
        "orbit_map": ("rev10_sobol_confirmatory.orbit_from_u reproduced with "
                      "the apolune prescribed instead of drawn; the "
                      "perilune-longitude solve, RAAN construction, "
                      "truth-degree rule and initial-state constructor are the "
                      "archived functions"),
        "flown_reference_orbits": FLOWN,
        "protocol_sha256": protocol_sha,
        "propagation_status": "frozen_pending_base_generation",
        "realized_perilune_km": [min(o["hp_km"] for o in orbits),
                                 max(o["hp_km"] for o in orbits)],
        "realized_apolune_km": [min(o["ha_km"] for o in orbits),
                                max(o["ha_km"] for o in orbits)],
        "realized_radial_span_km": [min(o["radial_span_km"] for o in orbits),
                                    max(o["radial_span_km"] for o in orbits)],
        "realized_inclination_range_deg": [min(o["incl_deg"] for o in orbits),
                                           max(o["incl_deg"] for o in orbits)],
        "realized_eccentricity": [min(o["eccentricity"] for o in orbits),
                                  max(o["eccentricity"] for o in orbits)],
        "inclination_regime_counts": counts,
        "adopted_truth_degrees": sorted({o["truth_degree"] for o in orbits}),
        "orbits": orbits,
    }
    payload["design_sha256"] = base.object_hash(payload)
    return payload


def main() -> int:
    protocol_sha = base.object_hash(
        {"rule": ("a paired apolune ladder on one draw of identities; the "
                  "archive's seed sequence continued past R31"),
         "seeds": {k: v["seed"] for k, v in BLOCKS.items()},
         "hp_km": list(HP_KM),
         "apolune_levels_km": list(APOLUNE_LEVELS_KM),
         "identities": IDENTITIES})
    model = base.load_model(300)

    frozen = {}
    for name, spec in BLOCKS.items():
        design = build(name, spec, model, protocol_sha)
        out = METRICS / f"r50_{name}_design_frozen.json"
        out.write_text(json.dumps(design, indent=2), encoding="utf-8")
        frozen[name] = {"file": out.name,
                        "design_sha256": design["design_sha256"],
                        "seed": design["seed"],
                        "design_key": design["design_key"]}
        print(f"[r50] {out.name}  seed={design['seed']}  "
              f"sha={design['design_sha256'][:16]}")
        print(f"      {design['identities']} identities x "
              f"{len(APOLUNE_LEVELS_KM)} apolunes = {design['sample_count']} "
              f"orbits")
        print(f"      hp {design['realized_perilune_km'][0]:.1f}-"
              f"{design['realized_perilune_km'][1]:.1f} km, "
              f"span {design['realized_radial_span_km'][0]:.0f}-"
              f"{design['realized_radial_span_km'][1]:.0f} km, "
              f"i {design['realized_inclination_range_deg'][0]:.1f}-"
              f"{design['realized_inclination_range_deg'][1]:.1f} deg, "
              f"e {design['realized_eccentricity'][0]:.3f}-"
              f"{design['realized_eccentricity'][1]:.3f}")
        print(f"      regimes {design['inclination_regime_counts']}, "
              f"adopted truth degrees {design['adopted_truth_degrees']}")

    prereg = {
        "schema": "r50_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R50: a paired radial-span ladder. Sixteen orbit "
                     "identities, each flown at four apolunes, propagated "
                     "through the budget ladder at four budgets, to measure "
                     "whether the budget crossing depends on radial span when "
                     "radial span is the only thing that varies."),
        "question": ("the paper reports that the winning allocation shifts with "
                     "budget and radial-span regime. The span half of that "
                     "sentence rests on one sub-box that restricts the span and "
                     "one step out to an elliptical population that moves "
                     "perilune, eccentricity, dwell and degree demand together. "
                     "Does the crossing move with radial span alone?"),
        "why_this_design": (
            "a sub-box drawn narrower would compare two populations and inherit "
            "the same confound in weaker form. Holding an identity fixed and "
            "changing only its apolune makes every orbit its own control, so a "
            "level-to-level difference cannot be a difference in perilune, "
            "inclination, argument of perilune or epoch geometry."),
        "levels": {
            "apolune_km": list(APOLUNE_LEVELS_KM),
            "inside_factor_box_km": list(IN_BOX_LEVELS_KM),
            "why_these": ("300 and 600 km are the narrowest and the widest arcs "
                          "the coverage box allows at this perilune; 1200 and "
                          "2400 km leave it, and 2400 km at a 100 km perilune "
                          "is the Okina geometry and the upper end of the R31 "
                          "population. The ladder therefore runs from the "
                          "geometry where the paper reports the constant degree "
                          "winning to the geometry where it reports the radial "
                          "endpoint winning."),
            "out_of_box_qualifier": ("levels 1200 and 2400 km are a scope "
                                     "extension exactly as R31 is, and every "
                                     "number from them carries that qualifier"),
        },
        "budgets": [1.00, 0.75, 0.62, 0.50],
        "budget_order": ("beta = 1.00 first, because it is the budget the "
                         "geometry sentence is written at, then 0.75, 0.62, "
                         "0.50 in that order. A budget the clock does not reach "
                         "is reported as not run."),
        "budget_note": ("0.62 is the amendment budget declared in "
                        "r28_calibration_amendment.json and computed by the "
                        "standard calibration; here it is pre-registered before "
                        "any orbit of this population propagates rather than "
                        "added afterwards."),
        "not_blind": ("the coverage-design, strata and R31 results are known. "
                      "The protections are the arithmetic seed, the single-draw "
                      "commitment, the level grid fixed here, and outcomes "
                      "fixed before propagation."),
        "single_draw_commitment": ("one draw of identities per block. If the "
                                   "numbers disagree with the sentence the "
                                   "paper prints, the disagreement is the "
                                   "result; no seed is changed, no level is "
                                   "added, moved or dropped after its numbers "
                                   "are known, and no identity is removed."),
        "protocol": ("identical to the R30 strata and R31: the same "
                     "truth-degree rule, tolerance levels, seven-day arc, "
                     "output grid, resolution rule, calibration and ladder, and "
                     "the same two-policy base scope declared in "
                     "r30_preregistration.json and inherited here"),
        "reference_degree": ("the adopted rule, which gives degree 300 at every "
                             "perilune in this band. The cap-lifted control R38 "
                             "ran for R31 is not run here and is not implied; "
                             "the cap binds at perilune, which the ladder holds "
                             "fixed, so it is common to all four levels."),
        "readouts": {
            "primary": ("per level and budget, the realized-work tally of the "
                        "budget-calibrated radial endpoint against its "
                        "work-matched constant degree. This is the comparison "
                        "the discussion's geometry sentence is about."),
            "secondary": ("per level and budget, the same tally for the "
                          "interior member of the span family, which is the "
                          "comparison the budget-axis sentence is about"),
            "paired": ("per identity, the level at which each verdict changes, "
                       "so the dependence is read within identities and not "
                       "only across level medians"),
        },
        "verdict_rule": ("per level, the R19 realized-work tally: radial or "
                         "interior if resolved wins exceed the constant "
                         "degree's, constant if the reverse, split if equal, "
                         "undecided if nothing resolves. A dependence is "
                         "claimed only if the level tallies order consistently "
                         "at a budget that decides."),
        "blocks": {
            "span_ladder_a": "the primary panel, sixteen identities",
            "span_ladder_b": ("a second sixteen identities at the same four "
                              "levels, run only if the clock allows after block "
                              "A carries every budget. It doubles the panel "
                              "behind each level. It is declared here so that "
                              "running it is a question of the clock and not a "
                              "decision taken after block A is known; if the "
                              "clock does not reach it, it is reported as "
                              "declared and not run."),
        },
        "pooling": ("blocks A and B are pooled with each other only level by "
                    "level, because they are the same design at the same "
                    "levels. Neither is pooled with the coverage designs, the "
                    "strata or R31, whose boxes they overlap."),
        "outcomes": OUTCOMES,
        "reporting_commitment": (
            "the outcome is written into the discussion whichever way it comes "
            "out. V_no_dependence contradicts a sentence the paper currently "
            "prints; if it is what returns, that sentence is corrected in the "
            "main text and the correction is not deferred to the supplement."),
        "completion_reporting": ("every number carries its level, its orbit "
                                 "count and its completion state; a level whose "
                                 "ladder is partial is reported as partial and "
                                 "is not read as agreement"),
        "prohibited": [
            "re-seeding a block or redrawing an identity",
            "adding, moving or dropping an apolune level after its numbers are "
            "known",
            "quoting the 1200 or 2400 km levels without the out-of-box "
            "qualifier",
            "pooling this population with the coverage designs, the strata or "
            "R31",
            "reading a level tally from a partial ladder",
        ],
        "strata": frozen,
        "sub_boxes": {k: (f"perilune {HP_KM[0]:.0f}-{HP_KM[1]:.0f} km, apolune "
                          f"prescribed at "
                          f"{', '.join(f'{h:.0f}' for h in APOLUNE_LEVELS_KM)} "
                          f"km, everything else the full range")
                      for k in BLOCKS},
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")
    print(f"[r50] {PREREG_OUT.name} "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
