"""R30: freeze four geometry strata, and the registration that governs them.

Designs A, B and C sample the whole factor box: perilune 30-150 km, apolune up
to 600 km, inclination 0-180 degrees. A population-level result on that box says
the interior member wins *on average over the box*. It does not say where in the
box it wins, and "where" is the first question a mission designer asks.

A stratum here is not a new experiment. It is the same scrambled Sobol
generator, the same map from the unit cube to an orbit, and the same protocol,
restricted to a sub-box of the cube:

    polar         inclination 80-100 deg
    equatorial    inclination 0-20 deg
    high_apolune  perilune 80-150 km with apolune 520-600 km
    frozen_like   inclination 84-90 deg, argument of perilune 85-95 deg

Only the sub-box changes. Every orbit is still produced by
rev10_sobol_confirmatory.orbit_from_u, which is the pinned map the three
coverage designs were built with, so a stratum point is a point the coverage
designs could have drawn -- just one the box makes likely instead of rare.

Two of these names deserve their exact meaning stated, because both could be
read as promising more than they deliver.

"high_apolune" is not a high-eccentricity orbit family. The frozen factor box
tops out at perilune 30 km with apolune 600 km, which is an eccentricity of
about 0.135: every orbit this paper samples is a near-circular low lunar orbit.
This stratum takes the widest arcs the box supports at a perilune high enough
not to confound the comparison with the low-perilune degree demand the paper
already attributes the effect to -- eccentricity 0.09-0.12 against a box-wide
0.01-0.135. It is the eccentricity contrast that exists here, and calling it a
high-eccentricity stratum would overstate it.

"frozen_like" is a geometry, not a solved frozen condition: near-polar with the
argument of perilune near 90 degrees, which is where lunar frozen orbits sit. No
stationarity condition is imposed, and the name carries the "_like" for that
reason.

Seeds continue the archive's arithmetic: A 20260723, B 20260724, C 20260725,
then 20260726-29 in the order the strata are listed above. Nothing is chosen by
inspection.

This script propagates nothing.

Usage:  python rev30_strata_freeze.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

import rev10_sobol_confirmatory as base

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

SEED_C = 20260725
PREREG_OUT = METRICS / "r30_preregistration.json"

# Each stratum maps the unit draw u -> a sub-box of the same unit cube, and the
# pinned orbit map then turns that into an orbit. The lambdas are written as
# fractions of the full factor range so the box is readable next to the map:
#   u0 -> hp = 30 + 120*u0 km,  u1 -> ha within [max(180,hp+30), 600] km,
#   u2 -> inclination = 180*u2 deg,  u3 -> argp = 360*u3 deg,  u4 -> longitude.
STRATA = {
    "polar": {
        "seed": SEED_C + 1, "design_key": "SP",
        "box": {"u2": (80.0 / 180.0, 100.0 / 180.0)},
        "means": "inclination 80-100 deg, everything else the full box",
    },
    "equatorial": {
        "seed": SEED_C + 2, "design_key": "SE",
        "box": {"u2": (0.0, 20.0 / 180.0)},
        "means": "inclination 0-20 deg, everything else the full box",
    },
    "high_apolune": {
        "seed": SEED_C + 3, "design_key": "SH",
        "box": {},
        "means": ("perilune 80-150 km with apolune 520-600 km: the widest arcs "
                  "the box supports, deliberately not low-perilune"),
    },
    "frozen_like": {
        "seed": SEED_C + 4, "design_key": "SF",
        "box": {"u2": (84.0 / 180.0, 90.0 / 180.0),
                "u3": (85.0 / 360.0, 95.0 / 360.0)},
        "means": ("near-polar with the argument of perilune near 90 deg; a "
                  "geometry, not a solved frozen condition"),
    },
    "low_perilune": {
        "seed": SEED_C + 5, "design_key": "SL",
        "box": {"u0": (0.0, 20.0 / 120.0)},
        "means": ("perilune 30-50 km, everything else the full box: the regime "
                  "the paper attributes the effect to"),
    },
}

OUTCOMES = {
    "K_holds_in_every_stratum": (
        "The interior member carries the realized-work tally at the declared "
        "budget in all four strata. The allocation result is then reported as "
        "holding across the geometries sampled, not only on average over the "
        "coverage box."),
    "L_geometry_dependent": (
        "The tally favours the interior member in some strata and the constant "
        "degree in others. The strata are then named individually, the "
        "population-level statement is qualified to the coverage design it was "
        "measured on, and no stratum is dropped or averaged into the others."),
    "M_stratum_undecided": (
        "A stratum resolves too few comparisons to produce a tally. It is "
        "reported as undecided with its resolved and unresolved counts, and it "
        "is not read as agreement with the others."),
}


def remap(row: np.ndarray, box: dict) -> np.ndarray:
    """Squeeze the unit draw into the stratum's sub-box, coordinate by
    coordinate. A coordinate the box does not mention is left alone."""
    out = np.array(row, dtype=float)
    for name, span in box.items():
        i = int(name[1:])
        if span is None:
            continue
        lo, hi = span
        out[i] = lo + (hi - lo) * out[i]
    return out


def high_ecc_row(row: np.ndarray) -> np.ndarray:
    """Perilune 80-150 km, apolune 520-600 km.

    The apolune coordinate has to be solved rather than scaled, because the map
    reads it relative to a lower limit that depends on the perilune it just
    produced: ha = ha_min + u1 * (600 - ha_min).
    """
    out = np.array(row, dtype=float)
    out[0] = (80.0 - 30.0) / 120.0 + out[0] * (150.0 - 80.0) / 120.0
    hp_km = 30.0 + 120.0 * out[0]
    ha_min = max(180.0, hp_km + 30.0)
    lo = (520.0 - ha_min) / (600.0 - ha_min)
    out[1] = lo + out[1] * (1.0 - lo)
    return out


def build(name: str, spec: dict, model, protocol_sha: str) -> dict:
    points = qmc.Sobol(d=5, scramble=True, seed=spec["seed"]).random_base2(m=6)
    family = f"sobol_{name}"
    orbits = []
    for index, row in enumerate(points):
        r = (high_ecc_row(row) if name == "high_apolune"
             else remap(row, spec["box"]))
        orbit = base.orbit_from_u(index, r, family, model)
        orbit["u_raw"] = [float(x) for x in row]
        orbits.append(orbit)
    counts = {k: 0 for k in ("prograde", "high_inclination", "retrograde")}
    for o in orbits:
        counts[base.inclination_regime(o["incl_deg"])] += 1
    payload = {
        "schema": "r30_stratum_design_v1",
        "family": family,
        "stratum": name,
        "design_key": spec["design_key"],
        "role": "geometry stratum of the frozen factor box",
        "seed": spec["seed"],
        "seed_rule": (f"the archive's seed sequence continued: A {SEED_C-2}, "
                      f"B {SEED_C-1}, C {SEED_C}, then one per stratum in the "
                      f"order they are declared"),
        "sample_count": len(orbits),
        "dimension": 5,
        "generator": "scrambled Sobol random_base2(m=6)",
        "sub_box": spec["means"],
        "orbit_map": ("rev10_sobol_confirmatory.orbit_from_u, unchanged; the "
                      "stratum changes only which part of the unit cube is "
                      "sampled"),
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
        "low_perilune_orbits": sum(1 for o in orbits if o["hp_km"] < 50.0),
        "orbits": orbits,
    }
    payload["design_sha256"] = base.object_hash(payload)
    return payload


def main() -> int:
    protocol_sha = base.object_hash(
        {"rule": "seed sequence continued past design C, one per stratum",
         "seed_c": SEED_C,
         "strata": {k: v["seed"] for k, v in STRATA.items()}})
    model = base.load_model(300)

    frozen = {}
    for name, spec in STRATA.items():
        design = build(name, spec, model, protocol_sha)
        out = METRICS / f"r30_{name}_design_frozen.json"
        out.write_text(json.dumps(design, indent=2), encoding="utf-8")
        frozen[name] = {"file": out.name,
                        "design_sha256": design["design_sha256"],
                        "seed": design["seed"],
                        "design_key": design["design_key"]}
        print(f"[r30] {out.name}  seed={design['seed']}  "
              f"sha={design['design_sha256'][:16]}")
        print(f"      hp {design['realized_perilune_km'][0]:.0f}-"
              f"{design['realized_perilune_km'][1]:.0f} km, "
              f"ha {design['realized_apolune_km'][0]:.0f}-"
              f"{design['realized_apolune_km'][1]:.0f} km, "
              f"i {design['realized_inclination_range_deg'][0]:.1f}-"
              f"{design['realized_inclination_range_deg'][1]:.1f} deg, "
              f"e {design['realized_eccentricity'][0]:.3f}-"
              f"{design['realized_eccentricity'][1]:.3f}, "
              f"hp<50km on {design['low_perilune_orbits']}")

    prereg = {
        "schema": "r30_preregistration_v1",
        "created_utc": base.utc_now(),
        "campaign": ("R30: four geometry strata of the frozen factor box, run "
                     "to locate the allocation result rather than to restate "
                     "it"),
        "question": ("designs A, B and C measure the allocation result averaged "
                     "over the whole factor box. These four strata ask where in "
                     "the box it holds."),
        "not_blind": ("the coverage-design results are known. The protections "
                      "are the arithmetic seed rule, the single-draw "
                      "commitment, and outcomes fixed before any stratum "
                      "propagates."),
        "single_draw_commitment": ("one draw per stratum. A stratum whose "
                                   "numbers disagree with the coverage designs "
                                   "is reported as it came out; no stratum is "
                                   "re-seeded, redefined, or dropped after its "
                                   "numbers are known, and no orbit is removed "
                                   "from a stratum after propagation."),
        "protocol": ("identical to designs A, B and C in every respect except "
                     "which part of the unit cube is sampled: the same orbit "
                     "map, truth-degree rule, tolerance levels, seven-day arc, "
                     "output grid, resolution rule, calibration and ladder"),
        "base_scope": (
            "one deliberate reduction, declared before any stratum propagates. "
            "The coverage designs carry six policies at two tolerance levels, "
            "because they also answer the degree-schedule questions of the "
            "convergence study. A stratum's base carries the two the budget "
            "ladder actually reads -- the truth and the critical-degree "
            "comparator -- and not the four schedule policies "
            "(schedule_empirical, schedule_up, schedule_down, fixed_work). "
            "This is a third of the propagation cost, which is what buys the "
            "fourth stratum. No stratum result is quoted from a schedule "
            "policy, and the empirical prepass is unchanged, so every quantity "
            "the ladder and the resolution rule use is present at full "
            "protocol fidelity."),
        "scope": ("strata locate the existing claim inside the box that was "
                  "already sampled. They do not extend it to a new body, arc "
                  "length, force model, or factor range, and a stratum is not "
                  "a replication of a coverage design: it is a sub-population "
                  "and is reported as one."),
        "budget": ("the declared budget beta = 1 is propagated on every "
                   "stratum first, because that is the budget the headline "
                   "claim is made at. Further budgets are propagated only if "
                   "every stratum already has beta = 1."),
        "order_and_conditionality": (
            "the strata run cheapest first, and low_perilune runs last because "
            "every one of its orbits carries the highest adopted truth degree "
            "and it is the most expensive population in the campaign. It is "
            "declared here, before any stratum has propagated a single orbit, "
            "precisely so that running it is a question of the clock and not a "
            "decision taken after the other four are known. If the clock does "
            "not reach it, it is reported as declared and not run; it is not "
            "quoted partially and its absence is not read as a null result."),
        "pooling": ("strata are never pooled with each other or with the "
                    "coverage designs into a single tally. Their sub-boxes "
                    "overlap the coverage box and each other, so a pooled "
                    "count would double-count the geometry it oversamples."),
        "strata": frozen,
        "sub_boxes": {k: v["means"] for k, v in STRATA.items()},
        "outcomes": OUTCOMES,
        "verdict_rule": ("per stratum, the R19 realized-work tally: interior if "
                         "resolved_interior_wins > resolved_fixed_wins, "
                         "constant if the reverse, split if equal, undecided if "
                         "nothing resolves. The panel verdict is K only if "
                         "every propagated stratum returns interior."),
        "completion_reporting": ("every stratum number carries its orbit count "
                                 "and its completion state; a stratum whose "
                                 "base or ladder is partial is reported as "
                                 "partial and is not read as agreement"),
        "prohibited": [
            "re-seeding a stratum",
            "redefining a sub-box after seeing its numbers",
            "dropping an orbit or a stratum after propagation",
            "pooling strata into one tally",
            "quoting a stratum whose base is incomplete",
        ],
    }
    prereg["preregistration_sha256"] = base.object_hash(prereg)
    PREREG_OUT.write_text(json.dumps(prereg, indent=2), encoding="utf-8")
    print(f"[r30] {PREREG_OUT.name} "
          f"sha256={prereg['preregistration_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
