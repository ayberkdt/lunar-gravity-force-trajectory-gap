"""Sidecars for the R41 reference trajectories.

Every other propagating campaign writes, beside each raw state array, a JSON
sidecar carrying the configuration that produced it and the digest of the array
itself. The manuscript's data-availability statement promises exactly that. The
R41 driver wrote the arrays and not the sidecars, which is a defect in the
driver rather than in the run; this repairs it and the driver is patched so a
re-run does it inline.

What the sidecars carry is derived, not invented:

  config      the parameters that fully determine the trajectory: the raised
              reference degree, the tolerance level, the arc, the output grid,
              the maximum step, and the orbit's frozen initial state. These are
              constants of the campaign and the design record, not
              reconstructions.

  status      derived from the array. A seven-day arc on a 120-s grid has 5041
              output epochs; a run that ended early would have fewer, and an
              impact would have terminated it. A file with the full grid ran to
              completion, and that is asserted only where the count agrees.

  telemetry   NOT retained. The driver computed step counts and call counts and
              discarded them. The sidecar says so rather than leaving a reader
              to assume they were never produced, and the field is null.

Usage:  python rev41_sidecars.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as bt
import rev41_reference_degree_control as r41

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
RAW = METRICS / "r41_raw"
CASES = METRICS / "r41_cases"

EXPECTED_EPOCHS = int(round(bt.DURATION / bt.OUTPUT_STEP)) + 1


def main() -> int:
    if not RAW.exists():
        print(f"[abort] {RAW} not found")
        return 2
    prereg = json.loads((METRICS / "r41_preregistration.json").read_text(
        encoding="utf-8"))
    declared = {(o["design"], int(o["sobol_index"])): o
                for o in prereg["selection_rule"]["orbits"]}

    src = {d: {int(r["sobol_index"]): r for r in
               json.loads(bt.DESIGNS[d]["rows"].read_text(encoding="utf-8"))["rows"]}
           for d in ("A", "B")}

    written, mismatched = 0, []
    for npz in sorted(RAW.rglob("reference_*.npz")):
        level = npz.stem.split("_", 1)[1]
        design = npz.parent.parent.name.split("_")[0]
        index = int(npz.parent.name.replace("sobolA_", ""))
        if (design, index) not in declared:
            mismatched.append(f"{design}{index:03d} not in the registration")
            continue
        row = src[design][index]
        with np.load(npz) as z:
            n_epochs = int(len(z["t_s"]))
            last_t = float(z["t_s"][-1])

        complete = (n_epochs == EXPECTED_EPOCHS
                    and abs(last_t - bt.DURATION) < 0.5 * bt.OUTPUT_STEP)
        tol = bt.LEVELS[level]
        config = {
            "campaign": "R41",
            "purpose": ("reference trajectory at the raised degree, for the "
                        "reference-degree control (O41)"),
            "design": design, "sobol_index": index,
            "reference_degree": r41.NEW_REFERENCE,
            "base_reference_degree": r41.BASE_REFERENCE,
            "tolerance_level": level,
            "rtol": tol["rtol"],
            "atol_position_m": tol["atol_position_m"],
            "atol_velocity_m_s": tol["atol_velocity_m_s"],
            "max_step_s": bt.MAX_STEP,
            "duration_s": bt.DURATION,
            "output_step_s": bt.OUTPUT_STEP,
            "initial_state_si": list(row["design_point"]["initial_state_si"]),
            "frame": ("inertial, Moon rotating uniformly at its sidereal rate "
                      "about the polar axis, gravity only"),
            "integrator": "DOP853",
            "propagator_source": ("rev14_budget_trajectory._propagate, called "
                                  "through rev41_reference_degree_control"),
        }
        sidecar = CASES / npz.relative_to(RAW).with_suffix(".json")
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        base.atomic_json(sidecar, {
            "schema": "r41_reference_trajectory_v1",
            "created_utc": base.utc_now(),
            "config": config,
            "config_sha256": base.object_hash(config),
            "status": "complete" if complete else "incomplete",
            "event": None if complete else "arc did not reach the full grid",
            "telemetry": None,
            "telemetry_note": (
                "not retained: the R41 driver computed the step and call "
                "counts and did not write them. The driver is patched to "
                "write them inline; this sidecar was generated after the run "
                "and does not invent them."),
            "raw_path": str(npz.relative_to(ROOT)).replace("\\", "/"),
            "raw_sha256": base.file_hash(npz),
            "n_output_epochs": n_epochs,
            "last_output_epoch_s": last_t,
        })
        written += 1
        if not complete:
            mismatched.append(f"{design}{index:03d}/{level}: {n_epochs} epochs, "
                              f"last t={last_t:.0f}s")

    digests = sorted(json.loads(p.read_text(encoding="utf-8"))["raw_sha256"]
                     for p in CASES.rglob("reference_*.json"))
    rollup = base.object_hash(digests)
    (METRICS / "r41_trajectory_index.json").write_text(json.dumps({
        "schema": "r41_trajectory_index_v1",
        "created_utc": base.utc_now(),
        "trajectories": written,
        "expected_output_epochs": EXPECTED_EPOCHS,
        "rolled_up_raw_digest": rollup,
        "note": ("one sidecar per raw state array, each carrying the config "
                 "that determines the trajectory and the digest of the array. "
                 "The arrays themselves are excluded from the public "
                 "repository as regenerable, which is the policy applied to "
                 "every campaign; the digests travel with the manifest."),
        "files": [str(p.relative_to(METRICS)).replace("\\", "/")
                  for p in sorted(CASES.rglob("reference_*.json"))],
    }, indent=2), encoding="utf-8")

    print(f"[written] {written} sidecars under {CASES.relative_to(ROOT)}")
    print(f"[written] r41_trajectory_index.json  rolled-up digest "
          f"{rollup[:16]}")
    if mismatched:
        print("[note] " + "; ".join(mismatched))
        return 1
    print(f"  every trajectory carries the full {EXPECTED_EPOCHS}-epoch grid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
