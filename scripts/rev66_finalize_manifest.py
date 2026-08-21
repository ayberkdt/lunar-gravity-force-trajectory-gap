"""SHA-256 integrity manifest for R66: the allocation-anatomy figure.

R66 is a derived product, not a campaign. It propagates nothing, rescores
nothing and adds no number the archived records do not already carry: the
editorial round that moved the qualification and control summaries out of the
main text used part of the space it freed for a figure showing what the
equal-budget construction actually does along one orbit.

Everything the figure draws is read from frozen records. The design point and
the verdict come from the (O14) design-A record at beta = 1; the binned degree
schedule comes from that orbit's own archived case configuration; the three
state arrays are the ones those campaigns wrote, pinned by the raw_sha256
fields their sidecars already carry, and stay owned by the R11 and R14
manifests. They are therefore not re-indexed here, following the convention
R59 and R60 use for reused inputs.

The orbit is selected by geometry alone --- the design-A member whose perilune
and apolune are jointly closest to the design medians in log distance --- so
the figure's subject cannot have been chosen on its outcome. The maker refuses
to draw unless the position RMS it recomputes reproduces the record's own
atallah_error_m and fixed_error_m, and unless the schedule it reconstructs
still meets the campaign's one-percent work-match target.

Usage:  python rev66_finalize_manifest.py
"""

from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
FIGURES = ROOT / "figures"
SCRIPTS = ROOT / "python_codes"

ORBIT = "sobolA_032"
REUSED = [
    "r14_trajectory_A_beta_1.00.json",
    "e1_band_shares_60_100_nmax300.json",
]
# rev66_claims_reanchor.py writes into the claims ledger, so it carries a
# number into the manuscript audit and needs an owner, as rev60_claims.py
# does for R60.
SCRIPTS_OWNED = ["make_figures_allocation_anatomy.py",
                 "rev66_claims_reanchor.py"]
FIGURE = "fig_allocation_anatomy.pdf"


def entry(path: Path) -> dict:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"sha256": h, "bytes": path.stat().st_size}


def main() -> int:
    payload = {
        "schema": "r66_final_experiment_manifest_v1",
        "created_utc": datetime.datetime.now(datetime.timezone.utc)
                        .isoformat().replace("+00:00", "Z"),
        "scope": ("R66: the main-text allocation-anatomy figure, drawn from "
                  "the archived (O14) design-A beta = 1 record, that orbit's "
                  "case configuration, and the R11/R14 state arrays those "
                  "campaigns wrote."),
        "why": ("The endpoint comparison was reported as a population verdict "
                "and explained dynamically, with nothing between the two "
                "showing the allocation itself. The figure supplies the "
                "degree history the budget buys and the in-track displacement "
                "it produces, on one orbit chosen by geometry alone."),
        "propagation": ("none. No orbit is propagated and none is re-scored. "
                        "The maker refuses to draw unless the recomputed "
                        "position RMS reproduces the record's own error "
                        "fields and the reconstructed schedule still holds "
                        "the declared budget."),
        "orbit_selection": ("design-A member minimising |log(hp/median hp)| + "
                            "|log(ha/median ha)|; outcome-blind"),
        "state_arrays": {
            "owner_manifests": ["r11", "r14"],
            "note": ("truth_tight.npz and fixed_critical_tight.npz under the "
                     "R11 convergence tree, atallah_budget_tight.npz and its "
                     "case configuration under the R14 beta = 1 tree; each "
                     "pinned by the raw_sha256 of its own sidecar and indexed "
                     "by its own manifest. Trajectory records are a partition, "
                     "so they are named here and re-indexed nowhere."),
        },
        "reused_inputs": {n: entry(METRICS / n) for n in REUSED},
        "scripts": {s: entry(SCRIPTS / s) for s in SCRIPTS_OWNED},
        "generated_figures": {f"figures/{FIGURE}": entry(FIGURES / FIGURE)},
    }
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["manifest_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    out = METRICS / "r66_final_experiment_manifest.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] metrics/{out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
