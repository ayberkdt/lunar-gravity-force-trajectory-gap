"""SHA-256 integrity manifest for R60: the editorial round's derivatives.

R60 is a derived product, not a campaign: it propagates nothing, rescores
nothing, and adds no number the archived records do not already carry. It
holds what the editorial round produced, and the round produced three things:
the claims-ledger entries pinning the realized-work excess, the main-text
allocation-interior figure drawn from the (O28) family records, and the beta=1
headline table split out of the full budget-grid table when that table moved
to the supplement.

rev60_claims.py is indexed here for the reason rev58_claims.py is indexed by
the R58 manifest: a script that writes entries into the claims ledger carries
numbers into the manuscript's audit, and it was the one artifact of this round
that no manifest owned.

Source records stay owned by their own campaign manifests; they are recorded
here as reused inputs by digest, the same convention R59 uses.

Usage:  python rev60_finalize_manifest.py
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

REUSED = ["r18_span_sweep_A_beta_1.00.json",
          "r18_span_sweep_B_beta_1.00.json"]
SCRIPTS_OWNED = ["make_figures_r60_interior.py", "rev60_claims.py"]
FIGURE = "fig_allocation_interior.pdf"
# The beta = 1 headline table is emitted by R14's own sealed generator
# (rev14_tables.py) from R14's frozen records; R60 owns the file because it
# exists for the same editorial pass and R14's finalizer cannot be re-run
# without sweeping later campaigns into its inventory.
TABLE = "r14_trajectory_pareto_headline_table.tex"


def entry(path: Path) -> dict:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"sha256": h, "bytes": path.stat().st_size}


def main() -> int:
    payload = {
        "schema": "r60_final_experiment_manifest_v1",
        "created_utc": datetime.datetime.now(datetime.timezone.utc)
                        .isoformat().replace("+00:00", "Z"),
        "scope": ("R60: the editorial round's derivatives. The claims-ledger "
                  "entries pinning the realized-work excess, the main-text "
                  "allocation-interior figure, drawn from the archived (O28) "
                  "five-member family records at beta=1 on designs A and B, "
                  "and the beta=1 headline table split out of the full "
                  "budget-grid table when that table moved to the "
                  "supplement."),
        "why": ("An editorial revision moved the interior-existence result "
                "into a main-text figure. The figure reads the frozen (O28) "
                "records only; its per-k counts are the records' own archived "
                "best_k_counts, cross-checked against the raw argmin before "
                "anything is drawn."),
        "propagation": ("none. No orbit is propagated and none is re-scored: "
                        "the errors plotted and the argmin counted are the "
                        "ones the (O28) campaign wrote."),
        "reused_inputs": {n: entry(METRICS / n) for n in REUSED},
        "scripts": {s: entry(SCRIPTS / s) for s in SCRIPTS_OWNED},
        "generated_figures": {f"figures/{FIGURE}": entry(FIGURES / FIGURE)},
        "generated_tables": {TABLE: entry(METRICS / TABLE)},
    }
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["manifest_sha256"] = hashlib.sha256(body.encode()).hexdigest()
    out = METRICS / "r60_final_experiment_manifest.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[written] metrics/{out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
