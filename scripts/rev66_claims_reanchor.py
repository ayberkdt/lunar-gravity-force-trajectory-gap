"""Re-anchor the band-share claim after the Introduction was made precise.

The claims ledger matches each entry against a fragment of the manuscript, so
a rewritten sentence breaks the check even when the number behind it has not
moved. The Introduction used to say the 61--100 band "carries 1.9% of the total
perturbing acceleration", which reads as an additive share; the band quantities
are vector RMS amplitudes and the supplement says so, and the sentence now says
so too. The number, its source record and its printing site are unchanged; only
the phrase the checker matches on is.

This script edits nothing else: it loads the ledger, replaces the one
`appears_as` fragment, and writes the file back.

Usage:  python rev66_claims_reanchor.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "metrics" / "claims_ledger.json"

CLAIM_ID = "band.share_61_100_at_100km"
OLD = "$1.9\\%$ of the total perturbing acceleration at\n100~km"
NEW = "RMS amplitude equal to $1.9\\%$ of the total\nperturbation RMS at 100~km"


def main() -> int:
    payload = json.loads(LEDGER.read_text(encoding="utf-8"))
    entries = payload["claims"] if isinstance(payload, dict) else payload
    seq = entries if isinstance(entries, list) else list(entries.values())
    hit = [e for e in seq if e.get("id") == CLAIM_ID]
    if len(hit) != 1:
        raise SystemExit(f"[refuse] expected one {CLAIM_ID}, found {len(hit)}")
    entry = hit[0]
    if entry.get("appears_as") == NEW:
        print("[ok] already re-anchored")
        return 0
    if entry.get("appears_as") != OLD:
        raise SystemExit("[refuse] fragment is not the one this script expects: "
                         f"{entry.get('appears_as')!r}")
    entry["appears_as"] = NEW
    LEDGER.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[re-anchored] {CLAIM_ID}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
