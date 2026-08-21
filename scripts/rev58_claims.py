"""Adds the R58 endpoint work-matched numbers to the claims ledger.

Every number the manuscript prints is supposed to be pinned to a sealed record,
so the three quantities the control puts in the main text and the supplement
get entries here rather than standing on prose alone. The control was run after
the registered outcomes were scored, so the entries are exploratory.

Idempotent: re-running replaces the entries it owns instead of appending
duplicates, and the digests are recomputed from the records on disk.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
LEDGER = METRICS / "claims_ledger.json"

OEU1 = "r58_endpoint_equal_work_OEU_beta_1.00.json"
OEU075 = "r58_endpoint_equal_work_OEU_beta_0.75.json"

# appears_as must be a literal fragment of the typeset corpus; these two are
# the sentence in Section VIII and the sentence in the supplement subsection.
MAIN_FRAGMENT = "wins 26--15 at a median ratio of $1.52$"
SUPP_FRAGMENT = "20--41 at $\\beta=0.50$, $0.62$ and $0.75$"

ENTRIES = [
    {
        "id": "workmatched.oeu.tally_b100_radial",
        "claim": "at equal realized total quadratic work the ceiling-free "
                 "radial endpoint wins 26 resolved comparisons at the "
                 "declared budget",
        "status": "exploratory",
        "appears_in": ["chapters/08_budget.tex",
                       "chapters/supp_operational.tex"],
        "source": OEU1,
        "check": {"kind": "path", "path": ["summary", "resolved_radial_wins"]},
        "expect": 26,
        "tol": 0,
        "appears_as": MAIN_FRAGMENT,
    },
    {
        "id": "workmatched.oeu.tally_b100_constant",
        "claim": "at equal realized total quadratic work the constant degree "
                 "takes 15 resolved comparisons at the declared budget",
        "status": "exploratory",
        "appears_in": ["chapters/08_budget.tex",
                       "chapters/supp_operational.tex"],
        "source": OEU1,
        "check": {"kind": "path",
                  "path": ["summary", "resolved_constant_wins"]},
        "expect": 15,
        "tol": 0,
        "appears_as": MAIN_FRAGMENT,
    },
    {
        "id": "workmatched.oeu.tally_b075_radial",
        "claim": "the beta = 0.75 cell changes side under work matching: the "
                 "radial endpoint takes 20 of the resolved comparisons there",
        "status": "exploratory",
        "appears_in": ["chapters/supp_operational.tex"],
        "source": OEU075,
        "check": {"kind": "path", "path": ["summary", "resolved_radial_wins"]},
        "expect": 20,
        "tol": 0,
        "appears_as": SUPP_FRAGMENT,
    },
]

OWNED = {e["id"] for e in ENTRIES}


def digest(name: str) -> str:
    h = hashlib.sha256()
    h.update((METRICS / name).read_bytes())
    return h.hexdigest()


def main() -> int:
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    claims = [c for c in led["claims"] if c.get("id") not in OWNED]
    for e in ENTRIES:
        e = dict(e)
        e["source_sha256"] = digest(e["source"])
        claims.append(e)
    led["claims"] = claims
    LEDGER.write_text(json.dumps(led, indent=2), encoding="utf-8")
    print(f"[r58-claims] ledger now holds {len(claims)} claims "
          f"({len(ENTRIES)} owned here)")
    for e in ENTRIES:
        print(f"  {e['id']}: expect {e['expect']} from {e['source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
