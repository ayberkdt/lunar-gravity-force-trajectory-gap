"""Pins the realized-work excess.

The median 29% figure once appeared in four places; the compression rounds of
August 2026 left it in Section VI alone, with the abstract keeping the
qualitative claim ("at greater realized cost") and the Discussion and
Conclusion dropping it. Until this entry it carried no ledger record: the
checker verified 68 other numbers and never asked about this one. That is the
blind spot of a ledger -- it audits the claims it holds, not the claims it
lacks -- so the number is pinned here on both designs.

Idempotent: re-running replaces the entries it owns.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
LEDGER = METRICS / "claims_ledger.json"

# Re-anchored twice. 2026-08-17: the editorial pass dropped "of the two" and
# the median-of-products aside. 2026-08-20: the sentence was rewritten to say
# the cost once and directly, so the phrase moved again. The number and its
# four printing sites are unchanged; only the phrase the checker matches on is.
FRAGMENT = ("median $29\\%$ more realized quadratic work than the "
            "constant degree on both")

ENTRIES = [
    {
        "id": "budget.realized_work_excess.A",
        "claim": "on design A the radial endpoint's realized total quadratic "
                 "work is a median 1.29 times the constant degree's at the "
                 "declared budget",
        "status": "confirmatory",
        "appears_in": ["chapters/08_budget.tex"],
        "source": "r14_trajectory_A_beta_1.00.json",
        "check": {"kind": "path",
                  "path": ["summary", "total_work_ratio", "median"]},
        "expect": 1.2892193847594222,
        "tol": 0.0005,
        "appears_as": FRAGMENT,
    },
    {
        "id": "budget.realized_work_excess.B",
        "claim": "on design B the same ratio has median 1.29",
        "status": "confirmatory",
        "appears_in": ["chapters/08_budget.tex"],
        "source": "r14_trajectory_B_beta_1.00.json",
        "check": {"kind": "path",
                  "path": ["summary", "total_work_ratio", "median"]},
        "expect": 1.2895691167034848,
        "tol": 0.0005,
        "appears_as": FRAGMENT,
    },
]

OWNED = {e["id"] for e in ENTRIES}


def digest(name: str) -> str:
    return hashlib.sha256((METRICS / name).read_bytes()).hexdigest()


def main() -> int:
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    claims = [c for c in led["claims"] if c.get("id") not in OWNED]
    for e in ENTRIES:
        e = dict(e)
        e["source_sha256"] = digest(e["source"])
        claims.append(e)
    led["claims"] = claims
    LEDGER.write_text(json.dumps(led, indent=2), encoding="utf-8")
    print(f"[r60-claims] ledger now holds {len(claims)} claims")
    for e in ENTRIES:
        print(f"  {e['id']}: expect {e['expect']:.6f} from {e['source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
