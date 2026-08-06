"""Which campaign owns which budget subtree, in one place.

Several campaigns reuse an earlier campaign's driver with a budget argument, so
their trajectories land under that earlier campaign's file-name prefix. R23 did
this for beta = 0.50 and R25 for the budgets it swept. The manifests must still
partition the trajectory records: a record indexed under two manifests turns
the stated partition into overlapping inventories, and the integrity check
fails on it.

Ownership is declared here as an allow-list per campaign rather than as a
deny-list, because a deny-list has to be edited every time a later campaign
adds a budget, and forgetting to edit it produces a silent overlap rather than
an error. An allow-list forgets in the safe direction: a budget nobody claims
shows up as an unindexed record, which the integrity check also reports.

The budgets below are the ones each campaign propagated itself, taken from the
records the campaign sealed. Later additions under the same prefix belong to
the campaign that ran them and are indexed only there.
"""

from __future__ import annotations

# (design, beta tag) pairs each campaign propagated under its own prefix.
R14_OWN = {("A", "0.50"), ("A", "0.75"), ("A", "1.00"), ("A", "1.50"),
           ("A", "3.00"), ("B", "0.50"), ("B", "1.00")}
R18_OWN = {("A", "0.50"), ("A", "1.00"), ("A", "1.50"),
           ("B", "0.50"), ("B", "1.00")}

# R19 is the beta = 1 campaign and beta = 1 is the unsuffixed case, so it owns
# exactly the subtrees carrying no budget suffix.
R19_SUFFIXLESS = True

# Everything else under those prefixes belongs to a later campaign:
#   R23  r19 at beta = 0.50 on both designs
#   R25  r14 at B 0.75, A 1.25, B 1.50
#        r18 at A 0.75, B 0.75, A 1.25, B 1.50
#        r19 at A 0.75, B 0.75, A 1.25, B 1.50


def _budget_of(part: str) -> tuple[str, str] | None:
    """Parse ``A_beta_0.75`` or ``A_beta_0.75_k_0.50`` into ``("A", "0.75")``."""
    if "_beta_" not in part:
        return None
    design, rest = part.split("_beta_", 1)
    design = design.split("_")[-1] or design
    return design, rest.split("_k_")[0]


def owned(path_parts, own: set[tuple[str, str]]) -> bool:
    """True when no part of the path names a budget outside ``own``."""
    for part in path_parts:
        budget = _budget_of(part)
        if budget is not None and budget not in own:
            return False
    return True


def owned_by_r14(path) -> bool:
    return owned(path.parts, R14_OWN)


def owned_by_r18(path) -> bool:
    return owned(path.parts, R18_OWN)
