"""R62 (O54): the controlled apolune ladder's interior panel, re-matched at
the scoring tolerance.

O53 established that the level the realized-work match is made at changes
which populations share the crossing bracket, and that two of those changes
survive every resolution cut. That leaves one asymmetry in the paper: the
geometry strata now carry the level-consistent match, while the *controlled*
test of the radial-span direction, the paired apolune ladder of (O49), still
carries the older tight-level convention on its interior panel. The
uncontrolled evidence is held to a stricter accounting than the controlled
evidence, which is the wrong way round.

This campaign runs (O42)'s construction, unchanged, on the ladder's interior
member, on the two independent identity blocks, at the two apolune levels
where the reference-degree ceiling does not bind.

Scope, and why it stops where it does. Only the 300 and 600~km levels are
run. (O49) reports the ceiling beginning to bind at 1200~km and binding on
every orbit at 2400~km, so a rematch there would confound the accounting
change with the ceiling, and a censored comparator is dropped rather than
clamped in any case. The question this campaign asks -- whether a controlled
apolune change still moves the crossing once the match is made at the scoring
tolerance -- is answered by the 300 to 600 step alone.

Partition: r62_* records, its own case and raw trees. R44's and R61's records
are untouched, and r18_cases is read read-only.

Usage:
    R62_POP=span_ladder_a python rev62_ladder_interior_rematch.py plan --beta 0.75
    R62_POP=span_ladder_b python rev62_ladder_interior_rematch.py run --beta 0.75 --workers 10
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"

POP = os.environ.get("R62_POP")
if POP not in ("span_ladder_a", "span_ladder_b"):
    raise SystemExit("R62_POP must be span_ladder_a or span_ladder_b")

# The levels the ceiling leaves alone. Stated here as data, not as a comment,
# because the campaign's scope claim is exactly this set.
APOLUNE_KM = (300.0, 600.0)

import population_registry as registry                      # noqa: E402
import rev14_budget_pareto as pareto                        # noqa: E402
import rev14_budget_trajectory as r14                       # noqa: E402
import rev18_span_sweep as r18                              # noqa: E402
import rev44_equal_work_tighter as r44                      # noqa: E402

KEY = registry.spec("r50", POP)["design_key"]
ROWS = METRICS / f"r50_{POP}_rows.json"
REUSE_CASE = METRICS / "r11_cases" / f"stratum_{POP}_convergence"
REUSE_RAW = METRICS / "r11_raw" / f"stratum_{POP}_convergence"
CALIBRATION = METRICS / f"r50_budget_pareto_{POP}.json"

# --- identity, registered at import so spawned children share it ------------
r14.DESIGNS[KEY] = {"rows": ROWS,
                    "reuse_case": REUSE_CASE,
                    "reuse_raw": REUSE_RAW}
pareto.DESIGNS[KEY] = {"rows": ROWS,
                       "r12_case": METRICS / "r50_cases" / f"atallah_{POP}",
                       "r11_raw": REUSE_RAW}
r14.PARETO = CALIBRATION

# --- partition: r62 trees, never r44's or r61's -----------------------------
r44.CASE_ROOT = METRICS / "r62_cases"
r44.RAW_ROOT = METRICS / "r62_raw"


def out_path(design: str, beta: float) -> Path:
    return METRICS / (f"r62_ladder_interior_{design}_"
                      f"{r18.beta_tag(beta)}.json")


r44.out_path = out_path


def kept_indices(design: str, beta: float) -> set[int]:
    """Sobol indices at the two unconfounded apolune levels."""
    span = json.loads(
        (METRICS / f"r18_span_sweep_{design}_{r18.beta_tag(beta)}.json"
         ).read_text(encoding="utf-8"))
    return {int(r["sobol_index"]) for r in span["rows"]
            if float(r.get("ha_km", -1)) in APOLUNE_KM}


_orig_build = r44.build_tasks


def build_tasks(design: str, beta: float, verbose: bool = False):
    """R44's task construction, restricted to the two levels this campaign
    declares. The restriction is applied to the built tasks rather than to
    the span record, so the degree estimate each kept orbit receives is
    byte-for-byte the one R44's own code would have produced for it.

    R44 writes its censored list under an r44_ name. Any censoring here comes
    from the levels this campaign excludes anyway, so that file is rewritten
    under the r62_ prefix and the r44_ one removed rather than left behind to
    look like part of the sealed R44 campaign.
    """
    tag = r18.beta_tag(beta)
    stray = METRICS / f"r44_censored_{design}_{tag}.json"
    pre_existing = stray.exists()
    tasks, censored, degenerate, missing = _orig_build(design, beta, verbose)
    keep = kept_indices(design, beta)
    if censored and not pre_existing and stray.exists():
        kept_censored = [c for c in censored
                         if int(c["sobol_index"]) in keep]
        stray.unlink()
        if kept_censored:
            (METRICS / f"r62_censored_{design}_{tag}.json").write_text(
                json.dumps(kept_censored, indent=2), encoding="utf-8")
    return ([t for t in tasks if t["index"] in keep],
            [c for c in censored if int(c["sobol_index"]) in keep],
            [d for d in degenerate if d in keep],
            [m for m in missing if m in keep])


r44.build_tasks = build_tasks


def main() -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "run", "summarize"):
        s = sub.add_parser(name)
        s.add_argument("--beta", type=float, default=1.0)
        if name == "run":
            s.add_argument("--workers", type=int, default=10)
            s.add_argument("--deadline-min", type=float, default=180.0)
        s.set_defaults(cmd_name=name)
    a = p.parse_args()
    a.design = KEY
    n = len(kept_indices(KEY, float(a.beta)))
    print(f"[r62] population={POP} key={KEY} cell=beta_{a.beta:.2f} "
          f"apolune={'/'.join(f'{x:.0f}' for x in APOLUNE_KM)} km "
          f"({n} of 64 identities kept)")
    return {"plan": r44.plan, "run": r44.run,
            "summarize": r44.summarize}[a.cmd_name](a)


if __name__ == "__main__":
    raise SystemExit(main())
