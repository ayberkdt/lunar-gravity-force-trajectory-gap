"""R63 (O55): the ceiling-free apolune ladder, interior panel, re-matched at
the scoring tolerance.

(O54) carried the level-consistent match to the controlled apolune ladder but
stopped at 300 and 600 km, because on the capped blocks the reference-degree
ceiling begins to bind at 1200 km and binds everywhere at 2400 km, which would
have confounded the accounting change with the ceiling. The ceiling-free
blocks of (O50) and (O51) exist precisely to remove that confound, so on them
the two widest levels can be run.

This campaign therefore does two things at once. On 300 and 600 km it
replicates (O54)'s result on independent ceiling-free identities. On 1200 and
2400 km it asks the question (O54) could not: whether the interior member's
advantage keeps widening with apolune once both the ceiling and the
accounting objection are removed.

Partition: r63_* records and trees. R44, R61, R62 and the ladder campaigns
themselves are untouched; r18_cases is read read-only.

Usage:
    R63_POP=span_ladder_a_uncapped python rev63_ladder_uncapped_rematch.py plan --beta 1.00
    R63_POP=span_ladder_b_uncapped python rev63_ladder_uncapped_rematch.py run --beta 1.00 --workers 10
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

# population -> the registration that declares it. The two ceiling-free blocks
# were filed under different registrations (block A with R51, block B with
# R52), which is why the registry is looked up per population rather than
# assumed.
REGISTRY_OF = {"span_ladder_a_uncapped": ("r51", "r51"),
               "span_ladder_b_uncapped": ("r52", "r52")}

POP = os.environ.get("R63_POP")
if POP not in REGISTRY_OF:
    raise SystemExit(f"R63_POP must be one of {sorted(REGISTRY_OF)}")

import population_registry as registry                      # noqa: E402
import rev14_budget_pareto as pareto                        # noqa: E402
import rev14_budget_trajectory as r14                       # noqa: E402
import rev18_span_sweep as r18                              # noqa: E402
import rev44_equal_work_tighter as r44                      # noqa: E402

REG, PREFIX = REGISTRY_OF[POP]
KEY = registry.spec(REG, POP)["design_key"]
ROWS = METRICS / f"{PREFIX}_{POP}_rows.json"
REUSE_CASE = METRICS / "r11_cases" / f"stratum_{POP}_convergence"
REUSE_RAW = METRICS / "r11_raw" / f"stratum_{POP}_convergence"
CALIBRATION = METRICS / f"{PREFIX}_budget_pareto_{POP}.json"

# --- identity, registered at import so spawned children share it ------------
r14.DESIGNS[KEY] = {"rows": ROWS,
                    "reuse_case": REUSE_CASE,
                    "reuse_raw": REUSE_RAW}
pareto.DESIGNS[KEY] = {"rows": ROWS,
                       "r12_case": METRICS / f"{PREFIX}_cases" / f"atallah_{POP}",
                       "r11_raw": REUSE_RAW}
r14.PARETO = CALIBRATION

# --- partition: r63 trees ---------------------------------------------------
r44.CASE_ROOT = METRICS / "r63_cases"
r44.RAW_ROOT = METRICS / "r63_raw"


def out_path(design: str, beta: float) -> Path:
    return METRICS / (f"r63_ladder_uncapped_{design}_"
                      f"{r18.beta_tag(beta)}.json")


r44.out_path = out_path

_orig_build = r44.build_tasks


def build_tasks(design: str, beta: float, verbose: bool = False):
    """R44's construction unchanged, with its censored-list filename moved
    into this campaign's prefix. R44 hard-codes that name, and leaving an
    r44_ file behind would grow the file family of a sealed campaign that
    never ran these populations."""
    tag = r18.beta_tag(beta)
    stray = METRICS / f"r44_censored_{design}_{tag}.json"
    pre_existing = stray.exists()
    result = _orig_build(design, beta, verbose)
    if not pre_existing and stray.exists():
        stray.replace(METRICS / f"r63_censored_{design}_{tag}.json")
    return result


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
    print(f"[r63] population={POP} registry={REG} key={KEY} "
          f"cell=beta_{a.beta:.2f} (all four apolune levels)")
    return {"plan": r44.plan, "run": r44.run,
            "summarize": r44.summarize}[a.cmd_name](a)


if __name__ == "__main__":
    raise SystemExit(main())
