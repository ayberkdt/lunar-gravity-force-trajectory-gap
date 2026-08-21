"""R61 (O42-ext): the scoring-tolerance rematch carried to the populations
R44 never covered.

R44 re-established the equal-realized-work match at the level the errors are
scored at, but only on the two confirmatory coverage designs. The manuscript
nevertheless reads the crossing bracket across seven populations and says so
in Section IX.B: "the scoring-tolerance rematch exists for the two
confirmatory designs only, and on design B it moves the tally crossing above
0.75". One of two populations moved, and a convention claim about seven rests
on that. Council round 22, item 3(c), named the same gap from the other side.

This campaign runs R44's comparator construction, unchanged, on the third
coverage design and the five geometry strata. Nothing about the method is new:
the driver here only supplies the population identity that R44's `--design`
switch assumes, exactly as rev30_stratum_ops.py supplies it to R18/R19.

Partition: this campaign writes r61_* records, its own case and raw trees. The
sealed R44 manifest covers eight A/B cells and keeps covering exactly those;
no r44_* file is written, read-only reuse of r18_cases excepted.

The population is passed through the environment, not argv, so that the
ProcessPoolExecutor children -- which re-import this module under spawn with
an argv of their own -- register the same identity the parent did. Registering
at import is what rev29 and rev30 do, and for the same reason.

Usage:
    R61_POP=low_perilune python rev61_equal_work_tighter_ext.py plan --beta 0.75
    R61_POP=C python rev61_equal_work_tighter_ext.py run --beta 0.75 --workers 10
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

ROOT = HERE.parents[0]
METRICS = ROOT / "metrics"

POP = os.environ.get("R61_POP")
if not POP:
    raise SystemExit("R61_POP must name the population "
                     "(C | polar | equatorial | frozen_like | low_perilune | "
                     "high_apolune)")

import population_registry as registry                      # noqa: E402
import rev14_budget_pareto as pareto                        # noqa: E402
import rev14_budget_trajectory as r14                       # noqa: E402
import rev18_span_sweep as r18                              # noqa: E402
import rev44_equal_work_tighter as r44                      # noqa: E402


def _identity(pop: str):
    """(key, rows, reuse_case, reuse_raw, calibration) for one population.

    Design C is declared in the R29 registration and the five strata in R30;
    the two registrations have different shapes on purpose, which is why the
    keys come from population_registry rather than from a table here.
    """
    if pop == "C":
        return ("C",
                METRICS / "r26_designC_rows.json",
                METRICS / "r11_cases" / "designC_convergence",
                METRICS / "r11_raw" / "designC_convergence",
                METRICS / "r29_budget_pareto_designC.json")
    key = registry.spec("r30", pop)["design_key"]
    return (key,
            METRICS / f"r30_{pop}_rows.json",
            METRICS / "r11_cases" / f"stratum_{pop}_convergence",
            METRICS / "r11_raw" / f"stratum_{pop}_convergence",
            METRICS / f"r30_budget_pareto_{pop}.json")


KEY, ROWS, REUSE_CASE, REUSE_RAW, CALIBRATION = _identity(POP)

# --- identity, registered at import so spawned children share it ------------
r14.DESIGNS[KEY] = {"rows": ROWS,
                    "reuse_case": REUSE_CASE,
                    "reuse_raw": REUSE_RAW}
pareto.DESIGNS[KEY] = {"rows": ROWS,
                       "r12_case": METRICS / "r30_cases" / f"atallah_{POP}",
                       "r11_raw": REUSE_RAW}
r14.PARETO = CALIBRATION

# --- partition: r61 trees, never r44's --------------------------------------
r44.CASE_ROOT = METRICS / "r61_cases"
r44.RAW_ROOT = METRICS / "r61_raw"


def out_path(design: str, beta: float) -> Path:
    return METRICS / (f"r61_equal_work_tighter_{design}_"
                      f"{r18.beta_tag(beta)}.json")


r44.out_path = out_path


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
    print(f"[r61] population={POP} key={KEY} cell=beta_{a.beta:.2f}")
    return {"plan": r44.plan, "run": r44.run,
            "summarize": r44.summarize}[a.cmd_name](a)


if __name__ == "__main__":
    raise SystemExit(main())
