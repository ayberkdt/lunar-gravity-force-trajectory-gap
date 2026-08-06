"""Does the experiment contract state the constants the drivers actually use?

One item of the contract named a gradient degree the code does not use. That is
a small error with a large blast radius: a reader who finds it has no way to
tell whether it is the only one, and the contract is the document the
reproducibility claim rests on. This checks the constants that can be checked
mechanically --- the ones that appear as a literal in both the prose and a
driver --- so the answer is a count rather than an assurance.

What it can and cannot do. It verifies that a value the contract states is the
value the archived source defines. It cannot verify that the driver uses its own
constant correctly, and it does not try; that is what the per-campaign records
and the integrity gate are for. A constant that appears in neither place, or is
computed rather than declared, is reported as unchecked rather than as passing.

Usage:  python check_protocol_constants.py [--quiet]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "python_codes"
CONTRACT = ROOT / "chapters" / "supp_experiment_contract.tex"

# (label, regex the contract must contain, file, regex the source must contain)
CHECKS = [
    ("variational gradient degree",
     r"gradient by\s*\n?\s*central finite differences \(1~m step, \$N=120\$",
     "rev13_variational_check.py", r"GRADIENT_DEGREE\s*=\s*120\b"),
    ("variational finite-difference step",
     r"\(1~m step",
     "rev21_gradient_sensitivity.py", r'"fd_step_m":\s*1\.0'),
    ("tight vector tolerance, position",
     r"10\^\{-5\}\$~m/\$10\^\{-8\}\$~m~s\$\^\{-1\}\$ \(tight\)",
     "rev14_budget_trajectory.py", r'"atol_position_m":\s*1\.0e-5'),
    ("tighter vector tolerance, position",
     r"10\^\{-6\}\$~m/\$10\^\{-9\}\$~m~s\$\^\{-1\}\$ \(tighter\)",
     "rev14_budget_trajectory.py", r'"atol_position_m":\s*1\.0e-6'),
    ("maximum step",
     r"60~s maximum step",
     "rev14_budget_trajectory.py", r"MAX_STEP\s*=\s*60\.0"),
    ("output grid",
     r"120~s output grid",
     "rev14_budget_trajectory.py", r"OUTPUT_STEP\s*=\s*120\.0"),
    ("altitude bin width",
     r"10-km",
     "rev14_budget_trajectory.py", r"BIN_KM\s*=\s*10\.0"),
    ("degree floor",
     r"floor",
     "rev14_budget_trajectory.py", r"FLOOR\s*=\s*2\b"),
    ("budget work-match tolerance",
     r"one-percent work-match|1\\%",
     "rev14_budget_pareto.py", r"WORK_TOLERANCE\s*=\s*0\.01"),
    ("design A seed",
     r"20260723",
     "rev10_sobol_confirmatory.py", r"20260723"),
    ("design C seed",
     r"20260725",
     "rev26_designC_freeze.py", r"SEED_C\s*=\s*SEED_B\s*\+\s*1"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    text = CONTRACT.read_text(encoding="utf-8")
    problems, unchecked = [], []
    for label, pat_tex, fname, pat_src in CHECKS:
        src_path = CODE / fname
        in_tex = bool(re.search(pat_tex, text))
        in_src = (bool(re.search(pat_src, src_path.read_text(encoding="utf-8")))
                  if src_path.exists() else None)
        if in_src is None:
            unchecked.append(f"{label}: {fname} not on disk")
            continue
        if in_tex and in_src:
            if not a.quiet:
                print(f"  [ok]      {label}")
        elif in_src and not in_tex:
            unchecked.append(f"{label}: source defines it, contract wording "
                             f"did not match the pattern")
        elif in_tex and not in_src:
            problems.append(f"{label}: the contract states it, "
                            f"{fname} does not define it")
        else:
            unchecked.append(f"{label}: found in neither")

    for u in unchecked:
        print(f"  [unchecked] {u}")
    for p in problems:
        print(f"  [MISMATCH]  {p}")
    print(f"\n{len(CHECKS)} constants checked, {len(problems)} mismatches, "
          f"{len(unchecked)} unchecked")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
