"""R30: generate the propagated base of one geometry stratum.

The same two steps designs A, B and C needed, pointed at a stratum: the
empirical prepass that fixes each orbit's adopted truth degree, its critical
degree and its work degree, and the convergence tree that holds the truth
trajectories and the comparator policies at both tolerance levels.

No new machinery, and the same spawn hazard as design C. rev11_full_convergence
resolves its output tree from environment variables read at import, because a
ProcessPoolExecutor on Windows re-imports it in every child and a parent-only
patch would leave those children writing into design A's tree. The stratum is
therefore parsed off the command line before any of that is imported, the
environment is set, and the parent's globals are realigned to match.

The design is read from the frozen file and never regenerated.

Usage:
    python rev30_stratum_base.py --stratum polar prepass --workers 11
    python rev30_stratum_base.py --stratum polar run --workers 11 --deadline "<iso>"
    python rev30_stratum_base.py --stratum polar status
    python rev30_stratum_base.py --registry r31 --stratum operational_elliptical run
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"



def _from_argv(flag: str, default: str | None = None) -> str:
    """Read an argument before anything imports the convergence driver, because
    the driver decides its output tree at import time from the environment."""
    for i, arg in enumerate(sys.argv):
        if arg == flag and i + 1 < len(sys.argv):
            return sys.argv[i + 1]
        if arg.startswith(flag + "="):
            return arg.split("=", 1)[1]
    if default is None:
        raise SystemExit(f"{flag} is required")
    return default


import population_registry as registry                    # noqa: E402

REGISTRY = _from_argv("--registry", "r30")
STRATUM = _from_argv("--stratum")
PREREG = METRICS / f"{REGISTRY}_preregistration.json"
SPEC = registry.spec(REGISTRY, STRATUM)
DESIGN = METRICS / SPEC["file"]
ROWS = METRICS / f"{REGISTRY}_{STRATUM}_rows.json"
OUTPUT = METRICS / f"{REGISTRY}_{STRATUM}_convergence.json"
SMOKE_OUTPUT = METRICS / f"{REGISTRY}_{STRATUM}_convergence_smoke.json"
TREE = f"stratum_{STRATUM}_convergence"
CASE_ROOT = METRICS / "r11_cases" / TREE
RAW_ROOT = METRICS / "r11_raw" / TREE

os.environ["R11_TREE"] = TREE
os.environ["R11_CORRECTED"] = str(ROWS)
os.environ["R11_OUTPUT"] = str(OUTPUT)
os.environ["R11_SMOKE_OUTPUT"] = str(SMOKE_OUTPUT)

import rev11_designB_convergence as dB      # noqa: E402
import rev11_full_convergence as fc         # noqa: E402


# The two policies the budget ladder reads. The four schedule policies of the
# convergence study are not re-propagated per stratum; r30_preregistration.json
# declares that reduction and what it does not license. Patched at module level
# so that task building and the summary agree, and read from the task payload in
# every worker, so a child that never saw this patch still writes the right run.
fc.POLICIES = ("truth", "fixed_critical")
fc.COMPARED = ("fixed_critical",)


def install() -> None:
    dB.DESIGN_B = DESIGN
    dB.ROWS = ROWS
    fc.CORRECTED = ROWS
    fc.OUTPUT = OUTPUT
    fc.SMOKE_OUTPUT = SMOKE_OUTPUT
    fc.CASE_ROOT = CASE_ROOT
    fc.RAW_ROOT = RAW_ROOT
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)


def guard() -> None:
    if not DESIGN.exists():
        raise SystemExit(f"{DESIGN.name} missing")
    import rev10_sobol_confirmatory as base
    sha = base.file_hash(DESIGN)
    d = json.loads(DESIGN.read_text(encoding="utf-8"))
    if d["design_sha256"] != SPEC["design_sha256"]:
        raise SystemExit(f"{DESIGN.name} does not match the registered design "
                         f"hash; refusing to propagate an unfrozen design "
                         f"(file {sha[:16]})")
    if len(d.get("orbits", [])) != 64:
        raise SystemExit("a stratum is 64 orbits")
    if os.environ.get("R11_TREE") != TREE:
        raise SystemExit("output tree is not this stratum; refusing to write")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--registry", default="r30")
    ap.add_argument("command", choices=("prepass", "run", "status"))
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--deadline")
    a = ap.parse_args()

    guard()
    install()

    if a.command == "prepass":
        return dB.prepass(a.workers)
    if a.command == "status":
        if not OUTPUT.exists():
            print(f"{STRATUM}: no convergence output yet")
            if ROWS.exists():
                n = len(json.loads(ROWS.read_text(encoding="utf-8"))["rows"])
                print(f"{STRATUM}: prepass done for {n} orbits")
            return 0
        d = json.loads(OUTPUT.read_text(encoding="utf-8"))
        print(json.dumps({"stratum": STRATUM, "complete": d.get("complete"),
                          "orbits": len(d.get("rows", [])),
                          "failures": len(d.get("failures", []))}, indent=2))
        return 0

    if not ROWS.exists():
        print(f"{STRATUM}: prepass has not run; doing it first")
        rc = dB.prepass(a.workers)
        if rc != 0 or not ROWS.exists():
            return rc or 1
    return fc.run(False, fc.parse_deadline(a.deadline), a.workers)


if __name__ == "__main__":
    raise SystemExit(main())
