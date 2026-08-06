"""Generate the propagated base for design C, so a budget ladder can run on it.

A design is not a list of 64 points. Before any budget campaign can use it, the
population needs the same base designs A and B carry: an empirical prepass that
fixes each orbit's adopted truth degree, its critical degree and its work
degree, and then the convergence tree that holds the truth trajectories and the
comparator policies at both tolerance levels.

Nothing here is new machinery. The prepass and the convergence run are the
drivers that produced design B, pointed at the frozen design C population. The
one thing that has to be done exactly right is where the workers write:
rev11_full_convergence resolves its output tree from environment variables read
at import, precisely because a ProcessPoolExecutor on Windows re-imports it in
every child and a parent-only monkey-patch would leave those children writing
into design A's tree. The environment is therefore set here before the pool is
built, and the parent's already-imported globals are realigned to match, which
is the same two-step design B used.

The design is read from the frozen file and never regenerated, so the points
this propagates are the points that were hashed before any of it ran.

Usage:
    python rev26_designC_base.py prepass --workers 11
    python rev26_designC_base.py run --workers 11 --deadline "<iso>"
    python rev26_designC_base.py status
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

DESIGN_C = METRICS / "r26_sobolC_design_frozen.json"
ROWS = METRICS / "r26_designC_rows.json"
OUTPUT = METRICS / "r26_designC_convergence.json"
SMOKE_OUTPUT = METRICS / "r26_designC_convergence_smoke.json"
TREE = "designC_convergence"
CASE_ROOT = METRICS / "r11_cases" / TREE
RAW_ROOT = METRICS / "r11_raw" / TREE

# Must be set before rev11_full_convergence is imported anywhere in this
# process, so that both the parent and every spawned worker resolve the same
# tree. Importing the design-B driver below pulls fc in, so this comes first.
os.environ["R11_TREE"] = TREE
os.environ["R11_CORRECTED"] = str(ROWS)
os.environ["R11_OUTPUT"] = str(OUTPUT)
os.environ["R11_SMOKE_OUTPUT"] = str(SMOKE_OUTPUT)

import rev11_designB_convergence as dB      # noqa: E402
import rev11_full_convergence as fc         # noqa: E402


def install() -> None:
    """Point the validated machinery at design C.

    The prepass functions are pure: they take an orbit dict and return a
    record, reading only constants that are identical for every design, so
    rebinding the design-B module's paths in this process is enough for them.
    The convergence run is not pure, which is why its tree comes from the
    environment set above rather than from these rebindings.
    """
    dB.DESIGN_B = DESIGN_C
    dB.ROWS = ROWS
    fc.CORRECTED = ROWS
    fc.OUTPUT = OUTPUT
    fc.SMOKE_OUTPUT = SMOKE_OUTPUT
    fc.CASE_ROOT = CASE_ROOT
    fc.RAW_ROOT = RAW_ROOT
    CASE_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)


def guard() -> None:
    """Refuse to run against anything but the frozen, hashed design."""
    if not DESIGN_C.exists():
        raise SystemExit(f"{DESIGN_C.name} missing; run rev26_designC_freeze.py")
    d = json.loads(DESIGN_C.read_text(encoding="utf-8"))
    if d.get("family") != "sobolC" or len(d.get("orbits", [])) != 64:
        raise SystemExit("frozen design C is not a 64-orbit sobolC population")
    if os.environ.get("R11_TREE") != TREE:
        raise SystemExit("output tree is not design C; refusing to write")


def main() -> int:
    ap = argparse.ArgumentParser()
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
            print("design C: no convergence output yet")
            if ROWS.exists():
                n = len(json.loads(ROWS.read_text(encoding="utf-8"))["rows"])
                print(f"design C: prepass done for {n} orbits")
            return 0
        d = json.loads(OUTPUT.read_text(encoding="utf-8"))
        print(json.dumps({"complete": d.get("complete"),
                          "orbits": len(d.get("rows", [])),
                          "failures": len(d.get("failures", []))}, indent=2))
        return 0

    if not ROWS.exists():
        print("design C: prepass has not run; doing it first")
        rc = dB.prepass(a.workers)
        if rc != 0 or not ROWS.exists():
            return rc or 1
    return fc.run(False, fc.parse_deadline(a.deadline), a.workers)


if __name__ == "__main__":
    raise SystemExit(main())
