#!/usr/bin/env python3
"""One command that says whether this manuscript is fit to submit.

The pieces already existed; what did not was an order. That order matters and
was learned the hard way: the R10 manifest hashes main.pdf and supplement.pdf,
so any edit at all makes it stale, and re-sealing before recompiling seals the
previous build. Every step below therefore runs after the thing it depends on:

    compile  ->  seal  ->  mirror  ->  verify  ->  read

  compile   both documents, twice each, so cross-document references settle
  seal      re-stamp the manifest that owns the compiled documents
  mirror    refresh the curated archive copies and their checksums
  verify    the manifest integrity gate: digests, partition, ownership,
            twin agreement, derivative freshness
  read      the four manuscript checkers and the claims ledger, which ask a
            different question from the manifests: not whether the records are
            intact, but whether the sentences still match them

Exit status is zero only when every stage passes. Nothing here is a substitute
for reading the paper; it is a substitute for remembering nine commands and the
order they go in.

Usage:
    python submission_gate.py                full gate
    python submission_gate.py --no-compile   skip LaTeX (already built)
    python submission_gate.py --no-seal      verify only, change nothing
    python submission_gate.py --quick        skip the slow integrity check
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# (label, argv, cwd, writes)  -- `writes` marks a stage that changes the repo,
# so --no-seal can drop exactly those and leave a read-only audit.
def stages(a) -> list[tuple[str, list[str], Path, bool]]:
    py = sys.executable
    out: list[tuple[str, list[str], Path, bool]] = []

    if a.compile:
        for doc in ("supplement", "main"):
            for pas in (1, 2):
                out.append((f"compile {doc} (pass {pas})",
                            ["latexmk", "-pdf", "-interaction=nonstopmode",
                             f"{doc}.tex"], ROOT, True))

    if a.seal:
        out.append(("seal the manuscript manifest",
                    [py, "rev10_finalize_manifest.py"], HERE, True))
        out.append(("mirror the evidence archive",
                    [py, "populate_evidence_archive.py"], HERE, True))

    if not a.quick:
        out.append(("verify manifest integrity",
                    [py, "check_manifest_integrity.py"], HERE, False))

    out += [
        ("check consistency", [py, "check_consistency.py"], HERE, False),
        ("check labels", [py, "check_labels.py"], HERE, False),
        ("check assets", [py, "check_assets.py"], HERE, False),
        ("check protocol constants",
         [py, "check_protocol_constants.py"], HERE, False),
        ("check the claims ledger", [py, "claims_ledger.py", "--quiet"],
         HERE, False),
        ("self-test the claims ledger", [py, "test_claims_ledger.py"],
         HERE, False),
    ]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-compile", dest="compile", action="store_false")
    ap.add_argument("--no-seal", dest="seal", action="store_false")
    ap.add_argument("--quick", action="store_true",
                    help="skip the integrity check, which hashes the archive")
    a = ap.parse_args()

    todo = stages(a)
    results = []
    print(f"submission gate: {len(todo)} stages\n")
    for label, argv, cwd, _writes in todo:
        t0 = time.time()
        print(f"  ... {label}", flush=True)
        try:
            p = subprocess.run(argv, cwd=cwd, capture_output=True, text=True)
            rc, tail = p.returncode, (p.stdout or p.stderr).strip().splitlines()
        except FileNotFoundError as exc:
            rc, tail = 127, [str(exc)]
        results.append((label, rc, time.time() - t0,
                        tail[-1] if tail else ""))
        if rc != 0:
            # a failing stage is shown in full: the point of the gate is the
            # diagnosis, not the verdict
            print(f"      FAILED (rc={rc})")
            for line in (tail[-25:] if tail else []):
                print(f"      {line}")

    print("\n" + "-" * 72)
    worst = 0
    for label, rc, secs, last in results:
        mark = "ok  " if rc == 0 else "FAIL"
        worst = max(worst, rc)
        print(f"  [{mark}] {label:<38s} {secs:6.1f}s  {last[:60]}")
    print("-" * 72)
    if worst == 0:
        print("\nall stages passed: the records, the manifests and the "
              "sentences agree.")
        return 0
    failed = [l for l, rc, *_ in results if rc != 0]
    print(f"\n{len(failed)} stage(s) failed: {', '.join(failed)}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
