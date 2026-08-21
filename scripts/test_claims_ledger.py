#!/usr/bin/env python3
"""Prove the claims ledger fails when it should.

A checker nobody has seen fail is a checker nobody should trust. This injects
each of the three defects the ledger is built to catch into a scratch copy and
asserts that the state reported is the right one. It never touches the real
ledger.

Usage:  python test_claims_ledger.py
"""
from __future__ import annotations

import copy
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import claims_ledger as CL                                   # noqa: E402


def state_of(claim, text):
    res = CL.evaluate(claim)
    if res["state"] == "PASS" and CL.check_wording(claim, text):
        return "ABSENT"
    return res["state"]


def main() -> int:
    ledger = json.loads(CL.LEDGER.read_text(encoding="utf-8"))
    text = CL.manuscript_text()
    by_id = {c["id"]: c for c in ledger["claims"]}
    failures = []

    def expect(label, got, want):
        mark = "ok  " if got == want else "FAIL"
        print(f"  [{mark}] {label}: {got}")
        if got != want:
            failures.append(label)

    # 1. a claim whose record no longer yields the pinned value must FAIL
    c = copy.deepcopy(by_id["ladder.replication.agree"])
    c["expect"] = c["expect"] - 1
    expect("a wrong expected value is caught", state_of(c, text), "FAIL")

    # 2. a claim whose source has moved since it was pinned must go STALE,
    #    not pass quietly and not fail
    c = copy.deepcopy(by_id["ladder.uncapped.a.cells"])
    c["source_sha256"] = "0" * 64
    expect("a moved source is reported stale", state_of(c, text), "STALE")

    # 3. a claim whose words are no longer in the manuscript must be ABSENT
    c = copy.deepcopy(by_id["ladder.probe.b.demand"])
    c["appears_as"] = "a phrase this manuscript does not contain anywhere"
    expect("a vanished sentence is caught", state_of(c, text), "ABSENT")

    # 4. a claim naming a record that does not exist must FAIL, not crash
    c = copy.deepcopy(by_id["ladder.uncapped.a.cells"])
    c["source"] = "r99_does_not_exist.json"
    expect("a missing record is caught", state_of(c, text), "FAIL")

    # 5. an unknown derivation must FAIL rather than be skipped
    c = copy.deepcopy(by_id["ladder.uncapped.a.cells"])
    c["check"] = {"kind": "no_such_derivation"}
    expect("an unknown derivation is caught", state_of(c, text), "FAIL")

    # 6. and the real ledger must still pass, or the tests above are moot
    bad = [c["id"] for c in ledger["claims"] if state_of(c, text) != "PASS"]
    expect("the live ledger passes", bad or "all claims pass", "all claims pass")

    print()
    if failures:
        print(f"{len(failures)} self-test(s) failed: {', '.join(failures)}")
        return 1
    print(f"all self-tests passed over {len(ledger['claims'])} live claims")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
