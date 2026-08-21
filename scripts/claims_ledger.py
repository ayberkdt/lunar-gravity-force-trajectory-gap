#!/usr/bin/env python3
"""A checkable ledger of the manuscript's quantitative claims.

Every number this paper prints comes from a sealed record, but until now the
link between the two lived in whoever last checked it. This file makes the link
an object: each claim names the record it derives from, the derivation, the
value it expects, and the words it appears as in the manuscript. The checker
recomputes the value from the record and compares it against both.

Three failure modes are distinguished, because they need different responses:

  FAIL   the record no longer yields the claimed value. The manuscript is
         wrong, or the derivation is.
  STALE  the record's digest has changed since the claim was pinned. The claim
         is not refuted, it is unverified: nobody has confirmed it still holds.
         This is the state a new campaign puts its neighbours in, and it is why
         a claim is pinned by digest rather than by file name.
  ABSENT the words the claim says it appears as are not in the manuscript. The
         number was edited, moved or cut without the ledger being told.

A claim that passes all three is one where the record, the arithmetic and the
printed sentence agree. That is the whole guarantee; it is not a proof that the
claim is interesting or that the experiment was well designed.

Adding a claim is deliberately cheap: append an entry to the ledger with a
derivation drawn from the vocabulary below, run with --pin to stamp the current
digest, and commit both. Adding a *derivation kind* is deliberately less cheap,
because every kind is a new way for the ledger to be wrong.

Usage:
    python claims_ledger.py                  check every claim
    python claims_ledger.py --id ladder.*    check a subset, glob on the id
    python claims_ledger.py --pin            restamp digests after a rerun
    python claims_ledger.py --pin --id x.y   restamp one claim
    python claims_ledger.py --list           show the ledger without checking
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import re
import sys
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
LEDGER = METRICS / "claims_ledger.json"


# ---------------------------------------------------------------- utilities

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def load(name: str) -> dict:
    p = METRICS / name
    if not p.exists():
        raise FileNotFoundError(name)
    return json.loads(p.read_text(encoding="utf-8"))


def dig(obj, path):
    """Walk a record.

    A string path is dotted, with list indices in brackets: a.b[0].c. Several
    record keys contain dots of their own (a regime-map cell is keyed
    "OEU|1.00"), so a path may also be given as a list of literal segments,
    which is unambiguous and is what those claims use.
    """
    parts = (path if isinstance(path, list)
             else re.findall(r"[^.\[\]]+|\[\d+\]", path))
    cur = obj
    for part in parts:
        if isinstance(part, str) and part.startswith("[") and part.endswith("]"):
            cur = cur[int(part[1:-1])]
        elif isinstance(part, int):
            cur = cur[part]
        else:
            cur = cur[part]
    return cur


def compiled_chapters() -> list:
    """The chapter files the two documents actually \\input, in order.

    Globbing chapters/*.tex would also read the sections that were archived
    out of the supplement; a claim quoting one of those would pass the wording
    check while being absent from both compiled documents.
    """
    seen, out = set(), []
    for doc in ("main.tex", "supplement.tex"):
        src = (ROOT / doc).read_text(encoding="utf-8")
        src = re.sub(r"(?<!\\)%.*", "", src)
        for name in re.findall(r"\\input\{chapters/([A-Za-z0-9_\-]+)\}", src):
            p = ROOT / "chapters" / f"{name}.tex"
            if p.exists() and p not in seen:
                seen.add(p)
                out.append(p)
    return out


def manuscript_text() -> str:
    """Every compiled .tex, comments stripped and whitespace collapsed.

    Collapsing whitespace is what lets a claim quote a phrase that the source
    happens to wrap across lines, which is most of them.
    """
    out = [p.read_text(encoding="utf-8") for p in compiled_chapters()]
    for extra in ("main.tex", "supplement.tex"):
        out.append((ROOT / extra).read_text(encoding="utf-8"))
    s = "\n".join(out)
    s = re.sub(r"(?<!\\)%.*", "", s)
    return " ".join(s.split())


def macro_values(ledger: dict) -> dict:
    """What claims_ledger_macros.tex must contain for this ledger."""
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
             6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    n = len(ledger["claims"])
    sources = len({c["source"] for c in ledger["claims"]})
    return {"CLLclaims": str(n),
            "CLLsources": str(words.get(sources, sources)),
            "CLLstates": "three"}


def check_macros(ledger: dict) -> list:
    """Is the generated macro file still the one this ledger implies?

    The integrity checker excludes this pair from its modification-time rule on
    the grounds that this script checks it instead, so this script has to.
    """
    out = METRICS / "claims_ledger_macros.tex"
    want = macro_values(ledger)
    if not out.exists():
        return [f"{out.name} is missing"]
    have = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{([^}]*)\}",
                           out.read_text(encoding="utf-8")))
    return [f"{out.name}: \\{k} is {have.get(k, 'absent')!r}, ledger implies {v!r}"
            for k, v in want.items() if have.get(k) != v]


# ------------------------------------------------------------- derivations
#
# Each takes the loaded record and the claim's `check` block, and returns a
# number, a string or a list. Keep them total and boring: a derivation that
# raises on a shape it did not expect is better than one that guesses.

def _cells(rec, level_min=None, level_max=None, budgets=None):
    """Flatten a verdict record's by_budget/levels into (beta, level, cell)."""
    out = []
    for beta, blk in rec["by_budget"].items():
        if budgets and beta not in budgets:
            continue
        for lvl, cell in blk["levels"].items():
            f = float(lvl)
            if level_min is not None and f < level_min:
                continue
            if level_max is not None and f > level_max:
                continue
            out.append((beta, f, cell))
    return out


def d_path(rec, c):
    return dig(rec, c["path"])


def d_len(rec, c):
    return len(dig(rec, c["path"]))


def d_count_cells(rec, c):
    sel = _cells(rec, c.get("level_min"), c.get("level_max"), c.get("budgets"))
    where = c.get("where")
    if not where:
        return len(sel)
    return sum(1 for _, _, cell in sel
               if dig(cell, where["path"]) == where["equals"])


def d_cell_field(rec, c):
    for beta, lvl, cell in _cells(rec):
        if beta == c["beta"] and lvl == float(c["level"]):
            return dig(cell, c["path"])
    raise KeyError(f"no cell at beta={c['beta']} level={c['level']}")


def d_cell_tally(rec, c):
    for beta, lvl, cell in _cells(rec):
        if beta == c["beta"] and lvl == float(c["level"]):
            base = dig(cell, c.get("under", "uncapped"))
            return f"{base['radial']}-{base['fixed']}"
    raise KeyError(f"no cell at beta={c['beta']} level={c['level']}")


def d_aggregate_cells(rec, c):
    sel = _cells(rec, c.get("level_min"), c.get("level_max"), c.get("budgets"))
    vals = []
    for _, _, cell in sel:
        try:
            v = dig(cell, c["path"])
        except (KeyError, IndexError):
            continue
        if v is not None:
            vals.append(float(v))
    if not vals:
        raise ValueError("no values selected")
    return {"median": median, "max": max, "min": min,
            "sum": sum, "count": len}[c.get("agg", "median")](vals)


def d_aggregate_rows(rec, c):
    """Aggregate a field over a record's rows, optionally filtered."""
    rows = dig(rec, c.get("rows", "rows"))
    where = c.get("where")
    vals = []
    for r in rows:
        if where and dig(r, where["path"]) != where["equals"]:
            continue
        v = dig(r, c["path"])
        if v is not None:
            vals.append(float(v))
    if not vals:
        raise ValueError("no values selected")
    return {"median": median, "max": max, "min": min,
            "sum": sum, "count": len}[c.get("agg", "median")](vals)


def d_ratio_rows(rec, c):
    """Aggregate the ratio of two fields over a record's rows.

    Added for the R68 claim that at equal measured time the matched constant
    degree runs at a median multiple of the critical degree: both quantities
    are per-row fields and the ratio is not stored, so aggregate_rows cannot
    reach it without the record being rewritten.
    """
    rows = dig(rec, c.get("rows", "rows"))
    where = c.get("where")
    vals = []
    for r in rows:
        if where and dig(r, where["path"]) != where["equals"]:
            continue
        num, den = dig(r, c["numerator"]), dig(r, c["denominator"])
        if num is None or den in (None, 0):
            continue
        vals.append(float(num) / float(den))
    if not vals:
        raise ValueError("no values selected")
    return {"median": median, "max": max, "min": min,
            "count": len}[c.get("agg", "median")](vals)


DERIVATIONS = {
    "path": d_path,
    "ratio_rows": d_ratio_rows,
    "len": d_len,
    "count_cells": d_count_cells,
    "cell_field": d_cell_field,
    "cell_tally": d_cell_tally,
    "aggregate_cells": d_aggregate_cells,
    "aggregate_rows": d_aggregate_rows,
}


# ------------------------------------------------------------------ checker

def evaluate(claim: dict) -> dict:
    """Return {state, actual, detail} for one claim."""
    src = claim["source"]
    try:
        rec = load(src)
    except FileNotFoundError:
        return {"state": "FAIL", "actual": None,
                "detail": f"source record {src} is missing"}

    pinned = claim.get("source_sha256")
    now = sha256(METRICS / src)
    stale = bool(pinned) and pinned != now

    check = claim["check"]
    kind = check.get("kind", "path")
    if kind not in DERIVATIONS:
        return {"state": "FAIL", "actual": None,
                "detail": f"unknown derivation kind {kind!r}"}
    try:
        actual = DERIVATIONS[kind](rec, check)
    except Exception as exc:                                # noqa: BLE001
        return {"state": "FAIL", "actual": None,
                "detail": f"{type(exc).__name__}: {exc}"}

    expect = claim["expect"]
    tol = claim.get("tol", 0)
    if isinstance(expect, (int, float)) and isinstance(actual, (int, float)) \
            and not isinstance(expect, bool):
        agree = abs(actual - expect) <= tol * max(abs(actual), 1e-12) \
            if tol else actual == expect
    else:
        agree = actual == expect

    if not agree:
        return {"state": "FAIL", "actual": actual,
                "detail": f"expected {expect!r}, record gives {actual!r}"}
    if stale:
        return {"state": "STALE", "actual": actual,
                "detail": f"{src} changed since this claim was pinned"}
    return {"state": "PASS", "actual": actual, "detail": ""}


def check_wording(claim: dict, text: str) -> str | None:
    phrase = claim.get("appears_as")
    if not phrase:
        return None
    return None if " ".join(phrase.split()) in text else phrase


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", default="*", help="glob over claim ids")
    ap.add_argument("--pin", action="store_true",
                    help="restamp source digests for the selected claims")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--macros", action="store_true",
                    help="write the counts the supplement quotes")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    if not LEDGER.exists():
        print(f"[abort] {LEDGER.name} does not exist")
        return 2
    ledger = json.loads(LEDGER.read_text(encoding="utf-8"))
    claims = [c for c in ledger["claims"] if fnmatch.fnmatch(c["id"], a.id)]
    if not claims:
        print(f"no claim matches {a.id!r}")
        return 2

    if a.list:
        for c in claims:
            print(f"{c['id']:<44s} {c['source']:<34s} {c['claim'][:60]}")
        print(f"\n{len(claims)} claims")
        return 0

    if a.pin:
        for c in claims:
            p = METRICS / c["source"]
            if p.exists():
                c["source_sha256"] = sha256(p)
        LEDGER.write_text(json.dumps(ledger, indent=2), encoding="utf-8")
        print(f"[pinned] {len(claims)} claims restamped in {LEDGER.name}")
        return 0

    text = manuscript_text()
    tally = {"PASS": 0, "FAIL": 0, "STALE": 0, "ABSENT": 0}
    problems = []
    for c in claims:
        res = evaluate(c)
        missing = check_wording(c, text)
        state = res["state"]
        if state == "PASS" and missing:
            state = "ABSENT"
            res["detail"] = f"manuscript does not contain {missing!r}"
        tally[state] += 1
        if state != "PASS":
            problems.append((state, c, res))
        if not a.quiet:
            mark = {"PASS": "ok  ", "FAIL": "FAIL", "STALE": "STALE",
                    "ABSENT": "ABS "}[state]
            val = res["actual"]
            shown = f"{val:.6g}" if isinstance(val, float) else str(val)
            print(f"  [{mark}] {c['id']:<44s} {shown}")

    print(f"\n{tally['PASS']} pass, {tally['FAIL']} fail, "
          f"{tally['STALE']} stale, {tally['ABSENT']} absent "
          f"({len(claims)} claims)")

    if a.macros:
        # The supplement quotes how many claims the ledger carries. Typing that
        # into a sentence would make it the one number in this paper with no
        # record behind it, so it is generated like every other count.
        want = macro_values(ledger)
        n = len(ledger["claims"])
        sources = len({c["source"] for c in ledger["claims"]})
        body = ["% generated by claims_ledger.py --macros; do not edit"]
        body += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in want.items()]
        out = METRICS / "claims_ledger_macros.tex"
        out.write_text("\n".join(body) + "\n", encoding="utf-8")
        print(f"[written] {out.name}: {n} claims over {sources} records")
    for state, c, res in problems:
        print(f"\n[{state}] {c['id']}")
        print(f"    claim  : {c['claim']}")
        print(f"    source : {c['source']}")
        print(f"    detail : {res['detail']}")
        if c.get("appears_in"):
            print(f"    printed: {', '.join(c['appears_in'])}")
    stale_macros = [] if a.macros else check_macros(ledger)
    for line in stale_macros:
        print(f"\n[MACROS] {line}")
        print("    fix    : rerun claims_ledger.py --macros")
    return 1 if (tally["FAIL"] or tally["ABSENT"] or tally["STALE"]
                 or stale_macros) else 0


if __name__ == "__main__":
    raise SystemExit(main())
