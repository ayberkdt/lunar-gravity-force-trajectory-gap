"""Check the campaign machinery against its own invariants.

Written after a night in which three mechanical errors reached working code:
a shell heredoc turned LaTeX escapes into control characters, a regex reorder
duplicated supervisor stages, and a probe with an unused tag started a full
re-scoring run. Each was caught, but only after it ran. The invariants below
are the ones that would have caught them before.

Exit status is non-zero if any check fails, so this can gate a queue.

Usage:
    python revJ_selfcheck.py
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
METRICS = HERE.parent / "metrics"

MODULES = ["revJ_common.py", "revJ1_crossfield.py", "revJ2_fullforce.py",
           "revJ3_tolerance.py", "revJ_fidelity.py", "revJ_compare.py",
           "revJ_tables.py", "revJ_budget_grid.py", "revJ_ladder_primary.py",
           "revJ_supervisor.py"]

failures: list[str] = []
notes: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f" -- {detail}" if detail
                                                   else ""))
    if not ok:
        failures.append(f"{name}: {detail}")


def main() -> int:
    print("== source ==")
    for m in MODULES:
        p = HERE / m
        if not p.is_file():
            check(f"{m} present", False, "missing")
            continue
        text = p.read_text(encoding="utf-8")
        try:
            ast.parse(text)
            parsed = True
        except SyntaxError as exc:
            parsed = False
            check(f"{m} parses", False, str(exc))
        if parsed:
            check(f"{m} parses", True)
        # Control characters are the signature of a mangled shell edit.
        bad = {c for c in text if ord(c) < 32 and c not in "\n\t"}
        check(f"{m} has no stray control characters", not bad,
              f"found {[hex(ord(c)) for c in bad]}" if bad else "")

    print("\n== supervisor stage list ==")
    sys.path.insert(0, str(HERE))
    import revJ_supervisor as S

    names = [s[0] for s in S.STAGES]
    dupes = {n for n in names if names.count(n) > 1}
    check("stage names unique", not dupes, f"duplicated {sorted(dupes)}")

    # Two stages may legitimately share an artifact only if they are the
    # deliberately re-entrant final pair.
    arts: dict[str, list[str]] = {}
    for n, _argv, art, _need, _env in S.STAGES:
        arts.setdefault(art, []).append(n)
    shared = {a: v for a, v in arts.items() if len(v) > 1
              and a != "rJ_never_satisfied.tex"}
    check("artifacts identify one stage each", not shared, f"{shared}")

    for n, argv, _art, need, env in S.STAGES:
        if not (HERE / argv[0]).is_file():
            check(f"{n} script exists", False, argv[0])
        if need <= 0:
            check(f"{n} declares a positive duration", False, str(need))
        for k, v in env.items():
            if not k.startswith("JCAMP_"):
                check(f"{n} env key namespaced", False, k)
            if not isinstance(v, str):
                check(f"{n} env value is a string", False, f"{k}={v!r}")
    check("every stage has a script, a duration and namespaced env", True)

    print("\n== generated LaTeX ==")
    for tex in sorted(METRICS.glob("rJ_*.tex")):
        text = tex.read_text(encoding="utf-8")
        bad = {c for c in text if ord(c) < 32 and c != "\n"}
        check(f"{tex.name} clean", not bad,
              f"control chars {[hex(ord(c)) for c in bad]}" if bad else "")
        if "tabular" in text:
            opens = text.count(r"\begin{tabular}")
            closes = text.count(r"\end{tabular}")
            check(f"{tex.name} tabular balanced", opens == closes,
                  f"{opens} begin / {closes} end")

    print("\n== records ==")
    for rec in sorted(METRICS.glob("rJ*.json")):
        if rec.stat().st_size == 0:
            check(f"{rec.name} non-empty", False, "zero bytes")
            continue
        try:
            payload = json.loads(rec.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            check(f"{rec.name} parses", False, str(exc))
            continue
        if "complete" in payload and not payload["complete"]:
            notes.append(f"{rec.name}: complete=false (in progress or partial)")

    # A scored record must not claim more resolved comparisons than orbits.
    for score in sorted(METRICS.glob("rJ[123]_score*.json")):
        d = json.loads(score.read_text(encoding="utf-8"))
        c = d.get("counts", {})
        n = c.get("orbits")
        for key in ("resolved_tighter", "resolved_in_both_dynamics",
                    "reversal_orbits", "radial_wins_force",
                    "radial_loses_trajectory"):
            v = c.get(key)
            if isinstance(v, int) and isinstance(n, int) and v > n:
                check(f"{score.name} {key} <= orbits", False, f"{v} > {n}")
    check("scored counts within their populations", True)

    print("\n== manuscript macros ==")
    numbers = METRICS / "rJ_numbers.tex"
    if numbers.is_file():
        # \newcommand{\Name}{value} -- the captured token carries the leading
        # backslash and the usage scan does not, so it is stripped here. The
        # first version of this check compared the two forms and reported every
        # macro as undefined while the document compiled cleanly.
        defined = {line.split("{")[1].split("}")[0].lstrip("\\")
                   for line in numbers.read_text(encoding="utf-8").splitlines()
                   if line.startswith(r"\newcommand")}
        used: set[str] = set()
        for ch in (HERE.parent / "chapters").glob("*.tex"):
            t = ch.read_text(encoding="utf-8")
            for token in defined:
                if "\\" + token in t:
                    used.add(token)
        missing = set()
        for ch in (HERE.parent / "chapters").glob("*.tex"):
            t = ch.read_text(encoding="utf-8")
            for word in ("Jcross", "Jfull", "Jtol", "Jxover", "Jfid",
                         "Jladder"):
                idx = 0
                while True:
                    idx = t.find("\\" + word, idx)
                    if idx < 0:
                        break
                    end = idx + 1
                    while end < len(t) and t[end].isalpha():
                        end += 1
                    name = t[idx + 1:end]
                    if name not in defined:
                        missing.add(name)
                    idx = end
        check("every campaign macro used in the text is defined", not missing,
              f"undefined {sorted(missing)}")
        notes.append(f"{len(defined)} macros defined, {len(used)} used")

    print("\n== summary ==")
    for n in notes:
        print(f"note  {n}")
    if failures:
        print(f"\n{len(failures)} FAILED")
        for f in failures:
            print(f"  {f}")
        return 1
    print("\nall invariants hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
