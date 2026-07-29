"""Cross-document consistency sweep for stale text.

The manuscript has accumulated several generations of edits, and the failures
that survive are not LaTeX errors -- they compile fine -- but counts and
descriptions that were true one revision ago. This checks the classes that
have actually gone stale here at least once:

  * experiment-family counts (O1..On) against the highest O number defined
  * campaign-manifest counts against the manifests that exist on disk
  * "one/two population" and "design A only" phrasing against what ran
  * figure and table counts quoted in prose
  * spelling variant drift
  * references to metrics files that do not exist

Exit code is non-zero if any blocking problem is found, so it can gate a
build. Usage:  python check_consistency.py
"""

from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CH = ROOT / "chapters"
METRICS = ROOT / "metrics"

# Spelled-out counts the prose uses. Kept in one place because the previous
# per-check dictionaries stopped at "ten" and "twenty-eight", so a count that
# grew past them passed by not being recognised at all.
NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
    "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
    "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15,
    "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
    "twenty": 20, "thirty": 30, "forty": 40,
}
NUMBER_WORDS.update({
    f"{tens}-{unit}": NUMBER_WORDS[tens] + NUMBER_WORDS[unit]
    for tens in ("twenty", "thirty", "forty")
    for unit in ("one", "two", "three", "four", "five", "six", "seven",
                 "eight", "nine")
})

TEX = sorted(glob.glob(str(CH / "*.tex"))) + [str(ROOT / "main.tex"),
                                              str(ROOT / "supplement.tex")]


def read(p) -> str:
    return Path(p).read_text(encoding="utf-8")


def all_text() -> dict:
    return {Path(p).name: read(p) for p in TEX}


def check_experiment_counts(texts, problems, notes):
    defined = set()
    for s in texts.values():
        defined |= {int(m) for m in re.findall(r'\(O(\d+)\)~\\emph', s)}
    if not defined:
        return
    hi = max(defined)
    notes.append(f"orbit-level experiments defined: O1..O{hi} "
                 f"({len(defined)} distinct)")
    # every prose claim about the range must match
    for name, s in texts.items():
        for m in re.finditer(r'O1--O(\d+)', s):
            if int(m.group(1)) != hi:
                problems.append(f"{name}: claims O1--O{m.group(1)} but "
                                f"O{hi} is defined")
        for m in re.finditer(r'([a-z]+(?:-[a-z]+)?|\d+) orbit-level experiments',
                             s):
            word = m.group(1)
            got = int(word) if word.isdigit() else NUMBER_WORDS.get(word)
            if got is not None and got != hi:
                problems.append(f"{name}: says {word} orbit-level experiments "
                                f"but O1..O{hi} are defined")


def check_manifest_counts(texts, problems, notes):
    # r?? rather than r1?: the campaign numbering passed 19 and the narrower
    # glob silently stopped counting, which is how a stale count survives.
    paths = sorted(glob.glob(str(METRICS / "r[0-9]*_final_experiment_manifest.json")),
                   key=lambda p: int(re.match(r'r(\d+)', Path(p).name).group(1)))
    on_disk = [Path(p).name for p in paths]
    campaign = [x.split('_')[0] for x in on_disk]
    n = len(on_disk)
    notes.append(f"campaign manifests on disk: {n} ({', '.join(campaign)})")
    # the range the data-availability statement quotes must be the real one
    lo, hi = campaign[0], campaign[-1]
    for name, s in texts.items():
        for m in re.finditer(r'[Rr](\d+) through [Rr](\d+)', s):
            if (f"r{m.group(1)}", f"r{m.group(2)}") != (lo, hi):
                problems.append(f"{name}: quotes R{m.group(1)} through "
                                f"R{m.group(2)} but {lo}..{hi} exist")
        for pat in (r'(\w+) campaign manifests', r'all (\w+) manifests',
                    r'the (\w+) manifests', r'the (\w+) can be read as a'):
            for m in re.finditer(pat, s):
                got = NUMBER_WORDS.get(m.group(1).lower())
                if got is not None and got != n:
                    problems.append(f"{name}: says {m.group(1)} manifests "
                                    f"but {n} exist")
    # every manifest must parse and carry the digest that seals it
    for p in paths:
        try:
            d = json.loads(Path(p).read_text(encoding="utf-8"))
        except Exception as exc:                                   # noqa: BLE001
            problems.append(f"{Path(p).name}: unreadable ({exc})")
            continue
        # R10 predates the rolled-up self-digest and seals itself with a flat
        # per-entry index instead; either form counts as sealed.
        sealed = d.get("manifest_sha256") or all(
            e.get("sha256") for e in d.get("entries", [{}]))
        if not sealed:
            problems.append(f"{Path(p).name}: neither a manifest_sha256 nor a "
                            f"fully hashed entry index")


def check_population_language(texts, problems, notes):
    """Design B ran for the span sweep; nothing may still say one design."""
    if not (METRICS / "r18_span_sweep_B_beta_1.00.json").exists():
        return
    for name, s in texts.items():
        for pat in (r'on one population', r'on design A only',
                    r'one coverage design'):
            if re.search(pat, s, re.I):
                problems.append(f"{name}: '{pat}' but design B has run")


def check_missing_inputs(texts, problems, notes):
    used = set()
    for s in texts.values():
        used |= set(re.findall(r'\\input\{(metrics/[^}]+)\}', s))
    missing = []
    for u in sorted(used):
        p = ROOT / u
        if not (p.exists() or p.with_suffix('.tex').exists()):
            missing.append(u)
    if missing:
        problems.append("missing \\input targets: " + ", ".join(missing))
    else:
        notes.append(f"all {len(used)} \\input targets present")


def check_spelling(texts, problems, notes):
    hits = []
    for name, s in texts.items():
        for m in re.finditer(r'\b\w*favou\w*', s):
            hits.append(f"{name}:{m.group(0)}")
        for m in re.finditer(r'\b\w*behaviou\w*', s):
            hits.append(f"{name}:{m.group(0)}")
    if hits:
        problems.append("British spellings: " + ", ".join(hits[:8]))
    else:
        notes.append("spelling variants consistent (US)")


def check_control_chars(texts, problems, notes):
    """Catch LaTeX commands mangled by a Python escape.

    Writing "\\textbf" as "\textbf" in a non-raw string yields a TAB followed by
    "extbf"; the same trap turns "\\ref" into CR + "ef" and "\\varepsilon" into
    a vertical tab. Each of those has happened here at least once, so all three
    shapes are checked: stray control characters, tabs anywhere in the source,
    and line fragments that look like a command missing its backslash.
    """
    bad = []
    # Byte-level first: Python's universal-newline translation turns a lone CR
    # into LF on text read, so "\r"+"ef" is invisible to a text-mode search.
    # Every instance found in this manuscript was of exactly that shape.
    ctrl_stems = {b'\r': (b'ef', b'ho', b'ight'),
                  b'\t': (b'extbf', b'extit', b'imes', b'ext{'),
                  b'\x0b': (b'arepsilon', b'arphi'),
                  b'\x0c': (b'rac',), b'\x08': (b'egin', b'ar{'),
                  b'\x07': (b'lpha',)}
    for p in TEX:
        raw = Path(p).read_bytes()
        for ctrl, stems in ctrl_stems.items():
            for stem in stems:
                if ctrl + stem in raw:
                    bad.append(f"{Path(p).name}: {ctrl!r}+{stem.decode()} "
                               f"(backslash eaten by a Python escape)")
    for name, s in texts.items():
        for ch in ('\x07', '\x08', '\x0b', '\x0c'):
            if ch in s:
                bad.append(f"{name}:U+{ord(ch):04X}")
        if '\t' in s:
            line = s[:s.index('\t')].count('\n') + 1
            bad.append(f"{name}:{line}: TAB in source (lost backslash?)")
    frag = (r'ef\*?\{|ewcommand|extbf|extit|emph\{|imes\b|abla|ho\b|'
            r'oindent|ewline|au\b|ightarrow')
    for name, s in texts.items():
        for m in re.finditer(rf'(?m)(?:^|(?<=[\s{{]))({frag})', s):
            # a real word like "the" must not trigger: require the fragment to
            # sit where a command would, i.e. preceded by whitespace or brace
            ctx = s[max(0, m.start() - 1):m.start()]
            if ctx and ctx not in ' \t\n{':
                continue
            bad.append(f"{name}: '{m.group(1)}' without backslash")
    if bad:
        problems.append("control-character damage: " + ", ".join(bad))
    else:
        notes.append("no control-character damage")


def main() -> int:
    texts = all_text()
    problems: list[str] = []
    notes: list[str] = []
    for fn in (check_experiment_counts, check_manifest_counts,
               check_population_language, check_missing_inputs,
               check_spelling, check_control_chars):
        fn(texts, problems, notes)
    for n in notes:
        print(f"  [ok] {n}")
    for p in problems:
        print(f"  [!!] {p}")
    print(f"blocking consistency problems: {len(problems)}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
