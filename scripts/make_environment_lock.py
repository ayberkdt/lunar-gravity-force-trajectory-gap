"""Emit the environment specifications and lock files the package promises.

The supplement's source-snapshot statement lists "environment specifications
and lock files" among the package contents. There were none, so that promise
was unbacked in exactly the way the manifest exhaustiveness claim was.

Two interpreters produced the results, not one, and a single lock file would
misrepresent that. The split is by date, not by kind of work: the interpreter
changed partway through the revision sequence and propagations ran under both.

  * py310 -- Python 3.10.20 with SciPy 1.15.3, the virtualenv carrying the
             Lunaris package. R10 through R13. Their case sidecars record this
             interpreter, so the archive itself says what ran them.
  * py312 -- Python 3.12.1 with SciPy 1.14.1 and Numba 0.63.1. The kernel
             timing, the field-level work, and every propagation from R14
             onward. R14 through R17 record it in their own scenario blocks;
             R18 through R23 record no version fields at all, and the bytecode
             cache tags in this directory are what places them here.

Rather than trust either the text or today's machine, this reads the versions
each environment reports now, reads the versions the archive recorded when the
runs happened, and writes both into the lock header with the comparison stated.
Drift is reported, not smoothed over: a lock file that hides it is worse than
no lock file. Exit code is non-zero when a recorded field and a live field
disagree, so the drift has to be looked at rather than shipped silently.

Interpreter locations are arguments, not constants, because they are paths on
one machine and mean nothing on another. Only versions go into the artifacts.

Usage:
  python make_environment_lock.py
  python make_environment_lock.py --orbit-python <path> --timing-python <path>
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

DEFAULT_PY310 = r"D:\Masaustu\LUNAR_SIMULATION\.venv\Scripts\python.exe"
DEFAULT_PY312 = sys.executable

# What each environment produced, and the direct dependencies a reader has to
# install to re-run it. The lock files carry the full closure; these are the
# names that belong in a human-readable spec.
ROLES = {
    "py310": {
        "what": ("R10 through R13: the confirmatory design, the "
                 "vector-tolerance convergence trees, the published-rule "
                 "benchmark and the resolution diagnosis"),
        "direct": ["numpy", "scipy", "lunaris"],
        "checked_against": "the R10 and R11 case sidecars",
    },
    "py312": {
        "what": ("the kernel timing and the field-level work, and every "
                 "propagation from R14 onward: the fixed-budget allocation "
                 "campaign, the audit response, the sixty-day extensions, the "
                 "span sweep and the R23 controls"),
        "direct": ["numpy", "scipy", "numba"],
        "checked_against": "the R14, R15 and R17 scenario blocks",
    },
}

PROBE = (
    "import json,platform,sys\n"
    "from importlib.metadata import distributions\n"
    "d=sorted({(x.metadata['Name'],x.version) for x in distributions()},"
    "key=lambda t:t[0].lower())\n"
    "print(json.dumps({'python':sys.version.split()[0],"
    "'python_full':sys.version,'platform':platform.platform(),"
    "'distributions':d}))\n"
)


def probe(interpreter: str) -> dict | None:
    """Ask an interpreter what it is and what is installed in it."""
    try:
        out = subprocess.run([interpreter, "-c", PROBE], capture_output=True,
                             text=True, timeout=180)
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"[warn] cannot probe {interpreter}: {exc}")
        return None
    if out.returncode != 0:
        print(f"[warn] probe failed for {interpreter}: "
              f"{out.stderr.strip()[:200]}")
        return None
    return json.loads(out.stdout.strip().splitlines()[-1])


def find_recorded(paths, keys=("python", "numpy", "scipy", "numba")) -> dict:
    """Pull the environment fields the archive recorded, wherever they sit.

    The sidecar schema changed across campaigns, so this walks the JSON rather
    than assuming a fixed location. Only the first value found for each key is
    kept, and the file it came from is reported alongside it.
    """
    found: dict[str, tuple[str, str]] = {}

    def walk(node, origin):
        if isinstance(node, dict):
            for k, v in node.items():
                if k in keys and isinstance(v, str) and k not in found:
                    version = re.match(r"[\d.]+(?:rc\d+)?", v.strip())
                    if version:
                        found[k] = (version.group(0), origin)
                walk(v, origin)
        elif isinstance(node, list):
            for v in node:
                walk(v, origin)

    for p in paths:
        if not p.is_file():
            continue
        try:
            walk(json.loads(p.read_text(encoding="utf-8")),
                 p.relative_to(ROOT).as_posix())
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        if all(k in found for k in keys):
            break
    return found


def py310_record_sources():
    """R10 and R11 sidecars, which carry the interpreter that ran them."""
    out = []
    for tree in ("r10_cases", "r11_cases"):
        base = METRICS / tree
        if base.exists():
            out.extend(sorted(base.rglob("*.json"))[:4])
    return out


def py312_record_sources():
    """R14, R15 and R17 record python, numpy and scipy in their scenario
    blocks; R8 is the only place numba is recorded at all."""
    return [METRICS / "r14_timing_budget.json",
            METRICS / "r15_deployable_calibration_A.json",
            METRICS / "r17_longarc60.json",
            METRICS / "r8_alpha_margin.json"]


def compare(live: dict, recorded: dict) -> tuple[list[str], list[str]]:
    """Field-by-field agreement between what ran and what runs now."""
    agree, drift = [], []
    installed = {n.lower(): v for n, v in live["distributions"]}
    installed["python"] = live["python"]
    for key, (was, origin) in sorted(recorded.items()):
        now = installed.get(key)
        if now is None:
            drift.append(f"{key}: recorded {was} ({origin}), not installed now")
        elif now == was:
            agree.append(f"{key} {was}")
        else:
            drift.append(f"{key}: recorded {was} ({origin}), now {now}")
    return agree, drift


def write_lock(name: str, live: dict, agree, drift) -> Path:
    role = ROLES[name]
    out = ROOT / f"requirements-{name}.lock"
    head = [
        f"# {name} environment -- {role['what']}",
        "#",
        f"# python {live['python_full'].splitlines()[0].strip()}",
        f"# platform {live['platform']}",
        f"# captured {datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')}",
        "#",
        f"# Cross-check against {role['checked_against']}:",
    ]
    head += [f"#   agrees: {a}" for a in agree] or ["#   agrees: (nothing recorded)"]
    head += [f"#   DRIFT:  {d}" for d in drift]
    if not drift:
        head.append("#   no drift: every recorded version is the installed version")
    head += ["#",
             "# Full closure of the environment as captured, not a minimal set.",
             ""]
    body = [f"{n}=={v}" for n, v in live["distributions"]]
    out.write_text("\n".join(head + body) + "\n", encoding="utf-8")
    return out


def write_spec(name: str, live: dict) -> Path:
    role = ROLES[name]
    out = ROOT / f"environment-{name}.yml"
    installed = {n.lower(): v for n, v in live["distributions"]}
    lines = [
        f"# {name} environment -- {role['what']}",
        f"# Companion lock with the full closure: requirements-{name}.lock",
        f"name: lunar-truncation-{name}",
        "channels:",
        "  - conda-forge",
        "dependencies:",
        f"  - python={live['python']}",
    ]
    for dep in role["direct"]:
        version = installed.get(dep)
        lines.append(f"  - {dep}={version}" if version
                     else f"  # - {dep}  (not installed in this environment)")
    if "lunaris" in role["direct"]:
        lines += [
            "# lunaris is not distributed on a channel: install the pinned",
            "# source snapshot recorded in the campaign manifests",
            "# (tag paper-truncation-v1.0, commit 27e9ab86ed61d623f78c4"
            "53ea2054348f1044c23).",
        ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--py310", default=DEFAULT_PY310)
    ap.add_argument("--py312", default=DEFAULT_PY312)
    args = ap.parse_args()

    interpreters = {"py310": args.py310, "py312": args.py312}
    sources = {"py310": py310_record_sources(),
               "py312": py312_record_sources()}

    any_drift, probed = False, 0
    for name, interpreter in interpreters.items():
        live = probe(interpreter)
        if live is None:
            print(f"[skip] {name}: interpreter not reachable, no lock written")
            any_drift = True
            continue
        probed += 1
        recorded = find_recorded(sources[name])
        agree, drift = compare(live, recorded)
        lock = write_lock(name, live, agree, drift)
        spec = write_spec(name, live)
        print(f"[{name}] python {live['python']}, "
              f"{len(live['distributions'])} distributions")
        print(f"  {lock.name}, {spec.name}")
        for a in agree:
            print(f"  [ok] {a} matches the archive")
        for d in drift:
            print(f"  [drift] {d}")
            any_drift = True

    if probed == 0:
        print("no interpreter could be probed; nothing written")
        return 1
    if any_drift:
        print("\nenvironment drift or an unreachable interpreter: look at the "
              "lines above before depositing")
        return 1
    print("\nboth environments reproduce every version the archive recorded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
