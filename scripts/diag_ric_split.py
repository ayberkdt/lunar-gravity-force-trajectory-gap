"""Does the beta = 1 endpoint ordering survive dropping the in-track component?

At the declared budget both policies' seven-day errors are almost entirely
in-track. A reader whose application estimates the orbit absorbs much of that
component, so the ranking has to be shown not to rest on it alone. This
recomputes the (O10/O14) beta = 1 comparison on the radial, cross-track and
radial+cross-track norms from the archived state arrays, under the same
resolution rule the paper uses everywhere (Eq. 9): the envelope is rebuilt from
the same level-to-level self-differences, projected onto the same components.

Diagnostic only; it decides nothing the manuscript reports as a verdict.

Usage:  python diag_ric_split.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import rev10_sobol_confirmatory as base
import rev14_budget_trajectory as r14

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "diag_ric_split.json"
LEVELS = ("tight", "tighter")


def load(p):
    d = np.load(p)
    return d["t_s"], d["state_si"]


def ric_rms(ref, other):
    """Component RMS of (other - ref) in the reference RIC frame."""
    t, x = ref
    _, y = other
    n = min(len(t), len(y))
    r = x[:n, :3]
    v = x[:n, 3:6]
    d = y[:n, :3] - r
    er = r / np.linalg.norm(r, axis=1, keepdims=True)
    h = np.cross(r, v)
    ec = h / np.linalg.norm(h, axis=1, keepdims=True)
    ei = np.cross(ec, er)
    comp = {"radial": (d * er).sum(1), "in_track": (d * ei).sum(1),
            "cross": (d * ec).sum(1)}
    out = {k: float(np.sqrt(np.mean(v_ ** 2))) for k, v_ in comp.items()}
    out["radial_cross"] = float(np.sqrt(out["radial"] ** 2 + out["cross"] ** 2))
    out["total"] = float(np.sqrt(np.mean((d ** 2).sum(1))))
    return out


def main() -> int:
    tally = {}
    rows = []
    for design in ("A", "B"):
        for idx in range(64):
            # At beta = 1 the constant comparator IS the critical degree, so
            # its arcs are the reused R11 ones; only the radial policy is
            # written under the R14 budget tree.
            try:
                arcs = {"truth": {}, "fixed_budget": {}, "atallah_budget": {}}
                for lv in LEVELS:
                    _, raw = r14.reuse_paths(design, idx, "truth", lv)
                    arcs["truth"][lv] = load(raw)
                    _, raw = r14.reuse_paths(design, idx, "fixed_critical", lv)
                    arcs["fixed_budget"][lv] = load(raw)
                    _, raw = r14.paths(design, 1.00, idx, "atallah_budget", lv)
                    arcs["atallah_budget"][lv] = load(raw)
            except Exception:
                continue
            ref = arcs["truth"]["tighter"]
            err = {w: ric_rms(ref, arcs[w]["tighter"])
                   for w in ("fixed_budget", "atallah_budget")}
            env = {w: ric_rms(arcs[w]["tight"], arcs[w]["tighter"])
                   for w in ("fixed_budget", "atallah_budget")}
            env_ref = ric_rms(arcs["truth"]["tight"], arcs["truth"]["tighter"])
            rec = {"design": design, "sobol_index": idx}
            for comp in ("radial", "in_track", "cross", "radial_cross",
                         "total"):
                gap = abs(err["fixed_budget"][comp] - err["atallah_budget"][comp])
                den = (env["fixed_budget"][comp] + env["atallah_budget"][comp]
                       + 2 * env_ref[comp])
                resolved = den > 0 and gap / den > 1.0
                win = ("fixed" if err["fixed_budget"][comp]
                       < err["atallah_budget"][comp] else "radial")
                den_r = err["atallah_budget"][comp]
                rec[comp] = {"resolved": bool(resolved), "winner": win,
                             "rho": (err["fixed_budget"][comp] / den_r
                                     if den_r > 0 else None)}
                if resolved:
                    key = (design, comp, win)
                    tally[key] = tally.get(key, 0) + 1
            rows.append(rec)

    print(f"{len(rows)} orbits scored\n")
    print(f"{'design':<8}{'component':<14}{'fixed':>7}{'radial':>8}")
    for design in ("A", "B"):
        for comp in ("total", "in_track", "radial", "cross", "radial_cross"):
            f = tally.get((design, comp, "fixed"), 0)
            r = tally.get((design, comp, "radial"), 0)
            print(f"{design:<8}{comp:<14}{f:>7}{r:>8}")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\n[written] {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
