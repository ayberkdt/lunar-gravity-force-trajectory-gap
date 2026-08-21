"""R30 tables: the geometry panel at the declared budget.

One table, one row per stratum, in the columns the coverage designs are already
reported in, with the coverage designs printed underneath for reference and
separated by a rule. The separation is not decoration: the strata are
sub-populations of the same box and the designs are draws from all of it, so a
reader who runs an eye down the column must be able to see where one kind of row
stops and the other begins.

Per-call and realized-work rows are both carried, because the whole point of the
campaign is that those two accountings can disagree, and a panel that showed
only one of them would hide exactly the effect the paper is about.

Nothing here is pooled or averaged. rev30_verdict.py holds the rule; this file
holds the presentation and reads the same records.

Usage:  python rev30_tables.py [--beta 1.00]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import rev19_tables as r19

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"
import population_registry as registry

LABEL = {"polar": "polar", "equatorial": "equatorial",
         "high_apolune": "high apolune", "frozen_like": "frozen-like",
         "low_perilune": "low perilune",
         "operational_elliptical": "operational elliptical$^{\\ast}$"}


def load(key: str, beta: float):
    """Every population outside designs A and B writes beta = 1 suffixed."""
    if key in ("A", "B"):
        return r19.load(key, beta)
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta:.2f}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def per_call(key: str, beta: float):
    return r19.per_call(key, beta)


BUDGETS = (1.00, 0.75, 0.50)


def rows_for(label: str, key: str, betas=BUDGETS) -> list[str]:
    """One block per population: a row pair for every budget it propagated."""
    out, first = [], True
    for beta in betas:
        d, pc = load(key, beta), per_call(key, beta)
        if not d or not pc:
            continue
        s = d["summary"]
        a = s["achieved_work_ratio"]
        head = label if first else ""
        first = False
        out.append(f"{head} & {beta:.2f} & nominal per call & {pc['interior']} "
                   f"& {pc['fixed']} & {pc['unresolved']} & "
                   f"{pc['median_of_ratios']:.2f} & target-matched \\\\")
        out.append(f" & & realized total & {s['resolved_interior_wins']} & "
                   f"{s['resolved_fixed_wins']} & {s['unresolved']} & "
                   f"{s['median_rho']:.2f} & "
                   f"${a['median']:.3f}$ $[{a['min']:.3f},{a['max']:.3f}]$ \\\\")
    return out


def table(beta: float, strata: dict, with_coverage: bool = False) -> str:
    lines = ["\\begin{tabular}{@{}l l l r r r r r@{}}", "\\toprule",
             "population & $\\beta$ & budget held equal & interior & fixed & "
             "unres. & median $\\rho$ & work match \\\\", "\\midrule"]
    wrote = False
    for name, spec in strata.items():
        rows = rows_for(LABEL.get(name, name), spec["design_key"])
        if not rows:
            continue
        if wrote:
            lines.append("\\addlinespace")
        lines += rows
        wrote = True
    if with_coverage:
        ref = []
        for key in ("A", "B", "C"):
            block = rows_for(f"design {key}", key)
            if block:
                if ref:
                    ref.append("\\addlinespace")
                ref += block
        if ref:
            lines.append("\\midrule")
            lines += ref
    lines += ["\\bottomrule", "\\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=1.00)
    ap.add_argument("--registry", default="r30")
    ap.add_argument("--with-coverage", action="store_true",
                    help="print the A/B/C coverage rows under the strata")
    a = ap.parse_args()
    prereg = registry.registration(a.registry)
    strata = registry.populations(a.registry)

    # The coverage designs used to be printed under the rule "for reference".
    # They are reported in full in the sections that own them (S11--S13), and
    # repeating eighteen of their rows here made the geometry table look like a
    # population table rather than a table about where a restriction moves the
    # crossing. Pass --with-coverage to restore them.
    coverage = [(f"design_{k}", {"design_key": k}) for k in ("A", "B", "C")] \
        if a.with_coverage else []
    out: dict = {}
    for name, spec in list(strata.items()) + coverage:
      key = spec["design_key"]
      for beta in BUDGETS:
        d, pc = load(key, beta), per_call(key, beta)
        if not d or not pc:
            continue
        s = d["summary"]
        out[f"{name}_beta_{beta:.2f}"] = {
            "population": name, "design_key": key, "beta": beta,
            "is_stratum": name in strata,
            "sub_box": (prereg.get("sub_boxes", {}).get(name)
                        or prereg.get("design", {}).get("factor_box")),
            "orbits": s["orbits"], "resolved": s["resolved"],
            "interior_wins": s["resolved_interior_wins"],
            "fixed_wins": s["resolved_fixed_wins"],
            "unresolved": s["unresolved"],
            "median_rho_realized": s["median_rho"],
            "achieved_work_ratio": s["achieved_work_ratio"],
            "per_call": pc,
        }

    tex = METRICS / f"{a.registry}_population_table.tex"
    desc = METRICS / f"{a.registry}_manuscript_descriptives.json"
    tex.write_text(table(a.beta, strata, a.with_coverage), encoding="utf-8")
    desc.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[written] {tex.name}, {desc.name}")
    for name, v in out.items():
        kind = "stratum " if v["is_stratum"] else "design  "
        print(f"  {kind}{name:<30} realized {v['interior_wins']}/"
              f"{v['fixed_wins']}/{v['unresolved']} of {v['orbits']}, rho "
              f"{v['median_rho_realized']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
