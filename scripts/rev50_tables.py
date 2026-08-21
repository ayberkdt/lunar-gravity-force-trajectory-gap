"""R50: the span-ladder table, generated from the records rather than typed.

Two panels in one table. The upper panel is the resolved tally per apolune level
and budget; the lower is the median error ratio in the same cells. Levels are
rows because the claim is about levels, and budgets ascend left to right because
the crossing moves along that axis.

Blocks A and B are pooled level by level, which the registration allows because
they are the same design at the same levels, and the pooled orbit count is
printed so the pooling is visible rather than implied. Nothing here is pooled
across levels: a level-pooled tally would average away the only quantity this
population varies.

Two marks carry provenance into the table itself. A budget added by the
amendment is marked, so it is never read as part of the registered grid, and a
level outside the frozen factor box is marked, so it is never read as a result
about the sampled box.

Usage:
    python rev50_tables.py
    python rev50_tables.py --readout secondary
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "metrics"

import population_registry as registry                        # noqa: E402

REG = "r50"
REGISTERED = [0.50, 0.62, 0.75, 1.00]
AMENDED = [1.25, 1.50]
AMEND_MARK = r"$^{\ddagger}$"
OUTBOX_MARK = r"$^{\ast}$"


def design_index(name: str, spec: dict) -> dict:
    d = json.loads((METRICS / spec["file"]).read_text(encoding="utf-8"))
    if d["design_sha256"] != spec["design_sha256"]:
        raise SystemExit(f"{spec['file']} does not match the registered hash")
    return {o["sobol_index"]: (o["apolune_level_km"],
                              o["apolune_level_inside_factor_box"])
            for o in d["orbits"]}


def read(key: str, beta: float, readout: str) -> list[dict] | None:
    if readout == "primary":
        p = METRICS / f"r14_trajectory_{key}_beta_{beta:.2f}.json"
        if not p.exists():
            return None
        d = json.loads(p.read_text(encoding="utf-8"))
        out = []
        for r in d["rows"]:
            c = r["comparison"]
            winner = c.get("resolved_winner") or c.get("raw_winner")
            out.append({"i": r["sobol_index"], "resolved": bool(c["resolved"]),
                        "varying": winner == "atallah",
                        "rho": c.get("rho_budget")})
        return out
    p = METRICS / f"r19_equal_total_work_{key}_beta_{beta:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    if "summary" not in d:
        return None
    return [{"i": r["sobol_index"], "resolved": bool(r["resolved"]),
             "varying": r["winner"] == "interior",
             "rho": r.get("rho_workmatched")} for r in d["rows"]]


def censored_by_budget(betas: list[float]) -> dict[float, int]:
    """Orbits dropped before comparison, pooled over blocks, per budget.

    Both readouts are built from the same propagated set, so the censoring the
    trajectory record books is the censoring both tables inherit; it is read
    from there once rather than recounted per readout.
    """
    out: dict[float, int] = {}
    for name, spec in registry.populations(REG).items():
        key = spec["design_key"]
        for b in betas:
            p = METRICS / f"r14_trajectory_{key}_beta_{b:.2f}.json"
            if not p.exists():
                continue
            d = json.loads(p.read_text(encoding="utf-8"))
            out[b] = out.get(b, 0) + len(d.get("censored", []))
    return out


def collect(readout: str) -> tuple[dict, list[float], dict]:
    pops = registry.populations(REG)
    index = {n: design_index(n, s) for n, s in pops.items()}
    cells: dict[tuple[float, float], dict] = {}
    present: list[float] = []
    inbox: dict[float, bool] = {}
    for beta in REGISTERED + AMENDED:
        any_block = False
        for name, spec in pops.items():
            rows = read(spec["design_key"], beta, readout)
            if rows is None:
                continue
            any_block = True
            for r in rows:
                level, in_box = index[name][r["i"]]
                inbox[level] = in_box
                c = cells.setdefault((level, beta),
                                     {"orbits": 0, "resolved": 0, "win": 0,
                                      "loss": 0, "rho": [], "blocks": set()})
                c["orbits"] += 1
                c["blocks"].add(name)
                if r["rho"] is not None:
                    c["rho"].append(float(r["rho"]))
                if not r["resolved"]:
                    continue
                c["resolved"] += 1
                if r["varying"]:
                    c["win"] += 1
                else:
                    c["loss"] += 1
        if any_block:
            present.append(beta)
    return cells, present, inbox


def fmt_rho(values: list[float]) -> str:
    """siunitx is not loaded in either document, so the exponent is written."""
    if not values:
        return "--"
    m = median(values)
    if m >= 1000.0:
        exponent = 0
        while m >= 10.0:
            m /= 10.0
            exponent += 1
        return f"${m:.1f}\\times10^{{{exponent}}}$"
    if m >= 10.0:
        return f"${m:.0f}$"
    if m >= 1.0:
        return f"${m:.2f}$"
    return f"${m:.3f}$"


WORDS = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
         6: "six"}


def verdict_of(win: int, loss: int) -> str:
    if win + loss == 0:
        return "undecided"
    return "varying" if win > loss else "constant" if loss > win else "split"


def block_agreement(readout: str) -> tuple[int, int]:
    """Cells both blocks carry, and how many of them return the same verdict.

    Pooling answers a question about weight. The blocks are independent draws
    of identities at the same levels, so read apart they answer a different one,
    and the count of cells where they agree is the only form of that answer
    that does not go stale when a budget lands.
    """
    pops = registry.populations(REG)
    if len(pops) < 2:
        return (0, 0)
    per_block: dict[str, dict] = {}
    for name, spec in pops.items():
        index = design_index(name, spec)
        for beta in REGISTERED + AMENDED:
            rows = read(spec["design_key"], beta, readout)
            if rows is None:
                continue
            for r in rows:
                level, _ = index[r["i"]]
                c = per_block.setdefault(name, {}).setdefault(
                    (level, beta), [0, 0])
                if r["resolved"]:
                    c[0 if r["varying"] else 1] += 1
    names = sorted(per_block)
    shared = set(per_block[names[0]]) & set(per_block[names[1]])
    agree = sum(1 for k in shared
                if verdict_of(*per_block[names[0]][k])
                == verdict_of(*per_block[names[1]][k]))
    return agree, len(shared)


def _agreement_phrase(agreement: tuple[int, int]) -> str:
    return f"{agreement[0]} of the {agreement[1]}" if agreement[1] else "no"


def write_macros(levels: list[float], betas: list[float],
                 cell_counts: list[int], blocks: list[str],
                 agreement: tuple[int, int], readout: str) -> None:
    """Counts the prose needs, so that no count is typed into a sentence.

    Every quantity here changes when a block or a budget lands. A sentence that
    carried one of them as a literal would be true when it was written and
    stale an hour later, which is how this campaign's own supplement acquired a
    wrong manifest count once already.
    """
    per = sorted(set(cell_counts))
    registered = [b for b in betas if b in REGISTERED]
    amended = [b for b in betas if b in AMENDED]
    # These emit noun phrases rather than counts. A macro that emits "two" and
    # a sentence that supplies "block" agree only until the count changes, and
    # the first compile after block B landed printed "two block" and "none
    # amended ones" for exactly that reason.
    amended_phrase = {0: "neither amended budget",
                      1: "one of the two amended budgets",
                      2: "both amended budgets"}
    macros = {
        "RSblockPhrase": ("one block" if len(blocks) == 1
                          else f"{WORDS.get(len(blocks), len(blocks))} blocks"),
        "RSamendedPhrase": amended_phrase.get(
            len(amended), f"{len(amended)} amended budgets"),
        "RSlevels": WORDS.get(len(levels), str(len(levels))),
        "RSorbitsPerLevel": (str(per[0]) if len(per) == 1
                             else f"{min(per)}--{max(per)}"),
        # cells, not distinct orbits: an orbit appears once per budget, so the
        # sum over cells is a count of comparisons and is named as one.
        "RScomparisonsTotal": str(sum(cell_counts)),
        "RSregisteredBudgets": WORDS.get(len(registered), str(len(registered))),
        "RSbudgetList": ", ".join(f"{b:.2f}" for b in betas),
        # Both readouts, always, whichever table this invocation is writing.
        # The single unsuffixed macro this replaced was whatever the last run
        # happened to be, so a sentence citing it changed meaning with the
        # order the two tables were regenerated in and said which readout
        # nowhere.
        "RSblockAgreementPrimary": _agreement_phrase(
            agreement if readout == "primary" else block_agreement("primary")),
        "RSblockAgreementSecondary": _agreement_phrase(
            agreement if readout == "secondary"
            else block_agreement("secondary")),
        "RSspanRange": (f"{min(levels):.0f}--{max(levels):.0f}"),
    }
    body = ["% generated by rev50_tables.py; do not edit"]
    body += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()]
    (METRICS / "r50_span_ladder_macros.tex").write_text(
        "\n".join(body) + "\n", encoding="utf-8")
    print(f"[written] r50_span_ladder_macros.tex ({len(macros)} macros)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--readout", choices=("primary", "secondary"),
                    default="primary")
    a = ap.parse_args()

    cells, betas, inbox = collect(a.readout)
    if not cells:
        print("no records for this population yet; no table written")
        return 0
    levels = sorted({lvl for lvl, _ in cells})
    blocks = sorted({b for c in cells.values() for b in c["blocks"]})
    per_level = {lvl: max(cells[(lvl, b)]["orbits"] for b in betas
                          if (lvl, b) in cells) for lvl in levels}

    what = ("the budget-calibrated radial endpoint against its equal-budget "
            "constant degree" if a.readout == "primary" else
            "the interior member of the span family against its work-matched "
            "constant degree")
    ratio = ("$E_{\\mathrm{fix}}/E_{\\mathrm{rad}}$"
             if a.readout == "primary"
             else "$E_{\\mathrm{fix}}/E_{\\mathrm{policy}}$")

    head = " & ".join(f"$\\beta={b:.2f}$" + (AMEND_MARK if b in AMENDED else "")
                      for b in betas)
    lines = [
        r"% generated by rev50_tables.py; do not edit",
        r"\begin{tabular}{@{}l" + "r" * len(betas) + r"@{}}",
        r"\toprule",
        r"Apolune & " + head + r" \\",
        r"\midrule",
    ]
    # Blocks land one budget at a time, so a budget can be pooled over both
    # blocks while its neighbour is not. The count per cell is printed rather
    # than implied, because a table whose columns rest on different numbers of
    # orbits and does not say so is a table that misleads by omission.
    per_cell = []
    for b in betas:
        counts = sorted({cells[(lvl, b)]["orbits"] for lvl in levels
                         if (lvl, b) in cells})
        # A range, not a fraction: "14/16" reads as fourteen of sixteen, which
        # is what the tally means two rows below and not what this is. The
        # dash stays outside math mode, because inside it a range dash sets as
        # two minus signs.
        per_cell.append(f"${counts[0]}$--${counts[-1]}$" if len(counts) > 1
                        else f"${counts[0]}$" if counts else "--")
    # Where the orbit count falls below the full cell, the missing orbits are
    # censored rather than unresolved, and the two are not the same fact: an
    # unresolved orbit was compared and the comparison did not separate, a
    # censored one never entered the comparison. The censoring here is also not
    # random. Its recorded reason is the reference-degree ceiling, so it lands
    # on the level where that ceiling binds hardest, and a reader who takes the
    # count drop for a resolution effect draws the wrong inference from it.
    censored = censored_by_budget(betas)
    lines += [
        r"Orbits per cell & " + " & ".join(per_cell) + r" \\",
        r"Censored at the ceiling & "
        + " & ".join(f"${censored.get(b, 0)}$" for b in betas) + r" \\",
        r"\midrule",
        r"\multicolumn{%d}{@{}l}{\emph{Resolved tally, varying degree--constant"
        r" degree}} \\" % (len(betas) + 1),
    ]
    for lvl in levels:
        mark = "" if inbox.get(lvl, True) else OUTBOX_MARK
        row = [f"${lvl:.0f}$~km{mark}"]
        for b in betas:
            c = cells.get((lvl, b))
            row.append("--" if c is None
                       else f"${c['win']}$--${c['loss']}$")
        lines.append(" & ".join(row) + r" \\")
    # Where the verdict turns over, per budget, derived rather than asserted.
    # The threshold moves with the budget, so a sentence naming its place is
    # stale as soon as another budget lands; a generated row is not.
    thresholds = []
    for b in betas:
        varying = []
        for i, lvl in enumerate(levels):
            c = cells.get((lvl, b))
            if c and c["resolved"] and c["win"] > c["loss"]:
                varying.append(i)
        if not varying:
            thresholds.append("none")
        elif varying != list(range(varying[0], len(levels))):
            thresholds.append("mixed")
        elif varying[0] == 0:
            thresholds.append(f"$<{levels[0]:.0f}$")
        else:
            thresholds.append(f"${levels[varying[0] - 1]:.0f}$--"
                              f"${levels[varying[0]]:.0f}$")
    lines += [r"\addlinespace",
              r"Turnover (km) & " + " & ".join(thresholds) + r" \\"]
    lines += [r"\midrule",
              r"\multicolumn{%d}{@{}l}{\emph{Median %s in the same cells}} \\"
              % (len(betas) + 1, ratio)]
    for lvl in levels:
        mark = "" if inbox.get(lvl, True) else OUTBOX_MARK
        row = [f"${lvl:.0f}$~km{mark}"]
        for b in betas:
            c = cells.get((lvl, b))
            row.append("--" if c is None else fmt_rho(c["rho"]))
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]

    out = METRICS / f"r50_span_ladder_{a.readout}_table.tex"
    counts = ", ".join(f"{lvl:.0f}~km on {per_level[lvl]}" for lvl in levels)
    write_macros(levels, betas, [c["orbits"] for c in cells.values()], blocks,
                 block_agreement(a.readout), a.readout)
    note = (f"% {a.readout}: {what}. Blocks pooled: {', '.join(blocks)}. "
            f"Orbits per level: {counts}. Budgets present: "
            f"{', '.join(f'{b:.2f}' for b in betas)}.\n")
    out.write_text(note + "\n".join(lines) + "\n", encoding="utf-8")
    print(f"[written] {out.name}")
    print(f"  {what}")
    print(f"  blocks pooled: {', '.join(blocks)}")
    print(f"  orbits per level: {counts}")
    print(f"  budgets: {', '.join(f'{b:.2f}' for b in betas)}")
    missing = [f"{b:.2f}" for b in REGISTERED + AMENDED if b not in betas]
    if missing:
        print(f"  not present, reported as not run: {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
