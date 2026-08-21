"""Unattended chain for the night of 15 August 2026: finish the R53 column.

Five of the seven cells of the post-hoc beta = 0.62 column are on disk. The two
that are not, the equatorial and frozen-like sub-boxes, are what stands between
the manuscript's "computed on four of the six" and "computed on all six", and
they are the last declared-and-not-run cells of that registration.

What this chain does, in order, and nothing else:

  1. rev53_campaign.py at a worker count that has been observed to hold, with
     the registered prior scaled to that worker count so the campaign's own
     "never start a cell it cannot finish" guard stays honest.
  2. if a cell is still missing and the clock allows, the same campaign again at
     a lower worker count, down a fixed ladder. Each pass is a fresh supervisor
     with all of its own guards; nothing here reimplements them.
  3. rev53_verdict.py, once, after the cells have stopped changing.
  4. a morning report naming exactly what is now stale.

What it deliberately does not do: it writes no .tex, seals no manifest, redraws
no figure and touches no claim. Those are judgement, they are quick, and they
are wrong to do while the numbers they quote may still move.

Why the worker ladder. Cells 1-5 ran at eleven workers on 14 August. The
equatorial cell then broke a pool worker abruptly at eleven, at eight and at six,
including inside the span stage alone with the trajectory record already
complete; at four workers the pool held and all four workers were doing arithmetic
ninety seconds in. The failures are therefore read as contention on this machine
rather than as a property of the cell, and the response is fewer workers rather
than more attempts at the same number.

Usage:
    python launch_detached.py ../output/r54_night.log rev54_night_chain.py
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
METRICS = ROOT / "metrics"
OUTPUT = ROOT / "output"

BETA = 0.62
KEYS = ("SE", "SF")

# The campaign must be finished, scored and reported long before the machine is
# wanted back. Everything after the campaign is seconds of work.
CAMPAIGN_STOP = "2026-08-15 09:30"

# Tried in order, and only while the clock still allows a whole cell. Four is
# where the pool was observed to hold; the two below it are there so that a
# machine that is busier at 04:00 than it is now still has somewhere to go.
WORKER_LADDER = (4, 3, 2)

# rev53's registered priors were measured at eleven workers. At four the same
# ladder takes roughly the ratio of the worker counts, and the guard that
# refuses a cell it cannot finish compares against the prior, so the prior has
# to be scaled or the guard will wave through a cell it cannot complete.
PRIOR_SCALE = {4: 4.0, 3: 5.0, 2: 7.0}

# A campaign pass is given its own stop time plus this much slack before it is
# killed outright. The supervisor stops itself at its stop time; this only
# catches a supervisor that has stopped stopping.
PASS_SLACK_MIN = 45.0

# Below this the chain does not start another pass at all. A pass that cannot
# fit one cell only burns clock and writes an incomplete record.
MIN_PASS_MIN = 70.0

DISK_FLOOR_MB = 700.0

DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200

LOG = OUTPUT / "r54_night_chain.log"
REPORT = OUTPUT / "r54_morning_report.md"


def log(msg: str) -> None:
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
    print(line, flush=True)
    try:
        with LOG.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except OSError:
        pass


def free_mb() -> float:
    return shutil.disk_usage(METRICS).free / (1024.0 * 1024.0)


def remaining_min(stop: datetime) -> float:
    return (stop - datetime.now()).total_seconds() / 60.0


def ladder_done(key: str) -> bool:
    """The same test rev53_campaign uses, repeated rather than imported.

    Importing the supervisor would run its argument parser, and a chain that
    cannot answer "is this cell finished" without importing the thing it
    supervises is a chain that stops working the moment the supervisor grows a
    required flag.
    """
    p = METRICS / f"r19_equal_total_work_{key}_beta_{BETA:.2f}.json"
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        return "summary" in json.loads(p.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return False


def missing() -> list[str]:
    return [k for k in KEYS if not ladder_done(k)]


def spawn(args: list[str], log_path: Path, wall_limit_min: float) -> int:
    """Run one child with no console of its own, and kill it if it overruns.

    The console-less start is not optional here: three earlier launches of other
    campaigns died with STATUS_CONTROL_C_EXIT because a child that owns a console
    receives that console's control events, and an agent shell finishing an
    unrelated command was enough to deliver one.
    """
    log(f"start: {' '.join(args)}  (wall limit {wall_limit_min:.0f} min)")
    t0 = time.time()
    # unbuffered, so that both this chain's wall limit and the supervisor's own
    # stall guard read progress off a log that is written as it happens rather
    # than eight kilobytes at a time
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    with log_path.open("a", encoding="utf-8") as fh:
        proc = subprocess.Popen(
            [sys.executable] + args, cwd=str(HERE),
            stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
            env=env,
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)
        while True:
            try:
                rc = proc.wait(timeout=60)
                break
            except subprocess.TimeoutExpired:
                if (time.time() - t0) / 60.0 > wall_limit_min:
                    log("wall limit passed; killing the pass and its workers")
                    # /T because the pool workers are its children and killing
                    # the supervisor alone would leave them on the CPU
                    subprocess.call(["taskkill", "/PID", str(proc.pid),
                                     "/T", "/F"],
                                    stdout=subprocess.DEVNULL,
                                    stderr=subprocess.DEVNULL)
                    rc = proc.wait()
                    break
    log(f"end:   rc={rc} in {(time.time()-t0)/60.0:.1f} min")
    return rc


def run_campaign(stop: datetime) -> None:
    for workers in WORKER_LADDER:
        left = missing()
        if not left:
            log("both cells are on disk; no further pass is needed")
            return
        avail = remaining_min(stop)
        if avail < MIN_PASS_MIN:
            log(f"{avail:.0f} min to the stop time, below the {MIN_PASS_MIN:.0f} "
                f"a pass needs; {', '.join(left)} left unrun")
            return
        disk = free_mb()
        if disk < DISK_FLOOR_MB:
            log(f"{disk:.0f} MB free, below the {DISK_FLOOR_MB:.0f} floor; "
                f"stopping rather than filling the volume mid-tree")
            return
        log(f"pass at {workers} workers: {', '.join(left)} missing, "
            f"{avail:.0f} min left, {disk:.0f} MB free")
        spawn(["rev53_campaign.py", "--stop-at", CAMPAIGN_STOP,
               "--workers", str(workers),
               "--prior-scale", f"{PRIOR_SCALE[workers]:g}"],
              OUTPUT / f"r54_campaign_w{workers}.log",
              avail + PASS_SLACK_MIN)
    if missing():
        log(f"worker ladder exhausted; {', '.join(missing())} still missing")


def rescore() -> str:
    """Re-score the verdict, keeping the version the manifest sealed.

    The old verdict is not deleted, because the manifest and the claims ledger
    both pin its digest and the morning has to be able to see what changed
    rather than only what it now says.
    """
    live = METRICS / "r53_verdict.json"
    if live.exists():
        keep = METRICS / "r53_verdict.sealed_20260814.json"
        if not keep.exists():
            shutil.copy2(live, keep)
            log(f"kept the sealed verdict as {keep.name}")
    rc = spawn(["rev53_verdict.py"], OUTPUT / "r54_verdict.log", 20.0)
    if rc != 0:
        return f"rev53_verdict.py returned {rc}; the verdict was not re-scored"
    try:
        d = json.loads(live.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        return f"the re-scored verdict could not be read back: {exc}"
    return (f"outcome {d.get('outcome')}: {len(d.get('cells_run', []))} cells "
            f"run, {len(d.get('cells_declared_and_not_run', []))} declared "
            f"and not run")


def cell_lines() -> list[str]:
    """The two panels' tallies for every cell, read from the verdict."""
    out: list[str] = []
    try:
        d = json.loads((METRICS / "r53_verdict.json").read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return ["the verdict could not be read; take the tallies from the "
                "r19 and r14 records directly"]
    for key, rep in d.get("report", {}).items():
        for panel in ("endpoint", "interior"):
            c = rep.get("panels", {}).get(panel, {}).get("cells", {}).get("0.62")
            if not c:
                out.append(f"  {key:<4} {panel:<9} 0.62 not computed")
                continue
            out.append(f"  {key:<4} {panel:<9} 0.62  "
                       f"{c['wins']}--{c['losses']} of {c['resolved']} resolved, "
                       f"median rho {c['median_rho']:.3f}")
    return out


def write_report(verdict_note: str, started: datetime) -> None:
    left = missing()
    done = [k for k in KEYS if k not in left]
    lines = [
        "# R54 night chain, 15 August 2026",
        "",
        f"Started {started.isoformat(timespec='seconds')}, "
        f"finished {datetime.now().isoformat(timespec='seconds')}.",
        f"Disk free at the end: {free_mb():.0f} MB.",
        "",
        "## Cells",
        "",
        f"- landed this run or already present: {', '.join(done) or 'none'}",
        f"- still missing: {', '.join(left) or 'none'}",
        "",
        "## Verdict",
        "",
        f"{verdict_note}",
        "",
        "## The 0.62 cells, both panels",
        "",
        *cell_lines(),
        "",
        "## What is now stale and needs a hand",
        "",
        "1. `metrics/r53_verdict.json` was re-scored, so its digest moved: the",
        "   claims `posthoc.column.cells_complete` and `posthoc.column.outcome`",
        "   are pinned to it and will read STALE until re-pinned, and",
        "   `r53_final_experiment_manifest.json` needs re-sealing.",
        "2. `posthoc.column.cells_complete` expects 5 and prints",
        "   'Five of the seven'; with both cells in it becomes seven of seven,",
        "   in `supp_budget_pareto.tex` and `supp_experiment_contract.tex`.",
        "3. The outcome-Z carriers say two cells were declared and not run:",
        "   `supp_budget_pareto.tex` (O52), `supp_experiment_contract.tex` (O52)",
        "   and `supp_reproducibility.tex` (the R53 entry).",
        "4. Main text: `08_discussion.tex` says the post-hoc point is 'computed",
        "   on four of the six'. With SE and SF in it is all six.",
        "5. Figure 4 must be redrawn with `make_figures_r36_regime.py` (never",
        "   `make_figures.py`, which overwrites three manuscript figures from",
        "   stale E-series data).",
        "6. The registration window was 14 August with a 10:45 cutoff and this",
        "   run extends it a second time, so `rev53_write_departure.py` has to",
        "   state the new window rather than leaving it to be noticed.",
        "",
        "Logs: `output/r54_campaign_w*.log`, `output/r54_verdict.log`,",
        "`python_codes/r53_campaign.log`.",
        "",
    ]
    tmp = REPORT.with_suffix(".md.tmp")
    tmp.write_text("\n".join(lines), encoding="utf-8")
    tmp.replace(REPORT)
    log(f"report written to {REPORT}")


def main() -> int:
    started = datetime.now()
    OUTPUT.mkdir(exist_ok=True)
    log("=== R54 night chain: finish the R53 post-hoc budget column ===")
    log(f"missing at start: {', '.join(missing()) or 'none'}; "
        f"{free_mb():.0f} MB free")
    verdict_note = "the verdict was not re-scored"
    try:
        stop = datetime.fromisoformat(CAMPAIGN_STOP)
        run_campaign(stop)
    except Exception as exc:                       # noqa: BLE001
        log(f"the campaign stage raised {type(exc).__name__}: {exc}")
    try:
        if len(missing()) < len(KEYS):
            verdict_note = rescore()
        else:
            verdict_note = ("no new cell landed, so the sealed verdict was "
                            "left exactly as it was")
        log(verdict_note)
    except Exception as exc:                       # noqa: BLE001
        verdict_note = f"the verdict stage raised {type(exc).__name__}: {exc}"
        log(verdict_note)
    try:
        write_report(verdict_note, started)
    except Exception as exc:                       # noqa: BLE001
        log(f"the report could not be written: {type(exc).__name__}: {exc}")
    log("=== R54 night chain done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
