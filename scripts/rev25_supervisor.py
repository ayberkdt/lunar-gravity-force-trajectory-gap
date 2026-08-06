"""R25 supervisor: keep the machine working until a wall clock stop, and pick
each next budget from what the finished runs actually say.

The fixed queues that preceded this one each knew their budget in advance. That
stops being the right shape once the question is "where does the sign change
lie", because the informative next budget depends on the last answer. This
supervisor reads the records on disk between chains and chooses accordingly.

Work is taken in priority order, and the order is by what each job buys:

  1. bisect the crossing on design A. The realized-work verdict is 'fixed' at
     beta = 0.50 and 'interior' at 0.75, so the crossing sits between them and
     each chain halves the remaining bracket. This is the only job that makes
     the paper's interval narrower rather than merely better attested.
  2. replicate the bracket endpoints on design B. A crossing located on one
     scrambled design is a result about that design; located on two it is a
     result about the method.
  3. fill design B at beta = 1.50, which the manuscript currently names as a
     gap: "design B has no archived endpoint record at beta = 1.5, so that
     budget rests on one design".

A chain is started only if the window can hold it whole. Truncating the budget
stage would leave the two downstream stages running on a partial panel, and a
bracket localized that way is worth less than an honest gap. When no job is
both useful and affordable the supervisor stops rather than inventing work.

Budgets are rounded to two decimals because every driver tags its records with
f"beta_{beta:.2f}"; a budget that does not survive that rounding would write a
file whose name contradicts its contents.

Usage:
    python rev25_supervisor.py --stop-at "2026-07-31 09:30"
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "python_codes"
METRICS = ROOT / "metrics"
LOG = CODE / "r25_queue.log"

RESERVE_MIN = 20.0
# Measured on this machine at eleven workers: the span sweep took 70 min and
# the equal-work stage 21 min at beta = 0.75, and the budget stage runs about
# 256 propagations at the same rate, so a full chain is roughly 140 min. The
# thresholds sit above the measurement because starting a chain that cannot
# finish is the one failure this supervisor exists to avoid.
FULL_CHAIN_MIN = 150.0      # budget + span + equal-work
SHORT_CHAIN_MIN = 100.0     # span + equal-work, when the budget record exists
# design C's base is resumable at orbit granularity, so unlike a budget chain
# it loses nothing when the clock cuts it; a short window is still useful
DESIGNC_BASE_MIN = 40.0
TAG_RESOLUTION = 0.01       # the finest budget the record naming can express

PRIOR_QUEUE_DONE = re.compile(
    r"stage B complete|stage B skipped|B2 and B3 cannot run|"
    r"no useful window for B1|window closed after B")


def say(msg: str) -> None:
    line = f"[r25S {datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def minutes_left(stop_at: datetime) -> float:
    return (stop_at - datetime.now()).total_seconds() / 60.0


def tag(beta: float) -> str:
    return f"beta_{beta:.2f}"


def suffix(beta: float) -> str:
    return "" if abs(beta - 1.0) < 1e-9 else f"_{tag(beta)}"


# ------------------------------------------------------------------ reading
def verdict(design: str, beta: float) -> str | None:
    """What the realized-work comparison says at this budget, if it ran."""
    p = METRICS / f"r19_equal_total_work_{design}{suffix(beta)}.json"
    if not p.exists():
        return None
    try:
        s = json.loads(p.read_text(encoding="utf-8"))["summary"]
    except (OSError, KeyError, json.JSONDecodeError):
        return None
    win, lose = s.get("resolved_interior_wins", 0), s.get("resolved_fixed_wins", 0)
    if win == lose:
        return "split"
    return "interior" if win > lose else "fixed"


def known_budgets(design: str) -> dict[float, str]:
    out = {}
    for p in METRICS.glob(f"r19_equal_total_work_{design}*.json"):
        m = re.search(r"beta_(\d+\.\d+)", p.name)
        beta = float(m.group(1)) if m else 1.00
        v = verdict(design, beta)
        if v:
            out[beta] = v
    return out


def bracket(design: str) -> tuple[float, float] | None:
    """Widest pair of adjacent budgets whose verdicts differ."""
    known = known_budgets(design)
    lows = sorted(b for b, v in known.items() if v == "fixed")
    highs = sorted(b for b, v in known.items() if v == "interior")
    if not lows or not highs:
        return None
    lo = max(lows)
    above = [h for h in highs if h > lo]
    if not above:
        return None
    return (lo, min(above))


# ------------------------------------------------------------------ running
def run(cmd: list[str], label: str) -> int:
    say(f"{label}: $ {' '.join(cmd[1:])}")
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    for stream in (proc.stdout, proc.stderr):
        if stream:
            with LOG.open("a", encoding="utf-8") as fh:
                fh.write(stream)
            for ln in [l for l in stream.strip().splitlines() if l.strip()][-5:]:
                print("   " + ln, flush=True)
    say(f"{label}: exit {proc.returncode} after {(time.time()-t0)/60:.1f} min")
    return proc.returncode


def chain(design: str, beta: float, stop_at: datetime, workers: int) -> int:
    py = sys.executable
    worst = 0
    have_budget = (METRICS / f"r14_trajectory_{design}_{tag(beta)}.json").exists()

    if not have_budget:
        left = minutes_left(stop_at) - RESERVE_MIN
        b1_stop = datetime.now() + timedelta(minutes=min(left * 0.45, left))
        # R14 reads a naive deadline as UTC, so the local stop is made explicit
        worst = max(worst, run(
            [py, str(CODE / "rev14_budget_trajectory.py"), "run",
             "--design", design, "--beta", str(beta), "--workers", str(workers),
             "--deadline", b1_stop.astimezone().isoformat()],
            f"{design}@{beta:.2f} budget"))
        if not (METRICS / f"r14_trajectory_{design}_{tag(beta)}.json").exists():
            say(f"{design}@{beta:.2f}: budget record not produced; chain stops")
            return max(worst, 1)

    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 15:
        say(f"{design}@{beta:.2f}: window closed before the span stage")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev18_span_sweep.py"), "run",
         "--design", design, "--beta", str(beta), "--workers", str(workers),
         "--deadline-min", f"{min(left * 0.65, left):.0f}"],
        f"{design}@{beta:.2f} span"))
    run([py, str(CODE / "rev18_span_sweep.py"), "summarize",
         "--design", design, "--beta", str(beta), "--from-disk"],
        f"{design}@{beta:.2f} span summary")

    left = minutes_left(stop_at) - RESERVE_MIN
    if left <= 10:
        say(f"{design}@{beta:.2f}: window closed before the equal-work stage")
        return worst
    worst = max(worst, run(
        [py, str(CODE / "rev19_equal_total_work.py"), "run",
         "--design", design, "--beta", str(beta), "--workers", str(workers),
         "--deadline-min", f"{left:.0f}"], f"{design}@{beta:.2f} equal-work"))
    run([py, str(CODE / "rev19_equal_total_work.py"), "summarize",
         "--design", design, "--beta", str(beta)],
        f"{design}@{beta:.2f} equal-work summary")
    return worst


# ------------------------------------------------------------- job selection
def grid_budgets() -> set[float]:
    """Budgets the frozen Phase-A calibration actually carries.

    rev14_budget_trajectory does not calibrate a budget on demand; it reads one
    from r14_budget_pareto.json, whose grid is pre-registered. A budget outside
    that grid fails with a KeyError, and extending the grid is a registration
    decision rather than a scheduling one, so the supervisor treats the grid as
    the set of budgets it is allowed to propose.
    """
    try:
        d = json.loads((METRICS / "r14_budget_pareto.json").read_text(
            encoding="utf-8"))
        keys = d["designs"]["A"]["rows"][0]["budgets"].keys()
    except (OSError, KeyError, IndexError, json.JSONDecodeError):
        return set()
    out = set()
    for k in keys:
        m = re.fullmatch(r"beta_(\d+\.\d+)", k)
        if m:
            out.add(float(m.group(1)))
    return out


def cost_of(design: str, beta: float) -> float:
    return (SHORT_CHAIN_MIN
            if (METRICS / f"r14_trajectory_{design}_{tag(beta)}.json").exists()
            else FULL_CHAIN_MIN)


def next_job(affordable: float) -> tuple[str, float, str] | None:
    """The most valuable job that fits in the minutes still available.

    Bisection goes to whichever design currently has the wider bracket, which
    keeps the two designs in step instead of localizing one to the last decimal
    while the other has no midpoint at all. A crossing pinned on one scrambled
    design says something about that design; pinned on both, it says something
    about the method, and the second is what the paper needs.
    """
    grid = grid_budgets()
    candidates = []
    for design in ("A", "B"):
        br = bracket(design)
        if not br:
            continue
        lo, hi = br
        if hi - lo <= TAG_RESOLUTION * 1.5:
            continue
        # the midpoint is only runnable if the frozen calibration carries it;
        # the nearest grid budget strictly inside the bracket is the best the
        # supervisor may propose on its own authority
        inside = sorted(b for b in grid if lo < b < hi
                        and b not in known_budgets(design))
        if not inside:
            continue
        mid = (lo + hi) / 2.0
        pick = min(inside, key=lambda b: abs(b - mid))
        candidates.append((hi - lo, design, pick, lo, hi))
    # widest bracket first; design A breaks a tie because it is the design the
    # bracket was established on
    candidates.sort(key=lambda c: (-c[0], c[1]))
    for width, design, mid, lo, hi in candidates:
        if cost_of(design, mid) <= affordable:
            return (design, mid,
                    f"halve design {design}'s bracket ({lo:.2f}, {hi:.2f}], "
                    f"width {width:.2f}")

    # once no bracket can be narrowed on the frozen grid, fill the grid
    # budgets that were calibrated but never propagated. Design B at 1.50
    # first, because the manuscript names it: "design B has no archived
    # endpoint record at beta = 1.5, so that budget rests on one design".
    for design, beta in (("B", 1.50), ("A", 1.25), ("B", 1.25),
                         ("A", 2.00), ("B", 2.00)):
        if (beta in grid and beta not in known_budgets(design)
                and cost_of(design, beta) <= affordable):
            why = ("a gap the manuscript names" if (design, beta) == ("B", 1.50)
                   else "a calibrated grid budget never propagated")
            return (design, beta, f"design {design} at beta = {beta:.2f}, {why}")

    # last resort: start design C's base rather than leave the machine idle.
    # This produces nothing quotable on its own, but it is the long pole of a
    # third-design replication and it does not depend on which budget is run
    # later, so no part of it is wasted by a decision made afterwards.
    if designC_base_wanted() and DESIGNC_BASE_MIN <= affordable:
        return ("C", 0.0, "design C base: prepass and convergence tree")

    return None


def designC_base_wanted() -> bool:
    """True while design C is frozen but its convergence tree is unfinished."""
    if not (METRICS / "r26_sobolC_design_frozen.json").exists():
        return False
    out = METRICS / "r26_designC_convergence.json"
    if not out.exists():
        return True
    try:
        return not json.loads(out.read_text(encoding="utf-8")).get("complete")
    except (OSError, json.JSONDecodeError):
        return True


def designC_base(stop_at: datetime, workers: int) -> int:
    left = minutes_left(stop_at) - RESERVE_MIN
    b_stop = datetime.now() + timedelta(minutes=left)
    return run(
        [sys.executable, str(CODE / "rev26_designC_base.py"), "run",
         "--workers", str(workers),
         "--deadline", b_stop.astimezone().isoformat()], "design C base")


def wait_for_prior_queue(stop_at: datetime) -> None:
    say("waiting for the design-B queue to finish before taking over")
    while True:
        try:
            if PRIOR_QUEUE_DONE.search(
                    LOG.read_text(encoding="utf-8", errors="replace")):
                say("prior queue has ended; supervisor takes over")
                return
        except OSError:
            pass
        if minutes_left(stop_at) - RESERVE_MIN <= SHORT_CHAIN_MIN:
            say("window too small to take over; supervisor exits")
            return
        time.sleep(60)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True)
    ap.add_argument("--workers", type=int, default=11)
    args = ap.parse_args()

    stop_at = datetime.strptime(args.stop_at, "%Y-%m-%d %H:%M")
    wait_for_prior_queue(stop_at)

    worst = 0
    # A job that fails without consuming its window would otherwise be picked
    # again immediately and spin.
    stalls = 0
    while True:
        t_loop = time.time()
        affordable = minutes_left(stop_at) - RESERVE_MIN
        job = next_job(affordable)
        if job is None:
            say(f"no job both useful and affordable in {affordable:.0f} min; "
                f"stopping with the machine idle rather than inventing work")
            break
        design, beta, why = job
        if design == "C":
            say(f"next: {why} ({affordable:.0f} min available)")
            worst = max(worst, designC_base(stop_at, args.workers))
            continue
        say(f"next: design {design} at beta {beta:.2f} -- {why} "
            f"({affordable:.0f} min available, "
            f"{cost_of(design, beta):.0f} min expected)")
        worst = max(worst, chain(design, beta, stop_at, args.workers))
        say(f"design A budgets now: "
            f"{ {f'{b:.2f}': v for b, v in sorted(known_budgets('A').items())} }")
        if (time.time() - t_loop) / 60.0 < 2.0:
            stalls += 1
            if stalls >= 3:
                say("three jobs returned without doing work; stopping rather "
                    "than spinning")
                break
        else:
            stalls = 0

    br = bracket("A")
    say(f"supervisor done (worst exit {worst}); design A crossing bracket "
        f"{'(%.2f, %.2f]' % br if br else 'not established'}; "
        f"{minutes_left(stop_at):.0f} min to spare")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
