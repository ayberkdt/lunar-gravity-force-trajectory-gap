"""Hold the R23 queue until the machine has been quiet, then run it.

R19 is a timing-insensitive campaign -- it scores position error, not wall
clock -- but it is still a 12-way propagation, and starting it on top of an
unrelated six-process job would both slow that job and stretch this one past
the window the operator allowed. ``wait_for_idle.py`` already guards the timed
panels, but it releases on the first quiet poll and its watch list predates the
SPICE orbit-family campaign. This guard is stricter in the two ways that matter
here: it requires the quiet to be *sustained*, and it treats any other Python
campaign process as competing.

Idle means, for every poll in a continuous window: no other Python process is
running a campaign script, and system-wide CPU is below a threshold.

The queue it releases holds the two controls the review made gate items, in the
order that puts the required work first:

    r19A    R19 at beta = 0.5 on design A -- realized-work matching at the
            budget scale where the constructive claim is strongest and has
            never been tested this way. The review asks for design A only.
    oracle  R23-B, the fixed oracle applied to the interior member.
    r19B    the design-B replication of r19A, run only if the window allows.

Every stage is resume-safe: a completed orbit writes a sidecar whose config
hash is checked on re-entry, so a stage cut short by the stop clock costs
nothing but the orbits it had not reached.

Usage:
    python rev23_run_when_idle.py --stop-clock 23:00
    python rev23_run_when_idle.py --stop-clock 23:00 --stages r19A,oracle
"""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time
from pathlib import Path

import rev10_sobol_confirmatory as base

HERE = Path(__file__).resolve().parent
LOG = HERE / "r23_queue.log"

# Any other Python process whose command line contains one of these is a
# competing campaign. Prefix-matched, like the wait_for_idle guard, so a new
# campaign script has to be added here or the guard will pass while it runs.
CAMPAIGN = (
    "run_spice_orbit_families", "run_orbit_experiments", "run_experiments",
    "rev2_", "rev3_", "rev4_", "rev5_", "rev6_", "rev7_", "rev8_", "rev9_",
    "rev10_", "rev11_", "rev12_", "rev13_", "rev14_", "rev15_", "rev16_",
    "rev17_", "rev18_", "rev19_", "rev20_", "rev21_", "rev22_",
    "robustness_", "supplemental_",
    # R23 campaign scripts are named individually rather than by the "rev23_"
    # prefix, so that a second observer waits for the campaign without seeing
    # the first observer as competing work and deadlocking on it.
    "rev23_ultratight_span", "rev23_oracle_vs_interior",
    "rev23_cost_curve_unified",
)


def say(msg: str) -> None:
    stamp = dt.datetime.now().strftime("%H:%M:%S")
    line = f"[r23 {stamp}] {msg}"
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def competing() -> list[str]:
    out = []
    for p in base.other_python_processes():
        cmd = " ".join(p.get("cmdline") or [])
        hit = next((w for w in CAMPAIGN if w in cmd), None)
        if hit:
            out.append(f"{p['pid']}:{hit}")
    return out


def cpu_percent(sample_s: float) -> float:
    try:
        import psutil
    except ImportError:
        return 0.0
    return float(psutil.cpu_percent(interval=sample_s))


def wait_for_sustained_idle(idle_s: float, poll_s: float, cpu_max: float,
                            stop_at: float) -> bool:
    """True once the machine has been quiet for a continuous idle_s."""
    quiet_since: float | None = None
    while True:
        if time.time() >= stop_at:
            say("stop clock reached before the machine ever went quiet; "
                "nothing was run")
            return False
        busy = competing()
        cpu = cpu_percent(min(poll_s, 5.0))
        if busy or cpu > cpu_max:
            if quiet_since is not None:
                say(f"quiet broken after "
                    f"{(time.time() - quiet_since) / 60:.1f} min")
            quiet_since = None
            why = f"{len(busy)} campaign proc(s) {busy}" if busy \
                else f"cpu {cpu:.0f}% > {cpu_max:.0f}%"
            say(f"busy: {why}")
        else:
            if quiet_since is None:
                quiet_since = time.time()
                say(f"quiet started (cpu {cpu:.0f}%)")
            held = time.time() - quiet_since
            if held >= idle_s:
                say(f"quiet held {held / 60:.1f} min >= "
                    f"{idle_s / 60:.0f} min; releasing the queue")
                return True
            say(f"quiet {held / 60:.1f}/{idle_s / 60:.0f} min (cpu {cpu:.0f}%)")
        time.sleep(max(poll_s - 5.0, 1.0))


COST_OPTS = {"target_s": "0.25", "sessions": "1", "out": None}


def stage_command(name: str, beta: float, workers: int, budget_min: float):
    """(argv, human label) for a stage name, or None if the name is unknown."""
    if name in ("r19A", "r19B"):
        design = name[-1]
        return ([sys.executable, "-u", "rev19_equal_total_work.py", "run",
                 "--design", design, "--beta", f"{beta}",
                 "--workers", str(workers),
                 "--deadline-min", f"{budget_min:.1f}"],
                f"R19 realized-work match, design {design}, beta={beta}")
    if name == "oracle":
        return ([sys.executable, "-u", "rev23_oracle_vs_interior.py", "run",
                 "--workers", str(workers),
                 "--deadline-min", f"{budget_min:.1f}"],
                "R23-B fixed oracle against the interior member, beta=1")
    if name == "costcurve":
        # single process by construction: a timing panel is only valid on a
        # machine with nothing else on it, which is why it is staged alone
        cmd = [sys.executable, "-u", "rev23_cost_curve_unified.py",
               "--target-s", COST_OPTS["target_s"],
               "--sessions", COST_OPTS["sessions"]]
        if COST_OPTS["out"]:
            cmd += ["--out", COST_OPTS["out"]]
        return (cmd, f"R23-D cost curve, {COST_OPTS['sessions']} session(s) "
                     f"at {COST_OPTS['target_s']}s per cell")
    if name == "ultraspan":
        return ([sys.executable, "-u", "rev23_ultratight_span.py", "run",
                 "--workers", str(workers),
                 "--deadline-min", f"{budget_min:.1f}"],
                "R23-C third tolerance level for the beta=1 comparison")
    return None


def run_stage(name: str, beta: float, workers: int,
              stop_at: float, min_useful_min: float) -> int:
    left_min = (stop_at - time.time()) / 60.0
    spec = stage_command(name, beta, workers, max(left_min - 2.0, 1.0))
    if spec is None:
        say(f"stage {name}: unknown, skipping")
        return 0
    cmd, label = spec
    if left_min < min_useful_min:
        say(f"stage {name}: only {left_min:.0f} min left, below the "
            f"{min_useful_min:.0f} min floor; skipping")
        return 0
    say(f"stage {name} -- {label}: starting, {left_min:.0f} min of window "
        f"left, workers={workers}")
    t0 = time.time()
    with LOG.open("a", encoding="utf-8") as fh:
        rc = subprocess.call(cmd, cwd=str(HERE), stdout=fh,
                             stderr=subprocess.STDOUT)
    say(f"stage {name}: exit {rc} after {(time.time() - t0) / 60:.1f} min")
    return rc


def parse_clock(text: str) -> float:
    hh, mm = (int(v) for v in text.split(":"))
    now = dt.datetime.now()
    stop = now.replace(hour=hh, minute=mm, second=0, microsecond=0)
    if stop <= now:
        stop += dt.timedelta(days=1)
    return stop.timestamp()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-clock", default="23:00",
                    help="local HH:MM after which no new work is started")
    ap.add_argument("--idle-min", type=float, default=5.0)
    ap.add_argument("--poll-s", type=float, default=30.0)
    ap.add_argument("--cpu-max", type=float, default=40.0)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--beta", type=float, default=0.5)
    ap.add_argument("--stages", default="r19A,oracle,r19B",
                    help="ordered stage names: r19A, oracle, r19B")
    ap.add_argument("--min-useful-min", type=float, default=20.0)
    ap.add_argument("--cost-target-s", default="0.25")
    ap.add_argument("--cost-sessions", default="1")
    ap.add_argument("--cost-out", default=None)
    a = ap.parse_args()

    COST_OPTS["target_s"] = a.cost_target_s
    COST_OPTS["sessions"] = a.cost_sessions
    COST_OPTS["out"] = a.cost_out

    stop_at = parse_clock(a.stop_clock)
    say(f"observer up: waiting for {a.idle_min:.0f} min of sustained quiet, "
        f"stop clock {a.stop_clock} "
        f"({(stop_at - time.time()) / 60:.0f} min away)")

    if not wait_for_sustained_idle(a.idle_min * 60.0, a.poll_s,
                                   a.cpu_max, stop_at):
        return 3

    worst = 0
    for name in [s.strip() for s in a.stages.split(",") if s.strip()]:
        rc = run_stage(name, a.beta, a.workers, stop_at, a.min_useful_min)
        worst = max(worst, rc)
    say(f"queue complete (worst exit {worst})")
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
