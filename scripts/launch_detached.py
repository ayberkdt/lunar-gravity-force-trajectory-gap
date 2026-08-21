"""Start a long campaign with no console at all, and exit.

Three launches of R50 died before this file existed, and the third one named the
cause. Start-Process from the agent shell was killed with the session
(DBG_TERMINATE_PROCESS). A Win32_Process.Create launch and then a scheduled task
both died with STATUS_CONTROL_C_EXIT within minutes -- at 11:05:09 and 11:06:42,
each within seconds of an unrelated shell command finishing in the parent
session. A process that owns a console receives that console's control events,
and the campaign was inheriting one whichever way it was started.

DETACHED_PROCESS gives the child no console, so no console control event can
reach it. CREATE_NEW_PROCESS_GROUP keeps it out of the parent's group as well,
and CREATE_BREAKAWAY_FROM_JOB leaves any job object the shell sits inside -- that
flag is attempted and dropped if the job forbids breakaway, which is not an
error worth failing on. Output goes to a file because a console-less process has
nowhere else to write, and stdin is the null device so nothing can block on it.

The multiprocessing workers the campaign spawns inherit the same absence of a
console, so the pool is covered by the same protection as the supervisor.

Usage:
    python launch_detached.py <logfile> <script> [args...]
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200
CREATE_BREAKAWAY_FROM_JOB = 0x01000000
CREATE_NO_WINDOW = 0x08000000


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    log_path = Path(sys.argv[1])
    cmd = [sys.executable, *sys.argv[2:]]

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log = open(log_path, "ab", buffering=0)

    # CREATE_NO_WINDOW first, DETACHED_PROCESS only as the fallback. Both keep
    # the child out of this shell's console so no console control event can
    # reach it, which is the whole point of this file. They differ in what the
    # child's own multiprocessing workers inherit: a DETACHED parent has no
    # console at all, so Windows hands each spawned worker a brand new one and
    # a row of empty console windows appears on the desktop for the length of
    # the run. CREATE_NO_WINDOW gives the parent a console that is never
    # displayed, the workers inherit that, and nothing is drawn.
    for flags in (CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP
                  | CREATE_BREAKAWAY_FROM_JOB,
                  CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP,
                  DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
                  | CREATE_BREAKAWAY_FROM_JOB,
                  DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP):
        try:
            p = subprocess.Popen(cmd, cwd=str(HERE), stdout=log,
                                 stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL,
                                 creationflags=flags, close_fds=True)
        except OSError as exc:                     # job forbids breakaway
            print(f"[launch] flags {flags:#x} refused: {exc}")
            continue
        print(f"[launch] pid={p.pid} flags={flags:#x} log={log_path}")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
