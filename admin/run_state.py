#!/usr/bin/env python3
"""Report whether a Codex /goal session has finished, from its session record.

    admin/run_state.py ROLLOUT.jsonl PID

ROLLOUT.jsonl is the session record, ~/.codex/sessions/<date>/rollout-*-<id>.jsonl,
whose first line (session_meta) names the run directory as its cwd. PID is the
`codex` process whose working directory is the run directory (`lsof -a -d cwd -p`).

Under /goal, Codex works in turns: each turn ends with a task_complete event, and
the goal starts the next turn at once (in the finished runs' records, within the
same second; no record was ever silent for more than about 3 minutes). The record
has no event of its own for a finished goal. So this script exits 0, printing one
line, when the run reaches a terminal state:

  DONE     the last event is task_complete and the record has been idle 10 minutes
  ERROR    the last event is an error or an abort, and the record has been idle
           10 minutes
  EXITED   the codex process has ended
  STALLED  the record has been idle 40 minutes, whatever its last event

and exits 1 while the run is still going. It only reads the end of the record and
tests that the process exists (signal 0, which sends nothing), so it can watch a
live run without touching it:

    until admin/run_state.py ROLLOUT.jsonl PID; do sleep 60; done

Needs only the Python standard library.
"""

import json
import os
import sys
import time

if len(sys.argv) != 3:
    sys.exit(__doc__)
path, pid = sys.argv[1], int(sys.argv[2])

try:
    os.kill(pid, 0)
    alive = True
except OSError:
    alive = False

st = os.stat(path)
idle = time.time() - st.st_mtime
with open(path, "rb") as f:
    f.seek(max(0, st.st_size - 4_000_000))
    lines = f.read().decode("utf-8", "replace").splitlines()

last = {}
for line in reversed(lines):
    try:
        last = json.loads(line)
        break
    except ValueError:
        continue
payload = last.get("payload") if isinstance(last.get("payload"), dict) else {}
kind = payload.get("type") or last.get("type", "?")
summary = (payload.get("last_agent_message") or payload.get("message") or "")[:400]
info = f"last={last.get('timestamp')} type={kind} idle={idle / 60:.1f}min :: {summary!r}"

state = None
if not alive:
    state = "EXITED"
elif idle >= 600 and kind == "task_complete":
    state = "DONE"
elif idle >= 600 and ("error" in kind or "abort" in kind):
    state = "ERROR"
elif idle >= 2400:
    state = "STALLED"

if state:
    print(state, info)
    sys.exit(0)
sys.exit(1)
