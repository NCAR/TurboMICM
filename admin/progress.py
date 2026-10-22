#!/usr/bin/env python3
"""Summarize a live Codex run, read-only: its session record and its repository.

    admin/progress.py ROLLOUT.jsonl [RUN_DIR [BASE_COMMIT]]

From the session record: when it started and last wrote, the turns ended so
far (task_complete events) with the start of each one's final message, the
context compactions, the tokens (the last token_count event: total = input -
cached + output), and Codex's usage window. From RUN_DIR, if given: the commits
after BASE_COMMIT (default: the root), newest first, and the files changed since
BASE_COMMIT. Git is run with --no-optional-locks, so nothing in the run
directory is written; the record is only read. Needs only the Python standard
library.
"""

import datetime
import json
import subprocess
import sys
from pathlib import Path

if len(sys.argv) not in (2, 3, 4):
    sys.exit(__doc__)
record = Path(sys.argv[1])
run_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None
base = sys.argv[3] if len(sys.argv) > 3 else None


def when(stamp):
    """An ISO timestamp (UTC) as local time."""
    moment = datetime.datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    return moment.astimezone().strftime("%H:%M:%S")


first = last = None
turns, compactions, usage, limits, messages = [], 0, None, None, 0
with open(record, errors="replace") as file:
    for line in file:
        try:
            item = json.loads(line)
        except ValueError:
            continue
        stamp = item.get("timestamp")
        first = first or stamp
        last = stamp or last
        payload = item.get("payload") if isinstance(item.get("payload"), dict) else {}
        kind = payload.get("type")
        if kind == "task_complete":
            turns.append((stamp, (payload.get("last_agent_message") or "").strip()))
        elif kind == "token_count" and payload.get("info"):
            usage = payload["info"].get("total_token_usage")
            limits = payload.get("rate_limits") or limits
        elif item.get("type") == "compacted" or kind == "context_compacted":
            compactions += 1
        elif kind == "message" and payload.get("role") == "assistant":
            messages += 1

idle = (datetime.datetime.now(datetime.timezone.utc)
        - datetime.datetime.fromisoformat(last.replace("Z", "+00:00"))).total_seconds() / 60
print(f"record: {record.name}")
print(f"started {when(first)}, last wrote {when(last)} ({idle:.1f} min ago); {messages} agent messages; "
      f"{len(turns)} turns ended; {compactions} compactions")
for stamp, message in turns[-3:]:
    print(f"  turn ended {when(stamp)}: {message[:300]!r}")
if usage:
    total = usage.get("input_tokens", 0) - usage.get("cached_input_tokens", 0) + usage.get("output_tokens", 0)
    print(f"tokens: total {total:,}; input {usage.get('input_tokens', 0):,} (cached {usage.get('cached_input_tokens', 0):,}); "
          f"output {usage.get('output_tokens', 0):,} (reasoning {usage.get('reasoning_output_tokens', 0):,})")
if limits and limits.get("primary"):
    primary = limits["primary"]
    resets = datetime.datetime.fromtimestamp(primary.get("resets_at", 0)).strftime("%Y-%m-%d %H:%M")
    print(f"usage window: {primary.get('used_percent')} % used (resets {resets})")
if run_dir:
    def git(*args):
        return subprocess.run(["git", "--no-optional-locks", "-C", str(run_dir), *args],
                               capture_output=True, text=True).stdout.strip()
    span = f"{base}..HEAD" if base else "HEAD"
    log = git("log", "--format=%h %ad %s", "--date=format:%H:%M", span)
    commits = log.splitlines() if log else []
    print(f"commits after {base or 'the root'}: {len(commits)}")
    for commit in commits[:8]:
        print(f"  {commit[:150]}")
    if base:
        print("changed since base: " + (git("diff", "--shortstat", base) or "nothing committed or in the working tree"))
