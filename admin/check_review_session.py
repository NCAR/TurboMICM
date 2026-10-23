#!/usr/bin/env python3
"""Check that a blind reviewer's Claude Code session stayed inside its folder.

    admin/check_review_session.py SESSION_ID [--folder PATH]

Finds the session's transcript under ~/.claude/projects/ and lists every tool
call by name, then any absolute path outside the review folder (default
~/NSF-NCAR/tier2-work/blind-review/REVIEW) and any relative `../` in a tool
call's input. Exits 1 if it finds either. Needs only the Python standard
library.
"""

import argparse
import collections
import json
import re
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("session")
parser.add_argument("--folder", default=str(Path.home() / "NSF-NCAR/tier2-work/blind-review/REVIEW"))
args = parser.parse_args()

matches = list((Path.home() / ".claude/projects").rglob(f"{args.session}.jsonl"))
if not matches:
    sys.exit(f"no transcript for session {args.session}")
calls = collections.Counter()
outside = []
for line in matches[0].open():
    try:
        item = json.loads(line)
    except json.JSONDecodeError:
        continue
    message = item.get("message", {})
    if item.get("type") != "assistant" or not isinstance(message.get("content"), list):
        continue
    for block in message["content"]:
        if block.get("type") != "tool_use":
            continue
        calls[block["name"]] += 1
        text = json.dumps(block["input"])
        for path in re.findall(r"(/(?:Users|private|tmp|etc|opt|var|Volumes)[^\"\s']*)", text):
            if not path.startswith(args.folder):
                outside.append((block["name"], path))
        if "../" in text:
            outside.append((block["name"], "relative ../: " + text[:160]))
print(f"transcript: {matches[0]}")
print("tool calls:", dict(calls))
for name, what in outside:
    print(f"OUTSIDE {name}: {what}")
print("outside the folder:", len(outside))
sys.exit(1 if outside else 0)
