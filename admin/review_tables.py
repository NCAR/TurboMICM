#!/usr/bin/env python3
"""Decode Tier 2's blind-review reports and print the scores and ranks.

    admin/review_tables.py [BLIND_REVIEW_DIR] [--labels opus,opus-2,fable]

BLIND_REVIEW_DIR (default ~/NSF-NCAR/tier2-work/blind-review) holds key.txt
(letter, run) and each reviewer's report-LABEL.json from
tier2/review/run-review.sh. For every report that exists and is not an error,
the script reads the JSON block that ends the report (REVIEW.md's format) and
prints, decoded with the key: each entry's three scores and their sum per
reviewer and in total, sorted by the total; each reviewer's rank and exercise
difficulties; and each session's turns, minutes, list-price cost, and output
tokens. Needs only the Python standard library.
"""

import argparse
import json
import re
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", nargs="?", type=Path, default=Path.home() / "NSF-NCAR/tier2-work/blind-review")
parser.add_argument("--labels", default="opus,opus-2,fable")
args = parser.parse_args()

key = dict(line.split() for line in (args.dir / "key.txt").read_text().splitlines() if line.strip())
reports = []
for label in args.labels.split(","):
    path = args.dir / f"report-{label}.json"
    if not path.exists():
        print(f"{label}: no report")
        continue
    session = json.loads(path.read_text())
    if session.get("is_error"):
        print(f"{label}: error {session.get('api_error')}")
        continue
    blocks = re.findall(r"```json\s*(\{.*?\})\s*```", session["result"], re.S)
    reports.append((label, json.loads(blocks[-1]), session))

total = {letter: sum(sum(r["scores"][letter]) for _, r, _ in reports) for letter in key}
print("| Entry | " + " | ".join(label for label, _, _ in reports) + " | Total |")
print("|---|" + "---:|" * (len(reports) + 1))
for letter in sorted(key, key=lambda letter: -total[letter]):
    cells = [f"{'/'.join(map(str, r['scores'][letter]))} = {sum(r['scores'][letter])}" for _, r, _ in reports]
    print(f"| {key[letter]} ({letter}) | " + " | ".join(cells) + f" | {total[letter]} |")
print()
for label, r, session in reports:
    print(f"{label}: rank {', '.join(key[letter] for letter in r['rank'])}")
    print(f"  exercises (a, b, c): " + "; ".join(f"{key[letter]} {'/'.join(v)}" for letter, v in r["exercises"].items()))
    print(f"  {session.get('num_turns')} turns, {session.get('duration_ms', 0) / 60000:.1f} min, "
          f"${session.get('total_cost_usd', 0):.2f} at list prices, {session['usage'].get('output_tokens'):,} output tokens, "
          f"session {session.get('session_id')}")
