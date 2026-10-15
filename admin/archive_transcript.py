#!/usr/bin/env python3
"""Archive an agent's session transcript: redact account identifiers and gzip it.

    admin/archive_transcript.py SESSION.jsonl OUTPUT.jsonl.gz [--redact KEY ...]

SESSION.jsonl is a session file, one JSON object per line, such as Codex's
~/.codex/sessions/.../rollout-*.jsonl. The values of the redacted keys (by
default Codex's creator_account_id and creator_user_id) are collected from the
whole file, and every occurrence of those values is replaced by "<redacted>".
Lines without them are kept byte for byte. Every line that parsed before must
still parse; the result is written with gzip -9. Needs only the Python standard
library.
"""

import argparse
import gzip
import json
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("session", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument(
    "--redact",
    action="append",
    metavar="KEY",
    help="a key whose values are redacted (default: creator_account_id, creator_user_id)",
)
args = parser.parse_args()
keys = set(args.redact or ["creator_account_id", "creator_user_id"])


def values_of(obj):
    """Every string value stored under a redacted key, at any depth."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key in keys and isinstance(value, str):
                yield value
            yield from values_of(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from values_of(value)


lines = args.session.read_bytes().split(b"\n")
parsed = []
secrets = set()
for line in lines:
    try:
        obj = json.loads(line)
    except ValueError:
        parsed.append(False)
        continue
    parsed.append(True)
    secrets.update(value for value in values_of(obj) if len(value) >= 8)

redacted, changed = [], 0
for line, was_json in zip(lines, parsed):
    new = line
    for secret in secrets:
        new = new.replace(secret.encode(), b"<redacted>")
    if new != line:
        changed += 1
        if was_json:
            json.loads(new)  # still valid JSON
    redacted.append(new)

data = b"\n".join(redacted)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_bytes(gzip.compress(data, 9))
print(
    f"{args.output}: {len(lines)} lines, {len(secrets)} values redacted in {changed} lines,"
    f" {len(data) / 1e6:.1f} MB -> {args.output.stat().st_size / 1e6:.1f} MB",
    file=sys.stderr,
)
