#!/usr/bin/env python3
"""Create the run directories, ROOT/<run>/, from prompts/*.txt.

    admin/setup_runs.py               every run with a prompt
    admin/setup_runs.py luna-high     only the named runs

Each run directory is its own git repository, holding PROMPT.md, the run's
prompt. A run whose repository already has commits has started and is left
alone. Otherwise PROMPT.md is written from the current prompt, but only if it
differs, so a running session's directory is not touched. A prompt with an
unfilled placeholder is refused.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLACEHOLDERS = ("<run>", "<model>", "<root>", "<micm_commit>", "<micm_ms>")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="*", metavar="RUN", help="runs to set up (default: every run with a prompt)")
args = parser.parse_args()

prompts = sorted((ROOT / "prompts").glob("*.txt"))
if args.runs:
    unknown = set(args.runs) - {p.stem for p in prompts}
    if unknown:
        parser.error(f"no prompt for: {' '.join(sorted(unknown))}")
    prompts = [p for p in prompts if p.stem in args.runs]
unfilled = {p.name: [h for h in PLACEHOLDERS if h in p.read_text()] for p in prompts}
unfilled = {name: holes for name, holes in unfilled.items() if holes}
if unfilled:
    for name, holes in unfilled.items():
        print(f"{name}: unfilled {' '.join(holes)}; run prompts/make_prompts.py first", file=sys.stderr)
    sys.exit(1)

for prompt in prompts:
    run = ROOT / prompt.stem
    if not re.fullmatch(r"[a-z]+-(high|xhigh|max)(-[0-9]+)?", run.name):
        continue
    run.mkdir(exist_ok=True)
    if not (run / ".git").exists():
        subprocess.run(["git", "init", "--quiet", str(run)], check=True)
    started = subprocess.run(["git", "-C", str(run), "rev-parse", "--verify", "--quiet", "HEAD"], capture_output=True)
    if started.returncode == 0:
        print(f"{run.name}: has commits, left alone")
        continue
    text = prompt.read_text()
    target = run / "PROMPT.md"
    if target.exists() and target.read_text() == text:
        print(f"{run.name}: ready, PROMPT.md unchanged")
        continue
    target.write_text(text)
    print(f"{run.name}: ready")
