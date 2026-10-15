#!/usr/bin/env python3
"""Create the run directories, ROOT/<model>-<effort>/, from prompts/*.txt.

    admin/setup_runs.py

Each run directory is its own git repository, holding PROMPT.md, the run's
prompt. A run whose repository already has commits has started and is left
alone; otherwise PROMPT.md is rewritten from the current prompt. A prompt with
an unfilled placeholder is refused.
"""

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLACEHOLDERS = ("<run>", "<model>", "<root>", "<micm_commit>", "<micm_ms>")

prompts = sorted((ROOT / "prompts").glob("*.txt"))
unfilled = {p.name: [h for h in PLACEHOLDERS if h in p.read_text()] for p in prompts}
unfilled = {name: holes for name, holes in unfilled.items() if holes}
if unfilled:
    for name, holes in unfilled.items():
        print(f"{name}: unfilled {' '.join(holes)}; run prompts/make_prompts.py first", file=sys.stderr)
    sys.exit(1)

for prompt in prompts:
    run = ROOT / prompt.stem
    if not re.fullmatch(r"[a-z]+-(high|xhigh|max)", run.name):
        continue
    run.mkdir(exist_ok=True)
    if not (run / ".git").exists():
        subprocess.run(["git", "init", "--quiet", str(run)], check=True)
    started = subprocess.run(["git", "-C", str(run), "rev-parse", "--verify", "--quiet", "HEAD"], capture_output=True)
    if started.returncode == 0:
        print(f"{run.name}: has commits, left alone")
        continue
    (run / "PROMPT.md").write_text(prompt.read_text())
    print(f"{run.name}: ready")
