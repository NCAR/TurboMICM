#!/usr/bin/env python3
"""Create the Tier 3 run directories, ROOT/<run>/, from prompts/*.txt.

    admin/setup_runs.py               every run with a prompt
    admin/setup_runs.py sol-max       only the named runs

ROOT is this Tier 3 root. First, ROOT/micm and ROOT/reference are made links to
the Tier 1 root's, so the entries' Makefiles, which use ../micm and ../reference,
still work. Then, for each prompt:
- the Tier 1 run, ROOT/../<run>, is cloned into ROOT/<run> with its history;
- the clone must be at the Tier 1 commit the prompt names;
- the clone's remote is removed, so nothing can be pushed back into Tier 1;
- the Tier 3 prompt is written to PROMPT.md.
A Tier 3 run with commits beyond its Tier 1 commit has started and is left alone.
Otherwise PROMPT.md is written only if it differs, so a running session's
directory is not touched. A prompt with an unfilled placeholder is refused.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PLACEHOLDERS = ("<run>", "<model>", "<root>", "<micm_commit>", "<tier1_commit>")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="*", metavar="RUN", help="runs to set up (default: every run with a prompt)")
args = parser.parse_args()


def git(*command, check=True):
    return subprocess.run(["git", *command], capture_output=True, text=True, check=check)


for name in ("micm", "reference"):
    link = ROOT / name
    if not link.exists():
        os.symlink(f"../{name}", link)
        print(f"{name}: linked to ../{name}")

prompts = sorted((ROOT / "prompts").glob("*.txt"))
if args.runs:
    unknown = set(args.runs) - {p.stem for p in prompts}
    if unknown:
        parser.error(f"no prompt for: {' '.join(sorted(unknown))}")
    prompts = [p for p in prompts if p.stem in args.runs]
failed = False
for prompt in prompts:
    text = prompt.read_text()
    holes = [h for h in PLACEHOLDERS if h in text]
    if holes:
        print(f"{prompt.name}: unfilled {' '.join(holes)}; run prompts/make_prompts.py first", file=sys.stderr)
        failed = True
        continue
    found = re.search(r"repository\s+at\s+commit\s+([0-9a-f]{7,40})\b", text)
    if not found:
        print(f"{prompt.name}: cannot find the Tier 1 commit in the prompt", file=sys.stderr)
        failed = True
        continue
    tier1_commit = found.group(1)
    run = ROOT / prompt.stem
    tier1 = ROOT.parent / prompt.stem
    if not run.exists():
        if not (tier1 / ".git").exists():
            print(f"{run.name}: the Tier 1 run {tier1} is not a git repository", file=sys.stderr)
            failed = True
            continue
        # --no-hardlinks: the clone shares no files with the Tier 1 repository
        git("clone", "--quiet", "--no-hardlinks", str(tier1), str(run))
        git("-C", str(run), "remote", "remove", "origin")
        print(f"{run.name}: cloned from the Tier 1 run")
    elif not (run / ".git").exists():
        print(f"{run.name}: exists but is not a git repository; check it", file=sys.stderr)
        failed = True
        continue
    head = git("-C", str(run), "rev-parse", "HEAD").stdout.strip()
    tier1_full = git("-C", str(run), "rev-parse", f"{tier1_commit}^{{commit}}", check=False).stdout.strip()
    if not tier1_full:
        print(f"{run.name}: Tier 1 commit {tier1_commit} is not in the clone", file=sys.stderr)
        failed = True
        continue
    if head != tier1_full:
        later = git("-C", str(run), "rev-list", "--count", f"{tier1_full}..HEAD").stdout.strip()
        if later != "0" and git("-C", str(run), "merge-base", "--is-ancestor", tier1_full, "HEAD", check=False).returncode == 0:
            print(f"{run.name}: has {later} Tier 3 commit(s), left alone")
            continue
        print(f"{run.name}: HEAD is not the Tier 1 commit {tier1_commit}; check the clone", file=sys.stderr)
        failed = True
        continue
    if git("-C", str(run), "remote").stdout.strip():
        print(f"{run.name}: has a remote; remove it before launch", file=sys.stderr)
        failed = True
        continue
    target = run / "PROMPT.md"
    if target.exists() and target.read_text() == text:
        print(f"{run.name}: ready at Tier 1 commit {tier1_commit}, PROMPT.md unchanged")
        continue
    target.write_text(text)
    print(f"{run.name}: ready at Tier 1 commit {tier1_commit}")
sys.exit(1 if failed else 0)
