#!/usr/bin/env python3
"""Generate the Tier 2 run prompts from the template ../prompt.v1.

    prompts/make_prompts.py                every planned Tier 2 run
    prompts/make_prompts.py luna-max       only the named runs

Each Tier 2 run works in <root>/<run>/, a clone of its Tier 1 run <root>/../<run>/
at that run's last commit, which fills <tier1_commit>. Generate a run's prompt only
after its Tier 1 run has finished. A run whose Tier 1 repository has no commits is
skipped. Re-run after editing the template, then run admin/setup_runs.py.
"""

import argparse
import re
import subprocess
from pathlib import Path

PLANNED = ["sol-xhigh", "sol-max", "astra-xhigh", "astra-max"]
MODELS = {"sonnet": "Sonnet", "opus": "Opus", "fable": "Fable", "luna": "Luna", "sol": "Sol", "astra": "Astra"}
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="*", metavar="RUN", help=f"<model>-<effort> runs (default: {' '.join(PLANNED)})")
parser.add_argument("--template", type=Path, default=ROOT / "prompt.v1")
parser.add_argument("--root", default=str(ROOT))
parser.add_argument("--micm-commit", default="e15898ce (v3.7.0-259)")
args = parser.parse_args()

template = args.template.read_text()
written = 0
for run in args.runs or PLANNED:
    match = re.fullmatch(r"([a-z]+)-(high|xhigh|max)(-[0-9]+)?", run)
    if not match or match.group(1) not in MODELS:
        parser.error(f"not a run name: {run}")
    tier1 = ROOT.parent / run
    if not (tier1 / ".git").exists():
        print(f"{run}: skipped, its Tier 1 run ({tier1}) is not a git repository")
        continue
    head = subprocess.run(["git", "-C", str(tier1), "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    if head.returncode != 0:
        print(f"{run}: skipped, its Tier 1 run ({tier1}) has no commits yet")
        continue
    # --no-optional-locks: never write the index of a repository a session may be using
    dirty = subprocess.run(["git", "--no-optional-locks", "-C", str(tier1), "status", "--porcelain"],
                           capture_output=True, text=True).stdout
    if dirty.strip():
        print(f"{run}: note, its Tier 1 working tree has uncommitted changes; Tier 2 starts from the last commit")
    text = template
    for placeholder, value in {
        "<run>": run,
        "<model>": MODELS[match.group(1)],
        "<root>": args.root,
        "<micm_commit>": args.micm_commit,
        "<tier1_commit>": head.stdout.strip(),
    }.items():
        text = text.replace(placeholder, value)
    (HERE / f"{run}.txt").write_text(text)
    written += 1
    commit = subprocess.run(["git", "-C", str(tier1), "log", "-1", "--format=%ad %s", "--date=format:%Y-%m-%d %H:%M"],
                            capture_output=True, text=True).stdout.strip()
    print(f"{run}: Tier 1 commit {head.stdout.strip()}, {commit}")
print(f"Wrote {written} prompt{'' if written == 1 else 's'} to {HERE}")
