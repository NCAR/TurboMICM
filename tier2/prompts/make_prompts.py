#!/usr/bin/env python3
"""Generate the Tier 2 run prompts from the template ../prompt.v2.

    prompts/make_prompts.py                every planned Tier 2 run
    prompts/make_prompts.py luna-max       only the named runs

Each Tier 2 run works in <root>/<run>/, a copy of its Tier 1 run at that run's
last commit, which fills <tier1_commit>. The Tier 1 runs ran on another machine;
here their repositories are the TurboMICM branches runs/<model>/<effort>, which
hold each run's history unchanged (--repo, default ~/NSF-NCAR/TurboMICM; its
local branch if it has one, otherwise origin's). Re-run after editing the
template, then run admin/setup_runs.py.
"""

import argparse
import re
import subprocess
from pathlib import Path

PLANNED = ["astra-max", "astra-xhigh", "astra-high", "sol-max", "sol-xhigh", "sol-high",
           "luna-max", "luna-xhigh", "luna-high"]
MODELS = {"sonnet": "Sonnet", "opus": "Opus", "fable": "Fable", "luna": "Luna", "sol": "Sol", "astra": "Astra"}
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="*", metavar="RUN", help=f"<model>-<effort> runs (default: {' '.join(PLANNED)})")
parser.add_argument("--template", type=Path, default=ROOT / "prompt.v2")
parser.add_argument("--root", default=str(ROOT))
parser.add_argument("--repo", type=Path, default=Path.home() / "NSF-NCAR" / "TurboMICM")
parser.add_argument("--micm-commit", default="e15898ce (v3.7.0-259)")
args = parser.parse_args()


def tier1_ref(run):
    """The ref of run's Tier 1 branch in --repo, or None."""
    model, effort = run.split("-", 1)
    for ref in (f"refs/heads/runs/{model}/{effort}", f"refs/remotes/origin/runs/{model}/{effort}"):
        found = subprocess.run(["git", "-C", str(args.repo), "rev-parse", "-q", "--verify", f"{ref}^{{commit}}"],
                               capture_output=True, text=True)
        if found.returncode == 0:
            return ref
    return None


template = args.template.read_text()
written = 0
for run in args.runs or PLANNED:
    match = re.fullmatch(r"([a-z]+)-(high|xhigh|max)(-[0-9]+)?", run)
    if not match or match.group(1) not in MODELS:
        parser.error(f"not a run name: {run}")
    ref = tier1_ref(run)
    if ref is None:
        print(f"{run}: skipped, {args.repo} has no Tier 1 branch for it")
        continue
    head = subprocess.run(["git", "-C", str(args.repo), "rev-parse", "--short", ref], capture_output=True, text=True)
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
    commit = subprocess.run(["git", "-C", str(args.repo), "log", "-1", "--format=%ad %s", "--date=format:%Y-%m-%d %H:%M", ref],
                            capture_output=True, text=True).stdout.strip()
    print(f"{run}: Tier 1 commit {head.stdout.strip()} ({ref}), {commit}")
print(f"Wrote {written} prompt{'' if written == 1 else 's'} to {HERE}")
