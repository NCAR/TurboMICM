#!/usr/bin/env python3
"""Generate the run prompts from the template ../prompt.v3.

    prompts/make_prompts.py                         every run, first instance
    prompts/make_prompts.py --instance 2 sol-high   a repeat: run sol-high-2

Each run works in <root>/<run>/, where <run> is <model>-<effort>, plus -<n> for a
repeat. Re-run after editing the template, then run admin/setup_runs.py to copy
the prompts into the run directories.
"""

import argparse
from pathlib import Path

MODELS = ["Sonnet", "Opus", "Fable", "Terra", "Sol", "Astra"]
EFFORTS = ["high", "xhigh", "max"]
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("runs", nargs="*", metavar="RUN", help="<model>-<effort> runs to write (default: all)")
parser.add_argument("--instance", type=int, default=1, help="repeat number; 2 and up add a -<n> suffix")
parser.add_argument("--template", type=Path, default=ROOT / "prompt.v3")
parser.add_argument("--root", default=str(ROOT))
parser.add_argument("--micm-commit", default="e15898ce (v3.7.0-259)")
# reference/micm_ts1_ref built with -O3 -DNDEBUG, one thread, 10,000 cells x 100 steps
# of 30 s, make_cells_csv.py seed 2026: median of 3 runs on the M5 run machine.
# (prompt.v2's 54,381 was measured on an M5 Pro with other inputs.)
parser.add_argument("--micm-ms", default="137,383")
args = parser.parse_args()

pairs = [(model, effort) for model in MODELS for effort in EFFORTS]
if args.runs:
    unknown = set(args.runs) - {f"{model.lower()}-{effort}" for model, effort in pairs}
    if unknown:
        parser.error(f"unknown runs: {' '.join(sorted(unknown))}")
    pairs = [(model, effort) for model, effort in pairs if f"{model.lower()}-{effort}" in args.runs]
if args.instance < 1:
    parser.error("--instance must be 1 or more")

template = args.template.read_text()
for model, effort in pairs:
    run = f"{model.lower()}-{effort}" + (f"-{args.instance}" if args.instance > 1 else "")
    text = template
    for placeholder, value in {
        "<run>": run,
        "<model>": model,
        "<root>": args.root,
        "<micm_commit>": args.micm_commit,
        "<micm_ms>": args.micm_ms,
    }.items():
        if value is not None:
            text = text.replace(placeholder, value)
    (HERE / f"{run}.txt").write_text(text)

print(f"Wrote {len(pairs)} prompt{'' if len(pairs) == 1 else 's'} to {HERE}")
for placeholder in ("<micm_commit>", "<micm_ms>"):
    if placeholder in text:
        print(f"Still unfilled: {placeholder}")
