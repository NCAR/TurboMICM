#!/usr/bin/env python3
"""Generate one prompt per model and effort level from the template ../prompt.v3.

    prompts/make_prompts.py

Each run works in <root>/<model>-<effort>/. Re-run after editing the template,
then run admin/setup_runs.py to copy the prompts into the run directories.
"""

import argparse
from pathlib import Path

MODELS = ["Sonnet", "Opus", "Fable", "Terra", "Sol", "Astra"]
EFFORTS = ["high", "xhigh", "max"]
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--template", type=Path, default=ROOT / "prompt.v3")
parser.add_argument("--root", default=str(ROOT))
parser.add_argument("--micm-commit", default="e15898ce (v3.7.0-259)")
# reference/micm_ts1_ref built with -O3 -DNDEBUG, one thread, 10,000 cells x 100 steps
# of 30 s, make_cells_csv.py seed 2026: median of 3 runs on the M5 run machine.
# (prompt.v2's 54,381 was measured on an M5 Pro with other inputs.)
parser.add_argument("--micm-ms", default="137,383")
args = parser.parse_args()

template = args.template.read_text()
for model in MODELS:
    for effort in EFFORTS:
        run = f"{model.lower()}-{effort}"
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

print(f"Wrote {len(MODELS) * len(EFFORTS)} prompts to {HERE}")
for placeholder in ("<micm_commit>", "<micm_ms>"):
    if placeholder in text:
        print(f"Still unfilled: {placeholder}")
