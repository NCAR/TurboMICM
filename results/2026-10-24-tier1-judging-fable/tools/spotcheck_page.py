#!/usr/bin/env python3
"""Sample citations from the judge's page and print each beside the cited code.

The page cites several entries; the entry is taken from the nearest preceding
bold run name (**run, N.**) or table row (| run |) or heading. Usage:
spotcheck_page.py PAGE N SEED
"""
import random, re, sys
from pathlib import Path
page, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
RUNS = ["sol-high","sol-xhigh","sol-max","astra-high","astra-xhigh","astra-max","luna-max","luna-xhigh","luna-high"]
text = Path(page).read_text()
cites = []
for m in re.finditer(r"`((?:src|include|tools|tests|examples)/[\w./-]+\.(?:c|h|py|inc)):(\d+)(?:-(\d+))?`", text):
    before = text[:m.start()]
    # entry: last run name mentioned on the same line or in a bold lead / table cell before
    line_start = before.rfind("\n") + 1
    line = text[line_start:m.start()]
    run = None
    for r in RUNS:
        if re.search(r"(\*\*|\| |^)" + re.escape(r) + r"(\*\*|,| \||\b)", line):
            run = r
    if run is None:
        prior = [(before.rfind("**" + r + ","), r) for r in RUNS] + [(before.rfind("\n| " + r + " |"), r) for r in RUNS] + [(before.rfind("#### " + r), r) for r in RUNS]
        run = max(prior)[1]
    cites.append((run, m.group(1), int(m.group(2)), int(m.group(3) or m.group(2)), line[-120:]))
random.seed(seed)
for run, f, a, b, ctx in random.sample(cites, min(n, len(cites))):
    p = Path.home() / f"NSF-NCAR/tier1-work/verify/{run}/clang/{run}" / f
    lines = p.read_text(errors="replace").split("\n") if p.exists() else []
    code = "\n      ".join(lines[a-1:min(b, a+4)]) if lines else "(missing file)"
    print(f"--- [{run}] {f}:{a}-{b}\n  CLAIM: ...{ctx}\n  CODE: {code}\n")
print(f"[{len(cites)} citations on the page]")
