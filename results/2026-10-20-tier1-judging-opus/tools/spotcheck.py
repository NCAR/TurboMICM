#!/usr/bin/env python3
"""Pick citations at random from a page and print each claim beside the cited code."""
import random, re, sys
from pathlib import Path
page, run, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
root = Path.home() / f"NSF-NCAR/tier1-work/verify/{run}/clang/{run}"
text = Path(page).read_text()
cites = []
for m in re.finditer(r"`((?:src|include|tools|tests|examples)/[\w./-]+\.(?:c|h|py)):(\d+)(?:-(\d+))?`", text):
    start = max(0, m.start() - 160)
    cites.append((m.group(1), int(m.group(2)), int(m.group(3) or m.group(2)), text[start:m.start()].replace("\n", " ")))
random.seed(seed)
for f, a, b, ctx in random.sample(cites, min(n, len(cites))):
    p = root / f
    lines = p.read_text(errors="replace").split("\n") if p.exists() else []
    code = "\n      ".join(lines[a - 1:min(b, a + 5)]) if lines else "(missing file)"
    print(f"--- {f}:{a}-{b}\n  CLAIM: ...{ctx[-160:]}\n  CODE: {code}\n")
print(f"[{len(cites)} citations on the page]")
