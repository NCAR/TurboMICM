#!/usr/bin/env python3
"""Make anonymized copies of Tier 2 entries for the blind review.

    admin/anonymize_tier2.py OUT RUN=BRANCH=LABEL [RUN=BRANCH=LABEL ...]

As admin/anonymize_tier3.py does for Tier 3, on the M5 Pro's paths. Each entry's
branch is exported with `git archive` (no history) into OUT/LABEL/, without
PROMPT.md (the prompt is the same for every entry, and names the run). In every
text file, local paths, run names, and model and harness names are replaced,
longest first:
  /Users/dfillmor/NSF-NCAR/ODE/benchmarks/tier2/<run>  ->  ROOT/entry
  /Users/dfillmor/NSF-NCAR/ODE/benchmarks/tier2        ->  ROOT
  /Users/dfillmor and /Users/fillmore                   ->  HOME
  <model>-<effort> run names                            ->  entry
  model and harness names                               ->  the agent
Then every file, text or not, is searched for anything left that could
identify an entry, and the matches are listed (exit status 1 if any).
"""
import io
import re
import subprocess
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ROOT = "/Users/dfillmor/NSF-NCAR/ODE/benchmarks/tier2"
TIER1_ROOT = "/Users/fillmore/NSF-NCAR/ODE/benchmarks"
RUNS = [f"{m}-{e}" for m in ("sol", "astra", "luna", "terra") for e in ("xhigh", "high", "max")]
MODELS = [r"gpt-6\.1-sol", r"gpt-6-astra", r"gpt-6-luna", r"gpt-5\.6-terra", r"gpt-6(\.1)?", r"codex(-cli)?", r"openai"]
LEFT = re.compile(rb"(?i)(fillmore|dfillmor|NSF-NCAR|astra|\bsol\b|\bluna\b|terra|gpt|codex|openai|xhigh|chatgpt)")

out = Path(sys.argv[1])
found = 0
for spec in sys.argv[2:]:
    run, branch, label = spec.split("=")
    dest = out / label
    dest.mkdir(parents=True, exist_ok=False)
    archive = subprocess.run(["git", "-C", str(REPO), "archive", branch], capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        members = [m for m in tar.getmembers() if m.name != "PROMPT.md"]
        tar.extractall(dest, members=members, filter="data")
    rules = [(re.escape(f"{ROOT}/{run}"), "ROOT/entry"), (re.escape(ROOT), "ROOT"),
             (re.escape(f"{TIER1_ROOT}/{run}"), "ROOT/entry"), (re.escape(TIER1_ROOT), "ROOT"),
             (re.escape("/Users/dfillmor"), "HOME"), (re.escape("/Users/fillmore"), "HOME")]
    rules += [(re.escape(r), "entry") for r in sorted(RUNS, key=len, reverse=True)]
    rules += [(m, "the agent") for m in MODELS]
    changed = 0
    for path in sorted(p for p in dest.rglob("*") if p.is_file()):
        data = path.read_bytes()
        if b"\0" not in data[:8192]:
            text = data.decode("utf-8", errors="surrogateescape")
            new = text
            for pattern, replacement in rules:
                new = re.sub(pattern, replacement, new, flags=re.IGNORECASE)
            if new != text:
                path.write_bytes(new.encode("utf-8", errors="surrogateescape"))
                changed += 1
                data = path.read_bytes()
        for m in LEFT.finditer(data):
            found += 1
            context = data[max(0, m.start() - 50):m.end() + 50].decode("utf-8", errors="replace").replace("\n", " ")
            print(f"{label}/{path.relative_to(dest)}: ...{context}...")
    print(f"{label}: {run} ({branch}): {sum(1 for p in dest.rglob('*') if p.is_file())} files, {changed} changed")
sys.exit(1 if found else 0)
