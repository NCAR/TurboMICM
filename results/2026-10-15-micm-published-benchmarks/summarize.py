#!/usr/bin/env python3
"""Print the latest entry of each chart series in data/, with each value's history.

    summarize.py data > latest.txt

The data/*.js files are MICM's published benchmark charts (github-action-benchmark
format), fetched from the gh-pages branch of NCAR/micm. Needs only the Python
standard library.
"""

import datetime
import json
import statistics
import sys
from pathlib import Path

for path in sorted(Path(sys.argv[1]).glob("*.js")):
    text = path.read_text()
    data = json.loads(text[text.index("{"):].rstrip().rstrip(";"))
    for series, entries in data["entries"].items():
        last = entries[-1]
        when = datetime.datetime.fromtimestamp(last["date"] / 1000).strftime("%Y-%m-%d")
        print(f"\n== {path.stem}: '{series}', {len(entries)} entries, latest {when} commit "
              f"{last['commit']['id'][:8]}: {last['commit']['message'].splitlines()[0][:70]}")
        for bench in last["benches"]:
            history = [x["value"] for e in entries for x in e["benches"] if x["name"] == bench["name"]]
            print(f"   {bench['name']:<22} {bench['value']:>16,.2f} {bench['unit']:<12} (history n={len(history)}, "
                  f"median {statistics.median(history):,.2f}, min {min(history):,.2f}, max {max(history):,.2f})")
