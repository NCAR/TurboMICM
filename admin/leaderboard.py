#!/usr/bin/env python3
"""Rank every timed program, MICM configurations and entries alike, from timings.csv files.

    admin/leaderboard.py [--threads 10] [--out DIR] TIMINGS_CSV [TIMINGS_CSV ...]

Each TIMINGS_CSV is admin/time_entries.py output. Runs of the same program are
pooled across files; a program's time is the median of its runs at --threads.
MICM programs are named micm-<layout>-<storage>-<algorithm>[-g32], as
admin/micm_configs.sh names them; any other program is an entry.

Writes, to --out (default: the current directory):
- leaderboard.md: one row per program, fastest first, with its median and range,
  the number of runs, its speed relative to MICM's default (standard ordering,
  in-place Mozart LU, one group per thread) and to the fastest MICM
  configuration, the gate results, and the largest relative difference from
  MICM's default result. A program with fewer than three runs is marked
  preliminary.
- matrices.md: MICM's configurations as matrix tables, layout by LU variant, one
  per grouping.
Needs only the Python standard library.
"""

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

LAYOUTS = ["standard", "vector1", "vector2", "vector4", "vector8", "vector128"]
LUS = [("in-place", "mozart"), ("in-place", "doolittle"), ("separate", "mozart"), ("separate", "doolittle")]
DEFAULT = "micm-standard-in-place-mozart"
# Each entry's own cell grouping, from its DESIGN.md and benchmark source.
ENTRY_GROUPING = {"sol-high": "32 cells", "sol-xhigh": "32 cells", "sol-max": "32 cells", "astra-high": "32 cells"}
PATTERN = re.compile(r"micm-(standard|vector\d+)-(in-place|separate)-(mozart|doolittle)(-g32)?$")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("timings", nargs="+", type=Path)
parser.add_argument("--threads", type=int, default=10)
parser.add_argument("--out", type=Path, default=Path("."))
args = parser.parse_args()

runs = defaultdict(list)
for path in args.timings:
    with open(path, newline="") as file:
        for row in csv.DictReader(file):
            if int(row["threads"]) == args.threads and row["elapsed_ms"]:
                runs[row["program"]].append(row)


def median(program):
    return statistics.median(float(row["elapsed_ms"]) for row in runs[program])


def describe(program):
    match = PATTERN.match(program)
    if not match:
        return "entry", program, ""
    layout, storage, algorithm, g32 = match.groups()
    grouping = "32-cell groups" if g32 else "one group per thread"
    return "MICM", f"{layout}, {storage} {algorithm.capitalize()}", grouping


default_ms = median(DEFAULT)
micm_programs = [p for p in runs if PATTERN.match(p)]
fastest_micm = min(micm_programs, key=median)
fastest_ms = median(fastest_micm)

lines = [
    f"| Rank | Entrant | Configuration | Cell groups | Median ms (range) | Runs | vs MICM default | vs fastest MICM | Gate | Largest relative difference |",
    "|---:|---|---|---|---:|---:|---:|---:|---|---:|",
]
for rank, program in enumerate(sorted(runs, key=median), 1):
    kind, configuration, grouping = describe(program)
    times = [float(row["elapsed_ms"]) for row in runs[program]]
    gates = sorted({row["gate"] for row in runs[program]})
    largest = [float(row["max_relative_difference"]) for row in runs[program] if row["max_relative_difference"]]
    name = program if kind == "entry" else "MICM"
    if kind == "entry":
        configuration = "the agent's C port"
        grouping = ENTRY_GROUPING.get(program, "its own")
    else:
        grouping = "32 cells" if grouping == "32-cell groups" else "one per thread"
    flag = "" if len(times) >= 3 else " (preliminary)"
    largest_text = f"{max(largest):.1g}" if largest else ""
    lines.append(
        f"| {rank} | {'**' + name + '**' if kind == 'entry' else name} | {configuration} | {grouping} |"
        f" {median(program):,.0f} ({min(times):,.0f}–{max(times):,.0f}){flag} | {len(times)} |"
        f" {default_ms / median(program):.2f}× | {fastest_ms / median(program):.2f}× | {', '.join(gates)} |"
        f" {largest_text} |"
    )
(args.out / "leaderboard.md").write_text("\n".join(lines) + "\n")

matrix = []
header = "| Matrix layout | " + " | ".join(
    f"{storage} {algorithm.capitalize()}" + (" (default)" if (storage, algorithm) == ("in-place", "mozart") else "")
    for storage, algorithm in LUS) + " |"
for suffix, title in (("", "one cell group per thread"), ("-g32", "32-cell groups")):
    matrix += [f"Median ms at {args.threads} threads (speed relative to MICM's default), {title}:", "", header,
               "|---|" + "---:|" * len(LUS)]
    for layout in LAYOUTS:
        cells = []
        for storage, algorithm in LUS:
            program = f"micm-{layout}-{storage}-{algorithm}{suffix}"
            if program not in runs:
                cells.append("")
                continue
            flag = "" if len(runs[program]) >= 3 else "*"
            bold = "**" if program == fastest_micm else ""
            cells.append(f"{bold}{median(program):,.0f} ({default_ms / median(program):.2f}×){bold}{flag}")
        matrix.append(f"| {layout} | " + " | ".join(cells) + " |")
    matrix.append("")
(args.out / "matrices.md").write_text("\n".join(matrix))
print("\n".join(lines))
print()
print("\n".join(matrix))
