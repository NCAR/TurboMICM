#!/usr/bin/env python3
"""Rank Tier 3 operating points from a timing session of admin/time_entries.py --scenario.

    rank_tier3.py OUT [--threads 10] [--levels 2,3]

OUT holds OUT/<set>/<scenario>/timings.csv, one folder per scenario, as
time_tier3.sh writes them. Programs are named ENTRY/POINT, for example
sol-max/fast or micm/r1e-3-a1e-15. A point reaches level L on a scenario when
every run there exits normally and passes the gate, and its lowest overall
digits are at least L and its lowest worst species at least L - 1, the rule of
score.py --level. time_entries.py records digits to three decimals, so a value
that rounds to a threshold exactly is flagged, to be rescored with score.py
--level. An entry's rank at level L is that of its fastest point that reaches L
on every scenario, by the sum over all scenarios of its median elapsed_ms at
the given thread count; the sums over each set are shown too.
"""

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("out", type=Path)
parser.add_argument("--threads", type=int, default=10)
parser.add_argument("--levels", default="2,3", help="comma-separated levels (default: 2,3)")
args = parser.parse_args()
levels = [float(level) for level in args.levels.split(",")]

scenarios = []  # (set, scenario), in folder order
cells = defaultdict(dict)  # program -> (set, scenario) -> results there
programs = []  # in order of first appearance
for table in sorted(args.out.glob("*/*/timings.csv")):
    key = (table.parent.parent.name, table.parent.name)
    scenarios.append(key)
    runs = defaultdict(list)
    with open(table, newline="") as file:
        for row in csv.DictReader(file):
            if int(row["threads"]) == args.threads:
                runs[row["program"]].append(row)
    for program, rows in runs.items():
        if program not in cells:
            programs.append(program)
        times = [float(row["elapsed_ms"]) for row in rows if row["elapsed_ms"]]
        scored = [row for row in rows if row["digits"] != ""]
        worst = min(scored, key=lambda row: float(row["worst_species"].split()[-1]), default=None)
        cells[program][key] = {
            "runs": len(rows),
            "gate": all(row["gate"] == "PASS" for row in rows) and len(times) == len(rows),
            "digits": min((float(row["digits"]) for row in scored), default=math.nan),
            "worst": float(worst["worst_species"].split()[-1]) if worst else math.nan,
            "worst_species": worst["worst_species"].split()[0] if worst else "",
            "median_s": statistics.median(times) / 1000 if times else math.nan,
            "outputs": len({row["output_sha256"] for row in rows}),
        }
sets = sorted({name for name, _ in scenarios})


def complete(program):
    return all(key in cells[program] for key in scenarios)


def reaches(program, level):
    return complete(program) and all(
        cell["gate"] and cell["digits"] >= level and cell["worst"] >= level - 1 for cell in cells[program].values()
    )


def borderline(program, level):
    return [
        f"{key[0]}/{key[1]}"
        for key, cell in cells[program].items()
        if abs(cell["digits"] - level) < 0.0006 or abs(cell["worst"] - (level - 1)) < 0.0006
    ]


def seconds(program, which=None):
    return sum(cell["median_s"] for key, cell in cells[program].items() if which in (None, key[0]))


def entry(program):
    return program.split("/")[0]


lines = [f"# Tier 3 ranking: {len(scenarios)} scenarios, {args.threads} threads", ""]
lines += ["Time: the sum over scenarios of each point's median elapsed_ms, in seconds.", ""]
flags = []
for level in levels:
    best = {}
    for program in programs:
        if reaches(program, level) and (entry(program) not in best or seconds(program) < seconds(best[entry(program)])):
            best[entry(program)] = program
        if reaches(program, level) and borderline(program, level):
            flags.append(f"{program} at level {level:g}: rounds to a threshold on {', '.join(borderline(program, level))}")
    micm = best.get("micm")
    lines += [f"## Level {level:g}", ""]
    lines.append("| Rank | Entry | Fastest point reaching it everywhere | Time (s) | " + " | ".join(f"{s} (s)" for s in sets)
                 + " | Speedup over MICM |")
    lines.append("|---:|---|---|---:|" + "---:|" * len(sets) + "---:|")
    ranked = sorted(best.values(), key=seconds)
    for rank, program in enumerate(ranked, 1):
        speedup = f"{seconds(micm) / seconds(program):.2f}" if micm else ""
        lines.append(f"| {rank} | {entry(program)} | {program.split('/', 1)[-1]} | {seconds(program):,.2f} | "
                     + " | ".join(f"{seconds(program, s):,.2f}" for s in sets) + f" | {speedup} |")
    for name in dict.fromkeys(entry(program) for program in programs):
        if name not in best:
            lines.append(f"| | {name} | none | | " + " | " * len(sets) + " |")
    lines.append("")

lines += ["## Every point", ""]
lines.append("| Point | " + " | ".join(f"Level {level:g}" for level in levels)
             + " | Lowest overall | Lowest worst species | Time (s) | " + " | ".join(f"{s} (s)" for s in sets)
             + " | Runs | Scenarios with differing outputs |")
lines.append("|---|" + "---|" * len(levels) + "---:|---|---:|" + "---:|" * len(sets) + "---:|---:|")
for program in programs:
    mine = cells[program]
    low = min(mine.values(), key=lambda cell: cell["worst"])
    gate = "" if all(cell["gate"] for cell in mine.values()) else " (gate failed)"
    lines.append(
        f"| {program}{gate}{'' if complete(program) else ' (incomplete)'} | "
        + " | ".join("yes" if reaches(program, level) else "no" for level in levels)
        + f" | {min(cell['digits'] for cell in mine.values()):.3f} | {low['worst_species']} {low['worst']:.3f}"
        + f" | {seconds(program):,.2f} | " + " | ".join(f"{seconds(program, s):,.2f}" for s in sets)
        + f" | {sum(cell['runs'] for cell in mine.values())} | {sum(cell['outputs'] > 1 for cell in mine.values())} |"
    )
lines.append("")

lines += ["## Median elapsed time by scenario (s)", ""]
lines.append("| Point | " + " | ".join(f"{s}/{name}" for s, name in scenarios) + " |")
lines.append("|---|" + "---:|" * len(scenarios))
for program in programs:
    lines.append(f"| {program} | " + " | ".join(
        f"{cells[program][key]['median_s']:,.3f}" if key in cells[program] else "" for key in scenarios) + " |")
lines.append("")

if flags:
    lines += ["## Rescore exactly", ""] + [f"- {flag}" for flag in flags] + [""]
print("\n".join(lines))
