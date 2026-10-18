#!/usr/bin/env python3
"""Summarize a timing session: medians and ratios against reference programs.

    derecho/tools/ratios.py TIMINGS_DIR [--against NAME ...] [--markdown]

TIMINGS_DIR holds admin/time_entries.py's timings.csv. For each thread count
and program: runs passing the gate, median elapsed_ms, range, and the ratio of
each --against program's median to it (how many times faster than that program;
default: every program whose name starts with "micm"). Needs only the standard
library.
"""
import argparse
import csv
import statistics
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("timings_dir", type=Path)
parser.add_argument("--against", nargs="*", help="reference programs (default: those named micm*)")
parser.add_argument("--markdown", action="store_true")
args = parser.parse_args()

with open(args.timings_dir / "timings.csv", newline="") as file:
    rows = list(csv.DictReader(file))
programs = list(dict.fromkeys(row["program"] for row in rows))
threads = sorted({int(row["threads"]) for row in rows}, reverse=True)
against = args.against if args.against else [p for p in programs if p.startswith("micm")]


def median(program, t):
    runs = [row for row in rows if row["program"] == program and int(row["threads"]) == t]
    good = [float(row["elapsed_ms"]) for row in runs if row["gate"] == "PASS" and row["elapsed_ms"]]
    return (statistics.median(good), min(good), max(good), len(good), len(runs)) if good else (None, None, None, 0, len(runs))


for t in threads:
    refs = {name: median(name, t)[0] for name in against}
    header = ["program", "threads", "gate", "median ms", "range"] + [f"vs {name}" for name in against]
    lines = []
    for program in programs:
        med, lo, hi, good, total = median(program, t)
        if med is None:
            lines.append([program, str(t), f"{good}/{total}", "", ""] + [""] * len(against))
            continue
        ratios = [f"{refs[name] / med:.2f}x" if refs[name] else "" for name in against]
        lines.append([program, str(t), f"{good}/{total}", f"{med:,.0f}", f"{lo:,.0f}–{hi:,.0f}"] + ratios)
    if args.markdown:
        print("| " + " | ".join(header) + " |")
        print("|" + "|".join("---:" if i else "---" for i in range(len(header))) + "|")
        for line in lines:
            print("| " + " | ".join(line) + " |")
        print()
    else:
        widths = [max(len(h), *(len(line[i]) for line in lines)) for i, h in enumerate(header)]
        print("  ".join(h.ljust(w) for h, w in zip(header, widths)))
        for line in lines:
            print("  ".join(v.ljust(w) for v, w in zip(line, widths)))
        print()
