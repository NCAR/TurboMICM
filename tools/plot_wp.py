#!/usr/bin/env python3
"""Work-precision plots from scan CSVs (tools/scan.py's format).

    tools/plot_wp.py scans/S2-tableau-tolerance.csv --out results/wp-S2.png
                     [--x time|solves|steps] [--y digits|worst] [--aggregate sum|min]

One curve per method (the label up to the first '-'), one point per
configuration: x is the median elapsed time summed over the scenarios (or the
summed operation count), y is the lowest overall digits (or the lowest worst-
species digits) over the scenarios, so a point's y is what the ranking's
"every scenario" rule sees. Lines at 2 and 3 digits (1 and 2 for the worst
species) mark the levels.
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--x", choices=["time", "solves", "steps", "decompositions"], default="time")
    parser.add_argument("--y", choices=["digits", "worst"], default="digits")
    parser.add_argument("--title", default="")
    args = parser.parse_args()
    rows = []
    for path in args.csv:
        with open(path, newline="") as f:
            rows += [r for r in csv.DictReader(f) if r["status"] == "ok" and r["gate"] == "PASS" and r["digits"]]
    by_point = defaultdict(lambda: defaultdict(list))  # label -> scenario -> rows
    for r in rows:
        by_point[r["label"]][r["scenario"]].append(r)
    n_scenarios = len({r["scenario"] for r in rows})
    series = defaultdict(list)
    for label, scen in by_point.items():
        if len(scen) < n_scenarios:
            continue
        xs, ys = 0.0, []
        for name, rs in scen.items():
            if args.x == "time":
                xs += statistics.median(float(r["elapsed_ms"]) for r in rs) / 1000.0
            else:
                key = {"solves": "solves", "steps": "internal_steps", "decompositions": "decompositions"}[args.x]
                xs += statistics.median(float(r[key]) for r in rs)
            ys.append(min(float(r["digits"] if args.y == "digits" else r["worst_digits"]) for r in rs))
        series[label.split("-")[0]].append((xs, min(ys), label))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for method, pts in sorted(series.items()):
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", label=method)
    for level in ((2, 3) if args.y == "digits" else (1, 2)):
        ax.axhline(level, color="gray", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel({"time": "time, s (sum over scenarios of median elapsed)", "solves": "linear solves (sum)", "steps": "internal steps (sum)", "decompositions": "factorizations (sum)"}[args.x])
    ax.set_ylabel("lowest overall digits over the scenarios" if args.y == "digits" else "lowest worst-species digits")
    ax.set_title(args.title or ", ".join(p.stem for p in args.csv))
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}: {sum(len(v) for v in series.values())} points, {len(series)} methods")
    return 0


if __name__ == "__main__":
    sys.exit(main())
