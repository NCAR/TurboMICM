#!/usr/bin/env python3
"""Work-precision plots from scan CSVs (tools/scan.py's format).

    tools/plot_wp.py scans/S2-tableau-tolerance.csv --out plots/wp-S2.png
                     [--x time|solves|steps|decompositions] [--y digits|worst]
                     [--series REGEX] [--mark LABEL ...] [--title TEXT]

One curve per series, by default the method (the label up to the first
'-'); --series gives a regular expression whose capture groups, joined by
spaces, name the series (for example '^(\\w+)-.*-(micm|h211b)$' separates
the controllers). --mark draws the named configurations, the operating
points, with large markers and their labels. One point per configuration: x is the median elapsed time summed over the scenarios (or the
summed operation count), y is the lowest overall digits (or the lowest worst-
species digits) over the scenarios, so a point's y is what the ranking's
"every scenario" rule sees. Lines at 2 and 3 digits (1 and 2 for the worst
species) mark the levels.
"""
import argparse
import csv
import re
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
    parser.add_argument("--series", default=r"^([^-]+)", help="regular expression on the label; its groups name the series")
    parser.add_argument("--mark", nargs="*", default=[], help="labels to mark as operating points")
    args = parser.parse_args()
    pattern = re.compile(args.series)
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
        m = pattern.search(label)
        name = " ".join(g for g in m.groups() if g) if m else label
        series[name].append((xs, min(ys), label))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for method, pts in sorted(series.items()):
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", markersize=4, linewidth=1, label=method)
    marked = [(x, y, label) for pts in series.values() for (x, y, label) in pts if label in args.mark]
    for x, y, label in marked:
        ax.plot([x], [y], marker="*", markersize=14, color="black", zorder=5)
        ax.annotate(label, (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)
    missing = set(args.mark) - {m[2] for m in marked}
    if missing:
        print(f"not plotted (not in every scenario or not in the CSVs): {sorted(missing)}", file=sys.stderr)
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
