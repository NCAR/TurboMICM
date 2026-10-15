#!/usr/bin/env python3
"""Tabulate admin/micm_configs.sh's results: OUT/configs.md.

    admin/micm_configs_table.py OUT DEFAULT

One row per configuration, in OUT/check.csv's order: its median elapsed_ms and
range with one cell group per thread and with 32-cell groups, each relative to
DEFAULT with the same grouping; its internal steps; the gate results of its runs
against MICM_FINAL; and the largest relative difference from MICM_FINAL. Needs
only the Python standard library.
"""

import csv
import statistics
import sys
from pathlib import Path

out, default = Path(sys.argv[1]), sys.argv[2]
with open(out / "timing" / "timings.csv", newline="") as file:
    rows = list(csv.DictReader(file))
with open(out / "check.csv", newline="") as file:
    configs = [row["config"] for row in csv.DictReader(file)]
GROUPINGS = ("", "-g32")


def runs(program):
    return [row for row in rows if row["program"] == program]


def median_ms(program):
    times = [float(row["elapsed_ms"]) for row in runs(program) if row["elapsed_ms"]]
    return (statistics.median(times), min(times), max(times)) if times else None


base = {grouping: median_ms(f"micm-{default}{grouping}") for grouping in GROUPINGS}
lines = [
    "| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default"
    " | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |",
    "|---|---:|---:|---:|---:|---:|---|---:|",
]
for config in configs:
    cells, steps, gates, largest = [], [], set(), []
    for grouping in GROUPINGS:
        program = f"micm-{config}{grouping}"
        timing = median_ms(program)
        if timing is None or base[grouping] is None:
            cells += ["failed", ""]
        else:
            median, low, high = timing
            cells += [f"{median:,.0f} ({low:,.0f}–{high:,.0f})", f"{base[grouping][0] / median:.2f}×"]
        steps.append("/".join(sorted({row["internal_steps"] for row in runs(program)})) or "-")
        gates.update(row["gate"] for row in runs(program))
        largest += [float(row["max_relative_difference"]) for row in runs(program) if row["max_relative_difference"]]
    largest_text = f"{max(largest):.2g}" if largest else ""
    lines.append(f"| {config} | {' | '.join(cells)} | {'; '.join(steps)} | {', '.join(sorted(gates))} | {largest_text} |")
(out / "configs.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
