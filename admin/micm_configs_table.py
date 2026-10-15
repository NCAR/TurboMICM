#!/usr/bin/env python3
"""Tabulate admin/micm_configs.sh's results: OUT/configs.md.

    admin/micm_configs_table.py OUT DEFAULT [LABEL]

LABEL (default all) picks the timing session, OUT/timing-LABEL; the table goes
to OUT/configs-LABEL.md. A configuration that session did not time is left blank.

First come matrix tables, matrix layout by LU variant: the median elapsed_ms of
each configuration and its speed relative to DEFAULT, one table per grouping,
then the largest relative difference of its final states from MICM_FINAL. Then
one row per configuration, in OUT/check.csv's order: its median elapsed_ms and
range with one cell group per thread and with 32-cell groups, each relative to
DEFAULT with the same grouping; its internal steps; the gate results of its runs
against MICM_FINAL; and the largest relative difference from MICM_FINAL. Any other
programs timed in the same rounds follow, relative to both defaults. Needs only the
Python standard library.
"""

import csv
import statistics
import sys
from pathlib import Path

out, default = Path(sys.argv[1]), sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else "all"
with open(out / f"timing-{label}" / "timings.csv", newline="") as file:
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
LAYOUTS = ["standard", "vector1", "vector2", "vector4", "vector8", "vector128"]
LUS = [("in-place", "mozart"), ("in-place", "doolittle"), ("separate", "mozart"), ("separate", "doolittle")]
lu_header = "| Matrix layout | " + " | ".join(
    f"{storage} {algorithm.capitalize()}" + (" (default)" if f"standard-{storage}-{algorithm}" == default else "")
    for storage, algorithm in LUS) + " |"
lines = []
for grouping, title in (("", "one cell group per thread"), ("-g32", "32-cell groups")):
    lines += [f"Median elapsed_ms (speed relative to the default), {title}:", "", lu_header, "|---|" + "---:|" * len(LUS)]
    for layout in LAYOUTS:
        cells = []
        for storage, algorithm in LUS:
            timing = median_ms(f"micm-{layout}-{storage}-{algorithm}{grouping}")
            program = f"micm-{layout}-{storage}-{algorithm}{grouping}"
            if not runs(program):
                cells.append("")
            elif timing and base[grouping]:
                cells.append(f"{timing[0]:,.0f} ({base[grouping][0] / timing[0]:.2f}×)")
            else:
                cells.append("failed")
        lines.append(f"| {layout} | " + " | ".join(cells) + " |")
    lines.append("")
lines += ["Largest relative difference of the final states from MICM_FINAL, over both groupings:", "", lu_header,
          "|---|" + "---:|" * len(LUS)]
for layout in LAYOUTS:
    cells = []
    for storage, algorithm in LUS:
        values = [float(row["max_relative_difference"]) for grouping in GROUPINGS
                  for row in runs(f"micm-{layout}-{storage}-{algorithm}{grouping}") if row["max_relative_difference"]]
        cells.append(f"{max(values):.1g}" if values else "")
    lines.append(f"| {layout} | " + " | ".join(cells) + " |")
lines.append("")
lines += [
    "| Configuration | One group per thread: median ms (range) | vs default | 32-cell groups: median ms (range) | vs default"
    " | Internal steps (group per thread; 32-cell) | Gate | Largest relative difference |",
    "|---|---:|---:|---:|---:|---:|---|---:|",
]
for config in configs:
    if not any(runs(f"micm-{config}{grouping}") for grouping in GROUPINGS):
        continue
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
known = {f"micm-{config}{grouping}" for config in configs for grouping in GROUPINGS}
others = [name for name in dict.fromkeys(row["program"] for row in rows) if name not in known]
if others:
    lines += [
        "",
        "| Program | Median ms (range) | vs default, one group per thread | vs default, 32-cell groups | Gate | Largest relative difference |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for name in others:
        timing = median_ms(name)
        gates = ", ".join(sorted({row["gate"] for row in runs(name)}))
        largest = [float(row["max_relative_difference"]) for row in runs(name) if row["max_relative_difference"]]
        largest_text = f"{max(largest):.2g}" if largest else ""
        if timing is None:
            lines.append(f"| {name} | failed | | | {gates} | {largest_text} |")
            continue
        median, low, high = timing
        ratios = [f"{base[g][0] / median:.2f}×" if base[g] else "" for g in GROUPINGS]
        lines.append(f"| {name} | {median:,.0f} ({low:,.0f}–{high:,.0f}) | {ratios[0]} | {ratios[1]} | {gates} | {largest_text} |")
(out / f"configs-{label}.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
