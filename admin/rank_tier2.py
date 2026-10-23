#!/usr/bin/env python3
"""Rank the Tier 2 entries from a timing session.

    admin/rank_tier2.py TIMINGS.csv [--threads 18] [--out RANKING.md]

TIMINGS.csv is admin/time_entries.py's output for a session that timed, at the
given thread count, MICM's programs (micm, micm-g32, micm-vector8-g32) and, for
each entry RUN, RUN/specialized, RUN/generic (its Tier 2 commit's
bin/ts1_bench_cpu), and RUN/tier1 (its Tier 1 commit's, built on this machine).
Every run must pass the gate; a program with a failed run is flagged. Entries
are ranked by the median elapsed time of their specialized benchmark, the
primary ranking; the table gives, beside it, the spread of the five runs, the
entry's generic and Tier 1 medians, the gain from specialization (generic over
specialized), and the speedup over each of MICM's three programs. Needs only the
Python standard library.
"""

import argparse
import csv
import statistics
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("timings", type=Path)
parser.add_argument("--threads", type=int, default=18)
parser.add_argument("--out", type=Path)
args = parser.parse_args()

rows = [r for r in csv.DictReader(open(args.timings)) if int(r["threads"]) == args.threads]
programs = {}
for row in rows:
    programs.setdefault(row["program"], []).append(row)


def median(name):
    runs = programs.get(name, [])
    times = [float(r["elapsed_ms"]) for r in runs if r["gate"] == "PASS" and r["elapsed_ms"]]
    if not times or len(times) != len(runs):
        return None, None, len(times), len(runs)
    return statistics.median(times), (max(times) - min(times)) / statistics.median(times), len(times), len(runs)


micm = {name: median(name)[0] for name in ("micm", "micm-g32", "micm-vector8-g32")}
entries = sorted({name.split("/")[0] for name in programs if "/" in name})
table = []
for entry in entries:
    special, spread, passed, total = median(f"{entry}/specialized")
    generic = median(f"{entry}/generic")[0]
    tier1 = median(f"{entry}/tier1")[0]
    table.append((special, spread, entry, generic, tier1, passed, total))
table.sort(key=lambda row: (row[0] is None, row[0] or 0))

lines = [f"Ranking at {args.threads} thread{'s' if args.threads != 1 else ''} (median of the runs, ms; "
         f"MICM: default {micm['micm']:,.1f}, 32-cell groups {micm['micm-g32']:,.1f}, "
         f"vector8 with 32-cell groups {micm['micm-vector8-g32']:,.1f})\n",
         "| Rank | Entry | Specialized (ms) | Spread | Generic, Tier 2 (ms) | Tier 1 build (ms) | Gain from specialization "
         "| Over MICM's default | Over MICM, 32-cell groups | Over MICM's fastest |",
         "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
for rank, (special, spread, entry, generic, tier1, passed, total) in enumerate(table, 1):
    if special is None:
        lines.append(f"| — | {entry} | {passed} of {total} runs passed | | | | | | | |")
        continue
    def ms(value):
        return f"{value:,.1f}" if value is not None else "—"
    gain = f"{generic / special:.2f}" if generic is not None else "—"
    lines.append(f"| {rank} | {entry} | {special:,.1f} | {spread * 100:.1f} % | {ms(generic)} | {ms(tier1)} | "
                 f"{gain} | {micm['micm'] / special:.2f} | {micm['micm-g32'] / special:.2f} | "
                 f"{micm['micm-vector8-g32'] / special:.2f} |")
failed = [name for name, runs in programs.items() if any(r["gate"] != "PASS" for r in runs)]
lines.append("")
lines.append(f"Runs: {len(rows)}; programs with a run that failed the gate: {', '.join(failed) or 'none'}.")
text = "\n".join(lines) + "\n"
if args.out:
    args.out.write_text(text)
print(text)
