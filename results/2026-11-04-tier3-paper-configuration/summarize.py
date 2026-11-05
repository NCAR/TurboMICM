#!/usr/bin/env python3
"""Summarize a timing session: per-scenario medians over rounds for each program,
sums over the disclosed, held-out, and all scenarios, and the ratios the paper
quotes. Reads OUT/<set>/<scenario>/timings.csv as admin/time_entries.py writes them.

    python3 summarize.py OUT [--programs a,b,...] > summary.md
"""
import csv, statistics, sys
from pathlib import Path

ORDER = [("disclosed", s) for s in ("tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min", "polluted-10min", "clean-10min")] + \
        [("heldout", s) for s in ("summer-15min", "tropics-20min", "winter60-10min", "polluted-48h-30min", "clean-48h-5min")]
NAMES = {"tier1-30s": "Benchmark 30 s", "diurnal-30min": "Diurnal 30 min", "diurnal-10min": "Diurnal 10 min", "diurnal-5min": "Diurnal 5 min",
         "polluted-10min": "Polluted 10 min", "clean-10min": "Clean 10 min", "summer-15min": "Summer 15 min", "tropics-20min": "Tropics 20 min",
         "winter60-10min": "Winter 60N 10 min", "polluted-48h-30min": "Polluted 48 h 30 min", "clean-48h-5min": "Clean 48 h 5 min"}

out = Path(sys.argv[1])
rows = {}
for set_, s in ORDER:
    f = out / set_ / s / "timings.csv"
    if not f.exists():
        continue
    for r in csv.DictReader(open(f)):
        rows.setdefault((set_, s), {}).setdefault(r["program"], []).append(r)
programs = []
for k in rows:
    for p in rows[k]:
        if p not in programs:
            programs.append(p)
med = {}      # (set, scenario, program) -> median elapsed s
digits = {}
gate = {}
for k, byp in rows.items():
    for p, rs in byp.items():
        med[k + (p,)] = statistics.median(float(r["elapsed_ms"]) for r in rs) / 1000.0
        d = [float(r["digits"]) for r in rs if r.get("digits") not in ("", None)]
        digits[k + (p,)] = min(d) if d else float("nan")
        gate[k + (p,)] = all(r["gate"] == "PASS" for r in rs)
print("Median elapsed time per scenario (s), 10,000 cells, 10 threads; digits on the first 64 cells; every run passing the gate unless marked.\n")
print("| Scenario | " + " | ".join(programs) + " |")
print("|---|" + "---:|" * len(programs))
for set_, s in ORDER:
    if (set_, s) not in rows:
        continue
    cells = []
    for p in programs:
        k = (set_, s, p)
        if k in med:
            cells.append(f"{med[k]:.2f} ({digits[k]:.2f}){'' if gate[k] else ' GATE'}")
        else:
            cells.append("")
    print(f"| {NAMES[s]} | " + " | ".join(cells) + " |")
print()
def total(pred, p):
    ks = [(set_, s) for set_, s in ORDER if (set_, s) in rows and pred(set_)]
    if not ks or any((set_, s, p) not in med for set_, s in ks):
        return None
    return sum(med[(set_, s, p)] for set_, s in ks)
print("| Program | All scenarios (s) | Disclosed (s) | Held out (s) |")
print("|---|---:|---:|---:|")
sums = {}
for p in programs:
    a, d, h = total(lambda x: True, p), total(lambda x: x == "disclosed", p), total(lambda x: x == "heldout", p)
    sums[p] = a
    fmt = lambda x: "incomplete" if x is None else f"{x:.2f}"
    print(f"| {p} | {fmt(a)} | {fmt(d)} | {fmt(h)} |")
print()
def ratio(a, b):
    return sums[a] / sums[b] if sums.get(a) and sums.get(b) else float("nan")
pairs = [("loop/C2", "retaining/C2"), ("loop/C3", "retaining/C3"), ("micm/C2", "retaining/C2"), ("micm/C3", "retaining/C3"), ("micm/C2", "loop/C2"), ("micm/C3", "loop/C3")]
print("| Ratio of summed times | All | per-scenario range |")
print("|---|---:|---:|")
for a, b in pairs:
    if a in programs and b in programs:
        per = [med[(set_, s, a)] / med[(set_, s, b)] for set_, s in ORDER if (set_, s, a) in med and (set_, s, b) in med]
        print(f"| {a} / {b} | {ratio(a, b):.2f} | {min(per):.2f} to {max(per):.2f} |")
