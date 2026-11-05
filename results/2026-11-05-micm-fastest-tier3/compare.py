#!/usr/bin/env python3
"""Compare the 5 November session (MICM's fastest layout, vector ordering of 8,
with the retaining solver as anchor) with the 4 November paper session.

    python3 compare.py NEW_OUT OLD_OUT > compare.md

Both OUT directories hold <set>/<scenario>/timings.csv as admin/time_entries.py
writes them. Medians over rounds per scenario; sums over the eleven scenarios.
"""
import csv
import statistics
import sys
from pathlib import Path

ORDER = [("disclosed", s) for s in ("tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min",
                                     "polluted-10min", "clean-10min")] + \
        [("heldout", s) for s in ("summer-15min", "tropics-20min", "winter60-10min",
                                   "polluted-48h-30min", "clean-48h-5min")]


def load(out):
    med, digits, gate, rounds = {}, {}, {}, {}
    for set_, s in ORDER:
        f = Path(out) / set_ / s / "timings.csv"
        by = {}
        for r in csv.DictReader(open(f)):
            by.setdefault(r["program"], []).append(r)
        for p, rs in by.items():
            med[(s, p)] = statistics.median(float(r["elapsed_ms"]) for r in rs) / 1000.0
            d = [float(r["digits"]) for r in rs if r.get("digits") not in ("", None)]
            digits[(s, p)] = min(d) if d else float("nan")
            gate[(s, p)] = all(r["gate"] == "PASS" for r in rs)
            rounds[(s, p)] = len(rs)
    return med, digits, gate, rounds


new, new_digits, new_gate, new_rounds = load(sys.argv[1])
old, old_digits, _, _ = load(sys.argv[2])
names = [s for _, s in ORDER]


def total(med, p):
    return sum(med[(s, p)] for s in names)


print("# MICM's fastest layout on the eleven Tier 3 scenarios (5 November 2026)\n")
print("Per-scenario medians (s) over the rounds; digits on the first 64 cells; gate on all cells.\n")
print("| Scenario | micm8/C2 | micm8/C3 | retaining/C2 | retaining/C3 | micm8/ret C2 | micm8/ret C3 |")
print("|---|---:|---:|---:|---:|---:|---:|")
for s in names:
    row = [f"{new[(s, p)]:.2f} ({new_digits[(s, p)]:.2f}){'' if new_gate[(s, p)] else ' GATE'}"
           for p in ("micm8/C2", "micm8/C3", "retaining/C2", "retaining/C3")]
    row += [f"{new[(s, 'micm8/C2')] / new[(s, 'retaining/C2')]:.2f}",
            f"{new[(s, 'micm8/C3')] / new[(s, 'retaining/C3')]:.2f}"]
    print(f"| {s} | " + " | ".join(row) + " |")
print()
print(f"Rounds per program and scenario: {sorted(set(new_rounds.values()))}; "
      f"every run passing the gate: {all(new_gate.values())}\n")

print("| Sum over the eleven scenarios (s) | 5 Nov | 4 Nov | 5 Nov / 4 Nov |")
print("|---|---:|---:|---:|")
for p in ("retaining/C2", "retaining/C3"):
    a, b = total(new, p), total(old, p)
    print(f"| {p} (anchor) | {a:.2f} | {b:.2f} | {a / b:.3f} |")
for p in ("micm8/C2", "micm8/C3"):
    print(f"| {p} | {total(new, p):.2f} | | |")
for p in ("micm/C2", "micm/C3", "loop/C2", "loop/C3"):
    print(f"| {p} | | {total(old, p):.2f} | |")
print()

print("| Ratio | All | Per-scenario range |")
print("|---|---:|---:|")
for c in ("C2", "C3"):
    per = [new[(s, f"micm8/{c}")] / new[(s, f"retaining/{c}")] for s in names]
    print(f"| micm8/{c} / retaining/{c}, same session | "
          f"{total(new, f'micm8/{c}') / total(new, f'retaining/{c}'):.2f} | {min(per):.2f} to {max(per):.2f} |")
for c in ("C2", "C3"):
    # Each session's MICM time relative to its own retaining time, so session drift cancels.
    per = [(old[(s, f"micm/{c}")] / old[(s, f"retaining/{c}")]) / (new[(s, f"micm8/{c}")] / new[(s, f"retaining/{c}")])
           for s in names]
    allr = (total(old, f"micm/{c}") / total(old, f"retaining/{c}")) / \
           (total(new, f"micm8/{c}") / total(new, f"retaining/{c}"))
    print(f"| micm/{c} (default layout) / micm8/{c}, anchored | {allr:.2f} | {min(per):.2f} to {max(per):.2f} |")
for c in ("C2", "C3"):
    per = [(new[(s, f"micm8/{c}")] / new[(s, f"retaining/{c}")]) / (old[(s, f"loop/{c}")] / old[(s, f"retaining/{c}")])
           for s in names]
    allr = (total(new, f"micm8/{c}") / total(new, f"retaining/{c}")) / \
           (total(old, f"loop/{c}") / total(old, f"retaining/{c}"))
    print(f"| micm8/{c} / loop/{c} (library kernels), anchored | {allr:.2f} | {min(per):.2f} to {max(per):.2f} |")
