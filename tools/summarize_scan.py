#!/usr/bin/env python3
"""Summarize a scan CSV (tools/scan.py's format) point by point.

    tools/summarize_scan.py scans/S2-tableau-tolerance.csv [--sort time|digits|worst]
                            [--level 2|3] [--scenario NAME] [--markdown]

For each point (label): the number of scenarios run, whether every run exited
normally and passed the gate, the lowest overall digits and the lowest
worst-species digits over the scenarios (the ranking's "every scenario"
view), the summed median elapsed time, the summed counts per scenario step,
and the levels reached (L digits overall and L - 1 per species everywhere).
With --scenario, the rows of that scenario only, per round.
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, nargs="+")
    parser.add_argument("--sort", choices=["time", "digits", "worst", "label"], default="label")
    parser.add_argument("--scenario")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    rows = []
    for p in args.csv:
        with open(p, newline="") as f:
            rows += list(csv.DictReader(f))
    if args.scenario:
        rows = [r for r in rows if r["scenario"] == args.scenario]
    points = defaultdict(lambda: defaultdict(list))
    for r in rows:
        points[r["label"]][r["scenario"]].append(r)
    table = []
    for label, scen in points.items():
        ok = all(r["status"] == "ok" and r["gate"] == "PASS" for rs in scen.values() for r in rs)
        digits = [float(r["digits"]) for rs in scen.values() for r in rs if r["digits"]]
        worst = [float(r["worst_digits"]) for rs in scen.values() for r in rs if r["worst_digits"]]
        times = [statistics.median(float(r["elapsed_ms"]) for r in rs if r["elapsed_ms"]) / 1000 for rs in scen.values() if any(r["elapsed_ms"] for r in rs)]
        def per_step(key):
            vals = []
            for rs in scen.values():
                r = rs[0]
                if r.get(key):
                    steps = {"tier1-30s": 101, "diurnal-30min": 48, "diurnal-10min": 144, "diurnal-5min": 288, "polluted-10min": 144, "clean-10min": 144}.get(r["scenario"])
                    groups = (int(r["cells"]) + 31) // 32
                    if steps:
                        vals.append(float(r[key]) / (steps * groups))
            return statistics.mean(vals) if vals else float("nan")
        table.append({
            "label": label, "method": next(iter(scen.values()))[0]["method"], "n": len(scen), "ok": ok,
            "digits": min(digits) if digits else float("nan"), "worst": min(worst) if worst else float("nan"),
            "time": sum(times) if times else float("nan"),
            "steps": per_step("internal_steps"), "rej": per_step("rejected"), "lu": per_step("decompositions"),
            "solves": per_step("solves"), "f": per_step("function_calls"), "reuse": per_step("matrix_reuses"),
            "L2": ok and digits and min(digits) >= 2 and min(worst) >= 1, "L3": ok and digits and min(digits) >= 3 and min(worst) >= 2,
        })
    key = {"time": lambda t: t["time"], "digits": lambda t: -t["digits"], "worst": lambda t: -t["worst"], "label": lambda t: t["label"]}[args.sort]
    table.sort(key=key)
    sep = " | " if args.markdown else "  "
    head = ["label", "method", "scen", "ok", "digits", "worst", "time s", "steps/st", "rej/st", "LU/st", "solves/st", "F/st", "reuse/st", "L2", "L3"]
    if args.markdown:
        print("| " + " | ".join(head) + " |")
        print("|" + "|".join("---:" if h not in ("label", "method", "ok", "L2", "L3") else "---" for h in head) + "|")
    else:
        print(f"{head[0]:<30}{head[1]:<10}{head[2]:>5}{head[3]:>4}{head[4]:>8}{head[5]:>7}{head[6]:>8}{head[7]:>9}{head[8]:>7}{head[9]:>7}{head[10]:>10}{head[11]:>7}{head[12]:>9}{head[13]:>4}{head[14]:>4}")
    for t in table:
        cells = [t["label"], t["method"], str(t["n"]), "y" if t["ok"] else "N", f"{t['digits']:.2f}", f"{t['worst']:.2f}", f"{t['time']:.2f}",
                 f"{t['steps']:.2f}", f"{t['rej']:.2f}", f"{t['lu']:.2f}", f"{t['solves']:.2f}", f"{t['f']:.2f}", f"{t['reuse']:.2f}",
                 "y" if t["L2"] else "-", "y" if t["L3"] else "-"]
        if args.markdown:
            print("| " + " | ".join(cells) + " |")
        else:
            print(f"{cells[0]:<30}{cells[1]:<10}{cells[2]:>5}{cells[3]:>4}{cells[4]:>8}{cells[5]:>7}{cells[6]:>8}{cells[7]:>9}{cells[8]:>7}{cells[9]:>7}{cells[10]:>10}{cells[11]:>7}{cells[12]:>9}{cells[13]:>4}{cells[14]:>4}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
