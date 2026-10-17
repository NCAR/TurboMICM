#!/usr/bin/env python3
"""Run bin/ts1_bench_cpu in several settings, apply the experiment's gate to
each result, and write a table.

    check_benchmark.py --out DIR [--cells 10000 --seed 2026 --steps 100 --dt 30]
                       [--bench bin/ts1_bench_cpu] [--micm bin/micm_ts1_ref | --reference FINAL.csv]
                       [--variant NAME="ARGS" ...] [--repeat N]

Without --reference, MICM's final state comes from running --micm on the same
input. The default variants are the ranking configuration, one thread, other
group sizes, the generic kernels, and the other schedules. The gate is
reference/compare_final.py: |y - y_MICM| <= 1e-30 + 1e-10 |y_MICM|.
"""
import argparse
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
REFERENCE = Path("../reference")

DEFAULT_VARIANTS = [
    ("default (all threads, groups of 8)", ""),
    ("1 thread", "--threads 1"),
    ("groups of 4", "--group-size 4"),
    ("groups of 16", "--group-size 16"),
    ("groups of 32", "--group-size 32"),
    ("generic kernels", "--generic"),
    ("static schedule", "--schedule static"),
    ("guided schedule", "--schedule guided"),
    ("2 threads, groups of 3", "--threads 2 --group-size 3"),
]


def run(cmd, cwd):
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(cmd)}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cells", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--dt", type=float, default=30.0)
    ap.add_argument("--bench", type=Path, default=HERE / "bin/ts1_bench_cpu")
    ap.add_argument("--micm", type=Path, default=HERE / "bin/micm_ts1_ref")
    ap.add_argument("--reference", type=Path)
    ap.add_argument("--reference-dir", type=Path, default=HERE.parent / "reference")
    ap.add_argument("--variant", action="append", default=[], help='NAME="ARGS"; replaces the defaults')
    ap.add_argument("--repeat", type=int, default=1, help="runs per variant; the median is reported")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    python = sys.executable
    cells_csv = out / f"cells_{args.cells}_seed{args.seed}.csv"
    if not cells_csv.exists():
        run([python, str(args.reference_dir / "make_cells_csv.py"), str(args.cells), str(args.seed), "--output", str(cells_csv)], out)
    if args.reference:
        reference = args.reference.resolve()
    else:
        micm_dir = out / "micm"
        micm_dir.mkdir(exist_ok=True)
        reference = micm_dir / "ts1_final.csv"
        if not reference.exists():
            print("running MICM's reference ...", flush=True)
            text = run([str(args.micm.resolve()), str(args.cells), str(args.steps), str(args.dt), str(cells_csv)], micm_dir)
            (micm_dir / "stdout.txt").write_text(text)
            print(text.strip())
    variants = DEFAULT_VARIANTS
    if args.variant:
        variants = []
        for v in args.variant:
            name, _, arguments = v.partition("=")
            variants.append((name, arguments.strip('"')))
    rows = []
    for name, arguments in variants:
        slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
        run_dir = out / slug
        run_dir.mkdir(exist_ok=True)
        times, stats, gates = [], None, []
        for r in range(args.repeat):
            cmd = [str(args.bench.resolve())] + arguments.split() + [str(args.cells), str(args.steps), str(args.dt), str(cells_csv)]
            text = run(cmd, run_dir)
            (run_dir / f"stdout_{r}.txt").write_text(text)
            times.append(float(re.search(r"elapsed_ms=([0-9.]+)", text).group(1)))
            stats = re.search(r"internal_steps=.*", text).group(0)
            gate = subprocess.run([python, str(args.reference_dir / "compare_final.py"), str(run_dir / "ts1_final.csv"), str(reference)], capture_output=True, text=True)
            m = re.search(r"largest relative difference ([0-9.e+-]+)", gate.stdout)
            gates.append(("PASS" if gate.returncode == 0 else "FAIL", m.group(1) if m else "?"))
        median = statistics.median(times)
        worst = max(gates, key=lambda g: (g[0] == "FAIL", float(g[1]) if g[1] != "?" else 1.0))
        rows.append({"variant": name, "arguments": arguments, "elapsed_ms": times, "median_ms": median, "statistics": stats, "gate": worst[0], "max_relative_difference": worst[1]})
        print(f"{name:40s} {median:10.1f} ms  {worst[0]} {worst[1]}  {stats}", flush=True)
    (out / "results.json").write_text(json.dumps({"label": args.label, "cells": args.cells, "seed": args.seed, "steps": args.steps, "dt": args.dt, "rows": rows}, indent=1))
    lines = [f"# TS1 benchmark checks {args.label}".rstrip(), "", f"{args.cells} cells (seed {args.seed}), {args.steps} steps of {args.dt:g} s; the gate is rtol 1e-10, atol 1e-30 against MICM.", "",
             "| Variant | Arguments | Median ms | Runs | Gate | Largest relative difference | Statistics |", "|---|---|---:|---:|---|---:|---|"]
    for row in rows:
        lines.append(f"| {row['variant']} | `{row['arguments']}` | {row['median_ms']:,.0f} | {len(row['elapsed_ms'])} | {row['gate']} | {row['max_relative_difference']} | {row['statistics']} |")
    (out / "results.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out / 'results.md'}")


if __name__ == "__main__":
    main()
