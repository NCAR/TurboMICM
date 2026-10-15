"""Run the full TS1 gate and varied groups/schedules against the permitted MICM reference.

Build the desired C compiler and `make micm-ref` first. Python only orchestrates
executables; all timed chemistry remains in C. Results are kept under build/.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--compiler", required=True, choices=("clang", "gcc"))
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
base = root / "build/validation" / ("final-" + args.compiler)
base.mkdir(parents=True, exist_ok=True)
runs = []
for cells, seed, settings in (
        (10000, 2026, [(32, 10, "static")]),
        (127, 17, [(1, 1, "static"), (7, 4, "dynamic,1"), (32, 10, "guided,2"), (1000, 3, "static,2")]),
        (129, 73, [(3, 10, "static,1"), (64, 4, "dynamic,3"), (129, 1, "static")])):
    source = base / f"cells-{cells}-{seed}.csv"
    subprocess.run([sys.executable, root.parent / "reference/make_cells_csv.py", str(cells), str(seed),
                    "--output", source], check=True, capture_output=True, text=True)
    reference = base / f"reference-{cells}-{seed}"
    reference.mkdir(exist_ok=True)
    print(f"MICM: {cells} cells, seed {seed}", flush=True)
    run = subprocess.run([root / "bin/micm_ts1_ref", "--threads", "10" if cells == 10000 else "1",
                          str(cells), "100", "30", source], cwd=reference, capture_output=True, text=True, check=True)
    (reference / "stdout.txt").write_text(run.stdout)
    for group, threads, schedule in settings:
        directory = base / f"c-{cells}-{seed}-g{group}-t{threads}"
        directory.mkdir(exist_ok=True)
        env = dict(os.environ, CHEM_GROUP_SIZE=str(group), OMP_SCHEDULE=schedule)
        run = subprocess.run([root / "bin/ts1_bench_cpu", "--threads", str(threads), str(cells), "100", "30", source],
                             cwd=directory, env=env, capture_output=True, text=True, check=True)
        (directory / "stdout.txt").write_text(run.stdout)
        comparison = subprocess.run([sys.executable, root.parent / "reference/compare_final.py",
                                     directory / "ts1_final.csv", reference / "ts1_final.csv"],
                                    capture_output=True, text=True)
        (directory / "comparison.txt").write_text(comparison.stdout + comparison.stderr)
        result = dict(cells=cells, seed=seed, group=group, threads=threads, schedule=schedule,
                      stdout=run.stdout, comparison=comparison.stdout + comparison.stderr, passed=comparison.returncode == 0)
        runs.append(result)
        print(json.dumps(result), flush=True)
        comparison.check_returncode()
(base / "results.json").write_text(json.dumps(runs, indent=2) + "\n")
