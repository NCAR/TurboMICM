"""Make the small scenarios the agreement test and quick experiments use.

    python prototypes/make_mini.py [--source ../scenarios/diurnal-10min] [--cells 4] [--steps 12]
                                   [--start-step 36] [--out prototypes/mini] [--reference bin/micm_scenario]

Copies the first CELLS cells of a disclosed scenario and its photolysis
factors, shifted so that the mini scenario starts at step START_STEP of the
source (06:00 plus START_STEP steps, daylight by default), with STEPS steps and
a snapshot every STEPS / 2 steps, then makes a reference with MICM's six-stage
Rosenbrock at relative 1e-10 and absolute 1e-23 with its step limit raised
(make micm-scenario), as the disclosed references were made.
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, default=ENTRY.parent / "scenarios" / "diurnal-10min")
    parser.add_argument("--cells", type=int, default=4)
    parser.add_argument("--steps", type=int, default=12)
    parser.add_argument("--start-step", type=int, default=36)
    parser.add_argument("--out", type=Path, default=HERE / "mini")
    parser.add_argument("--reference", type=Path, default=ENTRY / "bin" / "micm_scenario")
    parser.add_argument("--no-reference", action="store_true")
    args = parser.parse_args()
    values = {}
    for line in (args.source / "scenario.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            k, _, v = line.partition("=")
            values[k.strip()] = v.strip()
    dt = float(values["dt"])
    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.source / values["cells"], newline="") as src, open(args.out / "cells.csv", "w", newline="") as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst)
        writer.writerow(next(reader))
        for i, row in enumerate(reader):
            if i >= args.cells:
                break
            writer.writerow(row)
    shift = args.start_step * dt
    with open(args.source / values["photolysis"], newline="") as src, open(args.out / "photolysis.csv", "w", newline="") as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst)
        writer.writerow(next(reader))
        for t, f in reader:
            t = float(t) - shift
            if t >= 0.0:
                writer.writerow([f"{t:.17g}", f])
    every = max(1, args.steps // 2)
    (args.out / "scenario.txt").write_text(
        f"# {args.cells} cells of {args.source.name} from its step {args.start_step}; {args.steps} steps of {dt:g} s\n"
        f"cells = cells.csv\ndt = {dt:g}\nsteps = {args.steps}\nphotolysis = photolysis.csv\nsnapshot_every = {every}\n")
    if args.no_reference:
        return 0
    if not args.reference.exists():
        sys.exit(f"{args.reference} is missing: make micm-scenario")
    for tolerances, name in ((("1e-10", "1e-23"), "reference_snapshots.csv"), (("1e-9", "1e-22"), "check_snapshots.csv")):
        run = subprocess.run([str(args.reference), "--threads", "1", "--method", "rosenbrock6dae", "--rtol", tolerances[0], "--atol", tolerances[1], "--max-steps", "10000000", "scenario.txt"],
                             cwd=args.out, capture_output=True, text=True)
        if run.returncode != 0:
            sys.exit(f"reference run failed: {run.stderr}")
        (args.out / "ts1_snapshots.csv").rename(args.out / name)
        (args.out / "ts1_final.csv").unlink(missing_ok=True)
    print(f"wrote {args.out}: {args.cells} cells, {args.steps} steps of {dt:g} s, reference and check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
