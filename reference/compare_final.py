#!/usr/bin/env python3
"""Check a TS1 final state against MICM's, cell by cell and species by species.

    compare_final.py YOURS MICM [--rtol 1e-10] [--atol 1e-30]

Both files are ts1_final.csv files: a header row of temperature, pressure, and
species names, then one row per cell, in the same cell order. Columns are
matched by name. A value passes when |y - y_MICM| <= atol + rtol |y_MICM|
(mol m-3); the defaults are the experiment's gate. Exits 0 when every value
passes and 1 otherwise. Needs only the Python standard library.
"""

import argparse
import csv
import math
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("yours")
parser.add_argument("micm")
parser.add_argument("--rtol", type=float, default=1e-10)
parser.add_argument("--atol", type=float, default=1e-30)
args = parser.parse_args()


def read(path):
    with open(path, newline="") as file:
        rows = list(csv.reader(file))
    header, data = rows[0], [row for row in rows[1:] if row]
    return {name: [float(row[i]) for row in data] for i, name in enumerate(header)}, len(data)


yours, n_yours = read(args.yours)
micm, n_micm = read(args.micm)
problems = []
if n_yours != n_micm:
    problems.append(f"{n_yours} cells, MICM has {n_micm}")
missing = sorted(set(micm) - set(yours))
extra = sorted(set(yours) - set(micm))
if missing:
    problems.append(f"missing columns: {' '.join(missing)}")
if extra:
    problems.append(f"unexpected columns: {' '.join(extra)}")
if problems:
    print("FAIL: " + "; ".join(problems))
    sys.exit(1)

failures, worst = 0, (0.0, None, None)
for name, reference in micm.items():
    for cell, (y, y_micm) in enumerate(zip(yours[name], reference)):
        difference = abs(y - y_micm)
        if not math.isfinite(y) or difference > args.atol + args.rtol * abs(y_micm):
            failures += 1
        if y_micm != 0.0:
            relative = difference / abs(y_micm) if math.isfinite(y) else math.inf
            if relative > worst[0]:
                worst = (relative, name, cell)

checked = n_micm * len(micm)
summary = f"{checked} values in {n_micm} cells x {len(micm)} columns; largest relative difference {worst[0]:.3g}"
if worst[1] is not None:
    summary += f" ({worst[1]}, cell {worst[2]})"
if failures:
    print(f"FAIL: {failures} of {summary} (rtol {args.rtol:g}, atol {args.atol:g})")
    sys.exit(1)
print(f"PASS: {summary} (rtol {args.rtol:g}, atol {args.atol:g})")
