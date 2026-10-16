#!/usr/bin/env python3
"""Score a TS1 result against a reference, in significant digits of accuracy.

    score.py RESULT REFERENCE [--threshold 1e-14] [--level L] [--json]

RESULT and REFERENCE are both ts1_snapshots.csv files (step, cell, temperature,
pressure, species...) or both ts1_final.csv files (temperature, pressure,
species...). Rows are matched by (step, cell), or by order for final states;
species are matched by name.

Following Sandu et al. (1997), each value counts when the reference value is at
least the threshold (default 1e-14 mol m-3). Its relative error is
|y - y_ref| / y_ref. Overall accuracy is -log10 of the root mean square of all
counted relative errors; a species' accuracy is the same over that species
alone. Accuracy is capped at 16 digits; an error of exactly zero counts as the cap.

The gate fails on any NaN or infinity, or any value below -1e-20 mol m-3.
--level L passes when the gate passes, overall accuracy is at least L, and every
species reaches at least L - 1. The exit status is 0 on a pass and 1 on a failure
(gate or level). Needs only the Python standard library.
"""

import argparse
import csv
import json
import math
import sys

CAP = 16.0
FLOOR = -1e-20

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("result")
parser.add_argument("reference")
parser.add_argument("--threshold", type=float, default=1e-14)
parser.add_argument("--level", type=float)
parser.add_argument("--json", action="store_true", help="print the results as JSON")
args = parser.parse_args()


def read(path):
    with open(path, newline="") as file:
        reader = csv.reader(file)
        header = next(reader)
        rows = [row for row in reader if row]
    keyed = header[:2] == ["step", "cell"]
    start = 4 if keyed else 2
    table = {}
    for index, row in enumerate(rows):
        key = (int(row[0]), int(row[1])) if keyed else index
        table[key] = [float(value) for value in row[start:]]
    return keyed, header[start:], table


def digits(mean_square):
    if mean_square <= 0.0:
        return CAP
    return min(CAP, -math.log10(math.sqrt(mean_square)))


result_keyed, result_species, result = read(args.result)
reference_keyed, reference_species, reference = read(args.reference)
problems = []
if result_keyed != reference_keyed:
    problems.append("one file has snapshots and the other does not")
missing = sorted(set(reference) - set(result))
if missing:
    problems.append(f"{len(missing)} reference rows missing from the result, e.g. {missing[0]}")
absent = sorted(set(reference_species) - set(result_species))
if absent:
    problems.append(f"species missing from the result: {' '.join(absent)}")
if problems:
    print("FAIL: " + "; ".join(problems))
    sys.exit(1)

column = {name: i for i, name in enumerate(result_species)}
nonfinite = negative = 0
worst_negative = 0.0
for values in result.values():
    for value in values:
        if not math.isfinite(value):
            nonfinite += 1
        elif value < FLOOR:
            negative += 1
            worst_negative = min(worst_negative, value)

sum_squares, count = 0.0, 0
species_sums = {}
for key, reference_values in reference.items():
    values = result[key]
    for j, name in enumerate(reference_species):
        y_ref = reference_values[j]
        if y_ref < args.threshold:
            continue
        y = values[column[name]]
        error = abs(y - y_ref) / y_ref if math.isfinite(y) else math.inf
        squared = error * error
        sum_squares += squared
        count += 1
        total, n = species_sums.get(name, (0.0, 0))
        species_sums[name] = (total + squared, n + 1)

overall = digits(sum_squares / count) if count else float("nan")
per_species = {name: digits(total / n) for name, (total, n) in species_sums.items()}
worst = sorted(per_species.items(), key=lambda item: item[1])[:5]
gate = nonfinite == 0 and negative == 0
summary = {
    "overall_digits": round(overall, 3),
    "worst_species": [[name, round(value, 3)] for name, value in worst],
    "values_counted": count,
    "species_counted": len(per_species),
    "rows": len(reference),
    "nonfinite": nonfinite,
    "below_floor": negative,
    "gate": "PASS" if gate else "FAIL",
}
passed = gate
if args.level is not None:
    level_ok = overall >= args.level and all(value >= args.level - 1 for value in per_species.values())
    summary["level"] = args.level
    summary["level_result"] = "PASS" if level_ok else "FAIL"
    passed = gate and level_ok

if args.json:
    print(json.dumps(summary))
else:
    print(
        f"overall {overall:.2f} digits over {count} values ({len(per_species)} species, {len(reference)} rows);"
        f" worst species: {', '.join(f'{name} {value:.2f}' for name, value in worst)}"
    )
    print(
        f"gate {summary['gate']}: {nonfinite} NaN or infinite values, {negative} below {FLOOR:g}"
        + (f" (lowest {worst_negative:.3g})" if negative else "")
    )
    if args.level is not None:
        print(f"level {args.level:g}: {summary['level_result']}")
sys.exit(0 if passed else 1)
