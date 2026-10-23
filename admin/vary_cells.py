#!/usr/bin/env python3
"""Write a cells file of distinct cells for a canonical mechanism.

    admin/vary_cells.py CANONICAL.json OUT.csv CELLS SEED

For exercise (c) of the Tier 2 review: cells.csv from admin/v0_to_canonical.py
holds identical cells, too few to fill a 32-cell group, so an entry whose
generated kernels run only on full groups falls back to its generic code. This
file has CELLS cells (use a multiple of 32), each perturbing the canonical
initial state, in the same column format:

    temperature        times a uniform factor in 0.9-1.1
    pressure           times a uniform factor in 0.8-1.2
    concentrations     times a lognormal factor (sigma 0.3); zeros stay zero
    PHOTO_* rates      times a sun factor, uniform in 0-1 and shared by the cell
    EMIS_* rates       times a lognormal factor (sigma 0.3)

The same arguments always give the same file. Needs only the Python standard
library.
"""

import argparse
import csv
import json
import math
import random
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("canonical", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("cells", type=int)
parser.add_argument("seed", type=int)
args = parser.parse_args()

mechanism = json.loads(args.canonical.read_text())
initial = mechanism["initial"]
species = [s["name"] for s in mechanism["species"] if not s["third_body"]]
columns = ["temperature", "pressure"] + species + mechanism["parameters"]
rng = random.Random(args.seed)
with open(args.out, "w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(columns)
    for _ in range(args.cells):
        sun = rng.uniform(0.0, 1.0)
        row = [initial["temperature"] * rng.uniform(0.9, 1.1), initial["pressure"] * rng.uniform(0.8, 1.2)]
        row += [initial["concentrations"].get(name, 0.0) * math.exp(rng.gauss(0.0, 0.3)) for name in species]
        for label in mechanism["parameters"]:
            value = initial["parameters"][label]
            row.append(value * sun if label.startswith("PHOTO_") else value * math.exp(rng.gauss(0.0, 0.3)))
        writer.writerow([repr(float(value)) for value in row])
