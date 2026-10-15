#!/usr/bin/env python3
"""Write a TS1 CELLS_CSV in which every cell differs.

    make_cells_csv.py CELLS SEED [--output FILE]

Each cell perturbs TS1's InitState (ts1_defaults.csv, written by
`micm_ts1_ref --write-defaults`). The same CELLS and SEED always give the same
file; it needs only the Python standard library.

    temperature        uniform in 260-310 K
    pressure           InitState pressure times a uniform factor in 0.6-1.0
    O2, N2             InitState mixing ratio times the cell's air density
    other species      InitState concentration, scaled with air density, times
                       a lognormal factor (sigma 0.5)
    photolysis (j*)    InitState rate times a sun factor, uniform in 0-1 and
                       shared by the cell's rates, times a lognormal factor
                       (sigma 0.1); jo2_a, jh2o_b, jh2o_c, jch4_a, jch4_b, and
                       jhf are 0, because InitState's 1 s-1 is unphysical
    het*, usr_*_aer    InitState value times a lognormal factor (sigma 0.3)

Columns whose InitState value is 0 (sink, usr_CO_OH, usr_DMS_OH) are left
out, so those take InitState's value.
"""

import argparse
import csv
import math
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1
UNPHYSICAL = {"jo2_a", "jh2o_b", "jh2o_c", "jch4_a", "jch4_b", "jhf"}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("cells", type=int)
parser.add_argument("seed", type=int)
parser.add_argument("--output", type=Path, help="default: ts1_cells_<CELLS>_seed<SEED>.csv")
parser.add_argument("--defaults", type=Path, default=HERE / "ts1_defaults.csv")
args = parser.parse_args()

with open(args.defaults, newline="") as file:
    names, values = list(csv.reader(file))
defaults = {name: float(value) for name, value in zip(names, values)}
columns = [name for name in names if defaults[name] != 0.0]
air_density_0 = defaults["pressure"] / (GAS_CONSTANT * defaults["temperature"])

rng = random.Random(args.seed)
output = args.output or Path(f"ts1_cells_{args.cells}_seed{args.seed}.csv")
with open(output, "w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(columns)
    for _ in range(args.cells):
        temperature = rng.uniform(260.0, 310.0)
        pressure = defaults["pressure"] * rng.uniform(0.6, 1.0)
        density_ratio = pressure / (GAS_CONSTANT * temperature) / air_density_0
        sun = rng.uniform(0.0, 1.0)
        row = []
        for name in columns:
            default = defaults[name]
            if name == "temperature":
                value = temperature
            elif name == "pressure":
                value = pressure
            elif name in ("O2", "N2"):
                value = default * density_ratio
            elif name in UNPHYSICAL:
                value = 0.0
            elif name.startswith("j"):
                value = default * sun * math.exp(rng.gauss(0.0, 0.1))
            elif name.startswith(("het", "usr_")):
                value = default * math.exp(rng.gauss(0.0, 0.3))
            else:
                value = default * density_ratio * math.exp(rng.gauss(0.0, 0.5))
            row.append(repr(value))
        writer.writerow(row)
print(f"wrote {output}: {args.cells} cells, {len(columns)} columns")
