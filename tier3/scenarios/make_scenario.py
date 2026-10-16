#!/usr/bin/env python3
"""Write a Tier 3 TS1 scenario: cells, a photolysis factor f(t), and the scenario file.

    make_scenario.py DIR --cells N --seed S --dt DT --steps K
                     [--diurnal --lat DEG --doy DAY --start HOUR | --constant]
                     [--nox-scale X] [--snapshot-every M]

DIR gets cells.csv, photolysis.csv, and scenario.txt (see README.md).

Cells perturb TS1's InitState (ROOT/reference/ts1_defaults.csv) as Tier 1's
make_cells_csv.py does: temperature uniform in 260-310 K; pressure 0.6-1.0 times
InitState's; O2 and N2 scaled with air density; other species scaled with air
density times a lognormal factor (sigma 0.5); heterogeneous and aerosol parameters
times a lognormal factor (sigma 0.3); the six photolysis rates that InitState sets
unphysically to 1 s-1 set to 0. Photolysis rates (labels beginning with "j") are
InitState's times a lognormal factor (sigma 0.1), and then:
- with --diurnal, they are noon values, and f(t) = max(0, cos(SZA(t))) / cos(SZA(noon))
  for latitude --lat, day of year --doy, starting at local solar time --start;
- with --constant, they are also multiplied by a per-cell sun factor, uniform in
  0-1, and f(t) = 1, as in Tier 1.
--nox-scale multiplies the initial NO, NO2, NO3, N2O5, and HNO3 (clean or polluted
cells). The same arguments always give the same files; needs only the Python
standard library.
"""

import argparse
import csv
import math
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1
UNPHYSICAL = {"jo2_a", "jh2o_b", "jh2o_c", "jch4_a", "jch4_b", "jhf"}
NOX = {"NO", "NO2", "NO3", "N2O5", "HNO3"}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("directory", type=Path)
parser.add_argument("--cells", type=int, required=True)
parser.add_argument("--seed", type=int, required=True)
parser.add_argument("--dt", type=float, required=True, help="seconds per step")
parser.add_argument("--steps", type=int, required=True)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument("--diurnal", action="store_true")
mode.add_argument("--constant", action="store_true")
parser.add_argument("--lat", type=float, default=40.0, help="latitude, degrees north (default 40)")
parser.add_argument("--doy", type=int, default=80, help="day of year (default 80, the March equinox)")
parser.add_argument("--start", type=float, default=6.0, help="local solar time at t = 0, hours (default 6)")
parser.add_argument("--nox-scale", type=float, default=1.0)
parser.add_argument("--snapshot-every", type=int, default=0, help="steps between snapshots (default: final only)")
parser.add_argument("--defaults", type=Path, default=HERE.parent / "reference" / "ts1_defaults.csv")
args = parser.parse_args()

with open(args.defaults, newline="") as file:
    names, values = list(csv.reader(file))
defaults = {name: float(value) for name, value in zip(names, values)}
columns = [name for name in names if defaults[name] != 0.0]
air_density_0 = defaults["pressure"] / (GAS_CONSTANT * defaults["temperature"])

args.directory.mkdir(parents=True, exist_ok=True)
rng = random.Random(args.seed)
with open(args.directory / "cells.csv", "w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(columns)
    for _ in range(args.cells):
        temperature = rng.uniform(260.0, 310.0)
        pressure = defaults["pressure"] * rng.uniform(0.6, 1.0)
        density_ratio = pressure / (GAS_CONSTANT * temperature) / air_density_0
        sun = 1.0 if args.diurnal else rng.uniform(0.0, 1.0)
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
                if name in NOX:
                    value *= args.nox_scale
            row.append(repr(value))
        writer.writerow(row)


def cos_zenith(hours):
    """Cosine of the solar zenith angle at local solar time `hours` (any real number)."""
    declination = math.radians(23.44) * math.sin(2.0 * math.pi * (284 + args.doy) / 365.0)
    latitude = math.radians(args.lat)
    hour_angle = 2.0 * math.pi * (hours - 12.0) / 24.0
    return math.sin(latitude) * math.sin(declination) + math.cos(latitude) * math.cos(declination) * math.cos(hour_angle)


duration = args.dt * args.steps
with open(args.directory / "photolysis.csv", "w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(["time_s", "factor"])
    if args.constant:
        writer.writerow([0, 1])
    else:
        noon = cos_zenith(12.0)
        if noon <= 0.0:
            raise SystemExit("the sun does not rise at this latitude and day")
        sample = 300.0  # seconds between tabulated points
        for i in range(int(math.ceil(duration / sample)) + 1):
            t = i * sample
            factor = max(0.0, cos_zenith(args.start + t / 3600.0)) / noon
            writer.writerow([repr(t), repr(factor)])

description = (
    f"{args.cells} cells (seed {args.seed}, NOx x{args.nox_scale:g}); "
    f"{args.steps} steps of {args.dt:g} s; "
    + (f"diurnal photolysis, {args.lat:g} N, day {args.doy}, from {args.start:g} h local time"
       if args.diurnal else "constant photolysis")
)
with open(args.directory / "scenario.txt", "w") as file:
    file.write(f"# {description}\n")
    file.write("cells = cells.csv\n")
    file.write(f"dt = {args.dt:g}\n")
    file.write(f"steps = {args.steps}\n")
    file.write("photolysis = photolysis.csv\n")
    if args.snapshot_every:
        file.write(f"snapshot_every = {args.snapshot_every}\n")
print(f"wrote {args.directory}: {description}")
