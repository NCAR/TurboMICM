#!/usr/bin/env python3
"""Resolve a MUSICA v0 (CAMP) mechanism into one canonical form for the Tier 2 exercises.

    admin/v0_to_canonical.py MECHANISM_DIR OUT_DIR [--cells 4] [--drop TYPE]

The review's third exercise runs each entry's generator on another mechanism:
Chapman and Carbon Bond 5 from MUSICA's configs/v0. So that every entry solves
the same problem, this script settles the configuration's meaning once, as
MICM's own v0 parser did before MICM removed it (MICM 6755647d^,
include/micm/configure/solver_config.hpp):
- rate constants in mol m-3 units: an Arrhenius A times C^(n-1), a Troe k0_A
  times C^n and kinf_A times C^(n-1), where n counts the reactants (each "qty"
  times, the third body included) and C = 1e-6 N_A converts molecule cm-3 to
  mol m-3; an Arrhenius Ea becomes C = -Ea / k_B; MICM's defaults fill the
  missing parameters;
- PHOTOLYSIS and EMISSION become user-defined rates, a label times a scaling
  factor; an emission has no reactants and makes one mole of its species;
- the third body: a species with tracer type THIRD_BODY, or named M (Carbon
  Bond 5 declares M without one; MICM's parser would have made it an ordinary
  species starting at zero, which switches off its five reactions). As in TS1,
  M holds no concentration: MICM parameterizes it on the air density.
Labels become C identifiers: MICM's PHOTO.<name> is PHOTO_<name> and
EMIS.<name> is EMIS_<name>, with "->" written "_to_" and any other character
outside [A-Za-z0-9_] written "_". Chapman's initial conditions name its
photolysis rates R1, R3, and R5, the numbers of its reactions; they are matched
to those reactions. A rate with no initial value is zero, as in MICM (Carbon
Bond 5 gives none for five photolysis reactions); canonical.json lists them.

OUT_DIR/canonical.json holds the species (in the configuration's order), the
reactions, and the initial state: temperature, pressure, concentrations, and
rate parameters. OUT_DIR/cells.csv holds the same state for CELLS identical
cells in the TS1 benchmarks' format (temperature, pressure, every species but
the third body, every rate parameter). --drop EMISSION leaves out that v0
reaction type and its rate parameters, for a variant without emissions. Needs
only the Python standard library.
"""

import argparse
import csv
import json
import re
from pathlib import Path

AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
CONVERSION = 1.0e-6 * AVOGADRO  # MOLES_M3_TO_MOLECULES_CM3 in MICM's parser
ARRHENIUS_DEFAULTS = {"A": 1.0, "B": 0.0, "C": 0.0, "D": 300.0, "E": 0.0}
TROE_DEFAULTS = {"k0_A": 1.0, "k0_B": 0.0, "k0_C": 0.0, "kinf_A": 1.0, "kinf_B": 0.0, "kinf_C": 0.0,
                 "Fc": 0.6, "N": 1.0}

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mechanism", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("--cells", type=int, default=4)
parser.add_argument("--drop", action="append", default=[], help="a v0 reaction type to leave out")
args = parser.parse_args()

records = []
for name in json.loads((args.mechanism / "config.json").read_text())["camp-files"]:
    document = json.loads((args.mechanism / name).read_text())
    records += document.get("camp-data", [])

species = []
for record in records:
    if record.get("type") == "CHEM_SPEC":
        third_body = record.get("tracer type") == "THIRD_BODY" or record["name"] == "M"
        species.append({"name": record["name"], "third_body": third_body})
names = {s["name"] for s in species}


def identifier(kind, name):
    return kind + "_" + re.sub(r"[^A-Za-z0-9_]", "_", name.replace("->", "_to_"))


def reactants_of(reaction):
    expanded = []
    for name, properties in (reaction.get("reactants") or {}).items():
        if name not in names:
            raise SystemExit(f"unknown reactant {name}")
        expanded += [name] * int(properties.get("qty", 1))
    return expanded


def products_of(reaction):
    products = []
    for name, properties in (reaction.get("products") or {}).items():
        if name not in names:
            raise SystemExit(f"unknown product {name}")
        products.append([name, float(properties.get("yield", 1.0))])
    return products


reactions = []
mechanisms = [record for record in records if record.get("type") == "MECHANISM"]
dropped = {}
for reaction in (r for mechanism in mechanisms for r in mechanism["reactions"]):
    kind = reaction["type"]
    if kind in args.drop:
        dropped[kind] = dropped.get(kind, 0) + 1
        continue
    if kind == "ARRHENIUS":
        reactants = reactants_of(reaction)
        parameters = {key: float(reaction.get(key, default)) for key, default in ARRHENIUS_DEFAULTS.items()}
        if "Ea" in reaction:
            if reaction.get("C", 0.0) != 0.0:
                raise SystemExit("an Arrhenius reaction gives both C and Ea")
            parameters["C"] = -float(reaction["Ea"]) / BOLTZMANN + 0.0  # no negative zero
        parameters["A"] *= CONVERSION ** (len(reactants) - 1)
        reactions.append({"type": "ARRHENIUS", "reactants": reactants, "products": products_of(reaction), **parameters})
    elif kind == "TROE":
        reactants = reactants_of(reaction)
        parameters = {key: float(reaction.get(key, default)) for key, default in TROE_DEFAULTS.items()}
        parameters["k0_A"] *= CONVERSION ** len(reactants)
        parameters["kinf_A"] *= CONVERSION ** (len(reactants) - 1)
        reactions.append({"type": "TROE", "reactants": reactants, "products": products_of(reaction), **parameters})
    elif kind == "PHOTOLYSIS":
        reactions.append({"type": "USER_DEFINED", "reactants": reactants_of(reaction), "products": products_of(reaction),
                          "label": identifier("PHOTO", reaction["MUSICA name"]),
                          "scaling_factor": float(reaction.get("scaling factor", 1.0)), "source": "PHOTOLYSIS"})
    elif kind == "EMISSION":
        if reaction["species"] not in names:
            raise SystemExit(f"unknown emitted species {reaction['species']}")
        reactions.append({"type": "USER_DEFINED", "reactants": [], "products": [[reaction["species"], 1.0]],
                          "label": identifier("EMIS", reaction["MUSICA name"]),
                          "scaling_factor": float(reaction.get("scaling factor", 1.0)), "source": "EMISSION"})
    else:
        raise SystemExit(f"reaction type {kind} is not handled")

initial = {"temperature": None, "pressure": None, "concentrations": {}, "parameters": {}}
labels = {r["label"] for r in reactions if r["type"] == "USER_DEFINED"}
with open(args.mechanism / "initial_conditions.csv", newline="") as file:
    for key, value in csv.reader(file):
        kind, _, name = key.partition(".")
        value = float(value)
        if kind == "ENV" and name in ("temperature", "pressure"):
            initial[name] = value
        elif kind == "ENV":
            initial.setdefault("environment", {})[name] = value
        elif kind == "CONC":
            if name not in names or name == "M":
                raise SystemExit(f"initial concentration for {name}, which is not a species with a concentration")
            initial["concentrations"][name] = value
        elif kind in ("PHOTO", "EMIS"):
            label = identifier(kind, name)
            numbered = re.fullmatch(r"R(\d+)", name)
            if label not in labels and kind == "PHOTO" and numbered:
                reaction = reactions[int(numbered.group(1)) - 1]
                if reaction["type"] != "USER_DEFINED":
                    raise SystemExit(f"PHOTO.{name}: reaction {numbered.group(1)} is not a photolysis")
                label = reaction["label"]
            if label not in labels and kind == "EMIS" and "EMISSION" in args.drop:
                continue
            if label not in labels:
                raise SystemExit(f"{key} matches no reaction")
            initial["parameters"][label] = value
        else:
            raise SystemExit(f"unexpected initial condition {key}")
if initial["temperature"] is None or initial["pressure"] is None:
    raise SystemExit("the initial conditions give no temperature or pressure")
zero = sorted(labels - set(initial["parameters"]))
for label in zero:
    initial["parameters"][label] = 0.0

counts = {}
for reaction in reactions:
    key = reaction.get("source", reaction["type"])
    counts[key] = counts.get(key, 0) + 1
canonical = {
    "name": args.out.name,
    "source": f"MUSICA configs/v0/{args.mechanism.name}, resolved by admin/v0_to_canonical.py",
    "species": species,
    "reactions": reactions,
    "parameters": sorted(labels, key=lambda label: [r.get("label") for r in reactions].index(label)),
    "initial": initial,
    "parameters_without_initial_value": zero,
    "dropped": dropped,
    "summary": {"species": len(species), "reactions": len(reactions), "by_source": counts},
}
args.out.mkdir(parents=True, exist_ok=True)
(args.out / "canonical.json").write_text(json.dumps(canonical, indent=1) + "\n")
columns = (["temperature", "pressure"] + [s["name"] for s in species if not s["third_body"]] + canonical["parameters"])
values = {"temperature": initial["temperature"], "pressure": initial["pressure"],
          **initial["concentrations"], **initial["parameters"]}
with open(args.out / "cells.csv", "w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(columns)
    for _ in range(args.cells):
        writer.writerow([repr(float(values.get(column, 0.0))) for column in columns])
print(json.dumps(canonical["summary"]))
