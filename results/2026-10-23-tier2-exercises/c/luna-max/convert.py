#!/usr/bin/env python3
"""Convert a canonical mechanism (canonical.json) to luna-max's mechanism JSON.

    convert.py CANONICAL_JSON OUTPUT_JSON

The output mirrors the entry's mechanisms/ts1.json, which its
tools/translate_ts1.py writes from MICM's TS1 header:

  {"name", "temperature", "pressure",
   "species": [{"name", "molecular_weight", "diffusion_coefficient",
                "third_body", "initial"}],
   "reactions": [{"type", "reactants", "products": [{"name", "yield"}],
                  "parameters": {...}, "label"?}],
   "custom_parameters": {label: value}}

* ARRHENIUS: parameters {A, B, C, D, E}, as in ts1.json.
* TROE: parameters {k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N}.
* USER_DEFINED (photolysis and emissions): parameters {scaling_factor} plus a
  top-level "label"; an emission keeps its empty reactant list.
* The third body M is a species with "third_body": true and initial 0.0; it
  stays in reactant and product lists by name, as in TS1's
  O2 + M + O -> O3 + M (the generator turns it into a factor of air density).
* Initial values: temperature, pressure, species concentrations (0.0 when
  unset), and every rate parameter in custom_parameters (0.0 when the
  canonical file gives none, as MICM's InitState sets them).
* Species have no molecular weight or diffusion coefficient in the canonical
  input; they are written as 0.0 and null, as ts1.json does for species
  without them (only SURFACE reactions read them).
"""

import json
import sys

ARRHENIUS_KEYS = ("A", "B", "C", "D", "E")
TROE_KEYS = ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")


def fail(message):
    raise SystemExit(f"convert.py: {message}")


def convert(canonical):
    initial = canonical["initial"]
    concentrations = initial.get("concentrations", {})
    names = [s["name"] for s in canonical["species"]]
    if len(set(names)) != len(names):
        fail("duplicate species names")
    third_bodies = {s["name"] for s in canonical["species"] if s.get("third_body")}
    for name in concentrations:
        if name not in names:
            fail(f"initial concentration for unknown species {name!r}")
        if name in third_bodies:
            fail(f"initial concentration for third body {name!r}")
    species = [
        {
            "name": s["name"],
            "molecular_weight": 0.0,
            "diffusion_coefficient": None,
            "third_body": bool(s.get("third_body", False)),
            "initial": float(concentrations.get(s["name"], 0.0)),
        }
        for s in canonical["species"]
    ]

    reactions = []
    labels = []
    for index, raw in enumerate(canonical["reactions"]):
        kind = raw["type"]
        for name in raw["reactants"]:
            if name not in names:
                fail(f"reaction {index}: unknown reactant {name!r}")
        products = []
        for pair in raw["products"]:
            name, value = pair
            if name not in names:
                fail(f"reaction {index}: unknown product {name!r}")
            products.append({"name": name, "yield": float(value)})
        reaction = {"type": kind, "reactants": list(raw["reactants"]), "products": products}
        if kind == "ARRHENIUS":
            reaction["parameters"] = {key: float(raw[key]) for key in ARRHENIUS_KEYS}
        elif kind == "TROE":
            reaction["parameters"] = {key: float(raw[key]) for key in TROE_KEYS}
        elif kind == "USER_DEFINED":
            reaction["parameters"] = {"scaling_factor": float(raw["scaling_factor"])}
            reaction["label"] = raw["label"]
            labels.append(raw["label"])
        else:
            fail(f"reaction {index}: unsupported type {kind!r}")
        reactions.append(reaction)

    parameter_names = list(canonical.get("parameters", []))
    for label in labels:
        if label not in parameter_names:
            parameter_names.append(label)
    values = initial.get("parameters", {})
    for name in values:
        if name not in parameter_names:
            fail(f"initial value for unknown rate parameter {name!r}")
    custom_parameters = {name: float(values.get(name, 0.0)) for name in parameter_names}

    return {
        "name": canonical["name"],
        "temperature": float(initial["temperature"]),
        "pressure": float(initial["pressure"]),
        "species": species,
        "reactions": reactions,
        "custom_parameters": custom_parameters,
    }


def main():
    if len(sys.argv) != 3:
        fail("usage: convert.py CANONICAL_JSON OUTPUT_JSON")
    with open(sys.argv[1]) as file:
        canonical = json.load(file)
    document = convert(canonical)
    with open(sys.argv[2], "w") as file:
        file.write(json.dumps(document, indent=2) + "\n")
    print(f"wrote {len(document['species'])} species, {len(document['reactions'])} reactions, "
          f"{len(document['custom_parameters'])} rate parameters to {sys.argv[2]}")


if __name__ == "__main__":
    main()
