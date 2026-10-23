#!/usr/bin/env python3
"""Convert a canonical mechanism (canonical.json) to luna-xhigh's descriptor JSON.

    python3 convert.py CANONICAL_JSON OUTPUT_JSON

The target is the schema of the entry's config/ts1.json (written by its
tools/translate_ts1.py), which tools/generate.py accepts directly when the file
has top-level "species" and "reactions":

  {"name", "temperature", "pressure",
   "species":   [{"name", "molecular_weight", "diffusion_coefficient", "third_body", "initial"}],
   "parameters": {label: initial value},
   "reactions": [{"type", "reactants", "products", "label", "parameters"}]}

Mapping, mirroring TS1:
  ARRHENIUS    -> ArrheniusRateConstant, parameters A B C D E (MICM's, mol m-3 units)
  TROE         -> TroeRateConstant, parameters k0_A k0_B k0_C kinf_A kinf_B kinf_C Fc N
  USER_DEFINED -> UserDefinedRateConstant, label = the rate parameter's name,
                  parameters {"scaling_factor"}; emissions keep an empty reactant list
  third body   -> species "M" with "third_body": true and initial 0.0; it stays in the
                  reactant and product lists by name, as in TS1 (the generator turns
                  each reactant M into a factor of air density and drops product M)
  reactants    -> names, repeated for stoichiometry; products -> {"species", "yield"}
  initial      -> top-level temperature and pressure, each species' "initial"
                  (0.0 when the canonical file gives none), and the "parameters"
                  values (0.0 when none is given)
The canonical file gives no molecular weights or diffusion coefficients; they are
used only by surface reactions, which these mechanisms lack, so they are written
as 0.0 and null.
"""

import json
import sys

ARRHENIUS_KEYS = ("A", "B", "C", "D", "E")
TROE_KEYS = ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")


def convert(canonical: dict) -> dict:
    initial = canonical.get("initial", {})
    concentrations = initial.get("concentrations", {})
    parameter_values = initial.get("parameters", {})
    names = {item["name"] for item in canonical["species"]}

    species = []
    for item in canonical["species"]:
        third_body = bool(item.get("third_body", False))
        species.append({
            "name": item["name"],
            "molecular_weight": 0.0,
            "diffusion_coefficient": None,
            "third_body": third_body,
            "initial": 0.0 if third_body else float(concentrations.get(item["name"], 0.0)),
        })

    parameters = {label: float(parameter_values.get(label, 0.0)) for label in canonical.get("parameters", [])}

    reactions = []
    for index, item in enumerate(canonical["reactions"]):
        kind = item["type"]
        for name in item["reactants"]:
            if name not in names:
                raise ValueError(f"reaction {index}: unknown reactant {name!r}")
        products = []
        for name, yield_ in item["products"]:
            if name not in names:
                raise ValueError(f"reaction {index}: unknown product {name!r}")
            products.append({"species": name, "yield": float(yield_)})
        if kind == "ARRHENIUS":
            out = {"type": "ArrheniusRateConstant", "label": None,
                   "parameters": {key: float(item[key]) for key in ARRHENIUS_KEYS}}
        elif kind == "TROE":
            out = {"type": "TroeRateConstant", "label": None,
                   "parameters": {key: float(item[key]) for key in TROE_KEYS}}
        elif kind == "USER_DEFINED":
            label = item["label"]
            if label not in parameters:
                parameters[label] = float(parameter_values.get(label, 0.0))
            out = {"type": "UserDefinedRateConstant", "label": label,
                   "parameters": {"scaling_factor": float(item.get("scaling_factor", 1.0))}}
        else:
            raise ValueError(f"reaction {index}: no mapping for reaction type {kind!r}")
        reactions.append({"type": out["type"], "reactants": list(item["reactants"]), "products": products,
                          "label": out["label"], "parameters": out["parameters"]})

    return {
        "name": canonical["name"],
        "temperature": float(initial["temperature"]),
        "pressure": float(initial["pressure"]),
        "species": species,
        "parameters": parameters,
        "reactions": reactions,
    }


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as stream:
        canonical = json.load(stream)
    mechanism = convert(canonical)
    with open(sys.argv[2], "w") as stream:
        stream.write(json.dumps(mechanism, indent=2) + "\n")
    print(f"wrote {sys.argv[2]}: {len(mechanism['species'])} species, {len(mechanism['reactions'])} reactions, "
          f"{len(mechanism['parameters'])} parameters")


if __name__ == "__main__":
    main()
