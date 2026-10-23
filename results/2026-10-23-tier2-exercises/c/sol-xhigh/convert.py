#!/usr/bin/env python3
"""Exercise C: write a canonical mechanism in sol-xhigh's input format.

    convert.py CANONICAL_JSON OUTPUT_JSON

The output mirrors the entry's TS1 configuration, mechanisms/ts1.json: one
CAMP JSON file with "concentration units", "initial state", and "camp-data"
(CHEM_SPEC entries, then one MECHANISM entry). Read by the entry's
tools/configure.py (generic descriptors) and tools/specialize.py (generator).

  - Units: "mol m-3", as in ts1.json; the canonical parameters are already in
    mol m-3 units, so they pass through unchanged (the reader's factor is 1).
  - Third body: the canonical "third_body": true species becomes, as M is in
    ts1.json, {"name": "M", "type": "CHEM_SPEC", "parameterization":
    "air-density"}; M stays in reactant and product lists, as in ts1.json
    (the reader multiplies k by the air density once per M reactant and drops
    M from the state).
  - ARRHENIUS: reactants (a list, repeats for stoichiometry), products as
    [{"species", "yield"}], then A, B, C, D, E.
  - TROE: k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N.
  - USER_DEFINED: "scaling_factor" and the label as "MUSICA name", as in
    ts1.json; emissions keep their empty reactant list.
  - Initial values: "initial state" with temperature, pressure,
    "concentrations" (the species the canonical file gives; others default
    to 0 in the reader) and "custom" (each label's value; labels without one
    are left out, so the library's zero default applies).
Only the standard library is used.
"""
import json
import math
import sys

ARRHENIUS = ("A", "B", "C", "D", "E")
TROE = ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")


def number(value, what):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise SystemExit(f"{what}: not a finite number: {value!r}")
    return float(value)


def convert(canonical):
    species = canonical["species"]
    names = [s["name"] for s in species]
    if len(set(names)) != len(names):
        raise SystemExit("duplicate species names")
    known = set(names)
    camp = []
    for s in species:
        entry = {"name": s["name"], "type": "CHEM_SPEC"}
        if s.get("third_body"):
            entry["parameterization"] = "air-density"
        camp.append(entry)

    reactions, labels = [], []
    for i, r in enumerate(canonical["reactions"]):
        kind = r["type"]
        for name in list(r["reactants"]) + [p[0] for p in r["products"]]:
            if name not in known:
                raise SystemExit(f"reaction {i}: unknown species {name}")
        out = {"type": kind,
               "reactants": list(r["reactants"]),
               "products": [{"species": name, "yield": number(y, f"reaction {i} yield")}
                            for name, y in r["products"]]}
        if kind == "ARRHENIUS":
            for key in ARRHENIUS:
                out[key] = number(r[key], f"reaction {i} {key}")
        elif kind == "TROE":
            for key in TROE:
                out[key] = number(r[key], f"reaction {i} {key}")
        elif kind == "USER_DEFINED":
            out["scaling_factor"] = number(r["scaling_factor"], f"reaction {i} scaling_factor")
            out["MUSICA name"] = r["label"]
            labels.append(r["label"])
        else:
            raise SystemExit(f"reaction {i}: no mapping for rate type {kind}")
        reactions.append(out)

    initial = canonical.get("initial", {})
    third = {s["name"] for s in species if s.get("third_body")}
    concentrations = {}
    for name, value in initial.get("concentrations", {}).items():
        if name not in known or name in third:
            raise SystemExit(f"initial concentration for {name}: not an integrated species")
        concentrations[name] = number(value, f"initial {name}")
    custom = {}
    for label, value in initial.get("parameters", {}).items():
        if label not in labels:
            raise SystemExit(f"initial value for {label}: no reaction uses this label")
        custom[label] = number(value, f"initial {label}")

    return {
        "concentration units": "mol m-3",
        "initial state": {
            "temperature": number(initial["temperature"], "temperature"),
            "pressure": number(initial["pressure"], "pressure"),
            "concentrations": concentrations,
            "custom": custom,
        },
        "camp-data": camp + [{"name": canonical.get("name", "mechanism"), "type": "MECHANISM",
                              "reactions": reactions}],
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    with open(sys.argv[1]) as stream:
        canonical = json.load(stream)
    result = convert(canonical)
    with open(sys.argv[2], "w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
