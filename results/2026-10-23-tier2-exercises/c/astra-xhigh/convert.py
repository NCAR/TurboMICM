#!/usr/bin/env python3
"""Convert a canonical mechanism (canonical.json) to astra-xhigh's configuration format.

    convert.py CANONICAL_JSON OUTPUT_JSON

The output mirrors the entry's configs/ts1.json, a CAMP-style document read by
its tools/configure.py:
  {"camp-data": [CHEM_SPEC entries..., {"type": "MECHANISM", "reactions": [...]}],
   "initial conditions": {"temperature", "pressure"},
   "initial rate parameters": {label: value}}
- a third body (canonical "third_body": true) is a CHEM_SPEC with
  "parameterization": "air density", as TS1's M is; it stays in reactant and
  product lists exactly where the canonical file has it, as in TS1;
- every other species gets "initial concentration [mol m-3]" (0.0 when the
  canonical file gives none, as MICM's InitState and cells.csv do);
- every reaction carries "concentration units": "mol m-3", as TS1's do (the
  canonical parameters are already in mol m-3; the reader's default,
  "molecule cm-3", would rescale A);
- ARRHENIUS keeps A, B, C, D, E; TROE keeps k0_A..kinf_C, Fc, N;
- USER_DEFINED keeps its label as "MUSICA name" and its factor as
  "scaling factor"; an emission (no reactants) stays USER_DEFINED with an
  empty reactant list, which the reader accepts;
- products are {"species", "yield"} objects in the canonical order; reactants
  keep their repeats;
- rate-parameter initial values go to "initial rate parameters"; labels
  without an initial value are left out (the state starts them at zero, and
  cells.csv sets every label).
Fails on any reaction type or field it does not know. Standard library only.
"""

import json
import sys

ARRHENIUS = ("A", "B", "C", "D", "E")
TROE = ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")
COMMON = {"type", "reactants", "products"}
KNOWN = {
    "ARRHENIUS": COMMON | set(ARRHENIUS),
    "TROE": COMMON | set(TROE),
    "USER_DEFINED": COMMON | {"label", "scaling_factor", "source"},
}


def convert(canonical):
    species = canonical["species"]
    names = [s["name"] for s in species]
    if len(set(names)) != len(names):
        raise ValueError("duplicate species names")
    initial = canonical["initial"]
    concentrations = initial.get("concentrations", {})
    unknown = set(concentrations) - set(names)
    if unknown:
        raise ValueError(f"initial concentrations for unknown species: {sorted(unknown)}")
    camp = []
    for s in species:
        if set(s) - {"name", "third_body"}:
            raise ValueError(f"unknown species fields: {s}")
        entry = {"name": s["name"], "type": "CHEM_SPEC"}
        if s.get("third_body"):
            if s["name"] in concentrations:
                raise ValueError(f"third body {s['name']} has a concentration")
            entry["parameterization"] = "air density"
        else:
            entry["initial concentration [mol m-3]"] = float(concentrations.get(s["name"], 0.0))
        camp.append(entry)

    reactions = []
    for index, r in enumerate(canonical["reactions"]):
        kind = r["type"]
        if kind not in KNOWN:
            raise ValueError(f"reaction {index}: unsupported type {kind}")
        extra = set(r) - KNOWN[kind]
        if extra:
            raise ValueError(f"reaction {index}: unknown fields {sorted(extra)}")
        for name in list(r["reactants"]) + [p[0] for p in r["products"]]:
            if name not in names:
                raise ValueError(f"reaction {index}: unknown species {name}")
        out = {
            "type": kind,
            "reactants": list(r["reactants"]),
            "products": [{"species": name, "yield": float(y)} for name, y in r["products"]],
            "concentration units": "mol m-3",
        }
        if kind == "ARRHENIUS":
            for key in ARRHENIUS:
                out[key] = float(r[key])
        elif kind == "TROE":
            for key in TROE:
                out[key] = float(r[key])
        else:
            out["MUSICA name"] = r["label"]
            out["scaling factor"] = float(r["scaling_factor"])
        reactions.append(out)
    camp.append({"name": canonical["name"], "type": "MECHANISM", "reactions": reactions})

    labels = {r["label"] for r in canonical["reactions"] if r["type"] == "USER_DEFINED"}
    parameters = initial.get("parameters", {})
    unknown = set(parameters) - labels
    if unknown:
        raise ValueError(f"initial values for unknown rate parameters: {sorted(unknown)}")
    return {
        "camp-data": camp,
        "initial conditions": {"temperature": float(initial["temperature"]),
                               "pressure": float(initial["pressure"])},
        "initial rate parameters": {k: float(v) for k, v in parameters.items()},
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        document = convert(json.load(f))
    with open(sys.argv[2], "w") as f:
        f.write(json.dumps(document, indent=2) + "\n")
    mechanism = document["camp-data"][-1]
    print(f"{sys.argv[2]}: {len(document['camp-data']) - 1} species, "
          f"{len(mechanism['reactions'])} reactions, "
          f"{len(document['initial rate parameters'])} rate-parameter initial values")
