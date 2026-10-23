#!/usr/bin/env python3
"""Write a canonical mechanism in astra-high's input format (CAMP JSON plus the
entry's documented extensions), mirroring how mechanisms/ts1.json expresses TS1.

    convert.py CANONICAL_JSON OUTPUT_JSON

Mapping, each item as in ts1.json:
- species -> {"name", "type": "CHEM_SPEC"}; the third body M gets
  "parameterization": "air density" (it holds no concentration).
- one MECHANISM record with "rate units": "mol m-3" (parameters are already in
  mol m-3, so the reader must not convert them again).
- reactants -> ordered list of names with repeats (M stays where it is, so the
  reader multiplies the rate by air density once per occurrence);
  products -> ordered list of {"species", "yield"} (M kept, as in TS1).
- ARRHENIUS -> A, B, C, D, E; TROE -> k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C,
  Fc, N; USER_DEFINED -> "scaling factor" and "MUSICA name" (the label, which is
  also the custom-parameter column name in cells.csv). Emissions are
  USER_DEFINED reactions with an empty reactant list.
- "initial state": the given concentrations, then the rate parameters, then
  temperature and pressure. Air density is not written: the entry computes it
  from the ideal gas law, as MICM's driver does.
"""
import json
import sys

PARAMETERS = {
    "ARRHENIUS": ("A", "B", "C", "D", "E"),
    "TROE": ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N"),
}


def convert(canonical):
    third_bodies = {s["name"] for s in canonical["species"] if s.get("third_body")}
    data = []
    for s in canonical["species"]:
        record = {"name": s["name"], "type": "CHEM_SPEC"}
        if s["name"] in third_bodies:
            record["parameterization"] = "air density"
        data.append(record)
    reactions = []
    for r in canonical["reactions"]:
        kind = r["type"]
        out = {"type": kind,
               "reactants": list(r["reactants"]),
               "products": [{"species": name, "yield": float(y)} for name, y in r["products"]]}
        if kind in PARAMETERS:
            for key in PARAMETERS[kind]:
                out[key] = r[key]
        elif kind == "USER_DEFINED":
            out["scaling factor"] = r["scaling_factor"]
            out["MUSICA name"] = r["label"]
        else:
            raise SystemExit(f"convert.py: no mapping for reaction type {kind}")
        reactions.append(out)
    initial = canonical["initial"]
    state = {}
    for name, value in initial["concentrations"].items():
        if name in third_bodies:
            continue  # parameterized: the reader rejects a value for it
        state[name] = value
    state.update(initial["parameters"])
    state["temperature"] = initial["temperature"]
    state["pressure"] = initial["pressure"]
    data.append({"name": canonical["name"], "type": "MECHANISM", "rate units": "mol m-3",
                 "reactions": reactions, "initial state": state})
    return {"camp-data": data}


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    with open(sys.argv[1]) as f:
        canonical = json.load(f)
    with open(sys.argv[2], "w") as f:
        json.dump(convert(canonical), f, indent=2)
        f.write("\n")
