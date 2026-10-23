#!/usr/bin/env python3
"""Convert a canonical mechanism (canonical.json) to luna-high's input format.

    convert.py CANONICAL_JSON OUTPUT_JSON

The output mirrors the entry's mechanisms/ts1.json, which its
tools/import_ts1.py writes from MICM's TS1 header ("c-mechanism-1"):

- species: {name, third_body, molecular_weight, diffusion}, in source order;
  M keeps "third_body": true, as TS1's M does (molecular weight and diffusion
  matter only to surface reactions, so they are 0.0, as for most TS1 species);
- reactions: lower-case type; reactants as one {name, coefficient: 1.0} per
  listed occurrence (repeats kept, M included, as in TS1); products as
  {name, coefficient: yield}; parameters keyed as in TS1, with
  "scaling_factor" on every reaction (1.0 unless the source gives one):
    ARRHENIUS    -> "arrhenius"     A, B, C, D, E
    TROE         -> "troe"          k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N
    USER_DEFINED -> "user_defined"  label, scaling_factor (photolysis and
                                    emissions; an emission has no reactants)
- defaults: every non-third-body species' initial concentration (0.0 where
  the source gives none, as TS1's InitState lists every species), every
  rate-parameter label's initial value, temperature, and pressure.
"""
import argparse
import json
from pathlib import Path

PARAMETER_KEYS = {
    "ARRHENIUS": ("arrhenius", ("A", "B", "C", "D", "E")),
    "TROE": ("troe", ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")),
}


def convert(canonical):
    species = canonical["species"]
    names = [s["name"] for s in species]
    if len(set(names)) != len(names):
        raise SystemExit("duplicate species names")
    third_body = {s["name"] for s in species if s.get("third_body", False)}
    out_species = [{"name": s["name"], "third_body": s["name"] in third_body,
                    "molecular_weight": 0.0, "diffusion": 0.0} for s in species]

    reactions, labels = [], []
    for i, reaction in enumerate(canonical["reactions"]):
        kind = reaction["type"]
        for name in reaction["reactants"] + [p[0] for p in reaction["products"]]:
            if name not in names:
                raise SystemExit(f"reaction {i}: unknown species {name!r}")
        if kind in PARAMETER_KEYS:
            out_type, keys = PARAMETER_KEYS[kind]
            missing = [k for k in keys if k not in reaction]
            if missing:
                raise SystemExit(f"reaction {i} ({kind}): missing {missing}")
            parameters = {k: float(reaction[k]) for k in keys}
            parameters["scaling_factor"] = float(reaction.get("scaling_factor", 1.0))
        elif kind == "USER_DEFINED":
            out_type = "user_defined"
            parameters = {"label": reaction["label"],
                          "scaling_factor": float(reaction.get("scaling_factor", 1.0))}
            if reaction["label"] not in labels:
                labels.append(reaction["label"])
        else:
            raise SystemExit(f"reaction {i}: no TS1 equivalent for type {kind!r}")
        reactions.append({
            "type": out_type,
            "reactants": [{"name": n, "coefficient": 1.0} for n in reaction["reactants"]],
            "products": [{"name": n, "coefficient": float(y)} for n, y in reaction["products"]],
            "parameters": parameters,
        })

    declared = canonical.get("parameters", labels)
    if sorted(declared) != sorted(labels):
        raise SystemExit(f"declared parameters {declared} differ from reaction labels {labels}")
    initial = canonical["initial"]
    concentrations = initial.get("concentrations", {})
    for name in concentrations:
        if name not in names or name in third_body:
            raise SystemExit(f"initial concentration for {name!r}, not a state species")
    values = initial.get("parameters", {})
    defaults = {
        "species": {n: float(concentrations.get(n, 0.0)) for n in names if n not in third_body},
        "rate_parameters": {label: float(values.get(label, 0.0)) for label in declared},
        "temperature": float(initial["temperature"]),
        "pressure": float(initial["pressure"]),
    }
    return {"format": "c-mechanism-1", "name": canonical["name"], "species": out_species,
            "reactions": reactions, "defaults": defaults}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("canonical", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    data = convert(json.loads(args.canonical.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
