#!/usr/bin/env python3
"""Translate the generated TS1 C++ data into the configuration dialect.

This is a data extractor, not a C++ evaluator. It deliberately rejects unknown
constructs and checks the number of extracted species and reactions. Run it
explicitly when updating the source mechanism; ordinary builds use the checked-in
configuration and do not need the source tree.
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def only(pattern, text, what):
    found = re.findall(pattern, text, re.S)
    if len(found) != 1:
        raise ValueError(f"expected one {what}, found {len(found)}")
    return found[0]


def extract(source):
    phase = source.split("inline micm::Phase CreateGasPhase()", 1)[1].split(
        "inline std::vector<micm::Process> CreateProcesses", 1
    )[0]
    species = []
    for body in re.findall(r'\{\s*auto species = micm::Species\("([^"\n]+)"\);(.*?)\n\s*\}', phase, re.S):
        name, content = body
        entry = {"name": name, "type": "CHEM_SPEC"}
        if "SetThirdBody()" in content:
            entry["parameterization"] = "air density"
        else:
            entry["molecular weight [kg mol-1]"] = float(
                only(r"MOLECULAR_WEIGHT, micm::Real\{\s*(" + NUMBER + r")\s*\}", content, name + " molecular weight")
            )
        diffusion = re.findall(r"phase_species\.emplace_back\(species,\s*(" + NUMBER + r")\)", content)
        if diffusion:
            entry["diffusion coefficient [m2 s-1]"] = float(diffusion[0])
        species.append(entry)
    if len(species) != 210 or len({s["name"] for s in species}) != len(species):
        raise ValueError(f"expected 210 distinct species, found {len(species)}")

    reaction_text = source.split("inline std::vector<micm::Process> CreateProcesses", 1)[1].split(
        "struct Ts1", 1
    )[0]
    chunks = reaction_text.split("processes.push_back(")[1:]
    reactions = []
    type_names = {"Arrhenius": "ARRHENIUS", "Troe": "TROE", "Surface": "SURFACE", "UserDefined": "USER_DEFINED"}
    for index, chunk in enumerate(chunks):
        chunk = chunk.split(".Build());", 1)[0]
        reactants_text = only(r"\.SetReactants\(\s*\{(.*?)\}\)", chunk, f"reaction {index} reactants")
        reactant_names = re.findall(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)', reactants_text)
        products_text = only(r"\.SetProducts\(\s*\{(.*?)\}\)", chunk, f"reaction {index} products")
        products = re.findall(r'StoichSpecies\{\s*SpeciesByName\(gas_phase,\s*"([^"]+)"\),\s*(' + NUMBER + r')\s*\}', products_text)
        if len({name for name, _ in products}) != len(products):
            raise ValueError(f"reaction {index} repeats a product")
        typ, fields = only(r"\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\s*\{(.*?)\}\)", chunk, f"reaction {index} rate")
        if typ not in type_names:
            raise ValueError(f"unknown rate type {typ}")
        counts = Counter(reactant_names)
        reaction = {
            "type": type_names[typ],
            "reactants": {name: ({"qty": count} if count != 1 else {}) for name, count in counts.items()},
            "products": {name: {"yield": float(value)} for name, value in products},
        }
        for key, value in re.findall(r'\.(\w+)_\s*=\s*("[^"]*"|' + NUMBER + r')', fields):
            reaction[{"label": "MUSICA name", "scaling_factor": "scaling factor", "reaction_probability": "reaction probability"}.get(key, key)] = (
                value[1:-1] if value.startswith('"') else float(value)
            )
        if typ == "Surface":
            reaction["gas-phase reactant"] = only(r'PhaseSpeciesByName\(gas_phase,\s*"([^"]+)"\)', fields, "surface species")
        reactions.append(reaction)
    if len(reactions) != 547:
        raise ValueError(f"expected 547 reactions, found {len(reactions)}")

    defaults = {name: float(value) for name, value in re.findall(r'concentration\("([^"]+)",\s*(' + NUMBER + r')\)', source)}
    for entry in species:
        entry["initial concentration [mol m-3]"] = defaults.get(entry["name"], 0.0)
    parameters = {name: float(value) for name, value in re.findall(r'parameter\("([^"]+)",\s*(' + NUMBER + r')\)', source)}
    temperature = float(only(r"conditions_\[c\]\.temperature_\s*=\s*(" + NUMBER + r")", source, "temperature"))
    pressure = float(only(r"conditions_\[c\]\.pressure_\s*=\s*(" + NUMBER + r")", source, "pressure"))
    return {"camp-data": species + [{"name": "TS1", "type": "MECHANISM", "reactions": reactions}],
            "initial conditions": {"temperature": temperature, "pressure": pressure, "rate parameters": parameters}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = extract(args.source.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"wrote {args.output}: {len(result['camp-data']) - 1} species, 547 reactions")


if __name__ == "__main__":
    main()
