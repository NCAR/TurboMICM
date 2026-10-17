#!/usr/bin/env python3
"""Translate MICM's TS1 benchmark mechanism (test/benchmark/ts1_mechanism.hpp,
generated C++) into the configuration format read by tools/mechanism.py.

    translate_ts1.py ts1_mechanism.hpp configs/ts1.json

The header is text-matched, never compiled. Every species, reaction, initial
concentration, and rate parameter must be recognized, or the translation stops.
Rate parameters in the header are already in MICM's internal units (mol m-3),
so the mechanism carries "rate constant units": "mol m-3" and no conversion is
applied when it is read.
"""
import json
import re
import sys
from pathlib import Path

RATE_TYPES = {"Arrhenius": "ARRHENIUS", "Troe": "TROE", "TernaryChemicalActivation": "TERNARY_CHEMICAL_ACTIVATION",
              "Tunneling": "TUNNELING", "UserDefined": "USER_DEFINED", "Surface": "SURFACE",
              "Branched": "BRANCHED", "TaylorSeries": "TAYLOR_SERIES", "Reversible": "REVERSIBLE"}
FIELD_NAMES = {"label": "MUSICA name", "scaling_factor": "scaling factor", "reaction_probability": "reaction probability"}


def number(text):
    text = text.strip().rstrip("f")
    return float(text)


def translate(header: Path):
    text = header.read_text()
    species = []
    species_re = re.compile(r'auto species = micm::Species\("([^"]+)"\);(.*?)phase_species\.emplace_back\(species(?:,\s*([^)]+))?\);', re.S)
    for name, body, diffusion in species_re.findall(text):
        entry = {"name": name, "type": "CHEM_SPEC"}
        if "SetThirdBody()" in body:
            entry["third body"] = True
        weight = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*([^}]+)\}', body)
        if weight:
            entry["molecular weight [kg mol-1]"] = number(weight[1])
        if diffusion:
            entry["diffusion coefficient [m2 s-1]"] = number(diffusion)
        species.append(entry)
    phase = re.search(r'return micm::Phase\{\s*"([^"]+)"', text)
    phase_name = phase[1] if phase else "gas"

    reactions = []
    for block in text.split("processes.push_back(")[1:]:
        block = block.split(".Build());", 1)[0]
        reactants = re.search(r'\.SetReactants\(\s*\{(.*?)\}\s*\)', block, re.S)
        products = re.search(r'\.SetProducts\(\s*\{(.*?)\}\s*\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\s*\)\s*\.SetPhase', block, re.S)
        if not (reactants and rate and rate[1] in RATE_TYPES):
            raise SystemExit("unrecognized reaction declaration:\n" + block[:300])
        rx = {"type": RATE_TYPES[rate[1]]}
        rx["reactants"] = re.findall(r'SpeciesByName\(gas_phase, "([^"]+)"\)', reactants[1])
        rx["products"] = []
        if products:
            for name, yield_ in re.findall(r'StoichSpecies\{\s*SpeciesByName\(gas_phase, "([^"]+)"\),\s*([^}]+)\}', products[1]):
                rx["products"].append({"species": name, "yield": number(yield_)})
        for key, value in re.findall(r'\.(\w+)_\s*=\s*("[^"]*"|[-+0-9.eE]+)', rate[2]):
            if key == "phase_species":
                continue
            key = FIELD_NAMES.get(key, key)
            rx[key] = json.loads(value) if value.startswith('"') else number(value)
        if rx["type"] == "SURFACE":
            ps = re.search(r'PhaseSpeciesByName\(gas_phase, "([^"]+)"\)', rate[2])
            rx["gas-phase reactant"] = ps[1] if ps else rx["reactants"][0]
            rx["gas-phase products"] = rx.pop("products")
        reactions.append(rx)

    initial = {}
    for name, value in re.findall(r'concentration\("([^"]+)",\s*([^)]+)\);', text):
        initial[name] = number(value)
    for entry in species:
        if entry["name"] in initial:
            entry["initial concentration [mol m-3]"] = initial.pop(entry["name"])
    if initial:
        raise SystemExit(f"initial concentrations for unknown species: {sorted(initial)}")
    parameters = {name: number(value) for name, value in re.findall(r'parameter\("([^"]+)",\s*([^)]+)\);', text)}
    conditions = {}
    for field in ("temperature", "pressure"):
        found = re.findall(r'conditions_\[c\]\.' + field + r'_\s*=\s*([^;]+);', text)
        if len(found) != 1:
            raise SystemExit(f"expected one default {field}")
        conditions[field] = number(found[0])
    expected = re.search(r'// (\d+) species and (\d+) reactions', text)
    if expected and (int(expected[1]) != len(species) or int(expected[2]) != len(reactions)):
        raise SystemExit(f"found {len(species)} species and {len(reactions)} reactions, the header says {expected[1]} and {expected[2]}")
    mechanism = {"name": "TS1", "type": "MECHANISM", "phase": phase_name, "rate constant units": "mol m-3", "reactions": reactions}
    return {"camp-data": species + [mechanism], "initial conditions": conditions, "initial rate parameters": parameters}


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    result = translate(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=1) + "\n")
    n_species = len(result["camp-data"]) - 1
    print(f"wrote {sys.argv[2]}: {n_species} species, {len(result['camp-data'][-1]['reactions'])} reactions, "
          f"{len(result['initial rate parameters'])} rate parameters")


if __name__ == "__main__":
    main()
