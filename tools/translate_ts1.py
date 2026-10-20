#!/usr/bin/env python3
"""Translate the read-only MICM TS1 mechanism table into the project JSON format.

The script parses only the ordinary species, reaction, and InitState expressions
in MICM's generated TS1 header. It does not import the MICM or MUSICA libraries.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"translate_ts1: {message}")


def parse_species(source: str) -> list[dict]:
    start = source.index("inline micm::Phase CreateGasPhase()")
    end = source.index("inline std::vector<micm::Process> CreateProcesses", start)
    block = source[start:end]
    species = []
    pattern = re.compile(
        r'auto species = micm::Species\("([^"\\]*(?:\\.[^"\\]*)*)"\);'
        r'(.*?)phase_species\.emplace_back\(species(?:,\s*([^)]*))?\);',
        re.S,
    )
    for match in pattern.finditer(block):
        name, body, diffusion = match.groups()
        mw_match = re.search(r"MOLECULAR_WEIGHT,\s*micm::Real\{\s*([^}]+)\}", body)
        mw = float(mw_match.group(1)) if mw_match else 0.0
        species.append(
            {
                "name": name,
                "molecular_weight": mw,
                "diffusion_coefficient": float(diffusion) if diffusion else None,
                "third_body": "SetThirdBody()" in body,
                "initial": 0.0,
            }
        )
    if len(species) != 210:
        fail(f"expected 210 TS1 species, parsed {len(species)}")
    return species


def parse_reactions(source: str) -> list[dict]:
    start = source.index("inline std::vector<micm::Process> CreateProcesses")
    end = source.index("static void InitState", start)
    block = source[start:end]
    species_rx = re.compile(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)')
    product_rx = re.compile(
        r"micm::StoichSpecies\s*\{\s*SpeciesByName\(gas_phase,\s*\"([^\"]+)\"\),\s*([^}]+)\}"
    )
    statement_rx = re.compile(
        r"processes\.push_back\(\s*micm::ChemicalReactionBuilder\(\)(.*?)\.Build\(\)\);", re.S
    )
    reactions = []
    for statement in statement_rx.findall(block):
        reactant_match = re.search(r"\.SetReactants\(\s*\{(.*?)\}\s*\)\s*\.SetProducts", statement, re.S)
        product_match = re.search(r"\.SetProducts\(\s*\{(.*?)\}\s*\)\s*\.SetRateConstant", statement, re.S)
        rate_match = re.search(r"\.SetRateConstant\(\s*micm::(\w+)\s*\{(.*?)\}\s*\)", statement, re.S)
        if not reactant_match or not product_match or not rate_match:
            fail("could not parse a process builder expression")
        reactants = species_rx.findall(reactant_match.group(1))
        products = [
            {"name": name, "yield": float(value.strip())}
            for name, value in product_rx.findall(product_match.group(1))
        ]
        cxx_type, rate_body = rate_match.groups()
        if cxx_type == "ArrheniusRateConstantParameters":
            kind = "ARRHENIUS"
        elif cxx_type == "TroeRateConstantParameters":
            kind = "TROE"
        elif cxx_type == "SurfaceRateConstantParameters":
            kind = "SURFACE"
        elif cxx_type == "UserDefinedRateConstantParameters":
            kind = "USER_DEFINED"
        else:
            fail(f"unsupported TS1 rate type {cxx_type}")
        params = {}
        for name, value in re.findall(r"\.(\w+)_\s*=\s*([^,\s}]+)", rate_body):
            try:
                params[name] = float(value)
            except ValueError:
                pass
        label_match = re.search(r'\.label_\s*=\s*"([^"]+)"', rate_body)
        reaction = {"type": kind, "reactants": reactants, "products": products, "parameters": params}
        if label_match:
            reaction["label"] = label_match.group(1)
        surface_species = re.search(
            r'\.phase_species_\s*=\s*PhaseSpeciesByName\(gas_phase,\s*"([^"]+)"\)', rate_body
        )
        if surface_species:
            reaction["surface_species"] = surface_species.group(1)
        if "probability" in params:
            # C++ member spelling maps to the runtime's surface probability slot.
            reaction["parameters"]["reaction_probability"] = params["probability"]
        reactions.append(reaction)
    if len(reactions) != 547:
        fail(f"expected 547 TS1 reactions, parsed {len(reactions)}")
    return reactions


def parse_initial_state(source: str, species: list[dict]) -> tuple[float, float, dict[str, float]]:
    start = source.index("static void InitState")
    block = source[start:]
    for name, value in re.findall(r'concentration\("([^"]+)",\s*([^)]*)\);', block):
        by_name = next((item for item in species if item["name"] == name), None)
        if by_name is not None:
            by_name["initial"] = float(value)
    custom = {
        name: float(value)
        for name, value in re.findall(r'parameter\("([^"]+)",\s*([^)]*)\);', block)
    }
    temperature = re.search(r"conditions_\[c\]\.temperature_\s*=\s*([^;]+);", block)
    pressure = re.search(r"conditions_\[c\]\.pressure_\s*=\s*([^;]+);", block)
    if not temperature or not pressure:
        fail("could not parse TS1 temperature and pressure defaults")
    return float(temperature.group(1)), float(pressure.group(1)), custom


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("../micm/test/benchmark/ts1_mechanism.hpp"),
        help="MICM's read-only generated TS1 mechanism header",
    )
    parser.add_argument("--output", type=Path, default=Path("mechanisms/ts1.json"))
    args = parser.parse_args()
    source = args.source.read_text()
    species = parse_species(source)
    reactions = parse_reactions(source)
    temperature, pressure, custom = parse_initial_state(source, species)
    document = {
        "name": "TS1",
        "temperature": temperature,
        "pressure": pressure,
        "species": species,
        "reactions": reactions,
        "custom_parameters": custom,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    print(f"wrote {len(species)} species, {len(reactions)} reactions to {args.output}")


if __name__ == "__main__":
    main()
