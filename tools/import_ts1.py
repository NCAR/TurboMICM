#!/usr/bin/env python3
"""Translate MICM's declarative TS1 builder calls into portable JSON."""
import argparse
import json
import re
from pathlib import Path


def blocks(source, marker):
    start = 0
    while True:
        pos = source.find(marker, start)
        if pos < 0:
            return
        depth = 0
        quote = False
        escape = False
        for i in range(pos, len(source)):
            c = source[i]
            if quote:
                if escape:
                    escape = False
                elif c == "\\":
                    escape = True
                elif c == '"':
                    quote = False
                continue
            if c == '"':
                quote = True
            elif c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    yield source[pos:i + 1]
                    start = i + 1
                    break


def parse_species(source):
    phase = source.split("inline micm::Phase CreateGasPhase()", 1)[1].split(
        "inline std::vector<micm::Process> CreateProcesses", 1)[0]
    result = []
    pat = re.compile(r'auto species = micm::Species\("([^"]+)"\);(.*?)(?=\n\s*\})', re.S)
    for name, body in pat.findall(phase):
        mw = re.search(r'MOLECULAR_WEIGHT, micm::Real\{\s*([^}]+)', body)
        diffusion = re.search(r'phase_species\.emplace_back\(species,\s*([^\)]+)', body)
        result.append({"name": name, "third_body": "SetThirdBody" in body,
                       "molecular_weight": float(mw.group(1)) if mw else 0.0,
                       "diffusion": float(diffusion.group(1)) if diffusion else 0.0})
    if len(result) != 210:
        raise ValueError(f"expected 210 species, parsed {len(result)}")
    return result


def parse_terms(text, products=False):
    if products:
        return [{"name": n, "coefficient": float(c)} for n, c in
                re.findall(r'StoichSpecies\s*\{?\s*SpeciesByName\(gas_phase,\s*"([^"]+)"\),\s*([^}]+)\}', text)]
    return [{"name": n, "coefficient": 1.0} for n in
            re.findall(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)', text)]


def parse_reactions(source):
    body = source.split("inline std::vector<micm::Process> CreateProcesses", 1)[1].split(
        "struct Ts1", 1)[0]
    result = []
    for block in blocks(body, "processes.push_back("):
        rx = re.search(r'\.SetReactants\(\s*\{(.*?)\}\s*\)', block, re.S)
        prod = re.search(r'\.SetProducts\(\s*\{(.*?)\}\s*\)', block, re.S)
        rate = re.search(r'\.SetRateConstant\(\s*micm::(\w+)RateConstantParameters\{(.*?)\}\)', block, re.S)
        if not rate or not rx or not prod:
            raise ValueError("unrecognized reaction builder")
        kind = rate.group(1).upper()
        key = {"ARRHENIUS": "arrhenius", "TROE": "troe", "USERDEFINED": "user_defined",
               "SURFACE": "surface"}.get(kind)
        if key is None:
            raise ValueError(f"unsupported TS1 type: {kind}")
        params = {k: float(v) for k, v in re.findall(r'\.(\w+)_\s*=\s*([-+0-9.eE]+)', rate.group(2))}
        label = re.search(r'\.label_\s*=\s*"([^"]+)"', rate.group(2))
        if label:
            params["label"] = label.group(1)
        params["scaling_factor"] = float(params.pop("scaling_factor", 1.0))
        result.append({"type": key,
                       "reactants": parse_terms(rx.group(1)),
                       "products": parse_terms(prod.group(1), True),
                       "parameters": params})
    if len(result) != 547:
        raise ValueError(f"expected 547 reactions, parsed {len(result)}")
    return result


def parse_defaults(source):
    body = source.split("static void InitState", 1)[1].split("\n    }\n  };", 1)[0]
    concentrations = {name: float(value) for name, value in
                      re.findall(r'concentration\("([^"]+)",\s*([-+0-9.eE]+)\);', body)}
    parameters = {name: float(value) for name, value in
                  re.findall(r'parameter\("([^"]+)",\s*([-+0-9.eE]+)\);', body)}
    temperature = re.search(r'temperature_\s*=\s*([-+0-9.eE]+)', body)
    pressure = re.search(r'pressure_\s*=\s*([-+0-9.eE]+)', body)
    if not temperature or not pressure:
        raise ValueError("could not parse TS1 default temperature and pressure")
    return {"species": concentrations, "rate_parameters": parameters,
            "temperature": float(temperature.group(1)), "pressure": float(pressure.group(1))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path, help="MICM test/benchmark/ts1_mechanism.hpp")
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    source = args.source.read_text()
    data = {"format": "c-mechanism-1", "name": "TS1", "species": parse_species(source),
            "reactions": parse_reactions(source), "defaults": parse_defaults(source)}
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
