#!/usr/bin/env python3
"""Translate the supplied generated TS1 builder into the public JSON schema.

The input is read only; the resulting JSON is committed, so normal builds do not
depend on MICM sources being installed.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "micm/test/benchmark/ts1_mechanism.hpp"


def balanced(text: str, opening: int, left: str, right: str) -> tuple[str, int]:
    depth = 0
    quote = ""
    escaped = False
    for pos in range(opening, len(text)):
        char = text[pos]
        if quote:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = ""
            continue
        if char in ('"', "'"):
            quote = char
        elif char == left:
            depth += 1
        elif char == right:
            depth -= 1
            if depth == 0:
                return text[opening + 1 : pos], pos + 1
    raise ValueError(f"unclosed {left} at offset {opening}")


def value(text: str, key: str, default: float = 0.0) -> float:
    match = re.search(rf"\.{re.escape(key)}\s*=\s*([^,}}\s]+)", text)
    return float(match.group(1)) if match else default


def string_value(text: str, key: str) -> str | None:
    match = re.search(rf"\.{re.escape(key)}\s*=\s*\"([^\"]*)\"", text)
    return match.group(1) if match else None


def parse_species(source: str) -> list[dict]:
    start = source.index("inline micm::Phase CreateGasPhase()")
    end = source.index("return micm::Phase", start)
    body = source[start:end]
    species = []
    for match in re.finditer(r"auto species = micm::Species\(\"([^\"]+)\"\);", body):
        tail = body[match.end() :]
        block_end = tail.find("\n      }")
        block = tail[:block_end] if block_end >= 0 else tail
        mw_match = re.search(r"MOLECULAR_WEIGHT,\s*micm::Real\{\s*([^}]+)\}", block)
        diff_match = re.search(r"phase_species\.emplace_back\(species,\s*([-+0-9.eE]+)\)", block)
        third_body = "SetThirdBody()" in block
        species.append(
            {
                "name": match.group(1),
                "molecular_weight": float(mw_match.group(1)) if mw_match else 0.0,
                "diffusion_coefficient": float(diff_match.group(1)) if diff_match else None,
                "third_body": third_body,
            }
        )
    if len(species) != 210:
        raise ValueError(f"expected 210 TS1 species, found {len(species)}")
    return species


def parse_reactions(source: str) -> list[dict]:
    start = source.index("inline std::vector<micm::Process> CreateProcesses")
    end = source.index("return processes;", start)
    body = source[start:end]
    reactions = []
    pos = 0
    while True:
        match = re.search(r"processes\.push_back\s*\(", body[pos:])
        if not match:
            break
        opening = pos + match.end() - 1
        block, pos = balanced(body, opening, "(", ")")

        def call_body(name: str) -> str:
            found = re.search(rf"\.{name}\s*\(", block)
            if not found:
                return ""
            opening = found.end() - 1
            return balanced(block, opening, "(", ")")[0]

        reactants = re.findall(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)', call_body("SetReactants"))
        product_text = call_body("SetProducts")
        products = []
        for p in re.finditer(r'SpeciesByName\(gas_phase,\s*"([^"]+)"\)', product_text):
            rest = product_text[p.end() :]
            coeff = re.match(r"\s*,\s*([-+0-9.eE]+)\s*\}", rest)
            products.append({"species": p.group(1), "yield": float(coeff.group(1)) if coeff else 1.0})

        rate_match = re.search(r"micm::([A-Za-z0-9_]+Parameters)\s*\{", block)
        if not rate_match:
            raise ValueError(f"could not identify rate type in reaction {len(reactions)}")
        kind = rate_match.group(1).removesuffix("Parameters")
        rate_open = block.index("{", rate_match.start())
        rate, _ = balanced(block, rate_open, "{", "}")
        label = string_value(rate, "label_")

        if kind == "ArrheniusRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("A", 1), ("B", 0), ("C", 0), ("D", 300), ("E", 0))}
        elif kind == "TroeRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("k0_A", 1), ("k0_B", 0), ("k0_C", 0), ("kinf_A", 1), ("kinf_B", 0), ("kinf_C", 0), ("Fc", .6), ("N", 1))}
        elif kind == "TernaryChemicalActivationRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("k0_A", 1), ("k0_B", 0), ("k0_C", 0), ("kinf_A", 1), ("kinf_B", 0), ("kinf_C", 0), ("Fc", .6), ("N", 1))}
        elif kind == "TernaryChemicalActivationJPL19":
            params = {k: value(rate, f"{k}_", default) for k, default in (("k0_A", 1), ("k0_B", 0), ("k0_C", 0), ("k0_D", 298), ("kinf_A", 1), ("kinf_B", 0), ("kinf_C", 0), ("kinf_D", 298), ("kint_A", 1), ("kint_B", 0), ("kint_C", 0), ("kint_D", 298), ("Fc", .6), ("N", 1))}
        elif kind == "BranchedRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("X", 0), ("Y", 0), ("a0", 0), ("n", 0))}
            params["branch"] = "alkoxy" if "Branch::Alkoxy" in rate else "nitrate"
        elif kind == "UserDefinedRateConstant":
            params = {"scaling_factor": value(rate, "scaling_factor_", 1)}
        elif kind == "SurfaceRateConstant":
            params = {
                "reaction_probability": value(rate, "reaction_probability_", 1),
                "surface_species": string_value(rate, "phase_species_") or reactants[0],
            }
        elif kind == "TunnelingRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("A", 1), ("B", 0), ("C", 0))}
        elif kind == "TaylorSeriesRateConstant":
            coeff_match = re.search(r"coefficients_\s*=\s*\{([^}]*)\}", rate)
            params = {k: value(rate, f"{k}_", default) for k, default in (("A", 1), ("B", 0), ("C", 0), ("D", 300), ("E", 0))}
            params["coefficients"] = [float(x.strip()) for x in coeff_match.group(1).split(",") if x.strip()] if coeff_match else []
        elif kind == "ReversibleRateConstant":
            params = {k: value(rate, f"{k}_", default) for k, default in (("A", 1), ("C", 0), ("k_r", 0))}
        else:
            raise ValueError(f"unsupported TS1 rate type: {kind}")

        reactions.append({"type": kind, "reactants": reactants, "products": products, "label": label, "parameters": params})
    if len(reactions) != 547:
        raise ValueError(f"expected 547 TS1 reactions, found {len(reactions)}")
    return reactions


def parse_initial_state(source: str) -> tuple[dict, dict, float, float]:
    start = source.index("static void InitState")
    body = source[start:]
    concentrations = {name: float(number) for name, number in re.findall(r'concentration\("([^\"]+)",\s*([-+0-9.eE]+)\)', body)}
    parameters = {name: float(number) for name, number in re.findall(r'parameter\("([^\"]+)",\s*([-+0-9.eE]+)\)', body)}
    temp = float(re.search(r"temperature_\s*=\s*([-+0-9.eE]+)", body).group(1))
    pressure = float(re.search(r"pressure_\s*=\s*([-+0-9.eE]+)", body).group(1))
    return concentrations, parameters, temp, pressure


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source = args.source.read_text()
    initial, custom, temperature, pressure = parse_initial_state(source)
    species = parse_species(source)
    for item in species:
        item["initial"] = initial.get(item["name"], 0.0)
    mechanism = {
        "name": "TS1",
        "temperature": temperature,
        "pressure": pressure,
        "species": species,
        "parameters": custom,
        "reactions": parse_reactions(source),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(mechanism, indent=2) + "\n")
    print(f"wrote {args.output}: {len(species)} species, {len(mechanism['reactions'])} reactions, {len(custom)} parameters")


if __name__ == "__main__":
    main()
