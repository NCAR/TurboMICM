#!/usr/bin/env python3
"""Read a MICM-format YAML/JSON mechanism and emit a C descriptor."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import yaml


RATE_ORDER = {
    "ArrheniusRateConstant": 0,
    "TroeRateConstant": 1,
    "TernaryChemicalActivationRateConstant": 2,
    "TernaryChemicalActivationJPL19": 3,
    "BranchedRateConstant": 4,
    "TunnelingRateConstant": 5,
    "TaylorSeriesRateConstant": 6,
    "ReversibleRateConstant": 7,
    "UserDefinedRateConstant": 8,
    "SurfaceRateConstant": 9,
}
RATE_ENUM = {
    "ArrheniusRateConstant": "CHEM_RATE_ARRHENIUS",
    "TroeRateConstant": "CHEM_RATE_TROE",
    "TernaryChemicalActivationRateConstant": "CHEM_RATE_TERNARY",
    "TernaryChemicalActivationJPL19": "CHEM_RATE_TERNARY_JPL19",
    "BranchedRateConstant": "CHEM_RATE_BRANCHED",
    "TunnelingRateConstant": "CHEM_RATE_TUNNELING",
    "TaylorSeriesRateConstant": "CHEM_RATE_TAYLOR",
    "ReversibleRateConstant": "CHEM_RATE_REVERSIBLE",
    "UserDefinedRateConstant": "CHEM_RATE_USER",
    "SurfaceRateConstant": "CHEM_RATE_SURFACE",
}
PARAM_INDEX = {
    "A": 0, "B": 1, "C": 2, "D": 3, "E": 4,
    "k0_A": 5, "k0_B": 6, "k0_C": 7,
    "kinf_A": 8, "kinf_B": 9, "kinf_C": 10,
    "Fc": 11, "N": 12, "X": 13, "Y": 14,
    "a0": 15, "n": 16, "k_r": 17,
    "k0_D": 18, "kinf_D": 19, "kint_A": 20,
    "kint_B": 21, "kint_C": 22, "kint_D": 23,
}


def load(path: Path):
    with path.open() as stream:
        return yaml.safe_load(stream)


def camp_records(config_path: Path) -> tuple[list[dict], list[dict]]:
    root = load(config_path)
    records = []
    for filename in root.get("camp-files", []):
        child = load(config_path.parent / filename)
        records.extend(child.get("camp-data", []))
    records.extend(root.get("camp-data", []))
    species = [r for r in records if r.get("type") == "CHEM_SPEC"]
    mechanisms = [r for r in records if r.get("type") == "MECHANISM"]
    if not species or not mechanisms:
        raise ValueError("expected camp-files/camp-data with CHEM_SPEC and MECHANISM records")
    return species, [rxn for mech in mechanisms for rxn in mech.get("reactions", [])]


def term_map(mapping, coefficient_key: str, default: float) -> list[dict]:
    out = []
    for name, options in (mapping or {}).items():
        out.append({"species": name, "yield": float(options.get(coefficient_key, default) if options else default)})
    return out


def from_camp(config_path: Path) -> dict:
    raw_species, raw_reactions = camp_records(config_path)
    species = []
    for item in raw_species:
        species.append({
            "name": item["name"],
            "molecular_weight": float(item.get("molecular weight [kg mol-1]", 0.0)),
            "diffusion_coefficient": item.get("diffusion coefficient [m2 s-1]"),
            "initial": 0.0,
            "third_body": item.get("parameterized") == "air_density",
        })
    type_map = {
        "ARRHENIUS": "ArrheniusRateConstant",
        "BRANCHED": "BranchedRateConstant",
        "SURFACE": "SurfaceRateConstant",
        "TERNARY_CHEMICAL_ACTIVATION": "TernaryChemicalActivationRateConstant",
        "TERNARY_CHEMICAL_ACTIVATION_JPL19": "TernaryChemicalActivationJPL19",
        "TROE": "TroeRateConstant",
        "TUNNELING": "TunnelingRateConstant",
        "TAYLOR_SERIES": "TaylorSeriesRateConstant",
        "REVERSIBLE": "ReversibleRateConstant",
        "PHOTOLYSIS": "UserDefinedRateConstant",
        "USER_DEFINED": "UserDefinedRateConstant",
        "FIRST_ORDER_LOSS": "UserDefinedRateConstant",
        "EMISSION": "UserDefinedRateConstant",
    }
    reactions = []
    for item in raw_reactions:
        raw_type = item["type"]
        if raw_type not in type_map:
            raise ValueError(f"unsupported reaction type {raw_type!r}")
        kind = type_map[raw_type]
        if raw_type == "SURFACE":
            reactants = [item["gas-phase reactant"]]
            products = term_map(item.get("gas-phase products"), "yield", 1.0)
            surface_species = item["gas-phase reactant"]
        elif raw_type == "FIRST_ORDER_LOSS":
            reactants, products = [item["species"]], []
            surface_species = None
        elif raw_type == "EMISSION":
            reactants, products = [], [{"species": item["species"], "yield": 1.0}]
            surface_species = None
        elif kind == "BranchedRateConstant":
            reactants = [{"species": name, "quantity": options.get("qty", 1.0)}
                         for name, options in item.get("reactants", {}).items()]
            products = []
            surface_species = None
        else:
            reactants = [{"species": name, "quantity": options.get("qty", 1.0)}
                         for name, options in item.get("reactants", {}).items()]
            products = term_map(item.get("products"), "yield", 1.0)
            surface_species = None
        fields = {
            "A": item.get("A", 1), "B": item.get("B", 0), "C": item.get("C", 0),
            "D": item.get("D", 300), "E": item.get("E", 0),
            "k0_A": item.get("k0_A", 1), "k0_B": item.get("k0_B", 0), "k0_C": item.get("k0_C", 0),
            "kinf_A": item.get("kinf_A", 1), "kinf_B": item.get("kinf_B", 0), "kinf_C": item.get("kinf_C", 0),
            "k0_D": item.get("k0_D", 298), "kinf_D": item.get("kinf_D", 298),
            "kint_A": item.get("kint_A", 1), "kint_B": item.get("kint_B", 0), "kint_C": item.get("kint_C", 0), "kint_D": item.get("kint_D", 298),
            "Fc": item.get("Fc", .6), "N": item.get("N", 1),
            "X": item.get("X", 0), "Y": item.get("Y", 0), "a0": item.get("a0", 0), "n": item.get("n", 0),
            "k_r": item.get("k_r", 0), "scaling_factor": item.get("scaling factor", 1),
            "reaction_probability": item.get("reaction probability", 1),
        }
        if kind == "BranchedRateConstant":
            for branch_key, branch_name in (("alkoxy products", "alkoxy"), ("nitrate products", "nitrate")):
                if branch_key in item:
                    reactions.append({"type": kind, "branch": branch_name, "reactants": reactants,
                                      "products": term_map(item[branch_key], "yield", 1.0),
                                      "label": item.get("MUSICA name"),
                                      "surface_species": None, "parameters": dict(fields)})
            continue
        if raw_type == "TAYLOR_SERIES":
            fields["coefficients"] = item.get("coefficients", [])
        reactions.append({
            "type": kind,
            "reactants": reactants,
            "products": products,
            "label": item.get("MUSICA name"),
            "surface_species": surface_species,
            "parameters": fields,
        })
    return {"name": config_path.parent.name, "temperature": 298.15, "pressure": 101325.0,
            "species": species, "parameters": {}, "reactions": reactions}


def normalized(path: Path) -> dict:
    raw = load(path)
    if isinstance(raw, dict) and "species" in raw and "reactions" in raw:
        return raw
    return from_camp(path)


def c_string(text: str) -> str:
    return json.dumps(text, ensure_ascii=True)


def emit(config: dict, output: Path) -> None:
    all_species = config["species"]
    integrated = [s for s in all_species if not s.get("third_body", False)]
    species_names = {s["name"]: i for i, s in enumerate(integrated)}
    third_bodies = {s["name"] for s in all_species if s.get("third_body", False)}
    reactions = sorted(config["reactions"], key=lambda r: RATE_ORDER[r["type"]])
    params = dict(config.get("parameters", {}))
    for reaction in reactions:
        if reaction["type"] == "UserDefinedRateConstant":
            label = reaction.get("label")
            if label:
                params.setdefault(label, 0.0)
        elif reaction["type"] == "SurfaceRateConstant":
            label = reaction.get("label") or f"surface_{reactions.index(reaction)}"
            params.setdefault(label + ".effective radius [m]", 0.0)
            params.setdefault(label + ".particle number concentration [# m-3]", 0.0)
    param_names = list(params)
    param_ids = {name: i for i, name in enumerate(param_names)}
    reaction_terms: list[tuple[list[tuple[int, float]], list[tuple[int, float]], dict]] = []
    flat_reactants, flat_products = [], []
    for reaction in reactions:
        reactants = []
        density_power = 0
        for term in reaction.get("reactants", []):
            name = term if isinstance(term, str) else term["species"]
            coeff = 1.0 if isinstance(term, str) else float(term.get("quantity", 1.0))
            if name in third_bodies:
                quantity = int(round(coeff))
                if quantity < 0 or abs(quantity - coeff) > 1e-12:
                    raise ValueError("parameterized reactant quantities must be non-negative integers")
                density_power += quantity
            else:
                if name not in species_names:
                    raise ValueError(f"unknown reactant {name!r}")
                quantity = int(round(coeff))
                if quantity < 0 or abs(quantity - coeff) > 1e-12:
                    raise ValueError("reactant quantities must be non-negative integers")
                reactants.extend((species_names[name], 1.0) for _ in range(quantity))
        products = []
        for term in reaction.get("products", []):
            name, coeff = term["species"], float(term.get("yield", 1.0))
            if name not in third_bodies:
                if name not in species_names:
                    raise ValueError(f"unknown product {name!r}")
                products.append((species_names[name], coeff))
        meta = dict(reaction)
        meta["air_density_power"] = density_power
        reaction_terms.append((reactants, products, meta))
        flat_reactants.extend(reactants)
        flat_products.extend(products)

    def number(x: float) -> str:
        x = float(x)
        if not math.isfinite(x):
            raise ValueError("non-finite value in mechanism")
        return format(x, ".17g")

    lines = ['#include "chem.h"', ""]
    lines.append(f"static const chem_species generated_species[{len(integrated)}] = {{")
    for species in integrated:
        diff = species.get("diffusion_coefficient")
        lines.append("  { %s, %s, %s, %s, %s, 0u }," % (
            c_string(species["name"]), number(species.get("molecular_weight", 0)),
            number(diff if diff is not None else 0), number(species.get("initial", 0)),
            number(species.get("absolute_tolerance", 1e-3))))
    lines.extend(["};", ""])
    lines.append(f"static const chem_term generated_reactants[{max(1, len(flat_reactants))}] = {{")
    lines.extend("  { %du, %s }," % (idx, number(coeff)) for idx, coeff in flat_reactants)
    if not flat_reactants:
        lines.append("  { 0u, 0.0 },")
    lines.extend(["};", ""])
    lines.append(f"static const chem_term generated_products[{max(1, len(flat_products))}] = {{")
    lines.extend("  { %du, %s }," % (idx, number(coeff)) for idx, coeff in flat_products)
    if not flat_products:
        lines.append("  { 0u, 0.0 },")
    lines.extend(["};", ""])
    lines.append(f"static const chem_parameter generated_parameters[{max(1, len(param_names))}] = {{")
    for name in param_names:
        lines.append("  { %s, %s }," % (c_string(name), number(params[name])))
    if not param_names:
        lines.append('  { "", 0.0 },')
    lines.extend(["};", ""])
    for index, (_, _, reaction) in enumerate(reaction_terms):
        coefficients = reaction.get("parameters", {}).get("coefficients", [])
        if coefficients:
            lines.append(f"static const double generated_coefficients_{index}[{len(coefficients)}] = {{")
            lines.append("  " + ", ".join(number(x) for x in coefficients))
            lines.extend(["};", ""])
    lines.append(f"static const chem_reaction generated_reactions[{len(reactions)}] = {{")
    roff = poff = 0
    for reaction_index, (reactants, products, reaction) in enumerate(reaction_terms):
        kind = reaction["type"]
        label = reaction.get("label") or ""
        p = [0.0] * 24
        for key, val in reaction.get("parameters", {}).items():
            if key in PARAM_INDEX:
                p[PARAM_INDEX[key]] = float(val)
        if kind == "BranchedRateConstant":
            pass
        if kind == "SurfaceRateConstant":
            p[0] = float(reaction.get("parameters", {}).get("reaction_probability", 1.0))
        custom_index = 0
        surface_species = 0
        if kind == "UserDefinedRateConstant" and label:
            custom_index = param_ids[label]
            p[1] = float(reaction.get("parameters", {}).get("scaling_factor", 1.0))
        elif kind == "SurfaceRateConstant":
            custom_index = param_ids[label + ".effective radius [m]"]
            name = reaction.get("surface_species", reaction.get("parameters", {}).get("surface_species", reactants[0][0] if reactants else ""))
            surface_species = species_names[name]
        rates = ", ".join(number(v) for v in p)
        coeffs = reaction.get("parameters", {}).get("coefficients", [])
        coeff_pointer = f"generated_coefficients_{reaction_index}" if coeffs else "NULL"
        branch = 1 if reaction.get("branch") == "alkoxy" else 0
        lines.append("  { %du, %du, %du, %du, { .kind = %s, .label = %s, .p = { %s }, .coefficients = %s, .coefficient_count = %du, .custom_parameter = %du, .surface_species = %du, .air_density_power = %du, .branch = %du } }," % (
            roff, len(reactants), poff, len(products), RATE_ENUM[kind], c_string(label), rates,
            coeff_pointer, len(coeffs), custom_index, surface_species, reaction["air_density_power"], branch))
        roff += len(reactants)
        poff += len(products)
    lines.extend(["};", ""])
    lines.extend([
        "static const chem_mechanism generated_mechanism = {",
        f"  {c_string(config.get('name', 'mechanism'))}, {len(integrated)}u, {len(reactions)}u,",
        f"  {len(flat_reactants)}u, {len(flat_products)}u, {len(param_names)}u,",
        "  generated_species, generated_reactions, generated_reactants, generated_products, generated_parameters,",
        f"  {number(config.get('temperature', 298.15))}, {number(config.get('pressure', 101325.0))}, {number(config.get('relative_tolerance', 1e-6))}",
        "};",
        "",
        "const chem_mechanism *chem_generated_mechanism(void)",
        "{",
        "  return &generated_mechanism;",
        "}",
        "",
    ])
    output.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    config = normalized(args.config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    emit(config, args.output)
    print(f"generated {args.output}: {config['name']}, {len(config['species'])} declared species, {len(config['reactions'])} reactions")


if __name__ == "__main__":
    main()
