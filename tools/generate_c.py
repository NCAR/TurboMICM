#!/usr/bin/env python3
"""Read a project mechanism JSON or MICM YAML/JSON config and emit C tables."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

KINDS = {
    "ARRHENIUS": ("CHEM_ARRHENIUS", ("A", "B", "C", "D", "E"), (1.0, 0, 0, 300, 0)),
    "TROE": ("CHEM_TROE", ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N"), (1, 0, 0, 1, 0, 0, 0.6, 1)),
    "TERNARY_CHEMICAL_ACTIVATION": ("CHEM_TERNARY", ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N"), (1, 0, 0, 1, 0, 0, 0.6, 1)),
    "TERNARY_CHEMICAL_ACTIVATION_JPL19": ("CHEM_TERNARY_JPL19", ("k0_A", "k0_B", "k0_C", "k0_D", "kinf_A", "kinf_B", "kinf_C", "kinf_D", "kint_A", "kint_B", "kint_C", "kint_D", "Fc", "N"), (1, 0, 0, 298, 1, 0, 0, 298, 1, 0, 0, 298, 0.6, 1)),
    "BRANCHED": ("CHEM_BRANCHED", ("Y", "X", "a0", "n"), (0, 1, 0.5, 0)),
    "TUNNELING": ("CHEM_TUNNELING", ("A", "B", "C"), (1, 0, 0)),
    "TAYLOR_SERIES": ("CHEM_TAYLOR", ("A", "B", "C", "D", "E"), (1, 0, 0, 300, 0)),
    "REVERSIBLE": ("CHEM_REVERSIBLE", ("A", "C", "k_r"), (1, 0, 0)),
    "USER_DEFINED": ("CHEM_USER_DEFINED", ("scaling_factor",), (1,)),
    "PHOTOLYSIS": ("CHEM_USER_DEFINED", ("scaling_factor",), (1,)),
    "FIRST_ORDER_LOSS": ("CHEM_USER_DEFINED", ("scaling_factor",), (1,)),
    "EMISSION": ("CHEM_USER_DEFINED", ("scaling_factor",), (1,)),
    "SURFACE": ("CHEM_SURFACE", ("reaction_probability",), (1,)),
}
MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * 6.02214076e23


def load_document(path: Path):
    text = path.read_text()
    if path.suffix.lower() == ".json":
        return json.loads(text)
    if yaml is None:
        raise SystemExit("YAML input needs PyYAML; run `make env` and use .conda/bin/python")
    return yaml.safe_load(text)


def read_config(path: Path) -> dict:
    doc = load_document(path)
    if isinstance(doc, dict) and "species" in doc and "reactions" in doc:
        return doc
    items = []
    files = doc.get("camp-files", []) if isinstance(doc, dict) else []
    for filename in files:
        loaded = load_document(path.parent / filename)
        items.extend(loaded.get("camp-data", []))
    if not items:
        items = doc.get("camp-data", []) if isinstance(doc, dict) else []
    species = []
    mechanisms = []
    for item in items:
        kind = str(item.get("type", "")).upper()
        if kind in {"CHEM_SPEC", "SPECIES", "THIRD_BODY"}:
            species.append(
                {
                    "name": item["name"],
                    "molecular_weight": item.get("molecular weight [kg mol-1]", 0.0),
                    "diffusion_coefficient": item.get("diffusion coefficient [m2 s-1]"),
                    "absolute_tolerance": item.get("absolute tolerance", 0.0),
                    "third_body": kind == "THIRD_BODY" or bool(item.get("third body", False)),
                    "initial": item.get("initial concentration", 0.0),
                }
            )
        elif kind == "MECHANISM":
            mechanisms.append(item)
    if not species or len(mechanisms) != 1:
        raise SystemExit(f"{path}: expected species and exactly one MECHANISM entry")
    reactions = []
    for raw in mechanisms[0].get("reactions", []):
        kind = str(raw.get("type", "")).upper()
        if kind not in KINDS:
            raise SystemExit(f"{path}: unsupported reaction type {kind}")
        params = {key: value for key, value in raw.items() if key not in {
            "type", "reactants", "products", "alkoxy products", "nitrate products",
            "gas-phase reactant", "gas-phase products", "species", "MUSICA name",
            "name", "label", "gas phase species",
        }}
        label = raw.get("MUSICA name", raw.get("label", raw.get("name")))
        if kind == "SURFACE":
            # Surface custom parameters are keyed by the gas species label in
            # the native C rate-constant structure.
            label = raw.get("gas-phase reactant", label)
        if kind == "TROE":
            # CAMP falloff inputs use molecule/cm3 units; the runtime uses
            # mol/m3 concentrations and therefore molar rate coefficients.
            params["k0_A"] = float(params.get("k0_A", 1.0)) * MOLES_M3_TO_MOLECULES_CM3 ** 2
            params["kinf_A"] = float(params.get("kinf_A", 1.0)) * MOLES_M3_TO_MOLECULES_CM3
        elif kind == "TERNARY_CHEMICAL_ACTIVATION":
            params["kinf_A"] = float(params.get("kinf_A", 1.0)) / MOLES_M3_TO_MOLECULES_CM3
        if kind in {"PHOTOLYSIS", "FIRST_ORDER_LOSS", "EMISSION", "USER_DEFINED", "SURFACE"} and label is None:
            label = f"rate_{len(reactions)}"
        reactants = raw.get("reactants", {})
        if isinstance(reactants, dict):
            expanded = []
            for name, data in reactants.items():
                qty = int(data.get("qty", 1)) if isinstance(data, dict) else int(data)
                if qty != (data.get("qty", 1) if isinstance(data, dict) else data):
                    raise SystemExit(f"{path}: reactant quantities must be integral")
                expanded.extend([name] * qty)
            reactants = expanded
        if kind == "SURFACE":
            gas_species = raw.get("gas-phase reactant")
            if not gas_species:
                raise SystemExit(f"{path}: SURFACE reaction needs gas-phase reactant")
            reactants = [gas_species]
            raw_products = raw.get("gas-phase products", raw.get("products", {}))
            products = raw_products
        elif kind == "FIRST_ORDER_LOSS":
            reactants = [raw["species"]]
            products = {}
        elif kind == "EMISSION":
            reactants = []
            products = {raw["species"]: {}}
        else:
            products = raw.get("products", {})

        def product_list(source):
            if isinstance(source, dict):
                return [
                    {"name": name, "yield": (data.get("yield", 1.0) if isinstance(data, dict) else data)}
                    for name, data in source.items()
                ]
            return source

        base = {"type": kind, "reactants": reactants, "parameters": params}
        if label is not None:
            base["label"] = label
        if kind == "BRANCHED":
            base["alkoxy_products"] = product_list(raw.get("alkoxy products", {}))
            base["nitrate_products"] = product_list(raw.get("nitrate products", {}))
            reactions.extend([dict(base, products=base["alkoxy_products"], branch=0),
                              dict(base, products=base["nitrate_products"], branch=1)])
        else:
            base["products"] = product_list(products)
            if kind == "SURFACE":
                base["surface_species"] = gas_species
                base["parameters"]["reaction_probability"] = raw.get("reaction probability", 1.0)
            reactions.append(base)
    return {
        "name": mechanisms[0].get("name", path.stem),
        "temperature": 298.15,
        "pressure": 101325.0,
        "species": species,
        "reactions": reactions,
        "custom_parameters": {},
    }


def c_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def number(value) -> str:
    result = float(value)
    if not math.isfinite(result):
        raise SystemExit(f"non-finite mechanism parameter {value!r}")
    return format(result, ".17g")


def markowitz_order(species: list[dict], reactions: list[dict]) -> list[int]:
    """Build a generic diagonal Markowitz order from the reaction Jacobian graph."""
    state_index = {item["name"]: i for i, item in enumerate(species)}
    n = len(species)
    rows = [set([i]) for i in range(n)]
    for reaction in reactions:
        reactants = [state_index[x] for x in reaction.get("reactants", []) if x in state_index]
        products = [state_index[x["name"]] for x in reaction.get("products", []) if x["name"] in state_index]
        for column in reactants:
            for row in reactants + products:
                rows[row].add(column)
    active = set(range(n))
    order = []
    while active:
        def score(i):
            row_count = sum(column in active for column in rows[i])
            column_count = sum(i in rows[row] for row in active)
            return (max(0, row_count - 1) * max(0, column_count - 1),
                    row_count + column_count, i)

        pivot = min(active, key=score)
        below = [row for row in active if row != pivot and pivot in rows[row]]
        right = [column for column in rows[pivot] if column in active and column != pivot]
        for row in below:
            rows[row].update(right)
        active.remove(pivot)
        order.append(pivot)
    return order


def compile_model(document: dict) -> str:
    all_species = document["species"]
    state_species = [s for s in all_species if not s.get("third_body", False)]
    permutation = markowitz_order(state_species, document["reactions"])
    state_species = [state_species[i] for i in permutation]
    species_index = {s["name"]: i for i, s in enumerate(state_species)}
    if len(species_index) != len(state_species):
        raise SystemExit("species names must be unique")
    custom_values = dict(document.get("custom_parameters", {}))
    reactions = []
    flat_reactants: list[int] = []
    flat_products: list[int] = []
    flat_yields: list[float] = []
    coefficients: list[float] = []

    for source in document["reactions"]:
        kind = str(source["type"]).upper()
        if kind not in KINDS:
            raise SystemExit(f"unsupported reaction type {kind}")
        c_kind, parameter_names, defaults = KINDS[kind]
        params = source.get("parameters", {})
        p = [0.0] * 16
        for i, (key, default) in enumerate(zip(parameter_names, defaults)):
            p[i] = float(params.get(key, default))
        if kind == "BRANCHED":
            p[4] = float(source.get("branch", 0))
        if kind == "USER_DEFINED" or kind == "PHOTOLYSIS" or kind in {"FIRST_ORDER_LOSS", "EMISSION"}:
            label = source.get("label")
            if label is None:
                raise SystemExit("user-defined reaction is missing a label")
            custom_values.setdefault(label, 0.0)
        else:
            label = None
        custom0 = custom1 = surface_index = -1
        if kind == "SURFACE":
            label = source.get("label", f"surface_{len(reactions)}")
            radius_name = f"{label}.effective radius [m]"
            number_name = f"{label}.particle number concentration [# m-3]"
            custom_values.setdefault(radius_name, 0.0)
            custom_values.setdefault(number_name, 0.0)
            custom0, custom1 = list(custom_values).index(radius_name), list(custom_values).index(number_name)
            surf_name = source.get("surface_species", source.get("reactants", [None])[0])
            if surf_name not in species_index:
                raise SystemExit(f"surface reactant {surf_name!r} must be a state species")
            surface_index = species_index[surf_name]
            surface = state_species[surface_index]
            molecular_weight = surface.get("molecular_weight", surface.get("molecular weight [kg mol-1]", 0.0))
            diffusion = surface.get("diffusion_coefficient", surface.get("diffusion coefficient [m2 s-1]"))
            if not molecular_weight or not diffusion:
                raise SystemExit(f"surface species {surf_name!r} needs molecular weight and diffusion coefficient")
        elif label is not None:
            custom0 = list(custom_values).index(label)
        if kind == "TAYLOR_SERIES":
            series = params.get("coefficients", [1.0])
            if len(series) > 16:
                raise SystemExit("Taylor series supports at most 16 coefficients")
            coeff_offset, coeff_count = len(coefficients), len(series)
            coefficients.extend(float(x) for x in series)
        else:
            coeff_offset, coeff_count = 0, 0
        reactant_offset = len(flat_reactants)
        body_order = 0
        for name in source.get("reactants", []):
            match = next((s for s in all_species if s["name"] == name), None)
            if match is None:
                raise SystemExit(f"reaction references unknown reactant {name!r}")
            if match.get("third_body", False):
                body_order += 1
            else:
                flat_reactants.append(species_index[name])
        product_offset = len(flat_products)
        for item in source.get("products", []):
            name = item["name"]
            match = next((s for s in all_species if s["name"] == name), None)
            if match is None:
                raise SystemExit(f"reaction references unknown product {name!r}")
            if match.get("third_body", False):
                continue
            flat_products.append(species_index[name])
            flat_yields.append(float(item.get("yield", 1.0)))
        r = {
            "kind": c_kind,
            "rbegin": reactant_offset,
            "nr": len(flat_reactants) - reactant_offset,
            "pbegin": product_offset,
            "np": len(flat_products) - product_offset,
            "body": body_order,
            "coeff_offset": coeff_offset,
            "coeff_count": coeff_count,
            "custom0": custom0,
            "custom1": custom1,
            "surface": surface_index,
            "params": p,
        }
        reactions.append(r)

    custom_names = list(custom_values)
    custom = [(name, float(custom_values[name])) for name in custom_names]
    output = [
        '#include "kinetics.h"',
        "",
        "static const chem_species chem_species_data[] = {",
    ]
    for s in state_species:
        mw = s.get("molecular_weight", s.get("molecular weight [kg mol-1]", 0.0)) or 0.0
        diffusion = s.get("diffusion_coefficient", s.get("diffusion coefficient [m2 s-1]")) or 0.0
        output.append(
            f"  {{.name={c_string(s['name'])}, .molecular_weight={number(mw)}, "
            f".diffusion_coefficient={number(diffusion)}, .initial_value={number(s.get('initial', 0.0))}, "
            f".absolute_tolerance={number(s.get('absolute_tolerance', s.get('absolute tolerance', 0.0)))}}},"
        )
    output += ["};", "", "static const chem_reaction chem_reaction_data[] = {"]
    for r in reactions:
        p_values = ", ".join(number(x) for x in r["params"])
        output.append(
            "  {.kind=%s, .reactant_offset=%d, .reactant_count=%d, .product_offset=%d, .product_count=%d, "
            ".body_order=%d, .coefficient_offset=%d, .coefficient_count=%d, .custom0=%d, .custom1=%d, "
            ".surface_species=%d, .p={%s}},"
            % (r["kind"], r["rbegin"], r["nr"], r["pbegin"], r["np"], r["body"],
               r["coeff_offset"], r["coeff_count"], r["custom0"], r["custom1"], r["surface"], p_values)
        )
    output += ["};", ""]
    output.append("static const uint32_t chem_reactant_data[] = {" + ",".join(map(str, flat_reactants or [0])) + "};")
    output.append("static const uint32_t chem_product_data[] = {" + ",".join(map(str, flat_products or [0])) + "};")
    output.append("static const double chem_yield_data[] = {" + ",".join(number(x) for x in (flat_yields or [0])) + "};")
    output.append("static const double chem_coefficient_data[] = {" + ",".join(number(x) for x in (coefficients or [0])) + "};")
    output += ["", "static const chem_custom_parameter chem_custom_data[] = {"]
    for name, value in custom:
        output.append(f"  {{{c_string(name)}, {number(value)}}},")
    if not custom:
        output.append('  {"", 0.0},')
    output += ["};", "", "const chem_mechanism chem_model = {"]
    output += [
        f"  .name={c_string(document.get('name', 'mechanism'))},",
        f"  .species_count={len(state_species)}, .reaction_count={len(reactions)},",
        f"  .reactant_count={len(flat_reactants)}, .product_count={len(flat_products)}, .custom_count={len(custom)},",
        f"  .default_temperature={number(document.get('temperature', 298.15))},",
        f"  .default_pressure={number(document.get('pressure', 101325.0))},",
        "  .species=chem_species_data, .reactions=chem_reaction_data,",
        "  .reactants=chem_reactant_data, .products=chem_product_data, .yields=chem_yield_data,",
        "  .coefficients=chem_coefficient_data, .custom_parameters=chem_custom_data",
        "};",
        "",
    ]
    return "\n".join(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    document = read_config(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(compile_model(document))


if __name__ == "__main__":
    main()
