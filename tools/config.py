"""Read and validate the supported CAMP-style mechanism configuration."""

import json
import math
from pathlib import Path

import yaml


def read_config(path):
    """Resolve camp-files relative to their parent, preserving declaration order."""
    entries, initial, units, constraints = [], {}, [], []
    visited, active = set(), set()

    def visit(file):
        file = Path(file).resolve()
        if file in active:
            raise ValueError(f"cyclic camp-files reference: {file}")
        if file in visited:
            return
        active.add(file)
        with file.open(encoding="utf-8") as stream:
            data = json.load(stream) if file.suffix.lower() == ".json" else yaml.safe_load(stream)
        if not isinstance(data, dict):
            raise ValueError(f"{file}: expected a mapping")
        for child in data.get("camp-files", []):
            if not isinstance(child, str):
                raise ValueError(f"{file}: camp-files entries must be paths")
            visit(file.parent / child)
        items = data.get("camp-data", [])
        if not isinstance(items, list):
            raise ValueError(f"{file}: camp-data must be a list")
        entries.extend(items)
        constraint_items = data.get("constraints", [])
        if not isinstance(constraint_items, list):
            raise ValueError(f"{file}: constraints must be a list")
        constraints.extend(constraint_items)
        if "initial conditions" in data:
            initial.update(data["initial conditions"])
        if "rate units" in data:
            units.append(data["rate units"])
        active.remove(file)
        visited.add(file)

    visit(path)
    if len(set(units)) > 1:
        raise ValueError("conflicting rate units")
    return {"camp-data": entries, "initial conditions": initial, "constraints": constraints,
            "rate units": units[0] if units else "molecules cm-3"}


def number(value, context):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: expected a number, got {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{context}: number must be finite")
    return result


def normalize(document):
    units = document.get("rate units", "molecules cm-3")
    if units not in ("molecules cm-3", "mol m-3"):
        raise ValueError(f"unsupported rate units {units!r}")
    species, by_name, reactions = [], {}, []
    for entry in document["camp-data"]:
        if not isinstance(entry, dict) or "type" not in entry:
            raise ValueError("each camp-data entry must have a type")
        if entry["type"] == "CHEM_SPEC":
            name = entry["name"]
            if not isinstance(name, str) or not name or name in by_name:
                raise ValueError(f"invalid or duplicate species name: {name!r}")
            param = entry.get("parameterization")
            if param is not None and param != "air density" and not isinstance(param, dict):
                raise ValueError(f"{name}: unknown parameterization {param!r}")
            if isinstance(param, dict):
                allowed = {"constant", "temperature", "pressure", "air density"}
                if set(param) - allowed:
                    raise ValueError(f"{name}: unknown parameterization coefficients {set(param) - allowed}")
                param = {key: number(param.get(key, 0), name) for key in allowed}
            item = {"name": name, "parameterized": param is not None,
                    "parameterization": param,
                    "initial": number(entry.get("initial concentration [mol m-3]", 0), name),
                    "absolute_tolerance": number(entry.get("absolute tolerance", 1e-3), name),
                    "molecular_weight": number(entry.get("molecular weight [kg mol-1]", 0), name),
                    "diffusion": number(entry.get("diffusion coefficient [m2 s-1]", 0), name)}
            by_name[name] = item
            if not item["parameterized"]:
                item["index"] = len(species)
                species.append(item)
        elif entry["type"] == "MECHANISM":
            reactions.extend(entry["reactions"])
        else:
            raise ValueError(f"unsupported camp-data type {entry['type']!r}")
    if not species:
        raise ValueError("mechanism has no state species")

    initial = document.get("initial conditions", {})
    defaults = initial.get("rate parameters", {})
    parameters, indices = [], {}

    def parameter(name):
        if not isinstance(name, str) or not name:
            raise ValueError("a named rate parameter is required")
        if name not in indices:
            indices[name] = len(parameters)
            parameters.append((name, number(defaults.get(name, 0), name)))
        return indices[name]

    def members(mapping, products=False):
        if not isinstance(mapping, dict):
            raise ValueError("reaction species must be a mapping")
        result, air_power, parameterizations = [], 0, []
        for name, properties in mapping.items():
            if name not in by_name:
                raise ValueError(f"unknown species {name!r}")
            if not isinstance(properties, dict):
                raise ValueError(f"{name}: expected a property mapping")
            spec = by_name[name]
            if products:
                if not spec["parameterized"]:
                    result.append((spec["index"], number(properties.get("yield", 1), name)))
            else:
                qty = number(properties.get("qty", 1), name)
                if qty < 0 or int(qty) != qty:
                    raise ValueError(f"{name}: reactant qty must be a nonnegative integer")
                if spec["parameterized"]:
                    if spec["parameterization"] == "air density":
                        air_power += int(qty)
                    else:
                        parameterizations.extend([spec["parameterization"]] * int(qty))
                else:
                    result.extend([spec["index"]] * int(qty))
        if air_power + len(parameterizations) > 4:
            raise ValueError("a reaction may have at most four parameterized reactants")
        return result, air_power, parameterizations

    output = []
    for raw in reactions:
        typ = raw["type"]
        if typ in ("FIRST_ORDER_LOSS", "EMISSION"):
            name = raw["species"]
            reactants, products = ({name: {}}, {}) if typ == "FIRST_ORDER_LOSS" else ({}, {name: {}})
        elif typ == "SURFACE" and "gas-phase products" in raw:
            reactants, products = {raw["gas-phase reactant"]: {}}, raw["gas-phase products"]
        else:
            reactants, products = raw.get("reactants", {}), raw.get("products", {})
        reactant_ids, air_power, parameterizations = members(reactants)
        rates = [(typ, products)]
        if typ == "BRANCHED":
            rates = [("BRANCHED_ALKOXY", raw["alkoxy products"]), ("BRANCHED_NITRATE", raw["nitrate products"])]
        for kind, product_map in rates:
            product_ids, _, _ = members(product_map, products=True)
            spec = {"kind": kind, "raw": raw, "reactants": reactant_ids,
                    "products": product_ids, "air_power": air_power, "parameterizations": parameterizations, "parameter": -1,
                    "molecule_units": units == "molecules cm-3"}
            if kind in ("PHOTOLYSIS", "USER_DEFINED", "FIRST_ORDER_LOSS", "EMISSION"):
                spec["parameter"] = parameter(raw["MUSICA name"])
            elif kind == "SURFACE":
                label = raw["MUSICA name"]
                spec["parameter"] = parameter(label + ".effective radius [m]")
                second = parameter(label + ".particle number concentration [# m-3]")
                if second != spec["parameter"] + 1:
                    raise ValueError(f"{label}: surface parameter names conflict")
                surface_species = by_name[raw["gas-phase reactant"]]
                if surface_species["diffusion"] <= 0 or surface_species["molecular_weight"] <= 0:
                    raise ValueError(f"{label}: surface species needs positive molecular weight and diffusion")
                spec["surface_species"] = surface_species
            output.append(spec)
    unused_defaults = set(defaults) - set(indices)
    if unused_defaults:
        raise ValueError(f"defaults for unknown parameters: {sorted(unused_defaults)}")
    constraints = []
    used_rows, constraint_names = set(), set()
    for raw in document.get("constraints", []):
        kind = raw["type"]
        if kind not in ("LINEAR", "EQUILIBRIUM"):
            raise ValueError(f"unsupported constraint type {kind!r}")
        name = raw["name"]
        if not isinstance(name, str) or not name or name in constraint_names:
            raise ValueError(f"invalid or duplicate constraint name {name!r}")
        constraint_names.add(name)
        row_name = raw["algebraic species"]
        if row_name not in by_name or by_name[row_name]["parameterized"]:
            raise ValueError(f"unknown algebraic species {row_name!r}")
        row = by_name[row_name]["index"]
        if row in used_rows:
            raise ValueError(f"multiple constraints for {row_name!r}")
        used_rows.add(row)

        def constraint_terms(mapping):
            if not isinstance(mapping, dict) or not mapping:
                raise ValueError(f"{name}: constraint terms must be a nonempty mapping")
            terms = []
            for species_name, coefficient in mapping.items():
                if species_name not in by_name or by_name[species_name]["parameterized"]:
                    raise ValueError(f"{name}: unknown constraint species {species_name!r}")
                value = number(coefficient, name)
                if kind == "EQUILIBRIUM" and value <= 0:
                    raise ValueError(f"{name}: equilibrium stoichiometry must be positive")
                terms.append((by_name[species_name]["index"], value))
            return terms

        left = constraint_terms(raw["terms"] if kind == "LINEAR" else raw["reactants"])
        right = constraint_terms(raw["products"]) if kind == "EQUILIBRIUM" else []
        k = number(raw.get("equilibrium constant", 0), name)
        if kind == "EQUILIBRIUM" and k <= 0:
            raise ValueError(f"{name}: equilibrium constant must be positive")
        constraints.append({"kind": kind, "name": name, "row": row, "reactants": left, "products": right,
                            "constant": number(raw.get("constant", 0), name), "equilibrium_constant": k,
                            "enthalpy": number(raw.get("enthalpy", 0), name),
                            "gas_constant": number(raw.get("gas constant", 0), name),
                            "reference_temperature": number(raw.get("reference temperature", 0), name)})
    return {"species": species, "reactions": output, "parameters": parameters, "constraints": constraints,
            "temperature": number(initial.get("temperature", 0), "temperature"),
            "pressure": number(initial.get("pressure", 0), "pressure")}
