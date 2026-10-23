#!/usr/bin/env python3
"""Read a MICM-format YAML/JSON mechanism and emit a C descriptor."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import yaml


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
        # The tutorial configuration examples express Troe low/high-pressure
        # coefficients in molecule/cm units. Normalize them to the mol/m3
        # convention used by the rate evaluator and the matching hand-built
        # example parameters.
        molecule_cm3 = 1.0e-6 * 6.02214076e23
        if raw_type == "TROE":
            fields["k0_A"] *= molecule_cm3 * molecule_cm3
            fields["kinf_A"] *= molecule_cm3
        elif raw_type == "TERNARY_CHEMICAL_ACTIVATION":
            fields["kinf_A"] /= molecule_cm3
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


def solver_mechanism(config: dict):
    """Return the fixed species, rates, reactions, and sparse LU plan."""
    all_species = config["species"]
    integrated = [item for item in all_species if not item.get("third_body", False)]
    species_ids = {item["name"]: index for index, item in enumerate(integrated)}
    third_bodies = {item["name"] for item in all_species if item.get("third_body", False)}
    parameters = dict(config.get("parameters", {}))
    reactions = []
    for index, source in enumerate(config["reactions"]):
        reaction = dict(source)
        rate_parameters = dict(source.get("parameters", {}))
        label = source.get("label")
        if source["type"] == "UserDefinedRateConstant" and label:
            parameters.setdefault(label, 0.0)
        elif source["type"] == "SurfaceRateConstant":
            label = label or f"surface_{index}"
            parameters.setdefault(label + ".effective radius [m]", 0.0)
            parameters.setdefault(label + ".particle number concentration [# m-3]", 0.0)
            reaction["label"] = label
        reactants, density_power = [], 0
        for term in source.get("reactants", []):
            name = term if isinstance(term, str) else term["species"]
            coefficient = 1.0 if isinstance(term, str) else float(term.get("quantity", 1.0))
            quantity = int(round(coefficient))
            if quantity < 0 or abs(quantity - coefficient) > 1.0e-12:
                raise ValueError("reactant quantities must be non-negative integers")
            if name in third_bodies:
                density_power += quantity
            elif name in species_ids:
                reactants.extend([species_ids[name]] * quantity)
            else:
                raise ValueError(f"unknown reactant {name!r}")
        products = []
        for term in source.get("products", []):
            name, coefficient = term["species"], float(term.get("yield", 1.0))
            if name in third_bodies:
                continue
            if name not in species_ids:
                raise ValueError(f"unknown product {name!r}")
            products.append((species_ids[name], coefficient))
        reaction["_reactants"] = reactants
        reaction["_products"] = products
        reaction["_density_power"] = density_power
        reaction["_parameters"] = rate_parameters
        reactions.append(reaction)

    n = len(integrated)
    pattern = [[i == j for j in range(n)] for i in range(n)]
    for reaction in reactions:
        for column in reaction["_reactants"]:
            for row in reaction["_reactants"]:
                pattern[row][column] = True
            for row, _ in reaction["_products"]:
                pattern[row][column] = True

    # Choose the pivot with the fewest actual missing fill entries, then use
    # Markowitz cost and species index as deterministic tie breakers. The
    # no-pivot Mozart fill is fixed at generation time.
    graph = [row[:] for row in pattern]
    alive = [True] * n
    row_degree = [sum(row) for row in graph]
    column_degree = [sum(graph[row][column] for row in range(n)) for column in range(n)]
    permutation = []
    for _ in range(n):
        candidates = []
        for item in range(n):
            if not alive[item]:
                continue
            columns = [column for column in range(n) if column != item and alive[column] and graph[item][column]]
            incoming = [row for row in range(n) if row != item and alive[row] and graph[row][item]]
            actual_fill = sum(not graph[row][column] for row in incoming for column in columns)
            markowitz_cost = (row_degree[item] - 1) * (column_degree[item] - 1)
            candidates.append(((actual_fill, markowitz_cost, row_degree[item] + column_degree[item], item), item))
        pivot = min(candidates)[1]
        permutation.append(pivot)
        alive[pivot] = False
        columns = [item for item in range(n) if item != pivot and alive[item] and graph[pivot][item]]
        incoming = [item for item in range(n) if item != pivot and alive[item] and graph[item][pivot]]
        for row in incoming:
            for column in columns:
                if not graph[row][column]:
                    graph[row][column] = True
                    row_degree[row] += 1
                    column_degree[column] += 1
        for column in columns:
            graph[pivot][column] = False
            column_degree[column] -= 1
        for row in incoming:
            graph[row][pivot] = False
            row_degree[row] -= 1

    inverse = [0] * n
    for position, original in enumerate(permutation):
        inverse[original] = position
    permuted = [[pattern[permutation[row]][permutation[column]] for column in range(n)] for row in range(n)]
    for pivot in range(n):
        for row in range(pivot + 1, n):
            if permuted[row][pivot]:
                for column in range(pivot + 1, n):
                    if permuted[pivot][column]:
                        permuted[row][column] = True
    positions = {}
    for row in range(n):
        for column in range(n):
            if permuted[row][column]:
                positions[row, column] = len(positions)
    lower = [[row for row in range(pivot + 1, n) if permuted[row][pivot]] for pivot in range(n)]
    upper = [[column for column in range(pivot + 1, n) if permuted[pivot][column]] for pivot in range(n)]
    row_lower = [[column for column in range(row) if permuted[row][column]] for row in range(n)]
    row_upper = [[column for column in range(row + 1, n) if permuted[row][column]] for row in range(n)]
    updates = []
    for pivot in range(n):
        for column in upper[pivot]:
            for row in lower[pivot]:
                if permuted[row][column]:
                    updates.append((positions[row, column], positions[row, pivot], positions[pivot, column]))

    jobs = []
    for independent, original_independent in enumerate(permutation):
        for reaction_index, reaction in enumerate(reactions):
            for column in reaction["_reactants"]:
                if column != original_independent:
                    continue
                dependent, found = [], False
                for original_row in reaction["_reactants"]:
                    if original_row == original_independent and not found:
                        found = True
                        continue
                    dependent.append(inverse[original_row])
                jobs.append((reaction_index, independent, dependent, reaction["_products"]))
    return {
        "species": integrated, "species_ids": species_ids, "parameters": parameters,
        "reactions": reactions, "permutation": permutation, "inverse": inverse,
        "positions": positions, "lower": lower, "upper": upper,
        "row_lower": row_lower, "row_upper": row_upper, "updates": updates, "jobs": jobs,
    }


def _rate_code(reaction: dict, reaction_index: int, species: list[dict], parameter_ids: dict[str, int]) -> list[str]:
    p = reaction["_parameters"]
    kind = reaction["type"]
    out = [f"  {{ double *row = rates + {reaction_index}u * cell_count;", "    for (size_t c = 0; c < cell_count; ++c) {",
           "      const double temperature_c = temperature[c];", "      const double pressure_c = pressure[c];",
           "      const double air_density = pressure_c / (1.380649e-23 * 6.02214076e23 * temperature_c);",
           "      (void)pressure_c; (void)air_density;",
           "      double value = 0.0;"]
    def v(key: str, default: float = 0.0) -> str:
        return format(float(p.get(key, default)), ".17g")
    if kind == "ArrheniusRateConstant":
        out.append(f"      value = {v('A', 1)} * exp({v('C')} / temperature_c) * pow(temperature_c / {v('D', 300)}, {v('B')}) * (1.0 + {v('E')} * pressure_c);")
    elif kind in {"TroeRateConstant", "TernaryChemicalActivationRateConstant"}:
        multiply_density = "air_density" if kind == "TroeRateConstant" else "1.0"
        out.extend([
            f"      const double k0 = {v('k0_A', 1)} * exp({v('k0_C')} / temperature_c) * pow(temperature_c / 300.0, {v('k0_B')});",
            f"      const double kinf = {v('kinf_A', 1)} * exp({v('kinf_C')} / temperature_c) * pow(temperature_c / 300.0, {v('kinf_B')});",
            "      const double ratio = k0 * air_density / kinf;",
            f"      value = k0 * {multiply_density} / (1.0 + ratio) * pow({v('Fc', .6)}, {v('N', 1)} / ({v('N', 1)} + pow(log10(ratio), 2.0)));",
        ])
    elif kind == "TernaryChemicalActivationJPL19":
        out.extend([
            f"      const double k0 = {v('k0_A', 1)} * exp({v('k0_C')} / temperature_c) * pow(temperature_c / {v('k0_D', 298)}, {v('k0_B')});",
            f"      const double kinf = {v('kinf_A', 1)} * exp({v('kinf_C')} / temperature_c) * pow(temperature_c / {v('kinf_D', 298)}, {v('kinf_B')});",
            f"      const double kint = {v('kint_A', 1)} * exp({v('kint_C')} / temperature_c) * pow(temperature_c / {v('kint_D', 298)}, {v('kint_B')});",
            "      const double ratio = k0 * air_density / kinf;",
            f"      const double kf = k0 * air_density / (1.0 + ratio) * pow({v('Fc', .6)}, {v('N', 1)} / ({v('N', 1)} + pow(log10(ratio), 2.0)));",
            "      value = kf + kint * (1.0 - kf / kinf);",
        ])
    elif kind == "BranchedRateConstant":
        branch = int(reaction.get("branch") == "alkoxy")
        out.extend([
            "      const double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(" + v("n") + ");",
            "      const double air_ref = 2.45e19 / 6.02214076e23 * 1.0e6;",
            "      const double a_ref = k0 * air_ref;",
            "      const double b_ref = 0.43 * pow(293.0 / 298.0, -8.0);",
            "      const double aval_ref = a_ref / (1.0 + a_ref / b_ref) * pow(0.41, 1.0 / (1.0 + pow(log10(a_ref / b_ref), 2.0)));",
            f"      const double z = aval_ref * (1.0 - {v('a0')}) / {v('a0')};",
            "      const double a = k0 * air_density;",
            "      const double b = 0.43 * pow(temperature_c / 298.0, -8.0);",
            "      const double aval = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));",
            f"      const double pre = {v('X')} * exp(-{v('Y')} / temperature_c);",
            f"      value = {('pre * (z / (z + aval))' if branch else 'pre * (aval / (aval + z))')};",
        ])
    elif kind == "TunnelingRateConstant":
        out.append(f"      value = {v('A', 1)} * exp(-{v('B')} / temperature_c + {v('C')} / (temperature_c * temperature_c * temperature_c));")
    elif kind == "TaylorSeriesRateConstant":
        coefficients = p.get("coefficients", [])
        out.extend(["      double polynomial = 0.0;", "      double power = 1.0;"])
        for coefficient in coefficients:
            out.append(f"      polynomial += {format(float(coefficient), '.17g')} * power;")
            out.append("      power *= temperature_c;")
        out.append(f"      value = polynomial * {v('A', 1)} * exp({v('C')} / temperature_c) * pow(temperature_c / {v('D', 300)}, {v('B')}) * (1.0 + {v('E')} * pressure_c);")
    elif kind == "ReversibleRateConstant":
        out.append(f"      value = {v('A', 1)} * exp({v('C')} / temperature_c) * {v('k_r')};")
    elif kind == "UserDefinedRateConstant":
        label = reaction.get("label", "")
        if label not in parameter_ids:
            raise ValueError(f"user-defined rate {label!r} has no runtime parameter")
        out.append(f"      value = parameters[{parameter_ids[label]}u * cell_count + c] * {v('scaling_factor', 1)};")
    elif kind == "SurfaceRateConstant":
        label = reaction["label"]
        first, second = parameter_ids[label + ".effective radius [m]"], parameter_ids[label + ".particle number concentration [# m-3]"]
        surface_name = reaction.get("surface_species", p.get("surface_species"))
        if surface_name is None:
            surface_name = species[reaction["_reactants"][0]]["name"]
        surface_species = species[next(i for i, item in enumerate(species) if item["name"] == surface_name)]
        mw = float(surface_species.get("molecular_weight", 0.0))
        diffusion = float(surface_species.get("diffusion_coefficient") or 0.0)
        out.extend([
            f"      const double radius = parameters[{first}u * cell_count + c];",
            f"      const double number = parameters[{second}u * cell_count + c];",
            f"      const double speed = sqrt(8.0 * 8.31446261815324 * temperature_c / (3.14159265358979323846 * {format(mw, '.17g')}));",
            f"      value = 4.0 * number * 3.14159265358979323846 * radius * radius / (radius / {format(diffusion, '.17g')} + 4.0 / (speed * {v('reaction_probability', 1)}));",
        ])
    else:
        raise ValueError(f"unsupported rate type {kind!r}")
    for _ in range(reaction["_density_power"]):
        out.append("      value *= air_density;")
    out.extend(["      row[c] = value;", "    }", "  }"])
    return out


def specialized_source(config: dict) -> list[str]:
    model = solver_mechanism(config)
    species, reactions = model["species"], model["reactions"]
    n, nr, nnz = len(species), len(reactions), len(model["positions"])
    parameter_ids = {name: index for index, name in enumerate(model["parameters"])}
    def flatten(rows, value):
        offsets, values = [0], []
        for row_index, row in enumerate(rows):
            values.extend(value(row_index, column) for column in row)
            offsets.append(len(values))
        return offsets, values
    row_lower_offsets, row_lower_positions = flatten(
        model["row_lower"], lambda row, column: model["positions"][row, column])
    row_lower_columns = [column for row in model["row_lower"] for column in row]
    row_upper_offsets, row_upper_positions = flatten(
        model["row_upper"], lambda row, column: model["positions"][row, column])
    row_upper_columns = [column for row in model["row_upper"] for column in row]

    def size_array(name: str, values: list[int]) -> list[str]:
        result = [f"static const size_t {name}[{max(1, len(values))}] = {{"]
        entries = values or [0]
        for offset in range(0, len(entries), 12):
            result.append("  " + ", ".join(f"{value}u" for value in entries[offset:offset + 12]) + ",")
        result.append("};")
        return result

    lines = [
        "#include \"chem_specialized.h\"", "#include <float.h>", "#include <math.h>",
        "#include <stdint.h>", "#include <stdio.h>", "#include <stdlib.h>", "#include <string.h>", "",
        f"enum {{ GENERATED_SPECIES = {n}, GENERATED_REACTIONS = {nr}, GENERATED_NNZ = {nnz}, GENERATED_PARAMETERS = {len(parameter_ids)} }};",
        "static const size_t generated_permutation[GENERATED_SPECIES] = {",
        "  " + ", ".join(f"{item}u" for item in model["permutation"]), "};",
        "static const double generated_sorted_atol[GENERATED_SPECIES] = {",
        "  " + ", ".join(format(float(species[item].get("absolute_tolerance", 1.0e-3)), ".17g") for item in model["permutation"]), "};",
        "static const double generated_parameter_defaults[GENERATED_PARAMETERS ? GENERATED_PARAMETERS : 1] = {",
        "  " + (", ".join(format(float(value), ".17g") for value in model["parameters"].values()) or "0.0"), "};", "",
        "static void generated_rates(size_t cell_count, const double *restrict parameters,",
        "                           const double *restrict temperature, const double *restrict pressure,",
        "                           double *restrict rates)", "{",
    ]
    for index, reaction in enumerate(reactions):
        lines.extend(_rate_code(reaction, index, species, parameter_ids))
    lines.extend(["}", "", "static void generated_forcing(size_t cell_count, const double *restrict y,",
                  "                              const double *restrict rates, double *restrict forcing)", "{",
                  "  memset(forcing, 0, GENERATED_SPECIES * cell_count * sizeof(double));"])
    for index, reaction in enumerate(reactions):
        flux = f"rates[{index}u * cell_count + c]"
        for species_index in reaction["_reactants"]:
            flux += f" * y[{model['inverse'][species_index]}u * cell_count + c]"
        lines.extend(["  #pragma omp simd", "  for (size_t c = 0; c < cell_count; ++c) {", f"    const double flux = {flux};"])
        for species_index in reaction["_reactants"]:
            row = model["inverse"][species_index]
            lines.append(f"    forcing[{row}u * cell_count + c] -= flux;")
        for species_index, coefficient in reaction["_products"]:
            row = model["inverse"][species_index]
            lines.append(f"    forcing[{row}u * cell_count + c] += {format(coefficient, '.17g')} * flux;")
        lines.append("  }")
    lines.extend(["}", "", "static void generated_jacobian(size_t cell_count, const double *restrict y,",
                  "                               const double *restrict rates, double *restrict matrix)", "{",
                  "  memset(matrix, 0, GENERATED_NNZ * cell_count * sizeof(double));"])
    for reaction_index, independent, dependent, products in model["jobs"]:
        reaction = reactions[reaction_index]
        expr = f"rates[{reaction_index}u * cell_count + c]"
        for dep in dependent:
            expr += f" * y[{dep}u * cell_count + c]"
        lines.extend(["  #pragma omp simd", "  for (size_t c = 0; c < cell_count; ++c) {", f"    const double derivative = {expr};"])
        for row in dependent:
            position = model["positions"][row, independent]
            lines.append(f"    matrix[{position}u * cell_count + c] += derivative;")
        diagonal = model["positions"][independent, independent]
        lines.append(f"    matrix[{diagonal}u * cell_count + c] += derivative;")
        for species_index, coefficient in products:
            row = model["inverse"][species_index]
            position = model["positions"][row, independent]
            lines.append(f"    matrix[{position}u * cell_count + c] -= {format(coefficient, '.17g')} * derivative;")
        lines.append("  }")
    lines.append("}")
    lines.extend(size_array("generated_diagonal", [model["positions"][i, i] for i in range(n)]))
    lines.extend(size_array("generated_row_lower_offsets", row_lower_offsets))
    lines.extend(size_array("generated_row_lower_positions", row_lower_positions))
    lines.extend(size_array("generated_row_lower_columns", row_lower_columns))
    lines.extend(size_array("generated_row_upper_offsets", row_upper_offsets))
    lines.extend(size_array("generated_row_upper_positions", row_upper_positions))
    lines.extend(size_array("generated_row_upper_columns", row_upper_columns))
    lines.append("")
    for pivot in range(n):
        diagonal = model["positions"][pivot, pivot]
        lines.extend([f"static int generated_factor_pivot_{pivot}(size_t cell_count, double *restrict matrix,",
                      "                                      double *restrict inverse, size_t *bad_cell)", "{",
                      "  for (size_t c = 0; c < cell_count; ++c) {",
                      f"    const double value = matrix[{diagonal}u * cell_count + c];",
                      "    if (value == 0.0 || !isfinite(value)) { *bad_cell = c; return -1; }",
                      "    inverse[c] = 1.0 / value;", "  }"])
        for row in model["lower"][pivot]:
            position = model["positions"][row, pivot]
            lines.extend(["  #pragma omp simd", "  for (size_t c = 0; c < cell_count; ++c)",
                          f"    matrix[{position}u * cell_count + c] *= inverse[c];"])
        for column in model["upper"][pivot]:
            for row in model["lower"][pivot]:
                if (row, column) in model["positions"]:
                    target = model["positions"][row, column]
                    lower_position = model["positions"][row, pivot]
                    upper_position = model["positions"][pivot, column]
                    lines.extend(["  #pragma omp simd", "  for (size_t c = 0; c < cell_count; ++c)",
                                  f"    matrix[{target}u * cell_count + c] -= matrix[{lower_position}u * cell_count + c] * matrix[{upper_position}u * cell_count + c];"])
        lines.extend(["  return 0;", "}", ""])
    lines.extend(["static int generated_factor(size_t cell_count, double *restrict matrix,",
                  "                           double *restrict inverse, size_t *bad_pivot, size_t *bad_cell)", "{"])
    for pivot in range(n):
        lines.extend([f"  if (generated_factor_pivot_{pivot}(cell_count, matrix, inverse, bad_cell) != 0) {{",
                      f"    *bad_pivot = {pivot}u; return -1;", "  }"])
    lines.extend(["  return 0;", "}", "", "static void generated_solve(size_t cell_count, const double *restrict matrix, double *restrict rhs)", "{"])
    lines.extend([
                  "  for (size_t row = 0; row < GENERATED_SPECIES; ++row)",
                  "    for (size_t q = generated_row_lower_offsets[row]; q < generated_row_lower_offsets[row + 1u]; ++q) {",
                  "      const double *entry = matrix + generated_row_lower_positions[q] * cell_count;",
                  "      double *target = rhs + row * cell_count;",
                  "      const double *source = rhs + generated_row_lower_columns[q] * cell_count;",
                  "      #pragma omp simd", "      for (size_t c = 0; c < cell_count; ++c) target[c] -= entry[c] * source[c];", "    }",
                  "  for (size_t row = GENERATED_SPECIES; row-- > 0;) {",
                  "    for (size_t q = generated_row_upper_offsets[row]; q < generated_row_upper_offsets[row + 1u]; ++q) {",
                  "      const double *entry = matrix + generated_row_upper_positions[q] * cell_count;",
                  "      double *target = rhs + row * cell_count;",
                  "      const double *source = rhs + generated_row_upper_columns[q] * cell_count;",
                  "      #pragma omp simd", "      for (size_t c = 0; c < cell_count; ++c) target[c] -= entry[c] * source[c];", "    }",
                  "    const double *diagonal = matrix + generated_diagonal[row] * cell_count;",
                  "    double *target = rhs + row * cell_count;",
                  "    #pragma omp simd", "    for (size_t c = 0; c < cell_count; ++c) target[c] /= diagonal[c];", "  }", "}"])
    lines.extend(["", "static double generated_error(size_t cell_count, const double *restrict y,",
                  "                             const double *restrict ynew, const double *restrict error,",
                  "                             double relative_tolerance, double absolute_tolerance)", "{",
                  "  double sum = 0.0;"])
    lines.extend(["  for (size_t species = 0; species < GENERATED_SPECIES; ++species) {",
                  "    const double atol = absolute_tolerance >= 0.0 ? absolute_tolerance : generated_sorted_atol[species];",
                  "    for (size_t c = 0; c < cell_count; ++c) {",
                  "      const size_t i = species * cell_count + c;",
                  "      const double ymax = fmax(fabs(y[i]), fabs(ynew[i]));",
                  "      const double scale = atol + relative_tolerance * ymax;",
                  "      if (!isfinite(y[i]) || !isfinite(ynew[i]) || !isfinite(error[i]) || !isfinite(atol) || atol < 0.0 ||",
                  "          !isfinite(scale) || scale < 0.0 || (scale == 0.0 && error[i] != 0.0)) return NAN;",
                  "      const double scaled = scale == 0.0 ? 0.0 : error[i] / scale;",
                  "      sum += scaled * scaled;", "    }", "  }"])
    lines.extend(["  return fmax(sqrt(sum / (double)(cell_count * GENERATED_SPECIES)), 1.0e-10);", "}", ""])
    lines.extend(_emit_integrator(
        n, nr, nnz, len(parameter_ids), float(config.get("relative_tolerance", 1.0e-6)),
        [model["positions"][species, species] for species in range(n)],
    ))
    return lines


def _emit_integrator(species_count: int, reaction_count: int, nonzero_count: int, parameter_count: int,
                     default_relative_tolerance: float, diagonal_positions: list[int]) -> list[str]:
    """Emit the fixed three-stage MICM Rosenbrock integration path."""
    n = species_count
    lines = [
        "int chem_generated_integrate_group_steps(size_t cell_count, double *restrict y_external,",
        "    const double *restrict parameters_external, const double *restrict temperature,",
        "    const double *restrict pressure, double interval, size_t interval_count,",
        "    double relative_tolerance, double absolute_tolerance, chem_statistics *statistics)",
        "{",
        "  if (!cell_count || !interval_count || !y_external || !temperature || !pressure ||",
        "      !(interval > 0.0) || !isfinite(interval)) return -1;",
        f"  if (relative_tolerance < 0.0) relative_tolerance = {format(default_relative_tolerance, '.17g')};",
        "  if (!isfinite(relative_tolerance) || relative_tolerance < 0.0 ||",
        "      !isfinite(absolute_tolerance)) return -1;",
        "  if (cell_count > SIZE_MAX / GENERATED_SPECIES / sizeof(double) ||",
        "      cell_count > SIZE_MAX / GENERATED_NNZ / sizeof(double) ||",
        "      cell_count > SIZE_MAX / GENERATED_REACTIONS / sizeof(double)) return -1;",
        "  const size_t state_size = GENERATED_SPECIES * cell_count;",
        "  const size_t parameter_size = (GENERATED_PARAMETERS ? GENERATED_PARAMETERS : 1u) * cell_count;",
        "  double *y = calloc(state_size, sizeof(double));",
        "  double *ynew = calloc(state_size, sizeof(double));",
        "  double *initial_forcing = calloc(state_size, sizeof(double));",
        "  double *stage_forcing = calloc(state_size, sizeof(double));",
        "  double *error = calloc(state_size, sizeof(double));",
        "  double *parameters = calloc(parameter_size, sizeof(double));",
        "  double *rates = calloc(GENERATED_REACTIONS * cell_count, sizeof(double));",
        "  double *matrix = calloc(GENERATED_NNZ * cell_count, sizeof(double));",
        "  double *inverse = calloc(cell_count, sizeof(double));",
        "  double *k0 = calloc(state_size, sizeof(double));",
        "  double *k1 = calloc(state_size, sizeof(double));",
        "  double *k2 = calloc(state_size, sizeof(double));",
        "  if (!y || !ynew || !initial_forcing || !stage_forcing || !error || !parameters || !rates ||",
        "      !matrix || !inverse || !k0 || !k1 || !k2) {",
        "    free(y); free(ynew); free(initial_forcing); free(stage_forcing); free(error);",
        "    free(parameters); free(rates); free(matrix); free(inverse); free(k0); free(k1); free(k2);",
        "    return -1;",
        "  }",
        "  for (size_t c = 0; c < cell_count; ++c) {",
        "    for (size_t s = 0; s < GENERATED_SPECIES; ++s)",
        "      y[s * cell_count + c] = y_external[c * GENERATED_SPECIES + generated_permutation[s]];",
        "    for (size_t p = 0; p < GENERATED_PARAMETERS; ++p)",
        "      parameters[p * cell_count + c] = parameters_external ? parameters_external[c * GENERATED_PARAMETERS + p] : generated_parameter_defaults[p];",
        "  }",
        "  generated_rates(cell_count, parameters, temperature, pressure, rates);",
        "  const double h_min = 1.0e-15 * interval;",
        "  const double h_max = interval;",
        "  const double h_start = fmin(1.0e-6 * interval, h_max);",
        "  int failed = 0;",
        "  for (size_t interval_index = 0; interval_index < interval_count; ++interval_index) {",
        "    double h = fmin(fmax(h_min, h_start), h_max);",
        "    double present_time = 0.0;",
        "    int reject_last_h = 0, reject_more_h = 0;",
        "    size_t step_attempts = 0, accepted_this_interval = 0;",
        "    while ((present_time - interval + DBL_EPSILON) <= 0.0) {",
        "      if (step_attempts > 1000u) { fputs(\"solver exceeded 1000 internal steps\\n\", stderr); failed = 1; goto cleanup; }",
        "      if ((present_time + 0.1 * h) == present_time || h <= fmin(DBL_EPSILON, h_min)) {",
        "        fprintf(stderr, \"solver step size too small: t=%.17g h=%.17g\\n\", present_time, h); failed = 1; goto cleanup;",
        "      }",
        "      h = fmin(h, fabs(interval - present_time));",
        "      generated_forcing(cell_count, y, rates, initial_forcing);",
        "      if (statistics) ++statistics->function_calls;",
        "      generated_jacobian(cell_count, y, rates, matrix);",
        "      if (statistics) ++statistics->jacobian_updates;",
        "      int accepted = 0;",
        "      while (!accepted) {",
        "        const double alpha = 1.0 / (h * 0.43586652150845899941601945119356);",
        "        for (size_t s = 0; s < GENERATED_SPECIES; ++s) {",
        "          double *diagonal = matrix + generated_diagonal[s] * cell_count;",
        "          #pragma omp simd",
        "          for (size_t c = 0; c < cell_count; ++c) diagonal[c] += alpha;",
        "        }",
    ]
    lines.extend([
        "        size_t bad_pivot = 0, bad_cell = 0;",
        "        if (generated_factor(cell_count, matrix, inverse, &bad_pivot, &bad_cell) != 0) {",
        "          fprintf(stderr, \"linear factorization failed at cell %zu time %.17g h %.17g pivot %zu value %.17g\\n\",",
        "                  bad_cell, present_time, h, bad_pivot, matrix[generated_diagonal[bad_pivot] * cell_count + bad_cell]);",
        "          failed = 1; goto cleanup;",
        "        }",
        "        if (statistics) ++statistics->decompositions;",
        "        memcpy(stage_forcing, initial_forcing, state_size * sizeof(double));",
        "        memcpy(k0, initial_forcing, state_size * sizeof(double));",
        "        generated_solve(cell_count, matrix, k0);",
        "        if (statistics) ++statistics->solves;",
        "        for (size_t i = 0; i < state_size; ++i) ynew[i] = y[i] + k0[i];",
        "        generated_forcing(cell_count, ynew, rates, stage_forcing);",
        "        if (statistics) ++statistics->function_calls;",
        "        memcpy(k1, stage_forcing, state_size * sizeof(double));",
        "        for (size_t i = 0; i < state_size; ++i) k1[i] += (-1.0156171083877702091975600115545 / h) * k0[i];",
        "        generated_solve(cell_count, matrix, k1);",
        "        if (statistics) ++statistics->solves;",
        "        memcpy(k2, stage_forcing, state_size * sizeof(double));",
        "        for (size_t i = 0; i < state_size; ++i) {",
        "          k2[i] += (4.0759956452537699824805835358067 / h) * k0[i];",
        "          k2[i] += (9.2076794298330791242156818474003 / h) * k1[i];",
        "        }",
        "        generated_solve(cell_count, matrix, k2);",
        "        if (statistics) ++statistics->solves;",
        "        for (size_t i = 0; i < state_size; ++i) {",
        "          ynew[i] = y[i] + k0[i];",
        "          ynew[i] += 6.1697947043828245592553615689730 * k1[i];",
        "          ynew[i] += -0.42772256543218573326238373806514 * k2[i];",
        "          error[i] = 0.5 * k0[i];",
        "          error[i] += -2.907955871680469821793236208017 * k1[i];",
        "          error[i] += 0.22354069897811569627360909276199 * k2[i];",
        "        }",
        "        const double normalized_error = generated_error(cell_count, y, ynew, error, relative_tolerance, absolute_tolerance);",
        "        const double factor = fmin(6.0, fmax(0.2, 0.9 / pow(normalized_error, 1.0 / 3.0)));",
        "        double hnew = h * factor;",
        "        if (statistics) ++statistics->internal_steps;",
        "        ++step_attempts;",
        "        if (!isfinite(normalized_error)) {",
        "          fprintf(stderr, \"non-finite normalized error %.17g at time %.17g h %.17g\\n\", normalized_error, present_time, h);",
        "          failed = 1; goto cleanup;",
        "        }",
        "        if (normalized_error < 1.0 || h < h_min) {",
        "          ++accepted_this_interval; if (statistics) ++statistics->accepted;",
        "          present_time += h; memcpy(y, ynew, state_size * sizeof(double));",
        "          hnew = fmax(h_min, fmin(hnew, h_max));",
        "          if (reject_last_h) hnew = fmin(hnew, h);",
        "          reject_last_h = 0; reject_more_h = 0; h = hnew; accepted = 1;",
        "        } else {",
        "          if (reject_more_h) hnew = h * 0.1;",
        "          reject_more_h = reject_last_h; reject_last_h = 1; h = hnew;",
        "          if (statistics && accepted_this_interval >= 1) ++statistics->rejected;",
        "          generated_jacobian(cell_count, y, rates, matrix);",
        "          if (statistics) ++statistics->jacobian_updates;",
        "        }",
        "      }",
        "    }",
        "    for (size_t s = 0; s < GENERATED_SPECIES; ++s) {",
        "      const double atol = absolute_tolerance >= 0.0 ? absolute_tolerance : generated_sorted_atol[s];",
        "      double *row = y + s * cell_count;",
        "      #pragma omp simd",
        "      for (size_t c = 0; c < cell_count; ++c)",
        "        if (row[c] < 0.0 && row[c] >= -atol) row[c] = 0.0;",
        "    }",
        "  }",
        "  for (size_t c = 0; c < cell_count; ++c)",
        "    for (size_t s = 0; s < GENERATED_SPECIES; ++s)",
        "      y_external[c * GENERATED_SPECIES + generated_permutation[s]] = y[s * cell_count + c];",
        "cleanup:",
        "  free(y); free(ynew); free(initial_forcing); free(stage_forcing); free(error);",
        "  free(parameters); free(rates); free(matrix); free(inverse); free(k0); free(k1); free(k2);",
        "  return failed ? -1 : 0;",
        "}",
        "",
    ])
    return lines


def emit(config: dict, output: Path) -> None:
    all_species = config["species"]
    integrated = [s for s in all_species if not s.get("third_body", False)]
    species_names = {s["name"]: i for i, s in enumerate(integrated)}
    third_bodies = {s["name"] for s in all_species if s.get("third_body", False)}
    reactions = list(config["reactions"])
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
    lines.extend(specialized_source(config))
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
