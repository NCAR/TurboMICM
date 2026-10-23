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


def sparse_program(species_count: int, reactions: list[dict],
                   reactants: list[int], products: list[int]):
    """Return the fixed CSR pattern, position map, and elimination work lists."""
    rows = [{i} for i in range(species_count)]
    for reaction in reactions:
        reactant_ids = reactants[reaction["rbegin"]:reaction["rbegin"] + reaction["nr"]]
        product_ids = products[reaction["pbegin"]:reaction["pbegin"] + reaction["np"]]
        for column in reactant_ids:
            for row in reactant_ids + product_ids:
                rows[row].add(column)
    for pivot in range(species_count):
        right = [column for column in rows[pivot] if column > pivot]
        for row in range(pivot + 1, species_count):
            if pivot in rows[row]:
                rows[row].update(right)
    columns = [sorted(row) for row in rows]
    offsets = []
    position = [[-1] * species_count for _ in range(species_count)]
    cursor = 0
    for row, row_columns in enumerate(columns):
        offsets.append(cursor)
        for column in row_columns:
            position[row][column] = cursor
            cursor += 1
    offsets.append(cursor)
    lower = [[row for row in range(k + 1, species_count) if k in rows[row]]
             for k in range(species_count)]
    upper = [[column for column in columns[k] if column > k]
             for k in range(species_count)]
    return columns, offsets, position, lower, upper


def emit_rate_body(reaction: dict, index: int, species: list[dict],
                   coefficients: list[float]) -> list[str]:
    """Emit one mechanism-specific rate formula, with fixed inputs literal."""
    kind = reaction["kind"]
    p = reaction["params"]
    values = [number(value) for value in p]
    body = reaction["body"]
    factor = f" * pow(air, {body}.0)" if body else ""
    out = ["  {"]
    rate = f"rate_{index}"
    if kind == "CHEM_ARRHENIUS":
        formula = (f"{values[0]} * exp({values[2]} / temperature) * "
                   f"pow(temperature / {values[3]}, {values[1]}) * "
                   f"(1.0 + {values[4]} * pressure)")
    elif kind == "CHEM_TROE":
        formula = f"generated_falloff({', '.join(values[:8])}, temperature, air, air)"
    elif kind == "CHEM_TERNARY":
        formula = f"generated_falloff({', '.join(values[:8])}, temperature, air, 1.0)"
    elif kind == "CHEM_TERNARY_JPL19":
        out += [
            f"    const double {rate}_k0 = {values[0]} * exp({values[2]} / temperature) * pow(temperature / {values[3]}, {values[1]});",
            f"    const double {rate}_kinf = {values[4]} * exp({values[6]} / temperature) * pow(temperature / {values[7]}, {values[5]});",
            f"    const double {rate}_kint = {values[8]} * exp({values[10]} / temperature) * pow(temperature / {values[11]}, {values[9]});",
            f"    const double {rate}_k0_air = {rate}_k0 * air;",
            f"    const double {rate}_ratio = {rate}_k0_air / {rate}_kinf;",
            f"    const double {rate}_kf = {rate}_k0_air / (1.0 + {rate}_ratio) * pow({values[12]}, {values[13]} / ({values[13]} + pow(log10({rate}_ratio), 2.0)));",
            f"    const double {rate} = {rate}_kf + {rate}_kint * (1.0 - {rate}_kf / {rate}_kinf);",
            f"    rate_constants[{index}] = {rate}{factor};",
            "  }",
        ]
        return out
    elif kind == "CHEM_BRANCHED":
        out += [
            f"    const double {rate}_k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp({values[3]});",
            f"    const double {rate}_air_ref = 2.45e19 / 6.02214076e23 * 1.0e6;",
            f"    const double {rate}_a0 = {rate}_k0 * {rate}_air_ref;",
            f"    const double {rate}_b0 = 0.43 * pow(293.0 / 298.0, -8.0);",
            f"    const double {rate}_aval0 = {rate}_a0 / (1.0 + {rate}_a0 / {rate}_b0) * pow(0.41, 1.0 / (1.0 + pow(log10({rate}_a0 / {rate}_b0), 2.0)));",
            f"    const double {rate}_z = {rate}_aval0 * (1.0 - {values[2]}) / {values[2]};",
            f"    const double {rate}_a = {rate}_k0 * air;",
            f"    const double {rate}_b = 0.43 * pow(temperature / 298.0, -8.0);",
            f"    const double {rate}_aval = {rate}_a / (1.0 + {rate}_a / {rate}_b) * pow(0.41, 1.0 / (1.0 + pow(log10({rate}_a / {rate}_b), 2.0)));",
            f"    const double {rate}_pre = {values[1]} * exp(-{values[0]} / temperature);",
        ]
        branch = "z" if p[4] == 0.0 else "aval"
        out += [
            f"    const double {rate} = {rate}_pre * {rate}_{branch} / ({rate}_z + {rate}_aval);",
            f"    rate_constants[{index}] = {rate}{factor};",
            "  }",
        ]
        return out
    elif kind == "CHEM_TUNNELING":
        formula = (f"{values[0]} * exp(-{values[1]} / temperature + "
                   f"{values[2]} / (temperature * temperature * temperature))")
    elif kind == "CHEM_TAYLOR":
        polynomial = f"{rate}_polynomial"
        t_power = f"{rate}_t_power"
        out.append(f"    double {polynomial} = 0.0, {t_power} = 1.0;")
        series = coefficients[reaction["coeff_offset"]:
                              reaction["coeff_offset"] + reaction["coeff_count"]]
        for coefficient in series:
            out += [f"    {polynomial} += {number(coefficient)} * {t_power};",
                    f"    {t_power} *= temperature;"]
        formula = (f"{polynomial} * {values[0]} * exp({values[2]} / temperature) * "
                   f"pow(temperature / {values[3]}, {values[1]}) * "
                   f"(1.0 + {values[4]} * pressure)")
    elif kind == "CHEM_REVERSIBLE":
        formula = f"{values[0]} * exp({values[1]} / temperature) * {values[2]}"
    elif kind == "CHEM_USER_DEFINED":
        formula = f"custom_parameters[{reaction['custom0']}] * {values[0]}"
    elif kind == "CHEM_SURFACE":
        surface = species[reaction["surface"]]
        mw = surface.get("molecular_weight", surface.get("molecular weight [kg mol-1]", 0.0)) or 0.0
        diffusion = surface.get("diffusion_coefficient", surface.get("diffusion coefficient [m2 s-1]")) or 0.0
        pi = "3.141592653589793238462643383279502884"
        speed_factor = f"8.0 * 1.380649e-23 * 6.02214076e23 / ({pi} * {number(mw)})"
        formula = (
            f"4.0 * custom_parameters[{reaction['custom1']}] * {pi} * "
            f"custom_parameters[{reaction['custom0']}] * custom_parameters[{reaction['custom0']}] / "
            f"(custom_parameters[{reaction['custom0']}] / {number(diffusion)} + "
            f"4.0 / (sqrt(({speed_factor}) * temperature) * {values[0]}))"
        )
    else:
        raise SystemExit(f"cannot specialize generated rate kind {kind}")
    out += [f"    rate_constants[{index}] = ({formula}){factor};", "  }"]
    return out


def emit_specialized_kernels(species: list[dict], reactions: list[dict],
                             reactants: list[int], products: list[int],
                             yields: list[float], coefficients: list[float],
                             enabled: set[str]) -> str:
    """Emit fixed rate, RHS/Jacobian, no-pivot LU, and triangular solves."""
    n = len(species)
    columns, offsets, position, lower, upper = sparse_program(
        n, reactions, reactants, products
    )
    nnz = offsets[-1]
    lines = []
    if any(reaction["kind"] in {"CHEM_TROE", "CHEM_TERNARY"} for reaction in reactions):
        lines += [
            "static double generated_falloff(double k0_a, double k0_b, double k0_c,",
            "                                  double kinf_a, double kinf_b, double kinf_c,",
            "                                  double fc, double exponent, double temperature,",
            "                                  double air_density, double numerator_scale)",
            "{",
            "  const double k0 = k0_a * exp(k0_c / temperature) * pow(temperature / 300.0, k0_b);",
            "  const double kinf = kinf_a * exp(kinf_c / temperature) * pow(temperature / 300.0, kinf_b);",
            "  const double ratio = k0 * air_density / kinf;",
            "  return k0 * numerator_scale / (1.0 + ratio) *",
            "         pow(fc, exponent / (exponent + pow(log10(ratio), 2.0)));",
            "}",
            "",
        ]
    lines += [
        "static void generated_rates(double temperature, double pressure,",
        "                            const double *custom_parameters, double *rate_constants)",
        "{",
        "  const double air = pressure / (1.380649e-23 * 6.02214076e23 * temperature);",
        "  (void)air; (void)custom_parameters;",
    ]
    for i, reaction in enumerate(reactions):
        lines.extend(emit_rate_body(reaction, i, species, coefficients))
    lines += [
        "}", "", "static void generated_custom_rates(double temperature, double pressure,",
        "                                   const double *custom_parameters, double *rate_constants)",
        "{",
        "  const double air = pressure / (1.380649e-23 * 6.02214076e23 * temperature);",
        "  (void)air; (void)custom_parameters;",
    ]
    for i, reaction in enumerate(reactions):
        if reaction["kind"] == "CHEM_USER_DEFINED":
            value = number(reaction["params"][0])
            factor = f" * pow(air, {reaction['body']}.0)" if reaction["body"] else ""
            lines.append(
                f"  rate_constants[{i}] = custom_parameters[{reaction['custom0']}] * {value}{factor};"
            )
    if not any(reaction["kind"] == "CHEM_USER_DEFINED" for reaction in reactions):
        lines.append("  (void)rate_constants;")
    lines += [
        "}", "", "static int generated_evaluate(const double *concentrations,",
        "                              const double *rate_constants, double *forcing,",
        "                              double *negative_jacobian,",
        "                              const chem_sparse_pattern *pattern)", "{",
        f"  if (!pattern || pattern->n != {n} || pattern->nnz != {nnz}) return -1;",
        f"  for (size_t i = 0; i < {n}; ++i) forcing[i] = 0.0;",
        f"  if (negative_jacobian) for (size_t i = 0; i < {nnz}; ++i) negative_jacobian[i] = 0.0;",
    ]
    for i, reaction in enumerate(reactions):
        rids = reactants[reaction["rbegin"]:reaction["rbegin"] + reaction["nr"]]
        pids = products[reaction["pbegin"]:reaction["pbegin"] + reaction["np"]]
        lines += ["  {", f"    double flux = rate_constants[{i}];"]
        for species_id in rids:
            lines.append(f"    flux *= concentrations[{species_id}];")
        for species_id in rids:
            lines.append(f"    forcing[{species_id}] -= flux;")
        for product_offset, species_id in enumerate(pids):
            y = number(yields[reaction["pbegin"] + product_offset])
            lines.append(f"    forcing[{species_id}] += {y} * flux;")
        if rids:
            lines.append("    if (negative_jacobian) {")
            for a, column in enumerate(rids):
                derivative = f"derivative_{i}_{a}"
                lines.append(f"      double {derivative} = rate_constants[{i}];")
                for b, other in enumerate(rids):
                    if b != a:
                        lines.append(f"      {derivative} *= concentrations[{other}];")
                for row in rids:
                    pos = position[row][column]
                    if pos < 0:
                        raise AssertionError("reaction dependency missing from generated fill pattern")
                    lines.append(f"      negative_jacobian[{pos}] += {derivative};")
                for product_offset, row in enumerate(pids):
                    pos = position[row][column]
                    if pos < 0:
                        raise AssertionError("reaction dependency missing from generated fill pattern")
                    y = number(yields[reaction["pbegin"] + product_offset])
                    lines.append(f"      negative_jacobian[{pos}] -= {y} * {derivative};")
            lines.append("    }")
        lines.append("  }")
    lines += ["  return 0;", "}", "", "static int generated_factor(const chem_sparse_pattern *pattern, double alpha,",
              "                            const double *negative_jacobian, double *lu, void *context)", "{",
              "  (void)context;",
              f"  if (!pattern || pattern->n != {n} || pattern->nnz != {nnz}) return -1;",
              f"  memcpy(lu, negative_jacobian, {nnz} * sizeof(*lu));"]
    for i in range(n):
        lines.append(f"  lu[{position[i][i]}] += alpha;")
    for k in range(n):
        pivot = f"pivot_{k}"
        lines += [f"  {{ const double {pivot} = lu[{position[k][k]}];",
                  f"    if ({pivot} == 0.0 || !isfinite({pivot})) return -1;"]
        for row in lower[k]:
            multiplier = f"multiplier_{k}_{row}"
            lines += [f"    const double {multiplier} = lu[{position[row][k]}] / {pivot};",
                      f"    lu[{position[row][k]}] = {multiplier};"]
            for column in upper[k]:
                lines.append(
                    f"    lu[{position[row][column]}] = fma(-{multiplier}, lu[{position[k][column]}], lu[{position[row][column]}]);"
                )
        lines.append("  }")
    lines += ["  return 0;", "}", "", "static void generated_solve(const chem_sparse_pattern *pattern,",
              "                           const double *lu, double *right_hand_side, void *context)", "{",
              "  (void)context;",
              f"  (void)pattern; /* fixed pattern: n={n}, nnz={nnz} */"]
    for row, row_columns in enumerate(columns):
        lower_columns = [column for column in row_columns if column < row]
        lines += ["  {", f"    double value = right_hand_side[{row}];"]
        for column in lower_columns:
            lines.append(
                f"    value = fma(-lu[{position[row][column]}], right_hand_side[{column}], value);"
            )
        lines += [f"    right_hand_side[{row}] = value;", "  }"]
    for row in range(n - 1, -1, -1):
        upper_columns = [column for column in columns[row] if column > row]
        lines += ["  {", f"    double value = right_hand_side[{row}];"]
        for column in upper_columns:
            lines.append(
                f"    value = fma(-lu[{position[row][column]}], right_hand_side[{column}], value);"
            )
        lines += [f"    right_hand_side[{row}] = value / lu[{position[row][row]}];", "  }"]
    lines += [
        "}", "", "static const chem_specialized_kernels chem_generated_kernels = {",
        f"  .species_count = {n}, .pattern_nonzeros = {nnz},",
        f"  .rates = {'generated_rates' if 'rates' in enabled else 'NULL'},",
        f"  .custom_rates = {'generated_custom_rates' if 'rates' in enabled else 'NULL'},",
        f"  .evaluate = {'generated_evaluate' if 'chemistry' in enabled else 'NULL'},",
        f"  .factor = {'generated_factor' if 'factor' in enabled else 'NULL'},",
        f"  .solve = {'generated_solve' if 'solve' in enabled else 'NULL'}",
        "};", "",
    ]
    return "\n".join(lines)


def compile_model(document: dict, specialized: bool = False,
                  disabled_kernels: set[str] | None = None) -> str:
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
    if specialized:
        output[1:1] = ["#include <math.h>", "#include <string.h>", ""]
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
    output += ["};", ""]
    if specialized:
        output.append(emit_specialized_kernels(
            state_species, reactions, flat_reactants, flat_products,
            flat_yields, coefficients,
            {"rates", "chemistry", "factor", "solve"} - (disabled_kernels or set())
        ))
    output.append("const chem_mechanism chem_model = {")
    output += [
        f"  .name={c_string(document.get('name', 'mechanism'))},",
        f"  .species_count={len(state_species)}, .reaction_count={len(reactions)},",
        f"  .reactant_count={len(flat_reactants)}, .product_count={len(flat_products)}, .custom_count={len(custom)},",
        f"  .default_temperature={number(document.get('temperature', 298.15))},",
        f"  .default_pressure={number(document.get('pressure', 101325.0))},",
        "  .species=chem_species_data, .reactions=chem_reaction_data,",
        "  .reactants=chem_reactant_data, .products=chem_product_data, .yields=chem_yield_data,",
        "  .coefficients=chem_coefficient_data, .custom_parameters=chem_custom_data,",
        f"  .kernels={'&chem_generated_kernels' if specialized else 'NULL'}",
        "};",
        "",
    ]
    return "\n".join(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--specialized", action="store_true",
                        help="emit fixed rate, chemistry, LU, and solve kernels")
    parser.add_argument("--without-rate-kernel", action="store_true",
                        help="ablation: use generic rate constants")
    parser.add_argument("--without-chemistry-kernel", action="store_true",
                        help="ablation: use the generic reaction and Jacobian loops")
    parser.add_argument("--without-factor-kernel", action="store_true",
                        help="ablation: use the generic sparse LU factorization")
    parser.add_argument("--without-solve-kernel", action="store_true",
                        help="ablation: use the generic sparse triangular solves")
    args = parser.parse_args()
    disabled = set()
    if args.without_rate_kernel:
        disabled.add("rates")
    if args.without_chemistry_kernel:
        disabled.add("chemistry")
    if args.without_factor_kernel:
        disabled.add("factor")
    if args.without_solve_kernel:
        disabled.add("solve")
    if disabled and not args.specialized:
        parser.error("--without-* ablations require --specialized")
    document = read_config(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(compile_model(document, specialized=args.specialized,
                                         disabled_kernels=disabled))


if __name__ == "__main__":
    main()
