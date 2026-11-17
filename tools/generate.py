#!/usr/bin/env python3
"""Turn a CAMP-style YAML or JSON mechanism into a C11 data definition."""

import argparse
import json
import re
from pathlib import Path

from config import normalize, number, read_config


def reorder_mechanism(mechanism):
    """Symmetric minimum-fill ordering, derived from any reaction graph.

    The diagonal is included even for product-only species. The returned
    ordering changes only generated indices; names and numerical data follow
    their species, and reactions are remapped consistently.
    """
    n = len(mechanism["species"])
    rows = [{i} for i in range(n)]
    for reaction in mechanism["reactions"]:
        for col in reaction["reactants"]:
            for row in reaction["reactants"]:
                rows[row].add(col)
            for row, _ in reaction["products"]:
                rows[row].add(col)
    alive = set(range(n))
    order = []
    while alive:
        columns = {j: set() for j in alive}
        for i in alive:
            for j in rows[i] & alive:
                columns[j].add(i)
        pivot = min(alive, key=lambda i: ((len(rows[i] & alive)-1)*(len(columns[i])-1), i))
        order.append(pivot)
        outgoing = (rows[pivot] & alive) - {pivot}
        incoming = columns[pivot] - {pivot}
        for i in incoming:
            rows[i].update(outgoing)
        alive.remove(pivot)
    inverse = {old: new for new, old in enumerate(order)}
    mechanism["species"] = [mechanism["species"][i] for i in order]
    for reaction in mechanism["reactions"]:
        reaction["reactants"] = [inverse[i] for i in reaction["reactants"]]
        reaction["products"] = [(inverse[i], amount) for i, amount in reaction["products"]]
    return mechanism


def c_number(value):
    text = format(number(value, "generated number"), ".17g")
    return text if any(char in text for char in ".eE") else text + ".0"


def c_string(value):
    return json.dumps(value, ensure_ascii=True)


def fields(reaction):
    kind, raw = reaction["kind"], reaction["raw"]
    mapped = {"PHOTOLYSIS": "USER", "USER_DEFINED": "USER", "FIRST_ORDER_LOSS": "USER", "EMISSION": "USER",
              "TERNARY_CHEMICAL_ACTIVATION": "TERNARY", "TERNARY_CHEMICAL_ACTIVATION_JPL19": "TERNARY_JPL19",
              "BRANCHED_ALKOXY": "BRANCHED", "BRANCHED_NITRATE": "BRANCHED",
              "TUNNELING": "TUNNELING", "TAYLOR_SERIES": "TAYLOR", "REVERSIBLE": "REVERSIBLE",
              "TROE": "TROE", "ARRHENIUS": "ARRHENIUS", "SURFACE": "SURFACE"}
    if kind not in mapped:
        raise ValueError(f"unsupported reaction type {kind!r}")
    output = [f".kind = CHEM_{mapped[kind]}", f".parameter = {reaction['parameter']}"]

    order = len(reaction["reactants"]) + reaction["air_power"]
    molecule_factor = 1.0e-6 * 6.02214076e23 if reaction["molecule_units"] else 1.0

    def add(field, key, default, conversion_power=0):
        value = number(raw.get(key, default), key) * molecule_factor**conversion_power
        output.append(f".{field} = {c_number(value)}")

    if kind in ("ARRHENIUS", "TAYLOR_SERIES"):
        for field, key, default in [("a", "A", 1), ("b", "B", 0), ("c", "C", 0), ("d", "D", 300), ("e", "E", 0)]:
            add(field, key, default, order - 1 if key == "A" else 0)
        if kind == "TAYLOR_SERIES":
            coeff = raw.get("coefficients", [1])
            if not 0 < len(coeff) <= 16:
                raise ValueError("Taylor series must have between 1 and 16 coefficients")
            output.extend([f".coefficient_count = {len(coeff)}", ".coefficients = {" + ", ".join(map(c_number, coeff)) + "}"])
    elif kind in ("TROE", "TERNARY_CHEMICAL_ACTIVATION", "TERNARY_CHEMICAL_ACTIVATION_JPL19"):
        for group in ("k0", "kinf", "kint"):
            if group == "kint" and kind != "TERNARY_CHEMICAL_ACTIVATION_JPL19":
                continue
            for letter, default in (("A", 1), ("B", 0), ("C", 0)):
                conversion_power = 0
                if letter == "A":
                    if kind == "TROE":
                        conversion_power = order if group == "k0" else order - 1
                    elif kind == "TERNARY_CHEMICAL_ACTIVATION":
                        conversion_power = order - 1 if group == "k0" else order - 2
                    else:
                        conversion_power = order if group == "k0" else order - 1
                add(f"{group}_{letter.lower()}", f"{group}_{letter}", default, conversion_power)
            if kind == "TERNARY_CHEMICAL_ACTIVATION_JPL19":
                add(f"{group}_d", f"{group}_D", 298)
        add("fc", "Fc", 0.6)
        add("n", "N", 1)
    elif kind.startswith("BRANCHED_"):
        for field, key in (("branch_x", "X"), ("branch_y", "Y"), ("branch_a0", "a0"), ("branch_n", "n")):
            add(field, key, 0, order - 1 if key == "X" else 0)
        output.append(f".nitrate_branch = {int(kind == 'BRANCHED_NITRATE')}")
    elif kind == "TUNNELING":
        for field, key, default in (("a", "A", 1), ("b", "B", 0), ("c", "C", 0)):
            add(field, key, default, order - 1 if key == "A" else 0)
    elif kind == "REVERSIBLE":
        for field, key, default in (("a", "A", 1), ("c", "C", 0), ("b", "k_r", 0)):
            add(field, key, default)
    elif mapped[kind] == "USER":
        add("a", "scaling factor", 1)
    elif kind == "SURFACE":
        surface = reaction["surface_species"]
        output.extend([f".diffusion = {c_number(surface['diffusion'])}",
                       f".molecular_weight = {c_number(surface['molecular_weight'])}"])
        add("probability", "reaction probability", 1)
    return "{ " + ", ".join(output) + " }"


def generate(mechanism, symbol):
    if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", symbol):
        raise ValueError("symbol must be a C identifier")
    lines = ['#include "chem.h"', "", "/* Generated mechanism data; do not edit. */"]
    lines.append("static const ChemSpecies species[] = {")
    for s in mechanism["species"]:
        lines.append("    { " + ", ".join([c_string(s["name"])] + [c_number(s[k]) for k in ("initial", "absolute_tolerance", "molecular_weight", "diffusion")]) + " },")
    lines.append("};")
    reactants, products = [], []
    lines.append("static const ChemReaction reactions[] = {")
    for r in mechanism["reactions"]:
        lines.append(f"    {{ {len(reactants)}, {len(r['reactants'])}, {len(products)}, {len(r['products'])}, {r['air_power']}, {fields(r)} }},")
        reactants.extend(r["reactants"])
        products.extend(r["products"])
    if not mechanism["reactions"]:
        lines.append("    {0},")
    lines.append("};")
    lines.append("static const int reactants[] = { " + (", ".join(map(str, reactants)) or "0") + " };")
    lines.append("static const ChemProduct products[] = { " + (", ".join("{ %d, %s }" % (i, c_number(y)) for i, y in products) or "{0, 0.0}") + " };")
    params = mechanism["parameters"]
    lines.append("static const char *const parameters[] = { " + (", ".join(c_string(n) for n, _ in params) or '""') + " };")
    lines.append("static const double defaults[] = { " + (", ".join(c_number(v) for _, v in params) or "0.0") + " };")
    lines.append(f"const ChemMechanism {symbol} = {{ {len(mechanism['species'])}, {len(mechanism['reactions'])}, {len(params)}, species, reactions, reactants, products, parameters, defaults, {c_number(mechanism['temperature'])}, {c_number(mechanism['pressure'])} }};")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--symbol", default="chem_generated_mechanism")
    parser.add_argument("--reorder", action="store_true", help="apply a symmetric minimum-fill species ordering")
    args = parser.parse_args()
    mechanism = normalize(read_config(args.config))
    if args.reorder:
        mechanism = reorder_mechanism(mechanism)
    result = generate(mechanism, args.symbol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
