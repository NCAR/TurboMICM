#!/usr/bin/env python3
"""Generate immutable C mechanism tables from normalized JSON."""
import argparse
import json
import math
from pathlib import Path

KINDS = {
    "ARRHENIUS": "CHEM_ARRHENIUS", "arrhenius": "CHEM_ARRHENIUS",
    "TROE": "CHEM_TROE", "troe": "CHEM_TROE",
    "PHOTOLYSIS": "CHEM_PHOTOLYSIS", "photolysis": "CHEM_PHOTOLYSIS",
    "USER_DEFINED": "CHEM_USER_DEFINED", "user_defined": "CHEM_USER_DEFINED",
    "SURFACE": "CHEM_SURFACE", "surface": "CHEM_SURFACE",
    "TERNARY_CHEMICAL_ACTIVATION": "CHEM_TERNARY_ACTIVATION",
    "ternary_chemical_activation": "CHEM_TERNARY_ACTIVATION",
    "TERNARY_CHEMICAL_ACTIVATION_JPL19": "CHEM_TERNARY_ACTIVATION_JPL19",
    "ternary_chemical_activation_jpl19": "CHEM_TERNARY_ACTIVATION_JPL19",
    "TAYLOR_SERIES": "CHEM_TAYLOR_SERIES", "taylor_series": "CHEM_TAYLOR_SERIES",
    "REVERSIBLE": "CHEM_REVERSIBLE", "reversible": "CHEM_REVERSIBLE",
    "TUNNELING": "CHEM_TUNNELING", "tunneling": "CHEM_TUNNELING",
    "BRANCHED_ALKOXY": "CHEM_BRANCHED_ALKOXY", "branched_alkoxy": "CHEM_BRANCHED_ALKOXY",
    "BRANCHED_NITRATE": "CHEM_BRANCHED_NITRATE", "branched_nitrate": "CHEM_BRANCHED_NITRATE",
}
TYPE_ORDER = {
    "arrhenius": 0, "ARRHENIUS": 0,
    "troe": 1, "TROE": 1,
    "ternary_chemical_activation": 2, "TERNARY_CHEMICAL_ACTIVATION": 2,
    "ternary_chemical_activation_jpl19": 3, "TERNARY_CHEMICAL_ACTIVATION_JPL19": 3,
    "branched_alkoxy": 4, "BRANCHED_ALKOXY": 4,
    "branched_nitrate": 4, "BRANCHED_NITRATE": 4,
    "tunneling": 5, "TUNNELING": 5,
    "taylor_series": 6, "TAYLOR_SERIES": 6,
    "reversible": 7, "REVERSIBLE": 7,
    "user_defined": 8, "USER_DEFINED": 8,
    "photolysis": 8, "PHOTOLYSIS": 8,
    "surface": 9, "SURFACE": 9,
}


def markowitz_order(species, reactions):
    """Return the solver's diagonal Markowitz order for dynamic species."""
    dynamic = [sp for sp in species if not sp.get("third_body", False)]
    names = {sp["name"]: i for i, sp in enumerate(dynamic)}
    n = len(dynamic)
    outgoing = [set() for _ in range(n)]
    incoming = [set() for _ in range(n)]
    for reaction in reactions:
        reactants = [names[t["name"]] for t in reaction.get("reactants", [])
                     if t["name"] in names]
        products = [names[t["name"]] for t in reaction.get("products", [])
                    if t["name"] in names]
        for col in reactants:
            for row in reactants + products:
                outgoing[row].add(col)
                incoming[col].add(row)

    row_degree = [len(x) for x in outgoing]
    col_degree = [len(x) for x in incoming]
    alive = [True] * n
    order = []
    for _ in range(n):
        def cost(i):
            # MICM's Index is size_t; preserve its unsigned underflow for
            # isolated or product-only variables in the sparsity graph.
            mask = (1 << 64) - 1
            return (((row_degree[i]-1) & mask) *
                    ((col_degree[i]-1) & mask)) & mask
        pivot = min((i for i in range(n) if alive[i]), key=cost)
        order.append(pivot)
        alive[pivot] = False
        cols = [c for c in outgoing[pivot] if c != pivot and alive[c]]
        rows = [r for r in incoming[pivot] if r != pivot and alive[r]]
        for row in rows:
            for col in cols:
                if col not in outgoing[row]:
                    outgoing[row].add(col)
                    row_degree[row] += 1
                    incoming[col].add(row)
                    col_degree[col] += 1
        for col in cols:
            incoming[col].discard(pivot)
            col_degree[col] -= 1
        for row in rows:
            outgoing[row].discard(pivot)
            row_degree[row] -= 1
    return [dynamic[i] for i in order] + [sp for sp in species if sp.get("third_body", False)]


def cstr(value):
    return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    data = json.loads(args.input.read_text())
    reactions = sorted(data["reactions"], key=lambda r: TYPE_ORDER[r["type"]])
    species = markowitz_order(data["species"], reactions)
    defaults = data.get("defaults", {})
    names = {s["name"]: i for i, s in enumerate(species)}
    if len(names) != len(species):
        raise SystemExit("duplicate species names")

    custom_labels = []
    for reaction in reactions:
        params = reaction.get("parameters", {})
        label = params.get("label", params.get("MUSICA name"))
        if label is not None:
            labels = [label + ".effective radius [m]", label + ".particle number concentration [# m-3]"] if (
                str(reaction.get("type", "")).lower() == "surface" or
                str(reaction.get("type", "")).upper() == "SURFACE") else [label]
            for item in labels:
                if item not in custom_labels:
                    custom_labels.append(item)
    label_index = {label: i for i, label in enumerate(custom_labels)}

    out = ['#include "chemistry.h"', '']
    for ri, reaction in enumerate(reactions):
        for key, label in (("reactants", "r"), ("products", "p")):
            terms = reaction.get(key, [])
            if terms:
                out.append(f"static const chem_term rxn_{ri}_{label}[] = {{")
                for term in terms:
                    index = term.get("index", names.get(term["name"]))
                    out.append(f"  {{{index}, {float(term.get('coefficient', 1.0)):.17g}}},")
                out.append("};")
    out.append("")
    out.append("static const chem_species species[] = {")
    for sp in species:
        out.append("  {%s, %d, %.17g, %.17g, %.17g}," % (
            cstr(sp["name"]), int(bool(sp.get("third_body", False))),
            float(sp.get("molecular_weight", 0.0)), float(sp.get("diffusion", 0.0)),
            float(defaults.get("species", {}).get(sp["name"], 0.0))))
    out.append("};")
    out.append("static const chem_reaction reactions[] = {")
    for i, rxn in enumerate(reactions):
        kind = KINDS.get(rxn.get("type"))
        if kind is None:
            raise SystemExit(f"unsupported rate kind {rxn.get('type')!r} in reaction {i}")
        params = rxn.get("parameters", {})
        if isinstance(params, dict) and "parameters" in params:
            params = params["parameters"]
        values = [0.0] * 24
        if kind == "CHEM_ARRHENIUS":
            for j, key in enumerate(("A", "B", "C", "D", "E")):
                values[j] = float(params.get(key, {"A": 1.0, "D": 300.0}.get(key, 0.0)))
        elif kind == "CHEM_TROE":
            for j, key in enumerate(("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")):
                values[j] = float(params.get(key, {"k0_A": 1.0, "kinf_A": 1.0,
                                                   "Fc": 0.6, "N": 1.0}.get(key, 0.0)))
        elif kind == "CHEM_TERNARY_ACTIVATION":
            for j, key in enumerate(("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")):
                values[j] = float(params.get(key, {"k0_A": 1.0, "kinf_A": 1.0,
                                                   "Fc": 0.6, "N": 1.0}.get(key, 0.0)))
        elif kind == "CHEM_TERNARY_ACTIVATION_JPL19":
            for j, key in enumerate(("k0_A", "k0_B", "k0_C", "k0_D",
                                     "kinf_A", "kinf_B", "kinf_C", "kinf_D",
                                     "kint_A", "kint_B", "kint_C", "kint_D", "Fc", "N")):
                values[j] = float(params.get(key, {
                    "k0_A": 1.0, "k0_D": 298.0, "kinf_A": 1.0,
                    "kinf_D": 298.0, "kint_A": 1.0, "kint_D": 298.0,
                    "Fc": 0.6, "N": 1.0}.get(key, 0.0)))
        elif kind == "CHEM_TAYLOR_SERIES":
            for j, key in enumerate(("A", "B", "C", "D", "E")):
                values[j] = float(params.get(key, {"A": 1.0, "D": 300.0}.get(key, 0.0)))
            coeffs = params.get("coefficients", [1.0])
            values[5] = min(len(coeffs), 16)
            for j, value in enumerate(coeffs[:16]): values[6+j] = float(value)
        elif kind == "CHEM_REVERSIBLE":
            values[0] = float(params.get("A", 1.0))
            values[1] = float(params.get("C", 0.0))
            values[2] = float(params.get("k_r", 0.0))
        elif kind == "CHEM_TUNNELING":
            for j, key in enumerate(("A", "B", "C")):
                values[j] = float(params.get(key, 1.0 if j == 0 else 0.0))
        elif kind in ("CHEM_BRANCHED_ALKOXY", "CHEM_BRANCHED_NITRATE"):
            values[0] = float(params.get("X", 0.0))
            values[1] = float(params.get("Y", 0.0))
            values[2] = float(params.get("a0", 0.0))
            values[3] = float(params.get("n", 0.0))
            k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * math.exp(values[3])
            a = k0 * (2.45e19 / 6.02214076e23 * 1.0e6)
            b = 0.43 * (293.0 / 298.0) ** -8.0
            av = a / (1.0 + a / b) * 0.41 ** (1.0 / (1.0 + math.log10(a / b) ** 2))
            values[4] = av * (1.0 - values[2]) / values[2]
        else:
            label = params.get("label", params.get("MUSICA name"))
            if label is not None:
                lookup = label + ".effective radius [m]" if kind == "CHEM_SURFACE" else label
                values[7] = float(label_index[lookup])
            if kind == "CHEM_SURFACE":
                values[0] = float(params.get("reaction_probability", params.get("reaction probability", 1.0)))
            else:
                values[0] = float(params.get("scaling_factor", 1.0))
        name = params.get("label", params.get("MUSICA name", ""))
        reactants = rxn.get("reactants", [])
        products = rxn.get("products", [])
        if reactants and isinstance(reactants, dict):
            reactants = [{"name": n, "coefficient": float(v.get("qty", 1))} for n, v in reactants.items()]
        if products and isinstance(products, dict):
            products = [{"name": n, "coefficient": float(v.get("yield", 1))} for n, v in products.items()]
        rptr = f"rxn_{i}_r" if reactants else "NULL"
        pptr = f"rxn_{i}_p" if products else "NULL"
        surface_name = rxn.get("surface_species")
        if surface_name is None and kind == "CHEM_SURFACE" and reactants:
            surface_name = reactants[0].get("name")
        surface_index = names[surface_name] if surface_name is not None else "(size_t)-1"
        out.append("  {%s, %s, %s, %d, %s, %d, %s, {%s}, NULL}," % (
            kind, cstr(str(name)), rptr, len(reactants), pptr, len(products), surface_index,
                ", ".join(f"{v:.17g}" for v in values)))
    out.extend(["};", "", "static const char *const custom_labels[] = {"])
    out.extend(f"  {cstr(label)}," for label in custom_labels)
    out.extend(["};", "static const double default_custom_parameters[] = {"])
    out.extend(f"  {float(defaults.get('rate_parameters', {}).get(label, 0.0)):.17g},"
               for label in custom_labels)
    out.extend(["};", "", "const chem_mechanism chem_generated_mechanism = {",
                "  species, sizeof(species) / sizeof(species[0]),",
                "  reactions, sizeof(reactions) / sizeof(reactions[0]),",
                "  custom_labels, default_custom_parameters,",
                "  sizeof(custom_labels) / sizeof(custom_labels[0]),",
                f"  {float(defaults.get('temperature', 298.15)):.17g},",
                f"  {float(defaults.get('pressure', 101325.0)):.17g}", "};", ""])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(out))


if __name__ == "__main__":
    main()
