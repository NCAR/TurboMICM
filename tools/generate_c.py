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
    "TUNNELING": "CHEM_TUNNELING", "tunneling": "CHEM_TUNNELING",
    "BRANCHED_ALKOXY": "CHEM_BRANCHED_ALKOXY",
    "BRANCHED_NITRATE": "CHEM_BRANCHED_NITRATE",
}


def cstr(value):
    return '"' + value.replace('\\', '\\\\').replace('"', '\\"') + '"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    data = json.loads(args.input.read_text())
    species = data["species"]
    reactions = data["reactions"]
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
        out.append("  {%s, %d, %.17g, %.17g}," % (
            cstr(sp["name"]), int(bool(sp.get("third_body", False))),
            float(sp.get("molecular_weight", 0.0)), float(sp.get("diffusion", 0.0))))
    out.append("};")
    out.append("static const chem_reaction reactions[] = {")
    for i, rxn in enumerate(reactions):
        kind = KINDS.get(rxn.get("type"))
        if kind is None:
            raise SystemExit(f"unsupported rate kind {rxn.get('type')!r} in reaction {i}")
        params = rxn.get("parameters", {})
        if isinstance(params, dict) and "parameters" in params:
            params = params["parameters"]
        values = [0.0] * 8
        if kind == "CHEM_ARRHENIUS":
            for j, key in enumerate(("A", "B", "C", "D", "E")):
                values[j] = float(params.get(key, {"A": 1.0, "D": 300.0}.get(key, 0.0)))
        elif kind == "CHEM_TROE":
            for j, key in enumerate(("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")):
                values[j] = float(params.get(key, 0.0))
        elif kind == "CHEM_TERNARY_ACTIVATION":
            for j, key in enumerate(("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N")):
                values[j] = float(params.get(key, 0.0))
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
            values[0] = float(params.get("scaling_factor", params.get("reaction probability", 1.0)))
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
        out.append("  {%s, %s, %s, %d, %s, %d, %s, {%s}}," % (
            kind, cstr(str(name)), rptr, len(reactants), pptr, len(products), surface_index,
            ", ".join(f"{v:.17g}" for v in values)))
    out.extend(["};", "", "const chem_mechanism chem_generated_mechanism = {",
                "  species, sizeof(species) / sizeof(species[0]),",
                "  reactions, sizeof(reactions) / sizeof(reactions[0])", "};", ""])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(out))


if __name__ == "__main__":
    main()
