#!/usr/bin/env python3
"""Read CAMP-style JSON or YAML mechanism configurations into one document."""
import argparse
import json
from pathlib import Path

MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * 6.02214076e23


def load(path):
    text = path.read_text()
    if path.suffix.lower() == ".json":
        return json.loads(text)
    try:
        import yaml
    except ImportError as exc:
        raise SystemExit("YAML input requires PyYAML; run `make env` and use .conda/bin/python") from exc
    return yaml.safe_load(text)


def read_config(path):
    path = Path(path).resolve()
    root = load(path)
    if not isinstance(root, dict) or not isinstance(root.get("camp-files"), list):
        raise ValueError(f"{path}: expected a camp-files list")
    documents = [load(path.parent / name) for name in root["camp-files"]]
    entries = [entry for doc in documents for entry in doc.get("camp-data", [])]
    species = []
    reactions = []
    for entry in entries:
        if entry.get("type") == "CHEM_SPEC":
            species.append({
                "name": entry["name"],
                "molecular_weight": float(entry.get("molecular weight [kg mol-1]", 0.0)),
                "diffusion": float(entry.get("diffusion coefficient [m2 s-1]", 0.0)),
                "third_body": bool(entry.get("third body", False)),
                "density_parameterized": bool(entry.get("density parameterized", False)),
                "density_factor": float(entry.get("density factor", 1.0)),
            })
        elif entry.get("type") == "MECHANISM":
            reactions.extend(entry.get("reactions", []))
    if not species:
        raise ValueError(f"{path}: no CHEM_SPEC records")
    index = {s["name"]: i for i, s in enumerate(species)}
    if len(index) != len(species):
        raise ValueError(f"{path}: duplicate species names")
    normalized = []
    for reaction in reactions:
        rtype = reaction.get("type", "").upper()
        reactants = reaction.get("reactants", {})
        products = reaction.get("products", {})
        if rtype == "SURFACE":
            reactants = {reaction["gas-phase reactant"]: {}}
            products = reaction.get("gas-phase products", {})
        if rtype in ("EMISSION", "FIRST_ORDER_LOSS"):
            name = reaction["species"]
            reactants, products = ({name: {}} if rtype == "FIRST_ORDER_LOSS" else {}), ({name: {}} if rtype == "EMISSION" else {})
        for terms in (reactants, products):
            for name in terms:
                if name not in index:
                    raise ValueError(f"{path}: reaction uses unknown species {name!r}")
        branches = []
        if rtype == "BRANCHED":
            for branch_name, key in (("alkoxy", "alkoxy products"), ("nitrate", "nitrate products")):
                branches.append((f"BRANCHED_{branch_name.upper()}", reaction.get(key, {}), branch_name))
        else:
            branches.append((rtype, products, None))
        for normalized_type, branch_products, branch in branches:
            for name in branch_products:
                if name not in index:
                    raise ValueError(f"{path}: reaction uses unknown species {name!r}")
            parameters = {k: v for k, v in reaction.items()
                          if k not in {"type", "reactants", "products", "alkoxy products", "nitrate products",
                                       "gas-phase reactant", "gas-phase products", "species"}}
            if normalized_type == "TROE":
                parameters["k0_A"] = float(parameters["k0_A"]) * MOLES_M3_TO_MOLECULES_CM3**2
                parameters["kinf_A"] = float(parameters["kinf_A"]) * MOLES_M3_TO_MOLECULES_CM3
            elif normalized_type == "TERNARY_CHEMICAL_ACTIVATION":
                parameters["kinf_A"] = float(parameters["kinf_A"]) / MOLES_M3_TO_MOLECULES_CM3
            normalized.append({
                "type": "USER_DEFINED" if normalized_type in ("PHOTOLYSIS", "FIRST_ORDER_LOSS", "EMISSION") else normalized_type,
                "surface_species": reaction.get("gas-phase reactant") if rtype == "SURFACE" else None,
                "reactants": [{"name": n, "index": index[n], "coefficient": float(v.get("qty", 1))}
                              for n, v in reactants.items()],
                "products": [{"name": n, "index": index[n], "coefficient": float(v.get("yield", 1))}
                             for n, v in branch_products.items()],
                "parameters": {**parameters, **({"branch": branch} if branch else {})},
            })
    return {"format": "c-mechanism-1", "name": path.parent.name,
            "species": species, "reactions": normalized}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("configuration", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    args.output.write_text(json.dumps(read_config(args.configuration), indent=2) + "\n")


if __name__ == "__main__":
    main()
