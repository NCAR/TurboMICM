"""Read CAMP YAML/JSON configurations and emit standalone C11 builders.

No chemistry package is used. Numerical runtime data is supplied exclusively
through the public C mechanism API. See DECISIONS.md for format extensions.
"""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import re
import yaml

AVOGADRO_CM = 1e-6 * 6.02214076e23
TYPES = {
    "ARRHENIUS": ("ARRHENIUS", "arrhenius", {"A": 1, "B": 0, "C": 0, "D": 300, "E": 0}),
    "TROE": ("TROE", "falloff", {"k0_A": 1, "k0_B": 0, "k0_C": 0, "kinf_A": 1, "kinf_B": 0, "kinf_C": 0, "Fc": .6, "N": 1}),
    "TERNARY_CHEMICAL_ACTIVATION": ("TERNARY", "falloff", {"k0_A": 1, "k0_B": 0, "k0_C": 0, "kinf_A": 1, "kinf_B": 0, "kinf_C": 0, "Fc": .6, "N": 1}),
    "TERNARY_CHEMICAL_ACTIVATION_JPL19": ("TERNARY_JPL19", "jpl19", {"k0_A": 1, "k0_B": 0, "k0_C": 0, "k0_D": 298, "kinf_A": 1, "kinf_B": 0, "kinf_C": 0, "kinf_D": 298, "kint_A": 1, "kint_B": 0, "kint_C": 0, "kint_D": 298, "Fc": .6, "N": 1}),
    "BRANCHED": ("BRANCHED", "branched", {"X": 0, "Y": 0, "a0": 0, "n": 0}),
    "TUNNELING": ("TUNNELING", "tunneling", {"A": 1, "B": 0, "C": 0}),
    "TAYLOR_SERIES": ("TAYLOR", "taylor.arrhenius", {"A": 1, "B": 0, "C": 0, "D": 300, "E": 0}),
    "REVERSIBLE": ("REVERSIBLE", "reversible", {"A": 1, "C": 0, "reverse_rate": 0}),
    "USER_DEFINED": ("USER", "user", {"scaling": 1}),
    "SURFACE": ("SURFACE", "surface", {"probability": 1}),
}

def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError("Expected a configuration number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite configuration number")
    return value

def load_document(path: Path):
    text = path.read_text()
    doc = json.loads(text) if path.suffix.lower() == ".json" else yaml.safe_load(text)
    if not isinstance(doc, dict):
        raise ValueError(f"Expected mapping in {path}")
    return doc

def records(path: Path, stack=()):
    path = path.resolve()
    if path in stack:
        raise ValueError(f"Configuration include cycle: {path}")
    doc = load_document(path)
    data, files = doc.get("camp-data", []), doc.get("camp-files", [])
    if not isinstance(data, list) or any(not isinstance(r, dict) or not isinstance(r.get("type"), str) for r in data):
        raise ValueError("camp-data must be a list of typed records")
    if not isinstance(files, list) or any(not isinstance(f, str) or not f for f in files):
        raise ValueError("camp-files must be a list of filenames")
    result = list(data)
    for filename in files:
        result.extend(records(path.parent / filename, (*stack, path)))
    return result

def reactants(value):
    if isinstance(value, list):
        if any(not isinstance(v, str) or not v for v in value):
            raise ValueError("Reactants must name species")
        return list(value)  # ordered occurrences extension
    result = []
    if not isinstance(value, dict):
        raise ValueError("Reactants must be a mapping or ordered list")
    for name, spec in value.items():
        if spec is not None and not isinstance(spec, dict):
            raise ValueError("Reactant attributes must be a mapping")
        qty = number((spec or {}).get("qty", 1))
        if qty < 1 or not qty.is_integer():
            raise ValueError("Reactant qty must be a positive integer")
        result.extend([name] * int(qty))
    return result

def products(value):
    if isinstance(value, list):
        return [(v["species"], number(v.get("yield", 1))) for v in value]
    if not isinstance(value, dict) or any(spec is not None and not isinstance(spec, dict) for spec in value.values()):
        raise ValueError("Products must be a mapping or ordered list")
    return [(name, number((spec or {}).get("yield", 1))) for name, spec in value.items()]

def read_mechanism(path):
    data = records(Path(path))
    species = [copy.deepcopy(r) for r in data if r["type"] == "CHEM_SPEC"]
    names = [s.get("name") for s in species]
    if any(not isinstance(n, str) or not n or "\0" in n for n in names):
        raise ValueError("Species names must be nonempty strings without NUL")
    if len(set(names)) != len(names):
        raise ValueError("Duplicate species name")
    by_name = {s["name"]: s for s in species}
    for s in species:
        for key in ("absolute tolerance", "molecular weight [kg mol-1]", "diffusion coefficient [m2 s-1]"):
            if key in s and number(s[key]) < 0:
                raise ValueError(f"Negative species property {key}")
        if s.get("parameterization") not in {None, "air density"}:
            raise ValueError("Unknown species parameterization")
    reactions = []
    initial = {}
    for record in data:
        if record["type"] == "CHEM_SPEC":
            continue
        if record["type"] != "MECHANISM":
            raise ValueError(f"Unknown CAMP record type {record['type']}")
        state, definitions = record.get("initial state", {}), record.get("reactions", [])
        if not isinstance(state, dict) or not isinstance(definitions, list):
            raise ValueError("Expected initial-state mapping and reactions list")
        initial.update(state)
        for source in definitions:
            if not isinstance(source, dict) or not isinstance(source.get("type"), str):
                raise ValueError("Expected a typed reaction mapping")
            r = copy.deepcopy(source)
            kind = r["type"]
            if kind in {"PHOTOLYSIS", "FIRST_ORDER_LOSS", "EMISSION"}:
                if kind == "FIRST_ORDER_LOSS":
                    r["reactants"], r["products"] = {r["species"]: {}}, {}
                elif kind == "EMISSION":
                    r["reactants"], r["products"] = {}, {r["species"]: {}}
                kind = "USER_DEFINED"
            if kind not in TYPES:
                raise ValueError(f"Unknown reaction type {kind}")
            if kind == "SURFACE":
                r["reactants"] = {r["gas-phase reactant"]: {}}
                r["products"] = r.get("gas-phase products", {})
            lhs = reactants(r.get("reactants", {}))
            rhs = products(r.get("products", {}))
            for name in lhs + [n for n, _ in rhs]:
                if name not in by_name:
                    raise ValueError(f"Unknown reaction species {name}")
            ctype, member, defaults = TYPES[kind]
            params = {key: number(r.get(key, default)) for key, default in defaults.items()}
            units = r.get("rate units", record.get("rate units", "molecule cm-3"))
            if units not in {"molecule cm-3", "mol m-3"}:
                raise ValueError(f"Unsupported rate units {units}")
            factor = AVOGADRO_CM if units == "molecule cm-3" else 1.0
            order = len(lhs)
            if kind in {"ARRHENIUS", "TUNNELING", "TAYLOR_SERIES"}:
                params["A"] *= factor ** (order - 1)
            elif kind in {"TROE", "TERNARY_CHEMICAL_ACTIVATION", "TERNARY_CHEMICAL_ACTIVATION_JPL19"}:
                extra = 0 if kind == "TERNARY_CHEMICAL_ACTIVATION" else 1
                params["k0_A"] *= factor ** (order - 1 + extra)
                params["kinf_A"] *= factor ** (order - 2 + extra)
                if "kint_A" in params:
                    params["kint_A"] *= factor ** (order - 1)
            elif kind == "BRANCHED":
                params["X"] *= factor ** (order - 1)
            elif kind == "USER_DEFINED":
                params["scaling"] = number(r.get("scaling factor", r.get("scaling", 1)))
            elif kind == "SURFACE":
                sp = by_name[r["gas-phase reactant"]]
                params["probability"] = number(r.get("reaction probability", 1))
                params["diffusion"] = number(sp["diffusion coefficient [m2 s-1]"])
                params["molecular_weight"] = number(sp["molecular weight [kg mol-1]"])
            normal = dict(type=ctype, member=member, params=params, reactants=lhs,
                          products=rhs, label=r.get("MUSICA name", r.get("label", "")))
            if kind == "TAYLOR_SERIES":
                normal["coefficients"] = [number(x) for x in r.get("coefficients", [1])]
                if not 1 <= len(normal["coefficients"]) <= 16:
                    raise ValueError("Taylor series requires 1 to 16 coefficients")
            if not isinstance(normal["label"], str) or "\0" in normal["label"]:
                raise ValueError("Rate labels must be strings without NUL")
            if kind == "BRANCHED":
                if params["n"] < 0 or not params["n"].is_integer() or params["n"] > 4294967295:
                    raise ValueError("Branched n must be a nonnegative integer")
                for branch, key in enumerate(("alkoxy products", "nitrate products")):
                    item = copy.deepcopy(normal)
                    item["params"]["branch"] = branch
                    item["products"] = products(r[key])
                    if any(n not in by_name for n, _ in item["products"]):
                        raise ValueError("Unknown branch product")
                    reactions.append(item)
            else:
                reactions.append(normal)
    labels = set()
    for r in reactions:
        if r["type"] == "USER":
            labels.add(r["label"])
        elif r["type"] == "SURFACE":
            labels.update((r["label"] + ".effective radius [m]", r["label"] + ".particle number concentration [# m-3]"))
    variables = {s["name"] for s in species if not s.get("parameterization")}
    for key, value in initial.items():
        if key not in variables | labels | {"temperature", "pressure"}:
            raise ValueError(f"Unknown initial-state field {key}")
        number(value)
    return dict(species=species, reactions=reactions, initial=initial)

def literal(v):
    return format(number(v), ".17g")

def emit_c(mechanism, symbol):
    if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", symbol):
        raise ValueError("Invalid C identifier")
    species, reactions = mechanism["species"], mechanism["reactions"]
    indices = {s["name"]: i for i, s in enumerate(species)}
    q = lambda s: json.dumps(s, ensure_ascii=True)
    out = ['/* Generated configuration data. SPDX-License-Identifier: Apache-2.0 */', '#include "chem.h"']
    out.append(f'ChemStatus {symbol}_create(ChemMechanism **out)\n{{')
    out.append(f'    const ChemSpecies species[{max(1,len(species))}] = {{')
    for sp in species:
        fields = [f'.name={q(sp["name"])}']
        for key, member in (("absolute tolerance", "absolute_tolerance"), ("molecular weight [kg mol-1]", "molecular_weight"), ("diffusion coefficient [m2 s-1]", "diffusion")):
            if key in sp:
                fields.append(f'.{member}={literal(sp[key])}')
        if sp.get("parameterization") == "air density":
            fields.append('.parameter=chem_parameter_air_density')
        out.append('        {' + ', '.join(fields) + '},')
    if not species:
        out.append('        {0},')
    out.append('    };')
    out.append(f'    ChemReaction reactions[{max(1,len(reactions))}] = {{0}};')
    for i, r in enumerate(reactions):
        prefix = f'reactions[{i}]'
        if r['reactants']:
            values = ', '.join(str(indices[n]) for n in r['reactants'])
            out.append(f'    const size_t reactants_{i}[] = {{{values}}};')
            out.append(f'    {prefix}.reactants = reactants_{i}; {prefix}.reactant_count = {len(r["reactants"])};')
        if r['products']:
            values = ', '.join('{' + f'{indices[n]}, {literal(v)}' + '}' for n, v in r['products'])
            out.append(f'    const ChemProduct products_{i}[] = {{{values}}};')
            out.append(f'    {prefix}.products = products_{i}; {prefix}.product_count = {len(r["products"])};')
        out.append(f'    {prefix}.rate = chem_rate_default(CHEM_{r["type"]});')
        out.append(f'    {prefix}.rate.label = {q(r["label"])};')
        for key, value in r['params'].items():
            out.append(f'    {prefix}.rate.data.{r["member"]}.{key} = {literal(value)};')
        if 'coefficients' in r:
            out.append(f'    {prefix}.rate.data.taylor.count = {len(r["coefficients"])};')
            for j, v in enumerate(r['coefficients']):
                out.append(f'    {prefix}.rate.data.taylor.coefficients[{j}] = {literal(v)};')
    out.append(f'    return chem_mechanism_create({len(species)}, species, {len(reactions)}, reactions, out);\n}}')
    out.append(f'ChemStatus {symbol}_initialize(ChemState *state)\n{{')
    out.append('    ChemStatus status = CHEM_OK; double value; (void)value;')
    for key, value in mechanism['initial'].items():
        if key in {'temperature', 'pressure'}:
            out.append(f'    for (size_t c=0; c<state->cells; ++c) state->conditions[c].{key}={literal(value)};')
        else:
            setter = 'concentration' if key in indices else 'custom'
            out.append(f'    value={literal(value)}; status=chem_state_fill_{setter}(state,{q(key)},value);')
            out.append('    if (status != CHEM_OK) return status;')
    out.append('    for (size_t c=0; c<state->cells; ++c) chem_ideal_air_density(&state->conditions[c]);')
    out.append('    return status;\n}')
    return '\n'.join(out) + '\n'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--symbol', default='configured')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(emit_c(read_mechanism(args.configuration), args.symbol))
