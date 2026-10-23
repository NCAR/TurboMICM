#!/usr/bin/env python3
"""Write a canonical mechanism in the astra-max entry's configuration format.

    convert.py CANONICAL_JSON OUTPUT_JSON

The target format is the one of the entry's configs/ts1.json (written by its
tools/translate_ts1.py, read by tools/configure.py and tools/specialize.py):

  {"camp-data": [CHEM_SPEC objects in species order...,
                 {"name": ..., "type": "MECHANISM", "rate units": "mol m-3",
                  "reactions": [...]}],
   "initial state": {"temperature", "pressure", "concentrations", "parameters"}}

Mapping, mirroring configs/ts1.json:
  - species: {"name", "type": "CHEM_SPEC"}; the third body M (canonical
    "third_body": true) gets "parameterization": "air density", as TS1's M does;
  - ARRHENIUS: reactants (names, repeats for stoichiometry), products
    [{"species", "yield"}], A, B, C, D, E;
  - TROE: the same terms, k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N;
  - USER_DEFINED: the same terms, "MUSICA name" (the canonical label) and
    "scaling factor"; an emission is a USER_DEFINED reaction with an empty
    reactant list, as the canonical file gives it;
  - M stays a reactant or product wherever the canonical reaction lists it, as
    in TS1 (the entry's reader multiplies the rate by the air density for each
    M reactant and drops M from the products);
  - "rate units": "mol m-3" (canonical parameters are already in mol m-3);
  - initial state: temperature, pressure, concentrations (never M), and every
    rate parameter's value. The canonical "environment" (air density, time
    step) has no counterpart: the entry recomputes the air density from
    pressure and temperature, and the step is a benchmark argument.
Any other reaction type or field stops the conversion with an error.
"""
import argparse
import json
import math
from pathlib import Path

FIELDS = {
    'ARRHENIUS': ['A', 'B', 'C', 'D', 'E'],
    'TROE': ['k0_A', 'k0_B', 'k0_C', 'kinf_A', 'kinf_B', 'kinf_C', 'Fc', 'N'],
}
COMMON = {'type', 'reactants', 'products'}
EXTRA = {'ARRHENIUS': set(), 'TROE': set(), 'USER_DEFINED': {'label', 'scaling_factor', 'source'}}


def number(value, what):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'non-finite {what}: {value}')
    return result


def convert(canon):
    names = [s['name'] for s in canon['species']]
    if len(set(names)) != len(names):
        raise ValueError('duplicate species names')
    third = {s['name'] for s in canon['species'] if s.get('third_body')}
    species = []
    for s in canon['species']:
        unknown = set(s) - {'name', 'third_body'}
        if unknown:
            raise ValueError(f'species {s["name"]}: unmapped fields {sorted(unknown)}')
        item = {'name': s['name'], 'type': 'CHEM_SPEC'}
        if s.get('third_body'):
            item['parameterization'] = 'air density'
        species.append(item)
    known = set(names)
    reactions = []
    for i, r in enumerate(canon['reactions']):
        kind = r['type']
        if kind not in EXTRA:
            raise ValueError(f'reaction {i}: no mapping for type {kind}')
        unknown = set(r) - COMMON - set(FIELDS.get(kind, [])) - EXTRA[kind]
        if unknown:
            raise ValueError(f'reaction {i} ({kind}): unmapped fields {sorted(unknown)}')
        for name in list(r['reactants']) + [p[0] for p in r['products']]:
            if name not in known:
                raise ValueError(f'reaction {i}: unknown species {name}')
        item = {'type': kind, 'reactants': list(r['reactants']),
                'products': [{'species': name, 'yield': number(y, 'yield')} for name, y in r['products']]}
        if kind == 'USER_DEFINED':
            item['MUSICA name'] = r['label']
            item['scaling factor'] = number(r['scaling_factor'], 'scaling factor')
        else:
            for key in FIELDS[kind]:
                item[key] = number(r[key], key)
        reactions.append(item)
    initial = canon['initial']
    concentrations = {}
    for name, value in initial['concentrations'].items():
        if name in third:
            continue  # TS1's initial state never lists M: it holds no concentration
        if name not in known:
            raise ValueError(f'initial concentration of unknown species {name}')
        concentrations[name] = number(value, f'concentration of {name}')
    labels = {r['label'] for r in canon['reactions'] if r['type'] == 'USER_DEFINED'}
    parameters = {}
    for name, value in initial['parameters'].items():
        if name not in labels:
            raise ValueError(f'initial value for a parameter no reaction uses: {name}')
        parameters[name] = number(value, f'parameter {name}')
    mechanism = {'name': canon['name'], 'type': 'MECHANISM', 'rate units': 'mol m-3',
                 'reactions': reactions}
    return {'camp-data': species + [mechanism],
            'initial state': {'temperature': number(initial['temperature'], 'temperature'),
                              'pressure': number(initial['pressure'], 'pressure'),
                              'concentrations': concentrations,
                              'parameters': parameters}}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('canonical', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = convert(json.loads(args.canonical.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    mechanism = result['camp-data'][-1]
    print(f'{args.output}: {len(result["camp-data"]) - 1} species, {len(mechanism["reactions"])} reactions')


if __name__ == '__main__':
    main()
