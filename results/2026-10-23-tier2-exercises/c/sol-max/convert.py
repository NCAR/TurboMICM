#!/usr/bin/env python3
"""Write a canonical mechanism in sol-max's input format, as configs/ts1.json does.

    convert.py canonical.json OUTPUT.json

sol-max reads CAMP-style JSON (tools/configure.py). configs/ts1.json has:
  - "concentration units": "mol m-3" at the top, so no unit conversion is applied;
  - "camp-data": one CHEM_SPEC record per species; M is
    {"name": "M", "type": "CHEM_SPEC", "parameterization": "third body"};
  - one MECHANISM record whose reactions use "reactant sequence" (repeats for
    stoichiometry, M listed like any reactant) and "product sequence"
    ([{"species", "yield"}]);
  - ARRHENIUS with A, B, C, D, E; TROE with k0_A..C, kinf_A..C, Fc, N;
    USER_DEFINED with "scaling factor" and "MUSICA name" (the per-cell custom
    parameter, which is also the CSV column name);
  - "initial state": temperature, pressure, concentrations, rate parameters.
Every canonical field is either mapped or rejected; nothing is dropped silently,
except a USER_DEFINED reaction's "source" (PHOTOLYSIS or EMISSION), which is a
note: like TS1's photolysis, both become USER_DEFINED, and an emission is one
with an empty reactant sequence.
"""
import json
import sys

KEYS = {
    'ARRHENIUS': ['A', 'B', 'C', 'D', 'E'],
    'TROE': ['k0_A', 'k0_B', 'k0_C', 'kinf_A', 'kinf_B', 'kinf_C', 'Fc', 'N'],
    'USER_DEFINED': [],
}
COMMON = {'type', 'reactants', 'products'}
USER = {'label', 'scaling_factor', 'source'}


def convert(canon):
    names = [s['name'] for s in canon['species']]
    if len(set(names)) != len(names):
        raise ValueError('duplicate species names')
    camp = []
    for s in canon['species']:
        if set(s) - {'name', 'third_body'}:
            raise ValueError(f"unmapped species fields: {sorted(set(s) - {'name', 'third_body'})}")
        record = {'name': s['name'], 'type': 'CHEM_SPEC'}
        if s['third_body']:
            record['parameterization'] = 'third body'
        camp.append(record)
    reactions = []
    for i, r in enumerate(canon['reactions']):
        kind = r['type']
        if kind not in KEYS:
            raise ValueError(f'reaction {i}: unsupported type {kind}')
        allowed = COMMON | set(KEYS[kind]) | (USER if kind == 'USER_DEFINED' else set())
        if set(r) - allowed:
            raise ValueError(f'reaction {i}: unmapped fields {sorted(set(r) - allowed)}')
        for name in r['reactants'] + [p[0] for p in r['products']]:
            if name not in names:
                raise ValueError(f'reaction {i}: unknown species {name}')
        out = {'type': kind,
               'reactant sequence': list(r['reactants']),
               'product sequence': [{'species': name, 'yield': float(y)} for name, y in r['products']]}
        for key in KEYS[kind]:
            out[key] = float(r[key])
        if kind == 'USER_DEFINED':
            out['scaling factor'] = float(r['scaling_factor'])
            out['MUSICA name'] = r['label']
        reactions.append(out)
    camp.append({'name': canon['name'], 'type': 'MECHANISM', 'reactions': reactions})
    initial = canon['initial']
    for name in initial['concentrations']:
        if name not in names:
            raise ValueError(f'initial concentration for unknown species {name}')
    return {
        'concentration units': 'mol m-3',
        'camp-data': camp,
        'initial state': {
            'temperature': float(initial['temperature']),
            'pressure': float(initial['pressure']),
            'concentrations': {k: float(v) for k, v in initial['concentrations'].items()},
            'rate parameters': {k: float(v) for k, v in initial['parameters'].items()},
        },
    }


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        canonical = json.load(f)
    with open(sys.argv[2], 'w') as f:
        json.dump(convert(canonical), f, indent=2)
        f.write('\n')
