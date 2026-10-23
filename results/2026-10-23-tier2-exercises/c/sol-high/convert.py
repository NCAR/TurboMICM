#!/usr/bin/env python3
"""Convert a canonical mechanism (canonical.json) to sol-high's TS1 input format.

sol-high's build reads configs/ts1.json: extended CAMP JSON, one "camp-data"
list holding a CHEM_SPEC entry per species and one MECHANISM entry. This writes
the same layout that the entry's tools/translate_ts1.py produced for TS1:
  - species: {"name", "type": "CHEM_SPEC", "initial concentration"}; the third
    body as {"initial concentration": 0, "third body": true}, as TS1's M;
  - reactants as a list of names with repeats; products as a list of
    {"species", "yield"}; M kept in reactants and products, as in TS1;
  - ARRHENIUS (A, B, C, D, E) and TROE (k0_*, kinf_*, Fc, N) parameters under
    their own keys; USER_DEFINED with "scaling_factor" and "MUSICA name" (the
    label), emissions included (an empty reactant list);
  - MECHANISM "rate units": "mol m-3", "initial conditions" (temperature,
    pressure), and "initial rate parameters" (label -> value).

Usage: convert.py CANONICAL_JSON OUTPUT_JSON
"""
import json
import sys

KEYS = {
    'ARRHENIUS': ['A', 'B', 'C', 'D', 'E'],
    'TROE': ['k0_A', 'k0_B', 'k0_C', 'kinf_A', 'kinf_B', 'kinf_C', 'Fc', 'N'],
    'USER_DEFINED': ['scaling_factor'],
}


def convert(canon):
    init = canon['initial']
    conc = init.get('concentrations', {})
    names = [s['name'] for s in canon['species']]
    unknown = set(conc) - set(names)
    if unknown:
        raise ValueError(f'initial concentrations for unknown species: {sorted(unknown)}')
    species = []
    for s in canon['species']:
        entry = {'name': s['name'], 'type': 'CHEM_SPEC'}
        if s.get('third_body', False):
            entry['initial concentration'] = 0
            entry['third body'] = True
        else:
            entry['initial concentration'] = float(conc.get(s['name'], 0.0))
        species.append(entry)
    reactions = []
    for i, r in enumerate(canon['reactions']):
        kind = r['type']
        if kind not in KEYS:
            raise ValueError(f'reaction {i}: no TS1-format equivalent for type {kind}')
        out = {'type': kind,
               'reactants': list(r['reactants']),
               'products': [{'species': n, 'yield': float(y)} for n, y in r['products']]}
        for k in KEYS[kind]:
            out[k] = float(r[k])
        if kind == 'USER_DEFINED':
            out['MUSICA name'] = r['label']
        for n in out['reactants'] + [p['species'] for p in out['products']]:
            if n not in names:
                raise ValueError(f'reaction {i}: unknown species {n}')
        reactions.append(out)
    labels = [r['MUSICA name'] for r in reactions if 'MUSICA name' in r]
    if len(labels) != len(set(labels)):
        raise ValueError('a label is shared by two reactions (the reader requires unique labels)')
    params = init.get('parameters', {})
    initial_rate_parameters = {label: float(params[label]) for label in labels if label in params}
    mechanism = {'name': canon['name'], 'type': 'MECHANISM', 'reactions': reactions,
                 'rate units': 'mol m-3',
                 'initial conditions': {'temperature': float(init['temperature']),
                                        'pressure': float(init['pressure'])},
                 'initial rate parameters': initial_rate_parameters}
    return {'camp-data': species + [mechanism]}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        canon = json.load(f)
    with open(sys.argv[2], 'w') as f:
        f.write(json.dumps(convert(canon), indent=2) + '\n')


if __name__ == '__main__':
    main()
