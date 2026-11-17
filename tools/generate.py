"""Generate C mechanism and rate tables from a CAMP configuration."""
import argparse
import json
import math
from config import read_mechanism

CONVERSION = 6.02214076e17  # molecules cm-3 per mol m-3


def generate(path):
    species, source_reactions = read_mechanism(path)
    indices = {s['name']: i for i, s in enumerate(species)}
    species_by_name = {s['name']: s for s in species}
    reactants, products, rows, rate_rows, labels = [], [], [], [], []

    def label(name):
        if name not in labels:
            labels.append(name)
        return labels.index(name)

    def emit(reaction, reactant_map, product_map, kind, params, custom=0):
        rb, pb = len(reactants), len(products)
        for name, properties in reactant_map.items():
            qty = properties.get('qty', 1)
            if not isinstance(qty, int) or isinstance(qty, bool) or qty < 1:
                raise ValueError('reactant qty must be a positive integer')
            reactants.extend([indices[name]] * qty)
        for name, properties in product_map.items():
            products.append((indices[name], float(properties.get('yield', 1))))
        rows.append((rb, len(reactants)-rb, pb, len(products)-pb))
        rate_rows.append((kind, params, custom))

    for reaction in source_reactions:
        typ = reaction['type']
        rs = reaction.get('reactants', {})
        ps = reaction.get('products', {})
        if typ == 'SURFACE':
            name = reaction['gas-phase reactant']
            rs, ps = {name: {}}, reaction.get('gas-phase products', {})
        elif typ == 'FIRST_ORDER_LOSS':
            rs, ps = {reaction['species']: {}}, {}
        elif typ == 'EMISSION':
            rs, ps = {}, {reaction['species']: {}}
        order = sum(properties.get('qty', 1) for properties in rs.values())
        get = reaction.get
        if typ == 'ARRHENIUS':
            p = [float(get('A', 1))*CONVERSION**(order-1), get('B', 0), get('C', 0), get('D', 300), get('E', 0)]
            emit(reaction, rs, ps, 'ARRHENIUS', p)
        elif typ in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION'):
            low_order = order if typ == 'TROE' else order-1
            p = [float(get('k0_A', 1))*CONVERSION**low_order, get('k0_B', 0), get('k0_C', 0),
                 float(get('kinf_A', 1))*CONVERSION**(low_order-1), get('kinf_B', 0), get('kinf_C', 0),
                 get('Fc', .6), get('N', 1)]
            emit(reaction, rs, ps, 'TROE' if typ == 'TROE' else 'TERNARY', p)
        elif typ == 'TUNNELING':
            emit(reaction, rs, ps, 'TUNNELING', [float(get('A', 1))*CONVERSION**(order-1), get('B', 0), get('C', 0)])
        elif typ == 'BRANCHED':
            p = [float(get('X'))*CONVERSION**(order-1), get('Y'), get('a0'), get('n')]
            emit(reaction, rs, reaction['alkoxy products'], 'ALKOXY', p)
            emit(reaction, rs, reaction['nitrate products'], 'NITRATE', p)
        elif typ == 'SURFACE':
            sp = species_by_name[name]
            custom = label(name + '.effective radius [m]')
            if label(name + '.particle number concentration [# m-3]') != custom+1:
                raise ValueError('surface parameters must be consecutive')
            emit(reaction, rs, ps, 'SURFACE', [sp['diffusion coefficient [m2 s-1]'], sp['molecular weight [kg mol-1]'], get('reaction probability', 1)], custom)
        elif typ in ('PHOTOLYSIS', 'FIRST_ORDER_LOSS', 'EMISSION'):
            emit(reaction, rs, ps, 'USER', [1], label(reaction['MUSICA name']))
        else:
            raise ValueError(f'unsupported reaction type: {typ}')

    def number(value):
        value = float(value)
        if not math.isfinite(value):
            raise ValueError('non-finite rate parameter')
        return format(value, '.17g')
    def strings(values):
        return ', '.join(json.dumps(s, ensure_ascii=True) for s in values) or '0'
    lines = ['#include "chem.h"', '#include "rates.h"', '']
    lines.append('const char *const generated_species[] = {' + strings(s['name'] for s in species) + '};')
    lines.append('const char *const generated_custom_labels[] = {' + strings(labels) + '};')
    lines.append(f'const size_t generated_custom_count = {len(labels)};')
    lines.append('static const size_t reactants[] = {' + ', '.join(map(str, reactants or [0])) + '};')
    lines.append('static const ChemProduct products[] = {' + ', '.join('{%d, %s}' % (i, number(y)) for i, y in (products or [(0, 0)])) + '};')
    lines.append('static const ChemReaction reactions[] = {' + ', '.join('{%d, %d, %d, %d}' % row for row in (rows or [(0, 0, 0, 0)])) + '};')
    lines.append('const ChemRate generated_rates[] = {' + ', '.join('{CHEM_RATE_%s, {%s}, %d}' % (kind, ', '.join(number(x) for x in p), custom) for kind, p, custom in rate_rows) + '};')
    lines.append('const ChemMechanism generated_mechanism = {%d, %d, reactions, reactants, products};' % (len(species), len(rows)))
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()
    with open(args.output, 'w', encoding='utf-8') as out:
        out.write(generate(args.config))
