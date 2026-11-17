"""Generate C mechanism and rate tables from a CAMP configuration."""
import argparse
import json
import math
from config import read_mechanism, read_document

CONVERSION = 6.02214076e17  # molecules cm-3 per mol m-3


def markowitz_order(count, rows, reactants, products):
    """Choose a generic symmetric permutation to reduce no-pivot LU fill."""
    edges = [{i} for i in range(count)]
    for rb, rn, pb, pn in rows:
        inputs = set(reactants[rb:rb+rn])
        outputs = inputs | {species for species, _ in products[pb:pb+pn]}
        for row in outputs:
            edges[row].update(inputs)
    return markowitz_pattern(edges)


def markowitz_pattern(pattern):
    """Return a symmetric permutation of a structural adjacency list."""
    count = len(pattern)
    edges = [set(row) | {i} for i, row in enumerate(pattern)]
    active = set(range(count))
    order = []
    while active:
        def cost(i):
            return ((len(edges[i] & active)-1) *
                    (sum(i in edges[j] for j in active)-1), i)
        pivot = min(active, key=cost)
        order.append(pivot)
        active.remove(pivot)
        cols = edges[pivot] & active
        for row in active:
            if pivot in edges[row]:
                edges[row].update(cols)
    return order


def generate(path, reorder=False):
    species, source_reactions = read_mechanism(path)
    root = read_document(path)
    third_bodies = {s['name'] for s in species if s.get('third body', False)}
    species = [s for s in species if s['name'] not in third_bodies]
    indices = {s['name']: i for i, s in enumerate(species)}
    species_by_name = {s['name']: s for s in species}
    reactants, products, rows, rate_rows, labels, air_orders = [], [], [], [], [], []

    def label(name):
        if name not in labels:
            labels.append(name)
        return labels.index(name)

    def emit(reaction, reactant_map, product_map, kind, params, custom=0):
        rb, pb = len(reactants), len(products)
        air_order = 0
        for name, properties in reactant_map.items():
            qty = properties.get('qty', 1)
            if not isinstance(qty, int) or isinstance(qty, bool) or qty < 1:
                raise ValueError('reactant qty must be a positive integer')
            if name in third_bodies:
                air_order += qty
            else:
                reactants.extend([indices[name]] * qty)
        for name, properties in product_map.items():
            if name in third_bodies:
                continue
            products.append((indices[name], float(properties.get('yield', 1))))
        rows.append((rb, len(reactants)-rb, pb, len(products)-pb))
        rate_rows.append((kind, params, custom))
        air_orders.append(air_order)

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
        units = get('rate coefficient units', 'molecules cm-3')
        if units not in ('SI', 'molecules cm-3'):
            raise ValueError(f'unsupported rate coefficient units: {units}')
        conversion = 1 if units == 'SI' else CONVERSION
        if typ == 'ARRHENIUS':
            p = [float(get('A', 1))*conversion**(order-1), get('B', 0), get('C', 0), get('D', 300), get('E', 0)]
            emit(reaction, rs, ps, 'ARRHENIUS', p)
        elif typ in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION'):
            low_order = order if typ == 'TROE' else order-1
            p = [float(get('k0_A', 1))*conversion**low_order, get('k0_B', 0), get('k0_C', 0),
                 float(get('kinf_A', 1))*conversion**(low_order-1), get('kinf_B', 0), get('kinf_C', 0),
                 get('Fc', .6), get('N', 1)]
            emit(reaction, rs, ps, 'TROE' if typ == 'TROE' else 'TERNARY', p)
        elif typ == 'REVERSIBLE':
            emit(reaction, rs, ps, 'REVERSIBLE', [float(get('A', 1))*conversion**(order-1), get('C', 0), get('k_r', 0)])
        elif typ == 'TAYLOR_SERIES':
            coefficients = get('coefficients', [1])
            if not isinstance(coefficients, list) or not 1 <= len(coefficients) <= 16:
                raise ValueError('Taylor series requires 1 to 16 coefficients')
            p = [float(get('A', 1))*conversion**(order-1), get('B', 0), get('C', 0), get('D', 300), get('E', 0)] + coefficients
            emit(reaction, rs, ps, 'TAYLOR', p, len(coefficients))
        elif typ == 'TERNARY_CHEMICAL_ACTIVATION_JPL19':
            p = []
            for prefix, scale in [('k0', order), ('kinf', order-1), ('kint', order-1)]:
                p += [float(get(prefix+'_A', 1))*conversion**scale, get(prefix+'_B', 0),
                      get(prefix+'_C', 0), get(prefix+'_D', 298)]
            p += [get('Fc', .6), get('N', 1)]
            emit(reaction, rs, ps, 'TERNARY_JPL19', p)
        elif typ == 'TUNNELING':
            emit(reaction, rs, ps, 'TUNNELING', [float(get('A', 1))*conversion**(order-1), get('B', 0), get('C', 0)])
        elif typ == 'BRANCHED':
            p = [float(get('X'))*conversion**(order-1), get('Y'), get('a0'), get('n')]
            emit(reaction, rs, reaction['alkoxy products'], 'ALKOXY', p)
            emit(reaction, rs, reaction['nitrate products'], 'NITRATE', p)
        elif typ == 'SURFACE':
            sp = species_by_name[name]
            surface_label = get('surface parameter label', name)
            custom = label(surface_label + '.effective radius [m]')
            if label(surface_label + '.particle number concentration [# m-3]') != custom+1:
                raise ValueError('surface parameters must be consecutive')
            emit(reaction, rs, ps, 'SURFACE', [sp['diffusion coefficient [m2 s-1]'], sp['molecular weight [kg mol-1]'], get('reaction probability', 1)], custom)
        elif typ in ('PHOTOLYSIS', 'FIRST_ORDER_LOSS', 'EMISSION', 'USER_DEFINED'):
            emit(reaction, rs, ps, 'USER', [get('scaling factor', 1)], label(reaction['MUSICA name']))
        else:
            raise ValueError(f'unsupported reaction type: {typ}')

    if reorder:
        order = markowitz_order(len(species), rows, reactants, products)
        inverse = {old: new for new, old in enumerate(order)}
        species = [species[old] for old in order]
        reactants = [inverse[i] for i in reactants]
        products = [(inverse[i], y) for i, y in products]

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
    custom_defaults = root.get('custom parameter defaults', {})
    lines.append('const double generated_custom_defaults[] = {' + ', '.join(number(custom_defaults.get(name, 0)) for name in labels or ['']) + '};')
    lines.append('const double generated_initial_concentrations[] = {' + ', '.join(number(s.get('initial concentration [mol m-3]', 0)) for s in species) + '};')
    lines.append('const unsigned generated_air_orders[] = {' + ', '.join(map(str, air_orders or [0])) + '};')
    lines.append('static const size_t reactants[] = {' + ', '.join(map(str, reactants or [0])) + '};')
    lines.append('static const ChemProduct products[] = {' + ', '.join('{%d, %s}' % (i, number(y)) for i, y in (products or [(0, 0)])) + '};')
    lines.append('static const ChemReaction reactions[] = {' + ', '.join('{%d, %d, %d, %d}' % row for row in (rows or [(0, 0, 0, 0)])) + '};')
    lines.append('const ChemRate generated_rates[] = {' + ', '.join('{CHEM_RATE_%s, {%s}, %d}' % (kind, ', '.join(number(x) for x in p), custom) for kind, p, custom in rate_rows) + '};')
    lines.append('const ChemMechanism generated_mechanism = {%d, %d, reactions, reactants, products};' % (len(species), len(rows)))
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reorder', action='store_true')
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()
    with open(args.output, 'w', encoding='utf-8') as out:
        out.write(generate(args.config, reorder=args.reorder))
