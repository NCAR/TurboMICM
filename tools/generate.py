#!/usr/bin/env python3
"""Generate immutable C mechanism data from a validated configuration."""
import argparse
import json
import math
from pathlib import Path
from mechanism import MechanismError, read_mechanism


def number(x):
    x = float(x)
    if not math.isfinite(x):
        raise MechanismError('non-finite numeric parameter')
    return x.hex()


def reorder_species(species, rs, ps, reactions):
    """Choose a diagonal minimum-degree ordering and remap table indices."""
    n = len(species)
    rows = [{i} for i in range(n)]
    for _, rb, rc, pb, pc, _, _ in reactions:
        cols = {s for s, _ in rs[rb:rb+rc] if s >= 0}
        targets = cols | {s for s, _ in ps[pb:pb+pc]}
        for row in targets:
            rows[row].update(cols)
    live = set(range(n))
    order = []
    while live:
        def cost(i):
            row = len(rows[i] & live) - 1
            col = sum(i in rows[j] for j in live) - 1
            return (row * col, i)
        pivot = min(live, key=cost)
        remaining = live - {pivot}
        cols = rows[pivot] & remaining
        for row in remaining:
            if pivot in rows[row]:
                rows[row].update(cols)
        live.remove(pivot)
        order.append(pivot)
    inverse = {old: new for new, old in enumerate(order)}
    species = [species[i] for i in order]
    rs = [(inverse[i] if i >= 0 else i, value) for i, value in rs]
    ps = [(inverse[i], value) for i, value in ps]
    return species, rs, ps


def compile_mechanism(mech):
    species = [s for s in mech['species'] if not s.get('third body') and 'parameterization' not in s]
    parameterized = [s for s in mech['species'] if 'parameterization' in s]
    indices = {s['name']: i for i, s in enumerate(species)}
    indices.update({s['name']: -1 for s in mech['species'] if s.get('third body')})
    indices.update({s['name']: -2-i for i, s in enumerate(parameterized)})
    properties = {s['name']: s for s in mech['species']}
    params, rs, ps, reactions = [], [], [], []
    def param(name):
        if not isinstance(name, str) or not name:
            raise MechanismError('custom rate must have a nonempty MUSICA name')
        if name not in params:
            params.append(name)
        return params.index(name)
    for original in mech['reactions']:
        r = dict(original)
        typ = r['type']
        if typ == 'BRANCHED':
            variants = [('BRANCH_ALKOXY', r.get('alkoxy products', {})), ('BRANCH_NITRATE', r.get('nitrate products', {}))]
        else:
            variants = [(typ, r.get('products', {}))]
        for kind, products in variants:
            reactants = r.get('reactants', {})
            if typ == 'SURFACE' and 'gas-phase reactant' in r and not reactants:
                reactants = {r['gas-phase reactant']: {}}
                products = r.get('gas-phase products', {})
            if typ in ('FIRST_ORDER_LOSS', 'EMISSION'):
                reactants = {r['species']: {}} if typ == 'FIRST_ORDER_LOSS' else {}
                products = {r['species']: {}} if typ == 'EMISSION' else {}
            rb, pb = len(rs), len(ps)
            for name, attrs in reactants.items():
                rs.extend([(indices[name], 1.0)] * attrs.get('qty', 1))
            for name, attrs in products.items():
                if indices[name] >= 0:
                    ps.append((indices[name], attrs.get('yield', 1.0)))
            idx = -1
            if kind == 'ARRHENIUS':
                enum, values = kind, [r.get(k, d) for k, d in zip('ABCDE', (1, 0, 0, 300, 0))]
            elif kind in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION'):
                enum = 'TERNARY' if kind.startswith('TERNARY') else kind
                values = [r.get(k, d) for k, d in zip(('k0_A','k0_B','k0_C','kinf_A','kinf_B','kinf_C','Fc','N'), (1,0,0,1,0,0,.6,1))]
            elif kind == 'TERNARY_JPL19':
                enum = kind
                keys = ('k0_A','k0_B','k0_C','k0_D','kinf_A','kinf_B','kinf_C','kinf_D','kint_A','kint_B','kint_C','kint_D','Fc','N')
                defaults = (1,0,0,298,1,0,0,298,1,0,0,298,.6,1)
                values = [r.get(k, d) for k, d in zip(keys, defaults)]
            elif kind == 'REVERSIBLE':
                enum, values = kind, [r.get('A', 1), r.get('C', 0), r.get('k_r', 0)]
            elif kind == 'TAYLOR':
                coefficients = r.get('coefficients', [1])
                if not isinstance(coefficients, list) or not 1 <= len(coefficients) <= 16:
                    raise MechanismError('Taylor rate needs 1 to 16 coefficients')
                enum, values = kind, [r.get(k, d) for k, d in zip('ABCDE', (1, 0, 0, 300, 0))] + [len(coefficients)] + coefficients
            elif kind == 'TUNNELING':
                enum, values = kind, [r.get(k, d) for k, d in zip('ABC', (1,0,0))]
            elif kind in ('USER_DEFINED', 'PHOTOLYSIS', 'FIRST_ORDER_LOSS', 'EMISSION'):
                enum, values = 'USER', [r.get('scaling factor', 1)]
                idx = param(r.get('MUSICA name'))
            elif kind == 'SURFACE':
                enum = kind
                name = r.get('MUSICA name')
                idx = param(name + '.effective radius [m]')
                if param(name + '.particle number concentration [# m-3]') != idx + 1:
                    raise MechanismError('surface parameter pair must be contiguous')
                sp = properties[r['gas-phase reactant']]
                for key in ('diffusion coefficient [m2 s-1]', 'molecular weight [kg mol-1]'):
                    if key not in sp or float(sp[key]) <= 0:
                        raise MechanismError(f"species {sp['name']}: missing or invalid {key}")
                values = [sp['diffusion coefficient [m2 s-1]'], sp['molecular weight [kg mol-1]'], r.get('reaction probability', 1)]
            elif kind.startswith('BRANCH_'):
                enum, values = kind, [r[k] for k in ('X','Y','a0','n')]
            else:
                raise MechanismError('unsupported reaction ' + kind)
            # Tutorial rate constants are in molecule cm-3 units. TS1's
            # generated C++ constants have already been converted to SI.
            if mech['rate units'] == 'molecule cm-3':
                conversion = 1e-6 * 6.02214076e23
                order = len(rs)-rb
                if enum in ('ARRHENIUS', 'TUNNELING', 'TAYLOR'):
                    values[0] *= conversion ** (order-1)
                elif enum == 'TROE':
                    values[0] *= conversion ** order
                    values[3] *= conversion ** (order-1)
                elif enum == 'TERNARY_JPL19':
                    values[0] *= conversion ** order
                    values[4] *= conversion ** (order-1)
                    values[8] *= conversion ** (order-1)
                elif enum == 'TERNARY':
                    values[0] *= conversion ** (order-1)
                    values[3] *= conversion ** (order-2)
                elif enum.startswith('BRANCH_'):
                    values[0] *= conversion ** (order-1)
            reactions.append((enum, rb, len(rs)-rb, pb, len(ps)-pb, idx, values))
    species, rs, ps = reorder_species(species, rs, ps, reactions)
    return species, params, rs, ps, reactions


def emit(mech, symbol):
    species, params, rs, ps, reactions = compile_mechanism(mech)
    out = ['/* Generated by tools/generate.py. */', '#include "chem.h"', '']
    def array(ctype, name, items):
        out.append(f'static const {ctype} {name}[] = {{')
        out.extend('    ' + item + ',' for item in (items or ['0']))
        out.append('};')
    array('char *const', 'species_names', [json.dumps(s['name']) for s in species])
    array('char *const', 'parameter_names', [json.dumps(n) for n in params])
    parameterized = [s for s in mech['species'] if 'parameterization' in s]
    entries = []
    for item in parameterized:
        coefficients = item['parameterization']
        if not isinstance(coefficients, dict) or set(coefficients) - {'constant', 'temperature', 'pressure', 'air density'}:
            raise MechanismError('invalid species parameterization')
        entries.append('{' + ', '.join(number(coefficients.get(k, 0)) for k in ('constant','temperature','pressure','air density')) + '}')
    array('ChemParameterization', 'parameterized_species', entries or ['{0, 0, 0, 0}'])
    initial = mech['initial conditions']
    array('double', 'initial_species', [number(initial.get(s['name'], 0)) for s in species])
    array('double', 'initial_parameters', [number(initial.get(n, 0)) for n in params])
    for name, terms in [('reactants', rs), ('products', ps)]:
        array('ChemTerm', name, ['{%d, %s}' % (i, number(v)) for i, v in terms] or ['{0, 0}'])
    rows = []
    for kind, rb, rc, pb, pc, idx, values in reactions:
        rows.append('{CHEM_%s, %d, %d, %d, %d, %d, {%s}}' % (kind, rb, rc, pb, pc, idx, ', '.join(map(number, values))))
    array('ChemReaction', 'reactions', rows or ['{0}'])
    out.append(f'const ChemMechanism {symbol} = {{{len(species)}, {len(reactions)}, {len(params)}, species_names, parameter_names, initial_species, initial_parameters, reactions, reactants, products, {len(parameterized)}, parameterized_species}};')
    return '\n'.join(out) + '\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--symbol', default='generated_mechanism')
    args = parser.parse_args()
    if not args.symbol.isidentifier() or not args.symbol.isascii():
        parser.error('symbol must be an ASCII C identifier')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(emit(read_mechanism(args.config), args.symbol))

if __name__ == '__main__':
    main()
