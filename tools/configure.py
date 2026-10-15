"""Read CAMP YAML/JSON configurations and generate immutable C descriptions.

Only serialization is delegated to PyYAML; schema, mechanism semantics, units,
validation, includes, and code generation are implemented here.
"""
from dataclasses import dataclass
import argparse
import json
import math
import pathlib
import re

import yaml

AVOGADRO_CM = 1e-6 * 6.02214076e23
DEFAULTS = dict(A=1, B=0, C=0, D=300, E=0, k0_A=1, k0_B=0, k0_C=0, k0_D=300,
                kinf_A=1, kinf_B=0, kinf_C=0, kinf_D=300, kint_A=0, kint_B=0, kint_C=0,
                kint_D=300, Fc=0.6, N=1, X=0, Y=0, a0=0, n=0, probability=1,
                diffusion=0, molecular_weight=0, scaling=1, reverse=0)
TYPES = {'ARRHENIUS': 'CK_ARRHENIUS', 'TROE': 'CK_TROE', 'TERNARY_CHEMICAL_ACTIVATION': 'CK_TERNARY',
         'TERNARY_CHEMICAL_ACTIVATION_JPL19': 'CK_TERNARY_JPL19', 'TUNNELING': 'CK_TUNNELING',
         'BRANCH_ALKOXY': 'CK_BRANCH_ALKOXY', 'BRANCH_NITRATE': 'CK_BRANCH_NITRATE',
         'TAYLOR_SERIES': 'CK_TAYLOR', 'REVERSIBLE': 'CK_REVERSIBLE', 'USER_DEFINED': 'CK_USER',
         'PHOTOLYSIS': 'CK_USER', 'FIRST_ORDER_LOSS': 'CK_USER', 'EMISSION': 'CK_USER', 'SURFACE': 'CK_SURFACE'}


class ConfigurationError(ValueError):
    pass


def finite(value):
    if isinstance(value, bool):
        raise ConfigurationError('Boolean used as number')
    result = float(value)
    if not math.isfinite(result):
        raise ConfigurationError('Non-finite number')
    return result


def load(path, stack=()):
    path = pathlib.Path(path).resolve()
    if path in stack:
        raise ConfigurationError('Cyclic camp-files include')
    with path.open() as stream:
        data = json.load(stream) if path.suffix.lower() == '.json' else yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ConfigurationError(f'{path}: expected mapping')
    entries = list(data.get('camp-data', []))
    for child in data.get('camp-files', []):
        included = load(path.parent / child, (*stack, path))
        entries.extend(included['camp-data'])
    return {**data, 'camp-data': entries}


@dataclass
class Mechanism:
    name: str
    species: list
    reactions: list
    initial: dict


def read(path):
    data = load(path)
    entries = data['camp-data']
    species = [dict(e) for e in entries if e.get('type') == 'CHEM_SPEC']
    indices = {s['name']: i for i, s in enumerate(species)}
    if not species or len(indices) != len(species):
        raise ConfigurationError('Species must be nonempty and have unique names')
    units = data.get('concentration units', 'molecule cm-3')
    if units not in ('mol m-3', 'molecule cm-3'):
        raise ConfigurationError('Unknown concentration units')
    factor = AVOGADRO_CM if units == 'molecule cm-3' else 1.0
    for s in species:
        if s.get('parameterization') not in (None, 'air-density'):
            raise ConfigurationError('Unknown species parameterization')
    def reactants(value):
        result = []
        if isinstance(value, dict):
            for name, options in value.items():
                qty = options.get('qty', 1)
                if not isinstance(qty, int) or qty < 1:
                    raise ConfigurationError('Reactant qty must be a positive integer')
                result.extend([indices[name]] * qty)
        elif isinstance(value, list):
            result = [indices[name] for name in value]
        else:
            raise ConfigurationError('Reactants must be mapping or ordered list')
        return result
    def products(value):
        if isinstance(value, dict):
            return [(indices[name], finite(options.get('yield', 1))) for name, options in value.items()]
        if isinstance(value, list):
            return [(indices[p['species']], finite(p.get('yield', 1))) for p in value]
        raise ConfigurationError('Products must be mapping or ordered list')
    result = []
    names = []
    for entry in entries:
        if entry.get('type') != 'MECHANISM':
            continue
        names.append(entry.get('name', 'mechanism'))
        for raw in entry.get('reactions', []):
            raw = dict(raw)
            kind = raw['type']
            if kind == 'BRANCHED':
                branches = [('BRANCH_ALKOXY', raw['alkoxy products']), ('BRANCH_NITRATE', raw['nitrate products'])]
            else:
                branches = [(kind, raw.get('products', {}))]
            for kind, prods in branches:
                if kind not in TYPES:
                    raise ConfigurationError('Unsupported rate type ' + kind)
                react = raw.get('reactants', {})
                if kind == 'SURFACE':
                    react = {raw['gas-phase reactant']: {}}
                    prods = raw.get('gas-phase products', prods)
                if kind == 'FIRST_ORDER_LOSS':
                    react = {raw['species']: {}}
                if kind == 'EMISSION':
                    prods = {raw['species']: {}}
                rr, pp = reactants(react), products(prods)
                order = len(rr)
                rate = DEFAULTS.copy()
                for key in DEFAULTS:
                    if key in raw:
                        rate[key] = finite(raw[key])
                rate['scaling'] = finite(raw.get('scaling factor', raw.get('scaling_factor', 1)))
                rate['reverse'] = finite(raw.get('k_r', raw.get('reverse', 0)))
                rate['probability'] = finite(raw.get('reaction probability', 1))
                if kind in ('ARRHENIUS', 'TUNNELING', 'TAYLOR_SERIES'):
                    rate['A'] *= factor ** (order - 1)
                if kind in ('BRANCH_ALKOXY', 'BRANCH_NITRATE'):
                    rate['X'] *= factor ** (order - 1)
                if kind in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION', 'TERNARY_CHEMICAL_ACTIVATION_JPL19'):
                    offset = 0 if kind == 'TROE' or kind.endswith('JPL19') else -1
                    rate['k0_A'] *= factor ** (order + offset)
                    rate['kinf_A'] *= factor ** (order - 1 + offset)
                    rate['kint_A'] *= factor ** (order - 1 + offset)
                if TYPES[kind] == 'CK_USER':
                    rate['scaling'] *= factor ** (order - 1)
                if kind == 'SURFACE':
                    s = species[rr[0]]
                    rate['diffusion'] = finite(s['diffusion coefficient [m2 s-1]'])
                    rate['molecular_weight'] = finite(s['molecular weight [kg mol-1]'])
                    if rate['diffusion'] <= 0 or rate['molecular_weight'] <= 0:
                        raise ConfigurationError('Surface reaction requires positive diffusion and molecular weight')
                rate['type'] = TYPES[kind]
                rate['label'] = raw.get('MUSICA name', raw.get('label'))
                if TYPES[kind] in ('CK_USER', 'CK_SURFACE') and not isinstance(rate['label'], str):
                    raise ConfigurationError('User-defined/surface reaction needs a label')
                rate['coefficients'] = [finite(v) for v in raw.get('coefficients', [])]
                if len(rate['coefficients']) > 16:
                    raise ConfigurationError('At most 16 Taylor coefficients')
                result.append(dict(reactants=rr, products=pp, rate=rate))
    if not result:
        raise ConfigurationError('No reactions')
    return Mechanism(' + '.join(names), species, result, data.get('initial state', {}))


def number(value):
    return repr(finite(value))


def emit(m, symbol):
    if not re.fullmatch('[A-Za-z_][A-Za-z_0-9]*', symbol):
        raise ConfigurationError('Invalid C symbol')
    lines = ['/* Generated by tools/configure.py. Do not edit. */', '#include "chem.h"']
    def add(text):
        lines.append(text)
    concentrations = m.initial.get('concentrations', {})
    add(f'static const ck_species {symbol}_species[] = {{')
    for s in m.species:
        add('  {.name=' + json.dumps(s['name']) + ', .molecular_weight=' + number(s.get('molecular weight [kg mol-1]', 0)) +
            ', .diffusion=' + number(s.get('diffusion coefficient [m2 s-1]', 0)) +
            ', .absolute_tolerance=' + number(s.get('absolute tolerance', 1e-3)) +
            ', .initial=' + number(concentrations.get(s['name'], 0)) +
            ', .parameterized=' + ('1' if s.get('parameterization') else '0') + '},')
    add('};')
    for i, r in enumerate(m.reactions):
        if r['reactants']:
            add(f'static const size_t {symbol}_reactants_{i}[] = {{' + ','.join(map(str, r['reactants'])) + '};')
        if r['products']:
            add(f'static const ck_product {symbol}_products_{i}[] = {{' + ','.join('{' + str(s) + ',' + number(y) + '}' for s, y in r['products']) + '};')
    add(f'static const ck_reaction {symbol}_reactions[] = {{')
    for i, r in enumerate(m.reactions):
        rate = r['rate']
        fields = [f'.{key}={number(rate[key])}' for key in DEFAULTS]
        fields += ['.type=' + rate['type']]
        if rate['label'] is not None:
            fields += ['.label=' + json.dumps(rate['label'])]
        if rate['coefficients']:
            fields += ['.ncoefficients=' + str(len(rate['coefficients'])), '.coefficients={' + ','.join(map(number, rate['coefficients'])) + '}']
        add('  {.nreactants=' + str(len(r['reactants'])) + ', .reactants=' + (f'{symbol}_reactants_{i}' if r['reactants'] else 'NULL') +
            ', .nproducts=' + str(len(r['products'])) + ', .products=' + (f'{symbol}_products_{i}' if r['products'] else 'NULL') +
            ', .rate={' + ','.join(fields) + '}},')
    add('};')
    custom = m.initial.get('custom', {})
    if custom:
        add(f'static const ck_named_value {symbol}_custom[] = {{' + ','.join('{' + json.dumps(k) + ',' + number(v) + '}' for k, v in custom.items()) + '};')
    add(f'const ck_mechanism {symbol} = {{.name=' + json.dumps(m.name) + f', .nspecies={len(m.species)}, .nreactions={len(m.reactions)},' +
        f' .species={symbol}_species, .reactions={symbol}_reactions, .temperature=' + number(m.initial.get('temperature', 0)) +
        ', .pressure=' + number(m.initial.get('pressure', 0)) +
        (f', .ncustom_defaults={len(custom)}, .custom_defaults={symbol}_custom' if custom else '') + '};')
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    parser.add_argument('--symbol', default='mechanism')
    args = parser.parse_args()
    mechanism = read(args.configuration)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(emit(mechanism, args.symbol))


if __name__ == '__main__':
    main()
