"""Read CAMP YAML/JSON and emit immutable C mechanism data (no runtime parser)."""
import argparse
import copy
import json
import math
from pathlib import Path
import yaml


class ConfigurationError(ValueError):
    pass


def read_documents(path, visiting=None):
    path = Path(path).resolve()
    visiting = set() if visiting is None else visiting
    if path in visiting:
        raise ConfigurationError(f'configuration include cycle at {path}')
    visiting.add(path)
    doc = json.loads(path.read_text()) if path.suffix == '.json' else yaml.safe_load(path.read_text())
    if not isinstance(doc, dict):
        raise ConfigurationError(f'{path}: expected a mapping')
    data = list(doc.get('camp-data', []))
    for child in doc.get('camp-files', []):
        data.extend(read_documents(path.parent / child, visiting))
    visiting.remove(path)
    return data


def reactants(value):
    if isinstance(value, list):
        return value.copy()
    result = []
    for name, options in value.items():
        qty = options.get('qty', 1)
        if not isinstance(qty, int) or qty < 1:
            raise ConfigurationError(f'{name}: reactant qty must be a positive integer')
        result.extend([name] * qty)
    return result


def products(value):
    if isinstance(value, list):
        return [(v['species'], float(v.get('yield', 1))) for v in value]
    return [(name, float(v.get('yield', 1))) for name, v in value.items()]


RATE_PARAMETERS = {
    'ARRHENIUS': [('A', 1), ('B', 0), ('C', 0), ('D', 300), ('E', 0)],
    'TROE': [('k0_A', 1), ('k0_B', 0), ('k0_C', 0), ('kinf_A', 1), ('kinf_B', 0), ('kinf_C', 0), ('Fc', .6), ('N', 1)],
    'TUNNELING': [('A', 1), ('B', 0), ('C', 0)],
    'BRANCHED_ALKOXY': [('X', 1), ('Y', 0), ('a0', .5), ('n', 0)],
    'REVERSIBLE': [('A', 1), ('C', 0), ('k_r', 0)],
    'USER': [('scaling_factor', 1)],
}
RATE_PARAMETERS['TERNARY'] = RATE_PARAMETERS['TROE']
RATE_PARAMETERS['BRANCHED_NITRATE'] = RATE_PARAMETERS['BRANCHED_ALKOXY']
RATE_PARAMETERS['TAYLOR'] = RATE_PARAMETERS['ARRHENIUS']
RATE_PARAMETERS['TERNARY_JPL19'] = [(f'{s}_{k}', d) for s in ('k0', 'kinf', 'kint')
                                    for k, d in (('A', 1), ('B', 0), ('C', 0), ('D', 298))] + [('Fc', .6), ('N', 1)]


def normalize(path):
    data = read_documents(path)
    species = [s for s in data if s.get('type') == 'CHEM_SPEC']
    names = [s['name'] for s in species]
    if len(set(names)) != len(names) or not names:
        raise ConfigurationError('species must have unique names and must not be empty')
    by_name = {s['name']: s for s in species}
    dynamic = [s for s in species if not s.get('third body', False)]
    ids = {s['name']: i for i, s in enumerate(dynamic)}
    parameters, initial_parameters, reactions = [], [], []
    conditions = {'temperature': 300, 'pressure': 101325}
    initial = {}
    raw_reactions = []
    for entry in data:
        if entry.get('type') == 'MECHANISM':
            conditions.update(entry.get('initial conditions', {}))
            initial.update(entry.get('initial rate parameters', {}))
            raw_reactions.extend(entry.get('reactions', []))
        elif entry.get('type') != 'CHEM_SPEC':
            raise ConfigurationError(f'unsupported entry type {entry.get("type")}')
    for raw in raw_reactions:
        r = copy.deepcopy(raw)
        if r['type'] == 'BRANCHED':
            for branch in ('alkoxy', 'nitrate'):
                b = dict(r, type=f'BRANCHED_{branch.upper()}', products=r[f'{branch} products'])
                raw_reactions.append(b)
            continue
        kind = {'PHOTOLYSIS': 'USER', 'USER_DEFINED': 'USER', 'EMISSION': 'USER',
                'FIRST_ORDER_LOSS': 'USER', 'TERNARY_CHEMICAL_ACTIVATION': 'TERNARY',
                'TERNARY_CHEMICAL_ACTIVATION_JPL19': 'TERNARY_JPL19', 'TAYLOR_SERIES': 'TAYLOR'}.get(r['type'], r['type'])
        # Emission is zeroth order; loss is first order.
        if r['type'] == 'EMISSION':
            r.update(reactants={}, products={r['species']: {}})
        if r['type'] == 'FIRST_ORDER_LOSS':
            r.update(reactants={r['species']: {}}, products={})
        if kind == 'SURFACE':
            r.update(reactants={r['gas-phase reactant']: {}}, products=r['gas-phase products'])
        rr, pp = reactants(r.get('reactants', {})), products(r.get('products', {}))
        if any(n not in by_name for n in rr + [n for n, _ in pp]):
            raise ConfigurationError('reaction refers to unknown species')
        third_body_count = sum(n not in ids for n in rr)
        param = -1
        coeffs = []
        if kind == 'SURFACE':
            s = by_name[r['gas-phase reactant']]
            try:
                values = [s['diffusion coefficient [m2 s-1]'], s['molecular weight [kg mol-1]'], r.get('reaction probability', 1)]
            except KeyError as exc:
                raise ConfigurationError(f'surface species missing property {exc}') from exc
            labels = [r['MUSICA name'] + '.effective radius [m]',
                      r['MUSICA name'] + '.particle number concentration [# m-3]']
        else:
            if kind not in RATE_PARAMETERS:
                raise ConfigurationError(f'unsupported rate type {kind}')
            values = [r.get(k, d) for k, d in RATE_PARAMETERS[kind]]
            labels = [r['MUSICA name']] if kind == 'USER' else []
            if kind == 'TAYLOR':
                coeffs = r.get('coefficients', [])
                if len(coeffs) > 19:
                    raise ConfigurationError('at most 19 Taylor coefficients supported')
                values += coeffs
        if labels:
            param = len(parameters)
            for label in labels:
                if label in parameters:
                    raise ConfigurationError(f'duplicate custom rate label {label}')
                parameters.append(label)
                initial_parameters.append(float(initial.get(label, 0)))
        reactions.append(dict(reactants=[ids[n] for n in rr if n in ids],
                              products=[(ids[n], v) for n, v in pp if n in ids],
                              type=kind, values=values, parameter=param,
                              coefficient_count=len(coeffs), third_body_count=third_body_count))
    return dict(species=dynamic, reactions=reactions, parameters=parameters,
                initial_parameters=initial_parameters, conditions=conditions)


def number(v):
    f = float(v)
    if not math.isfinite(f):
        raise ConfigurationError('nonfinite numeric input')
    return repr(f)


def emit_c(mech, symbol):
    lines = ['/* Generated mechanism data; regenerate with tools/configure.py. */', '#include "chem.h"']
    def array(ctype, name, values):
        lines.append(f'static const {ctype} {name}[] = {{' + ','.join(values or ['0']) + '};')
    array('char *const', 'species_names', [json.dumps(s['name']) for s in mech['species']])
    array('char *const', 'parameter_names', [json.dumps(n) for n in mech['parameters']])
    array('double', 'initial_y', [number(s.get('initial concentration', 0)) for s in mech['species']])
    array('double', 'atol', [number(s.get('absolute tolerance', 1e-3)) for s in mech['species']])
    array('double', 'initial_p', [number(v) for v in mech['initial_parameters']])
    for i, r in enumerate(mech['reactions']):
        array('int', f'reactants_{i}', [str(v) for v in r['reactants']])
        array('ChemProduct', f'products_{i}', ['{' + str(n) + ',' + number(v) + '}' for n, v in r['products']])
    lines.append('static const ChemReaction reactions[] = {')
    for i, r in enumerate(mech['reactions']):
        lines.append('{%d,%d,reactants_%d,products_%d,{.type=CHEM_%s,.p={%s},.parameter=%d,.coefficient_count=%d,.third_body_count=%d}},' % (
            len(r['reactants']), len(r['products']), i, i, r['type'], ','.join(number(v) for v in r['values']),
            r['parameter'], r['coefficient_count'], r['third_body_count']))
    lines.append('};')
    lines.append(f'const ChemMechanismData {symbol} = {{' + ','.join([
        str(len(mech['species'])), str(len(mech['reactions'])), str(len(mech['parameters'])),
        'species_names', 'parameter_names', 'initial_y', 'initial_p', 'atol',
        number(mech['conditions']['temperature']), number(mech['conditions']['pressure']), 'reactions']) + '};')
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('config')
    p.add_argument('output')
    p.add_argument('--symbol', default='mechanism_data')
    a = p.parse_args()
    if not a.symbol.isidentifier() or not a.symbol.isascii():
        p.error('symbol must be an ASCII C identifier')
    Path(a.output).write_text(emit_c(normalize(a.config), a.symbol))
