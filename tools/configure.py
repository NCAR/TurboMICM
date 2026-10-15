"""Read CAMP-style YAML or JSON and generate a mechanism for the C library."""
import argparse
import json
import math
import re
from pathlib import Path

import yaml

KINDS = ['ARRHENIUS', 'TROE', 'TERNARY_CHEMICAL_ACTIVATION', 'TERNARY_CHEMICAL_ACTIVATION_JPL19',
         'BRANCHED', 'TUNNELING', 'TAYLOR_SERIES', 'REVERSIBLE', 'USER_DEFINED', 'SURFACE', 'LAMBDA']
C_KINDS = ['ARRHENIUS', 'TROE', 'TERNARY', 'JPL19', 'BRANCHED', 'TUNNELING', 'TAYLOR', 'REVERSIBLE', 'USER', 'SURFACE', 'LAMBDA']


def read_configuration(path):
    records, initial = [], {}
    active = set()

    def visit(path):
        path = Path(path).resolve()
        if path in active:
            raise ValueError(f'cyclic camp-files inclusion: {path}')
        active.add(path)
        data = json.loads(path.read_text()) if path.suffix.lower() == '.json' else yaml.safe_load(path.read_text())
        if not isinstance(data, dict):
            raise ValueError(f'{path}: configuration must be an object')
        for child in data.get('camp-files', []):
            visit(path.parent / child)
        records.extend(data.get('camp-data', []))
        initial.update(data.get('initial state', {}))
        active.remove(path)

    visit(path)
    species = [dict(r) for r in records if r['type'] == 'CHEM_SPEC']
    known = {r['name'] for r in species}
    if len(known) != len(species) or not species:
        raise ValueError('species names must be unique and nonempty')
    for r in records:
        if r['type'] not in ('CHEM_SPEC', 'MECHANISM'):
            raise ValueError(f"unsupported record type: {r['type']}")
    reactions = []
    for mechanism in (r for r in records if r['type'] == 'MECHANISM'):
        for raw in mechanism['reactions']:
            r = dict(raw)
            kind = r['type']
            reactants = r.get('reactant sequence')
            if reactants is None:
                reactants = []
                for name, details in r.get('reactants', {}).items():
                    quantity = details.get('qty', 1)
                    if not isinstance(quantity, int) or quantity < 1:
                        raise ValueError('reactant qty must be a positive integer')
                    reactants.extend([name] * quantity)
            products = r.get('product sequence')
            if products is None:
                products = [{'species': name, 'yield': details.get('yield', 1.0)} for name, details in r.get('products', {}).items()]
            if kind == 'BRANCHED':
                for branch, key in enumerate(['alkoxy products', 'nitrate products']):
                    if key not in r:
                        raise ValueError('BRANCHED requires both product branches')
                    ps = [{'species': name, 'yield': d.get('yield', 1.0)} for name, d in r[key].items()]
                    reactions.append(dict(r, reactants=reactants, products=ps, branch=branch))
                continue
            if kind == 'SURFACE':
                reactants = [r['gas-phase reactant']]
                if 'gas-phase products' in r:
                    products = [{'species': name, 'yield': d.get('yield', 1.0)} for name, d in r['gas-phase products'].items()]
            if kind == 'EMISSION':
                reactants, products = [], [{'species': r['species'], 'yield': 1.0}]
            if kind == 'FIRST_ORDER_LOSS':
                reactants, products = [r['species']], []
            if kind in ('PHOTOLYSIS', 'EMISSION', 'FIRST_ORDER_LOSS'):
                r['type'] = 'USER_DEFINED'
            if r['type'] not in KINDS or r['type'] == 'LAMBDA':
                raise ValueError(f'unsupported configuration rate type: {kind}')
            reactions.append(dict(r, reactants=reactants, products=products))
    for r in reactions:
        for name in r['reactants'] + [p['species'] for p in r['products']]:
            if name not in known:
                raise ValueError(f'unknown species {name}')
    reactions.sort(key=lambda r: KINDS.index(r['type']))
    return {'species': species, 'reactions': reactions, 'initial': initial}


def literal(value):
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=True)
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('configuration values must be finite')
    return repr(value)


def generate(config, name):
    species, reactions, initial = config['species'], config['reactions'], config['initial']
    ids = {s['name']: i for i, s in enumerate(species)}
    vmap = {s['name']: i for i, s in enumerate(s for s in species if not s.get('parameterization'))}
    n = len(vmap)
    pattern = {(i, i) for i in range(n)}
    for r in reactions:
        for col in (vmap[s] for s in r['reactants'] if s in vmap):
            for row in (vmap[s] for s in r['reactants']+[p['species'] for p in r['products']] if s in vmap):
                pattern.add((row, col))
    for i in range(n):
        for k in range(i+1, n):
            if (i, k) in pattern:
                for j in range(i+1, n):
                    if (j, i) in pattern:
                        pattern.add((j, k))
    indices = {pair: i for i, pair in enumerate(sorted(pattern))}
    lines = ['/* Generated by tools/configure.py; edit the configuration instead. */', '#include "chem.h"', '#include <stdlib.h>']
    groups = []
    for family in ('force', 'jacobian'):
        blocks = []
        if family == 'force':
            for ri, r in enumerate(reactions):
                rs = [vmap[s] for s in r['reactants'] if s in vmap]
                body = [f'double rate = k[{ri}*s+c];'] + [f'rate *= y[{v}*s+c];' for v in rs]
                body += [f'out[{v}*s+c] -= rate;' for v in rs]
                body += [f'out[{vmap[p["species"]]}*s+c] += {literal(p["yield"])}*rate;' for p in r['products'] if p['species'] in vmap]
                blocks.append(body)
        else:
            for col in range(n):
                for ri, r in enumerate(reactions):
                    rs = [vmap[s] for s in r['reactants'] if s in vmap]
                    for omitted, variable in enumerate(rs):
                        if variable != col:
                            continue
                        others = rs[:omitted]+rs[omitted+1:]
                        body = [f'double rate = k[{ri}*s+c];'] + [f'rate *= y[{v}*s+c];' for v in others]
                        body += [f'out[{indices[v,col]}*s+c] += rate;' for v in others+[col]]
                        body += [f'out[{indices[vmap[p["species"]],col]}*s+c] -= {literal(p["yield"])}*rate;' for p in r['products'] if p['species'] in vmap]
                        blocks.append(body)
        funcs = []
        for start in range(0, len(blocks), 32):
            fname = f'{name}_{family}_{start//32}'
            funcs.append(fname)
            lines.append(f'static void {fname}(size_t nc, size_t s, const double *restrict y, const double *restrict k, double *restrict out) {{')
            for block in blocks[start:start+32]:
                lines += ['  #pragma omp simd', '  for (size_t c = 0; c < nc; ++c) {'] + ['    '+stmt for stmt in block] + ['  }']
            lines.append('}')
        lines += [f'static void {name}_{family}(const ChemModel *m, size_t nc, size_t s, const double *y, const double *k, double *out) {{', '  (void)m;']
        lines += [f'  {fname}(nc, s, y, k, out);' for fname in funcs] + ['}']
        groups.append(f'{name}_{family}')
    lines += [f'ChemModel *chem_load_{name}(int specialized) {{', '  ChemModel *m = chem_model_create();', '  if (!m) return NULL;']
    for s in species:
        param = s.get('parameterization')
        if param not in (None, 'third body'):
            raise ValueError(f'unknown parameterization: {param}')
        init = initial.get('concentrations', {}).get(s['name'], 0.0)
        lines.append(f'  if (chem_model_species(m, {literal(s["name"])}, {literal(init)}, {int(bool(param))}, NULL)) goto fail;')
        for source, target in [('molecular weight [kg mol-1]', 'molecular_weight'), ('diffusion coefficient [m2 s-1]', 'diffusion'), ('absolute tolerance', 'absolute_tolerance')]:
            if source in s:
                lines.append(f'  m->species[{ids[s["name"]]}].{target} = {literal(s[source])};')
    rename = {'MUSICA name': 'label', 'scaling factor': 'A', 'reaction probability': 'probability', 'k_r': 'B'}
    fields = {'A', 'B', 'C', 'D', 'E', 'X', 'Y', 'a0', 'n', 'Fc', 'N', 'branch', 'k_r', 'MUSICA name', 'scaling factor', 'reaction probability'}
    fields |= {f'{kind}_{letter}' for kind in ('k0', 'kinf', 'kint') for letter in 'ABCD'}
    for r in reactions:
        kind = C_KINDS[KINDS.index(r['type'])]
        lines += ['  {', f'    ChemRate rate = chem_rate(CHEM_{kind});']
        for key in sorted(fields & r.keys()):
            value = r[key]
            lines.append(f'    rate.{rename.get(key,key)} = {literal(value)};')
        if kind == 'SURFACE':
            sp = next(s for s in species if s['name'] == r['gas-phase reactant'])
            for source, target in [('molecular weight [kg mol-1]', 'molecular_weight'), ('diffusion coefficient [m2 s-1]', 'diffusion')]:
                if source not in sp:
                    raise ValueError(f'surface reactant missing {source}')
                lines.append(f'    rate.{target} = {literal(sp[source])};')
        rs = ', '.join(str(ids[s]) for s in r['reactants']) or '0'
        ps = ', '.join(str(ids[p['species']]) for p in r['products']) or '0'
        ys = ', '.join(literal(p['yield']) for p in r['products']) or '0'
        lines += [f'    const size_t reactants[] = {{{rs}}}, products[] = {{{ps}}};', f'    const double yields[] = {{{ys}}};',
                  f'    if (chem_model_reaction(m, reactants, {len(r["reactants"])}, products, yields, {len(r["products"])}, rate)) goto fail;', '  }']
    lines.append('  if (chem_model_finalize(m)) goto fail;')
    lines.append(f'  if (m->nnz != {len(indices)}) goto fail;')
    for key, field in [('temperature','default_temperature'), ('pressure','default_pressure')]:
        if key in initial:
            lines.append(f'  m->{field} = {literal(initial[key])};')
    for label, value in initial.get('rate parameters', {}).items():
        lines += ['  {', f'    size_t i = chem_custom_id(m, {literal(label)});', '    if (i == SIZE_MAX) goto fail;', f'    m->custom_defaults[i] = {literal(value)};', '  }']
    lines += ['  if (specialized) {', f'    m->force_kernel = {groups[0]}; m->jacobian_kernel = {groups[1]};', '  }', '  return m;', 'fail:', '  chem_model_free(m); return NULL;', '}']
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('output')
    parser.add_argument('--name', required=True)
    args = parser.parse_args()
    if not re.fullmatch('[a-zA-Z_][a-zA-Z_0-9]*', args.name):
        parser.error('name must be a C identifier')
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(generate(read_configuration(args.config), args.name))
