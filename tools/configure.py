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

    def visit(path, inherited_units='molecule cm-3'):
        path = Path(path).resolve()
        if path in active:
            raise ValueError(f'cyclic camp-files inclusion: {path}')
        active.add(path)
        data = json.loads(path.read_text()) if path.suffix.lower() == '.json' else yaml.safe_load(path.read_text())
        if not isinstance(data, dict):
            raise ValueError(f'{path}: configuration must be an object')
        units = data.get('concentration units',inherited_units)
        if units not in ('molecule cm-3','mol m-3'):
            raise ValueError(f'unsupported concentration units: {units}')
        for child in data.get('camp-files', []):
            visit(path.parent / child,units)
        for record in data.get('camp-data', []):
            records.append(dict(record, **({'concentration units':record.get('concentration units',units)} if record.get('type') == 'MECHANISM' else {})))
        initial.update(data.get('initial state', {}))
        active.remove(path)

    visit(path)
    species = [dict(r) for r in records if r['type'] == 'CHEM_SPEC']
    known = {r['name'] for r in species}
    if len(known) != len(species) or not species or any(not isinstance(name,str) or not name or '\0' in name for name in known):
        raise ValueError('species names must be unique and nonempty')
    for r in records:
        if r['type'] not in ('CHEM_SPEC', 'MECHANISM'):
            raise ValueError(f"unsupported record type: {r['type']}")
    reactions = []
    for mechanism in (r for r in records if r['type'] == 'MECHANISM'):
        for raw in mechanism['reactions']:
            r = dict(raw)
            r['concentration units'] = mechanism['concentration units']
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
            else:
                products = [dict(p,**{'yield':p.get('yield',1.0)}) for p in products]
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
    conversion = 1e-6*6.02214076e23
    for r in reactions:
        if r['concentration units'] == 'mol m-3':
            continue
        order,kind = len(r['reactants']),r['type']
        if kind in ('TROE','TERNARY_CHEMICAL_ACTIVATION','TERNARY_CHEMICAL_ACTIVATION_JPL19'):
            low_order = order-1 if kind == 'TERNARY_CHEMICAL_ACTIVATION' else order
            r['k0_A'] = float(r.get('k0_A',1))*conversion**low_order
            r['kinf_A'] = float(r.get('kinf_A',1))*conversion**(low_order-1)
            if kind == 'TERNARY_CHEMICAL_ACTIVATION_JPL19': r['kint_A'] = float(r.get('kint_A',1))*conversion**(order-1)
        elif kind not in ('SURFACE','USER_DEFINED'):
            key = 'X' if kind == 'BRANCHED' else 'A'
            r[key] = float(r.get(key,1))*conversion**(order-1)
    reactions.sort(key=lambda r: KINDS.index(r['type']))
    return {'species': species, 'reactions': reactions, 'initial': initial}


def literal(value, numeric=False):
    if isinstance(value, str) and not numeric:
        if '\0' in value:
            raise ValueError('C names cannot contain a NUL character')
        # Emit UTF-8 bytes using only the basic C source character set. Three-digit
        # octal escapes avoid both surrogate escapes and hex-escape continuation.
        return '"'+''.join('\\'+chr(b) if chr(b) in '\\"?' else chr(b) if 32 <= b < 127 else f'\\{b:03o}' for b in value.encode('utf-8'))+'"'
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('configuration values must be finite')
    return repr(value)


def markowitz(n, pattern):
    pattern = set(pattern)
    alive = set(range(n))
    order = []
    while alive:
        def cost(v):
            rows = sum((v,j) in pattern for j in alive)
            columns = sum((j,v) in pattern for j in alive)
            return (rows-1)*(columns-1), v
        pivot = min(alive, key=cost)
        order.append(pivot)
        alive.remove(pivot)
        columns = [j for j in alive if (pivot,j) in pattern]
        rows = [j for j in alive if (j,pivot) in pattern]
        pattern.update((i,j) for i in rows for j in columns)
    return order


def mechanism_layout(config):
    """The reference variable permutation and row-major symbolic LU fill."""
    species, reactions = config['species'], config['reactions']
    vmap = {s['name']: i for i, s in enumerate(s for s in species if not s.get('parameterization'))}
    n = len(vmap)
    pattern = {(i, i) for i in range(n)}
    for r in reactions:
        for col in (vmap[s] for s in r['reactants'] if s in vmap):
            for row in (vmap[s] for s in r['reactants']+[p['species'] for p in r['products']] if s in vmap):
                pattern.add((row, col))
    order = markowitz(n, pattern)
    inverse = {old: new for new, old in enumerate(order)}
    vmap = {key: inverse[value] for key, value in vmap.items()}
    pattern = {(inverse[row],inverse[col]) for row,col in pattern}
    for i in range(n):
        for k in range(i+1, n):
            if (i, k) in pattern:
                for j in range(i+1, n):
                    if (j, i) in pattern:
                        pattern.add((j, k))
    indices = {pair: i for i, pair in enumerate(sorted(pattern))}
    return vmap, indices


def generate(config, name, kernel_chunk=32, fuse_kernels=False, cell_stride=0,
             emit_loader=True, dispatch_stride=0, initialize=False):
    if kernel_chunk < 1:
        raise ValueError('kernel chunk must be positive')
    species, reactions, initial = config['species'], config['reactions'], config['initial']
    ids = {s['name']: i for i, s in enumerate(species)}
    vmap, indices = mechanism_layout(config)
    n = len(vmap)
    lines = ['/* Generated by tools/configure.py; edit the configuration instead. */', '#include "chem.h"', '#include <stdlib.h>']
    groups = []
    for family in ('force', 'jacobian'):
        blocks = []
        written = set()

        def update(index, sign, expression):
            target = f'out[{index}*s+c]'
            if initialize and index not in written:
                statement = f'{target} = 0.0 {sign} {expression};'
            else:
                statement = f'{target} {sign}= {expression};'
            written.add(index)
            return statement
        if family == 'force':
            for ri, r in enumerate(reactions):
                rs = [vmap[s] for s in r['reactants'] if s in vmap]
                body = [f'double rate = k[{ri}*s+c];'] + [f'rate *= y[{v}*s+c];' for v in rs]
                body += [update(v,'-','rate') for v in rs]
                body += [update(vmap[p['species']],'+',f'{literal(p["yield"],True)}*rate') for p in r['products'] if p['species'] in vmap]
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
                        body += [update(indices[v,col],'+','rate') for v in others+[col]]
                        body += [update(indices[vmap[p['species']],col],'-',f'{literal(p["yield"],True)}*rate') for p in r['products'] if p['species'] in vmap]
                        blocks.append(body)
        funcs = []
        for start in range(0, len(blocks), kernel_chunk):
            fname = f'{name}_{family}_{start//kernel_chunk}'
            funcs.append(fname)
            lines.append(f'static void {fname}(size_t nc, size_t s, const double *restrict y, const double *restrict k, double *restrict out) {{')
            if cell_stride:
                lines[-1] = lines[-1].replace('size_t nc, size_t s,', 'size_t unused_nc, size_t unused_s,')
                lines += [f'  const size_t nc = {cell_stride}, s = {cell_stride};', '  (void)unused_nc; (void)unused_s;']
            lines += ['  (void)y; (void)k; (void)out;']
            if fuse_kernels:
                lines += ['  #pragma omp simd', '  for (size_t c = 0; c < nc; ++c) {']
                for block in blocks[start:start+kernel_chunk]:
                    lines += ['    {'] + ['      '+stmt for stmt in block] + ['    }']
                lines += ['  }']
            else:
                for block in blocks[start:start+kernel_chunk]:
                    lines += ['  #pragma omp simd', '  for (size_t c = 0; c < nc; ++c) {'] + ['    '+stmt for stmt in block] + ['  }']
            lines.append('}')
        lines += [f'static void {name}_{family}(const ChemModel *m, size_t nc, size_t s, const double *y, const double *k, double *out) {{',
                  '  (void)m; (void)nc; (void)s; (void)y; (void)k; (void)out;']
        if dispatch_stride:
            lines += [f'  if (nc == {dispatch_stride} && s == {dispatch_stride}) {{',
                      f'    {name}_fixed_{family}(m,nc,s,y,k,out); return;', '  }']
        if initialize:
            untouched = sorted(set(range(n if family == 'force' else len(indices))) - written)
            if untouched:
                lines += [f'  static const size_t zero[] = {{{", ".join(map(str,untouched))}}};',
                          f'  for (size_t i = 0; i < {len(untouched)}; ++i) {{',
                          '    double *restrict target = out+zero[i]*s;',
                          '    #pragma omp simd', '    for (size_t c = 0; c < nc; ++c) target[c] = 0.0;', '  }']
        lines += [f'  {fname}(nc, s, y, k, out);' for fname in funcs] + ['}']
        groups.append(f'{name}_{family}')
    if not emit_loader:
        return '\n'.join(lines) + '\n'
    lines += [f'ChemModel *chem_load_{name}(int specialized) {{', '  ChemModel *m = chem_model_create();', '  if (!m) return NULL;']
    for s in species:
        param = s.get('parameterization')
        if param not in (None, 'third body'):
            raise ValueError(f'unknown parameterization: {param}')
        init = initial.get('concentrations', {}).get(s['name'], 0.0)
        lines.append(f'  if (chem_model_species(m, {literal(s["name"])}, {literal(init,True)}, {int(bool(param))}, NULL)) goto fail;')
        for source, target in [('molecular weight [kg mol-1]', 'molecular_weight'), ('diffusion coefficient [m2 s-1]', 'diffusion'), ('absolute tolerance', 'absolute_tolerance')]:
            if source in s:
                lines.append(f'  m->species[{ids[s["name"]]}].{target} = {literal(s[source],True)};')
    rename = {'MUSICA name': 'label', 'scaling factor': 'A', 'reaction probability': 'probability', 'k_r': 'B'}
    fields = {'A', 'B', 'C', 'D', 'E', 'X', 'Y', 'a0', 'n', 'Fc', 'N', 'branch', 'k_r', 'MUSICA name', 'scaling factor', 'reaction probability'}
    fields |= {f'{kind}_{letter}' for kind in ('k0', 'kinf', 'kint') for letter in 'ABCD'}
    for r in reactions:
        kind = C_KINDS[KINDS.index(r['type'])]
        lines += ['  {', f'    ChemRate rate = chem_rate(CHEM_{kind});']
        for key in sorted(fields & r.keys()):
            value = r[key]
            lines.append(f'    rate.{rename.get(key,key)} = {literal(value,key != "MUSICA name")};')
        if kind == 'TAYLOR':
            coefficients = r.get('coefficients',[1.0])
            if not isinstance(coefficients,list) or not 1 <= len(coefficients) <= 16:
                raise ValueError('TAYLOR_SERIES requires 1 to 16 coefficients')
            lines.append(f'    rate.ncoefficients = {len(coefficients)};')
            lines += [f'    rate.coefficients[{i}] = {literal(value,True)};' for i,value in enumerate(coefficients)]
        if kind == 'SURFACE':
            sp = next(s for s in species if s['name'] == r['gas-phase reactant'])
            for source, target in [('molecular weight [kg mol-1]', 'molecular_weight'), ('diffusion coefficient [m2 s-1]', 'diffusion')]:
                if source not in sp:
                    raise ValueError(f'surface reactant missing {source}')
                lines.append(f'    rate.{target} = {literal(sp[source],True)};')
        rs = ', '.join(str(ids[s]) for s in r['reactants']) or '0'
        ps = ', '.join(str(ids[p['species']]) for p in r['products']) or '0'
        ys = ', '.join(literal(p['yield'],True) for p in r['products']) or '0'
        lines += [f'    const size_t reactants[] = {{{rs}}}, products[] = {{{ps}}};', f'    const double yields[] = {{{ys}}};',
                  f'    if (chem_model_reaction(m, reactants, {len(r["reactants"])}, products, yields, {len(r["products"])}, rate)) goto fail;', '  }']
    lines.append('  if (chem_model_reorder(m) || chem_model_finalize(m)) goto fail;')
    lines.append(f'  if (m->nnz != {len(indices)}) goto fail;')
    for key, field in [('temperature','default_temperature'), ('pressure','default_pressure')]:
        if key in initial:
            lines.append(f'  m->{field} = {literal(initial[key],True)};')
    for label, value in initial.get('rate parameters', {}).items():
        lines += ['  {', f'    size_t i = chem_custom_id(m, {literal(label)});', '    if (i == SIZE_MAX) goto fail;', f'    m->custom_defaults[i] = {literal(value,True)};', '  }']
    lines += ['  if (specialized) {', f'    m->force_kernel = {groups[0]}; m->jacobian_kernel = {groups[1]};', '  }', '  return m;', 'fail:', '  chem_model_free(m); return NULL;', '}']
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('output')
    parser.add_argument('--name', required=True)
    parser.add_argument('--kernel-chunk',type=int,default=32)
    args = parser.parse_args()
    if not re.fullmatch('[a-zA-Z_][a-zA-Z_0-9]*', args.name):
        parser.error('name must be a C identifier')
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(generate(read_configuration(args.config), args.name,args.kernel_chunk))
